# Stan model code and cached compilation for Bayesian historical-control arms.
# Models are compiled once per session (saved under stan_cache/) and reused
# via sampling(), avoiding repeated temp-directory compilation per replicate.

.sim_stan_registry <- new.env(parent = emptyenv())

STAN_EXP_CODE <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  real lam0_mn;
  real lam0_t;
  int<lower=0,upper=1> status[N];
  vector[N] x;
}

parameters {
  real<lower=0> lambda0;
  real beta;
}

model {
  lambda0 ~ normal(lam0_mn, lam0_t);
  beta ~ normal(0, 10);

  for (i in 1:N) {
    real lambda_i;
    lambda_i = exp(-lambda0 + beta * x[i]);
    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }
}"

STAN_EXP_COMM_CODE <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  vector[N] x;
  real lam0_mn;
}

parameters {
  real lambda0;
  real beta;
  real<lower=0> tau;
}

model {
  lambda0 ~ normal(lam0_mn, 1 / sqrt(tau));
  beta ~ normal(0, 10);
  tau ~ gamma(1, 1);

  for (i in 1:N) {
    real lambda_i;
    lambda_i = exp(-lambda0 + beta * x[i]);
    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }
}"

STAN_EXP_CASEW_CODE <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  vector[N] x;
  int<lower=0> N0;
  vector<lower=0>[N0] timec;
  int<lower=0,upper=1> statusc[N0];
  vector[N0] a0;
}

parameters {
  real lambda0;
  real beta;
}

model {
  lambda0 ~ normal(0, 10);
  beta ~ normal(0, 10);

  for (i in 1:N) {
    real lambda_i;
    lambda_i = exp(-lambda0 + beta * x[i]);
    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }

  for (j in 1:N0) {
    real lambdac_j;
    lambdac_j = exp(-lambda0);
    if (statusc[j] == 1)
      target += a0[j] * exponential_lpdf(timec[j] | lambdac_j);
    else
      target += a0[j] * exponential_lccdf(timec[j] | lambdac_j);
  }
}"

STAN_PSC_FIXED_CODE <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] H0;
  vector<lower=0>[N] h0;
  vector[N] lp;
  array[N] int<lower=0, upper=1> status;
}

parameters {
  real beta;
}

model {
  beta ~ normal(0, 31.6228);

  for (i in 1:N) {
    real H;
    real h;
    real S;
    H = H0[i] * exp(lp[i] + beta);
    h = h0[i] * exp(lp[i] + beta);
    S = exp(-H);
    if (status[i] == 1) {
      target += log(S * h + 1e-16);
    } else {
      target += log(S + 1e-16);
    }
  }
}"

STAN_MODEL_SPECS <- list(
  exp = STAN_EXP_CODE,
  comm = STAN_EXP_COMM_CODE,
  casew = STAN_EXP_CASEW_CODE,
  psc_fixed = STAN_PSC_FIXED_CODE
)

sim_stan_cache_dir <- function(paths = sim_paths()) {
  file.path(paths$sim_dir, "stan_cache")
}

sim_parallel_cores <- function(defaults = SIM_DEFAULTS) {
  requested <- as.integer(defaults$parallel_cores %||% 1L)
  requested <- max(1L, requested)
  max_cores <- as.integer(defaults$max_parallel_cores %||% 4L)
  available <- parallel::detectCores(logical = TRUE)
  if (is.na(available) || available < 1L) {
    available <- 1L
  }
  min(requested, max_cores, available)
}

sim_init_stan <- function(paths = sim_paths(), recompile = FALSE, verbose = TRUE) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for Stan models.")
  }
  cache_dir <- sim_stan_cache_dir(paths)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  options(rstan_model_dir = cache_dir)
  rstan::rstan_options(auto_write = TRUE)
  options(rstan.pdf_plot = FALSE)

  for (key in names(STAN_MODEL_SPECS)) {
    code <- STAN_MODEL_SPECS[[key]]
    rds_path <- file.path(cache_dir, paste0(key, ".stanmodel.rds"))
    if (!recompile && file.exists(rds_path)) {
      mod <- readRDS(rds_path)
      if (verbose) {
        cat("Loaded cached Stan model:", key, "\n")
      }
    } else {
      if (verbose) {
        cat("Compiling Stan model:", key, " (one-time; cached in stan_cache/)\n")
      }
      mod <- rstan::stan_model(
        model_code = code,
        model_name = paste0("psc_sim_", key)
      )
      saveRDS(mod, rds_path)
    }
    assign(key, mod, envir = .sim_stan_registry)
  }

  invisible(TRUE)
}

sim_stan_is_initialized <- function() {
  all(vapply(names(STAN_MODEL_SPECS), exists, logical(1), envir = .sim_stan_registry))
}

sim_stan_fit <- function(model_key, data, defaults = SIM_DEFAULTS) {
  model_key <- match.arg(model_key, names(STAN_MODEL_SPECS))
  if (!sim_stan_is_initialized()) {
    stop("Stan models not initialised. Call sim_init_stan() after loading rstan.")
  }
  mod <- get(model_key, envir = .sim_stan_registry)
  rstan::sampling(
    mod,
    data = data,
    iter = defaults$stan_iter,
    chains = defaults$stan_chains,
    seed = defaults$stan_seed,
    refresh = 0
  )
}
