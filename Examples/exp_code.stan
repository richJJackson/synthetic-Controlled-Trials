//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector<lower=0>[N] time;
  int<lower=0,upper=1> status[N];
  vector[N] x;  // covariate: treatment
}

parameters {
  real<lower=0> lambda0;
  real beta;
}

model {
  // Priors
  lambda0 ~ gamma(1,1);
  beta ~ normal(0,10);

  // Likelihood with right censoring
  for (i in 1:N) {
    real lambda_i;
    lambda_i = lambda0 * exp(beta * x[i]);

    if (status[i] == 1)
      target += exponential_lpdf(time[i] | lambda_i);
    else
      target += exponential_lccdf(time[i] | lambda_i);
  }
}

