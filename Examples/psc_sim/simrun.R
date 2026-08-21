### Run simulations

###
setwd("~/Documents/GitHub/synthetic-Controlled-Trials/Examples/psc_sim")

# Start (or first time)
source("run_all.R")

# --- stop with Ctrl+C when you need to leave ---
#     safest: after "checkpoint saved" or between scenarios

# Check where you are (no computation)
source("show_progress.R")

# Next day: resume
source("run_resume.R")
