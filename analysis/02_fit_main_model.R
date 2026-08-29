suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(rstan)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "fitting.R"))

data_file <- here::here("data", "derived", "hmd_cohort_data.rds")
stan_file <- here::here("stan", "gamma_gompertz_rw.stan")
output_directory <- here::here("results", "fits", "main")

assert_file(data_file)
assert_file(stan_file)

cohort_data <- readRDS(data_file)
fit_all_series(
  data = cohort_data,
  stan_file = stan_file,
  output_directory = output_directory,
  model_name = "gamma_gompertz_rw"
)
