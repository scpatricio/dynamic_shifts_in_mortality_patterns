suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(rstan)
  library(HDInterval)
  library(modeest)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "model_functions.R"))

fit_directory <- here::here("results", "fits", "main")
summary_directory <- here::here("results", "summaries")

fits <- load_fit_objects(fit_directory)
landmark_summary <- lapply(fits, summarize_main_fit) |>
  dplyr::bind_rows() |>
  dplyr::arrange(country, sex, year)

write_summary_files(
  landmark_summary,
  stem = "mortality_landmarks",
  directory = summary_directory
)

message("Saved mortality landmark summaries.")
