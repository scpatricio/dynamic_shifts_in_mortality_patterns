country_codes <- c(
  "AUS", "CAN", "DNK", "GBRTENW", "NZL_NM", "FRATNP",
  "ITA", "NLD", "NOR", "SWE", "CHE", "JPN"
)

sexes <- c("Female", "Male")
age_min <- 65L
cohort_min <- 1840L
cohort_max <- 1910L
base_seed <- 42L

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  parsed
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

stan_chains <- env_integer("STAN_CHAINS", 4L)
stan_warmup <- env_integer("STAN_WARMUP", 4000L)
stan_sampling <- env_integer("STAN_SAMPLING", 2000L)
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.na(detected_cores) || detected_cores < 2L) {
  detected_cores <- 1L
} else {
  detected_cores <- detected_cores - 1L
}
stan_cores <- min(stan_chains, env_integer("STAN_CORES", detected_cores))
overwrite_fits <- env_flag("OVERWRITE", FALSE)

sex_colours <- c("Female" = "#e34a1c", "Male" = "#2f82b4")
