assert_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  invisible(path)
}

posterior_mode <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  estimate <- modeest::mlv(x, na.rm = TRUE)
  if (is.list(estimate)) estimate <- estimate$M
  as.numeric(estimate)[1L]
}

posterior_summary <- function(x, prefix) {
  x <- x[is.finite(x)]
  suffixes <- c("", "_sd", "_mode", "_HDI_low", "_HDI_hi")
  if (!length(x)) {
    return(stats::setNames(rep(NA_real_, length(suffixes)), paste0(prefix, suffixes)))
  }
  interval <- HDInterval::hdi(x, credMass = 0.95, na.rm = TRUE)
  values <- c(mean(x), stats::sd(x), posterior_mode(x), interval[1L], interval[2L])
  stats::setNames(values, paste0(prefix, suffixes))
}

standardize_country_name <- function(x) {
  dplyr::recode(
    x,
    "England and Wales (Total Population)" = "England and Wales",
    "New Zealand Non-Maori" = "New Zealand (Non-Maori)",
    .default = x
  )
}

country_figure_slug <- function(country) {
  country <- standardize_country_name(country)
  if (country == "New Zealand (Non-Maori)") return("New_Zealand")
  slug <- gsub("[^A-Za-z0-9]+", "_", country)
  gsub("^_|_$", "", slug)
}

series_seed <- function(code, sex) {
  code_index <- match(code, country_codes)
  sex_index <- match(sex, sexes)
  base_seed + 10L * code_index + sex_index
}

fit_filename <- function(code, sex) {
  paste0(code, "_", tolower(sex), ".rds")
}

load_fit_objects <- function(directory) {
  files <- sort(list.files(directory, pattern = "\\.rds$", full.names = TRUE))
  if (!length(files)) {
    stop("No fit files found in ", directory, call. = FALSE)
  }
  stats::setNames(lapply(files, readRDS), basename(files))
}

write_summary_files <- function(object, stem, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file.path(directory, paste0(stem, ".rds")))
  utils::write.csv(object, file.path(directory, paste0(stem, ".csv")), row.names = FALSE)
}
