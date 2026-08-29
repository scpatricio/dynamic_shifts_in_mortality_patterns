prepare_series <- function(data, code, sex, age_0 = age_min) {
  observations <- data |>
    dplyr::filter(
      CNTRY == code,
      Sex == sex,
      Year >= cohort_min,
      Year <= cohort_max,
      Age >= age_0,
      Ex > 0,
      is.finite(Ex),
      !is.na(Dx)
    ) |>
    dplyr::mutate(Dx = as.integer(round(Dx))) |>
    dplyr::arrange(Year, Age)

  if (!nrow(observations)) {
    stop("No observations for ", code, " / ", sex, call. = FALSE)
  }

  years <- sort(unique(observations$Year))
  stan_data <- list(
    Y = length(years),
    N = nrow(observations),
    yr = as.integer(match(observations$Year, years)),
    Dx = observations$Dx,
    Ex = as.numeric(observations$Ex),
    t = as.numeric(observations$Age - age_0)
  )

  list(observations = observations, years = years, stan_data = stan_data)
}

sample_series <- function(stan_model, prepared, code, sex, model_name) {
  message("Fitting ", model_name, ": ", code, " / ", sex)
  fit <- rstan::sampling(
    object = stan_model,
    data = prepared$stan_data,
    chains = stan_chains,
    iter = stan_warmup + stan_sampling,
    warmup = stan_warmup,
    thin = 1L,
    cores = stan_cores,
    seed = series_seed(code, sex),
    refresh = 100L,
    control = list(adapt_delta = 0.95, max_treedepth = 12L)
  )

  country <- unique(prepared$observations$Country)
  if (length(country) != 1L) country <- code

  list(
    fit = fit,
    model = model_name,
    code = code,
    country = standardize_country_name(country),
    sex = sex,
    years = prepared$years,
    age_0 = age_min,
    observations = prepared$observations,
    stan_data = prepared$stan_data
  )
}

fit_all_series <- function(data, stan_file, output_directory, model_name) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  rstan::rstan_options(auto_write = TRUE)
  compiled_model <- rstan::stan_model(file = stan_file)

  for (code in country_codes) {
    for (sex in sexes) {
      output_file <- file.path(output_directory, fit_filename(code, sex))
      if (file.exists(output_file) && !overwrite_fits) {
        message("Skipping existing fit: ", basename(output_file))
        next
      }
      prepared <- prepare_series(data, code, sex)
      fit_object <- sample_series(compiled_model, prepared, code, sex, model_name)
      saveRDS(fit_object, output_file)
    }
  }
  invisible(NULL)
}
