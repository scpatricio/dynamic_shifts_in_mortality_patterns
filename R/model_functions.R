gamma_gompertz_hazard <- function(age, a, b, gamma, age_0 = age_min) {
  elapsed <- age - age_0
  growth <- exp(b * elapsed)
  a * growth / (1 + gamma * a * (growth - 1) / b)
}

gamma_gompertz_survival <- function(age, a, b, gamma, age_0 = age_min) {
  elapsed <- age - age_0
  denominator <- 1 + gamma * a * (exp(b * elapsed) - 1) / b
  denominator^(-1 / gamma)
}

mortality_landmarks <- function(a, b, gamma, age_0 = age_min) {
  ratio <- (b - gamma * a) / (gamma * a)
  valid <- is.finite(ratio) & ratio > 0 & a > 0 & b > 0 & gamma > 0

  deceleration <- plateau <- rep(NA_real_, length(ratio))
  deceleration[valid] <- age_0 + log((2 - sqrt(3)) * ratio[valid]) / b[valid]
  plateau[valid] <- age_0 + log(ratio[valid]) / b[valid]

  data.frame(deceleration = deceleration, plateau = plateau)
}

modal_age_at_death <- function(a, b, gamma, age_0 = age_min) {
  numerator <- b - gamma * a
  valid <- is.finite(numerator) & numerator > 0 & a > 0 & b > 0
  mode_age <- rep(NA_real_, length(a))
  mode_age[valid] <- age_0 + log(numerator[valid] / a[valid]) / b[valid]
  pmax(mode_age, age_0)
}

equivalent_plateau_hazard <- function(b, gamma, series_terms = 80L) {
  valid <- is.finite(b) & is.finite(gamma) & b > 0 & gamma > 0
  result <- rep(NA_real_, length(b))
  if (!any(valid)) return(result)

  p <- 1 / gamma[valid]
  n <- 0:series_terms
  weights <- 0.5^n
  reciprocal_terms <- 1 / outer(p, n, `+`)
  integral_constant <- rowSums(sweep(reciprocal_terms, 2L, weights, `*`))
  result[valid] <- b[valid] / integral_constant
  result
}

summarize_main_fit <- function(fit_object) {
  draws <- rstan::extract(fit_object$fit, pars = c("a", "b", "gamma"), permuted = TRUE)
  years <- fit_object$years

  summaries <- lapply(seq_along(years), function(index) {
    a <- draws$a[, index]
    b <- draws$b[, index]
    gamma <- draws$gamma[, index]

    landmarks <- mortality_landmarks(a, b, gamma, fit_object$age_0)
    mode_age <- modal_age_at_death(a, b, gamma, fit_object$age_0)
    equivalent_hazard <- equivalent_plateau_hazard(b, gamma)
    survival_dec <- gamma_gompertz_survival(
      landmarks$deceleration, a, b, gamma, fit_object$age_0
    )
    survival_plateau <- gamma_gompertz_survival(
      landmarks$plateau, a, b, gamma, fit_object$age_0
    )

    values <- c(
      posterior_summary(mode_age, "mad"),
      posterior_summary(landmarks$deceleration, "dec"),
      posterior_summary(landmarks$plateau, "onset"),
      posterior_summary(equivalent_hazard, "ehaz"),
      posterior_summary(survival_dec, "prob_dec"),
      posterior_summary(survival_plateau, "prob_onset")
    )

    data.frame(
      country = fit_object$country,
      code = fit_object$code,
      sex = fit_object$sex,
      year = years[index],
      as.list(values),
      check.names = FALSE
    )
  })

  dplyr::bind_rows(summaries)
}
