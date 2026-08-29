suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(here)
  library(rstan)
  library(HDInterval)
  library(modeest)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))

fit_directory <- here::here("results", "fits", "main")
summary_directory <- here::here("results", "summaries")
figure_directory <- here::here("figures", "diagnostics")
dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

summarize_diagnostic_fit <- function(fit_object) {
  observations <- fit_object$observations
  draws <- rstan::extract(
    fit_object$fit,
    pars = c("Dx_rep", "a", "b", "gamma"),
    permuted = TRUE
  )

  thresholds <- sort(unique(observations$Dx))
  replicated_cdf <- t(vapply(
    seq_len(nrow(draws$Dx_rep)),
    function(index) {
      findInterval(thresholds, sort(draws$Dx_rep[index, ])) / nrow(observations)
    },
    numeric(length(thresholds))
  ))
  cdf_interval <- t(apply(replicated_cdf, 2L, HDInterval::hdi, credMass = 0.95))

  pp <- data.frame(
    country = fit_object$country,
    code = fit_object$code,
    sex = fit_object$sex,
    observed_cdf = stats::ecdf(observations$Dx)(thresholds),
    predictive_cdf = colMeans(replicated_cdf),
    predictive_cdf_low = cdf_interval[, 1L],
    predictive_cdf_hi = cdf_interval[, 2L]
  )

  parameters <- lapply(seq_along(fit_object$years), function(index) {
    values <- c(
      posterior_summary(draws$a[, index], "a"),
      posterior_summary(draws$b[, index], "b"),
      posterior_summary(draws$gamma[, index], "gamma")
    )
    data.frame(
      country = fit_object$country,
      code = fit_object$code,
      sex = fit_object$sex,
      year = fit_object$years[index],
      as.list(values),
      check.names = FALSE
    )
  }) |>
    dplyr::bind_rows()

  list(pp = pp, parameters = parameters)
}

fits <- load_fit_objects(fit_directory)
diagnostics <- lapply(fits, summarize_diagnostic_fit)
pp_all <- lapply(diagnostics, `[[`, "pp") |> dplyr::bind_rows()
parameters_all <- lapply(diagnostics, `[[`, "parameters") |> dplyr::bind_rows()

saveRDS(
  list(pp = pp_all, parameters = parameters_all),
  file.path(summary_directory, "model_diagnostics.rds")
)

parameter_plot <- function(data, parameter, subtitle) {
  ggplot(data, aes(x = year, colour = sex, fill = sex)) +
    geom_ribbon(
      aes(
        ymin = .data[[paste0(parameter, "_HDI_low")]],
        ymax = .data[[paste0(parameter, "_HDI_hi")]]
      ),
      alpha = 0.25,
      colour = NA
    ) +
    geom_line(aes(y = .data[[paste0(parameter, "_mode")]]), linewidth = 0.7) +
    scale_y_log10() +
    scale_colour_manual(values = sex_colours) +
    scale_fill_manual(values = sex_colours) +
    labs(x = "Cohort", y = "Value (log scale)", subtitle = subtitle) +
    theme_minimal()
}

for (country_name in sort(unique(parameters_all$country))) {
  parameter_data <- dplyr::filter(parameters_all, country == country_name)
  pp_data <- dplyr::filter(pp_all, country == country_name)

  pp_plot <- ggplot(
    pp_data,
    aes(x = observed_cdf, y = predictive_cdf, colour = sex, fill = sex)
  ) +
    geom_ribbon(
      aes(ymin = predictive_cdf_low, ymax = predictive_cdf_hi),
      alpha = 0.25,
      colour = NA
    ) +
    geom_line(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_colour_manual(values = sex_colours) +
    scale_fill_manual(values = sex_colours) +
    labs(x = "Empirical CDF", y = "Posterior predictive CDF", subtitle = "P--P plot") +
    theme_minimal()

  combined <-
    parameter_plot(parameter_data, "a", expression("Baseline mortality (" * italic(a) * ")")) +
    parameter_plot(parameter_data, "b", expression("Gompertz slope (" * italic(b) * ")")) +
    parameter_plot(parameter_data, "gamma", expression("Frailty variance (" * gamma * ")")) +
    pp_plot +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = country_name) &
    theme(legend.position = "bottom")

  ggsave(
    filename = file.path(
      figure_directory,
      paste0("fit_", country_figure_slug(country_name), ".pdf")
    ),
    plot = combined,
    width = 9,
    height = 6,
    device = cairo_pdf
  )
}

message("Saved model diagnostics.")
