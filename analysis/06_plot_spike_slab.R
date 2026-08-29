suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(rstan)
  library(HDInterval)
  library(modeest)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))

fit_directory <- here::here("results", "fits", "spike_slab")
summary_directory <- here::here("results", "summaries")
figure_directory <- here::here("figures", "manuscript")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

fits <- load_fit_objects(fit_directory)

summarize_spike_fit <- function(fit_object) {
  draws <- rstan::extract(
    fit_object$fit,
    pars = c("p_no_plateau", "w"),
    permuted = TRUE
  )

  summaries <- lapply(seq_along(fit_object$years), function(index) {
    probability <- draws$p_no_plateau[, index]
    interval <- HDInterval::hdi(probability, credMass = 0.95)
    data.frame(
      country = fit_object$country,
      code = fit_object$code,
      sex = fit_object$sex,
      cohort = fit_object$years[index],
      p_no_plateau = mean(probability),
      p_no_plateau_sd = stats::sd(probability),
      p_no_plateau_HDI_low = interval[1L],
      p_no_plateau_HDI_hi = interval[2L],
      prior_weight = mean(draws$w)
    )
  })
  dplyr::bind_rows(summaries)
}

spike_summary <- lapply(fits, summarize_spike_fit) |>
  dplyr::bind_rows() |>
  dplyr::arrange(country, sex, cohort)

write_summary_files(
  spike_summary,
  stem = "spike_slab_summary",
  directory = summary_directory
)

plot_data <- spike_summary |>
  dplyr::mutate(
    p_no_plateau_plot = pmax(p_no_plateau, 1e-8),
    prior_weight_plot = pmax(prior_weight, 1e-8)
  )

spike_plot <- ggplot(
  plot_data,
  aes(x = cohort, y = p_no_plateau_plot, colour = sex, group = sex)
) +
  geom_line(linewidth = 0.7) +
  geom_hline(
    aes(yintercept = prior_weight_plot, colour = sex),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  scale_y_log10() +
  labs(
    x = "Cohort (birth year)",
    y = expression(Pr(gamma == 0 ~ "|" ~ data)),
    colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(figure_directory, "spike_and_slab.pdf"),
  plot = spike_plot,
  width = 10,
  height = 7,
  device = cairo_pdf
)

message("Saved spike-and-slab summary and figure.")
