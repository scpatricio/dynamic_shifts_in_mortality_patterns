suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(here)
  library(rstan)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))

summary_file <- here::here("results", "summaries", "mortality_landmarks.rds")
data_file <- here::here("data", "derived", "hmd_cohort_data.rds")
fit_directory <- here::here("results", "fits", "main")
figure_directory <- here::here("figures", "manuscript")
summary_directory <- here::here("results", "summaries")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

assert_file(summary_file)
assert_file(data_file)

landmarks <- readRDS(summary_file)
cohort_data <- readRDS(data_file)

base_theme <- theme_minimal() +
  theme(legend.position = "bottom")

markers_long <- dplyr::bind_rows(
  dplyr::transmute(
    landmarks, country, sex, year, marker = "Modal age at death",
    age = mad_mode, low = mad_HDI_low, high = mad_HDI_hi
  ),
  dplyr::transmute(
    landmarks, country, sex, year, marker = "Mortality deceleration",
    age = dec_mode, low = dec_HDI_low, high = dec_HDI_hi
  ),
  dplyr::transmute(
    landmarks, country, sex, year, marker = "Plateau onset",
    age = onset_mode, low = onset_HDI_low, high = onset_HDI_hi
  )
) |>
  dplyr::mutate(
    marker = factor(
      marker,
      levels = c("Modal age at death", "Mortality deceleration", "Plateau onset")
    )
  )

schedule_plot <- ggplot(
  markers_long,
  aes(
    x = year,
    y = age,
    colour = sex,
    fill = sex,
    linetype = marker,
    group = interaction(sex, marker)
  )
) +
  geom_ribbon(
    aes(ymin = low, ymax = high),
    alpha = 0.10,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  scale_fill_manual(values = sex_colours) +
  scale_linetype_manual(values = c("solid", "dashed", "12")) +
  labs(x = "Cohort (birth year)", y = "Age (years)", colour = NULL, linetype = NULL) +
  base_theme

ggsave(
  file.path(figure_directory, "plot_1_new.pdf"),
  schedule_plot,
  width = 10,
  height = 7,
  device = cairo_pdf
)

gap_data <- landmarks |>
  dplyr::mutate(
    mode_to_deceleration = dec_mode - mad_mode,
    deceleration_to_plateau = onset_mode - dec_mode
  )

gap_one <- ggplot(
  gap_data,
  aes(x = year, y = mode_to_deceleration, colour = sex)
) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  labs(x = NULL, y = "Gap (years)", subtitle = "Mode to mortality deceleration") +
  base_theme

gap_two <- ggplot(
  gap_data,
  aes(x = year, y = deceleration_to_plateau, colour = sex)
) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  labs(x = "Cohort (birth year)", y = "Gap (years)", subtitle = "Deceleration to plateau onset") +
  base_theme

gap_plot <- gap_one / gap_two +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_directory, "gap_landmarks.pdf"),
  gap_plot,
  width = 10,
  height = 11,
  device = cairo_pdf
)

level_plot <- ggplot(
  landmarks,
  aes(x = year, y = ehaz_mode, colour = sex, fill = sex)
) +
  geom_ribbon(aes(ymin = ehaz_HDI_low, ymax = ehaz_HDI_hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  scale_fill_manual(values = sex_colours) +
  scale_y_log10() +
  labs(x = "Cohort (birth year)", y = "Equivalent plateau hazard", colour = NULL, fill = NULL) +
  base_theme

ggsave(
  file.path(figure_directory, "plot_2_new.pdf"),
  level_plot,
  width = 10,
  height = 7,
  device = cairo_pdf
)

survival_plot <- ggplot(
  landmarks,
  aes(x = year, y = prob_dec, colour = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = prob_dec_HDI_low, ymax = prob_dec_HDI_hi),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ country, ncol = 4) +
  scale_colour_manual(values = sex_colours) +
  scale_fill_manual(values = sex_colours) +
  scale_y_log10() +
  labs(x = "Cohort (birth year)", y = "Survival from age 65 to deceleration", colour = NULL, fill = NULL) +
  base_theme

ggsave(
  file.path(figure_directory, "survival.pdf"),
  survival_plot,
  width = 10,
  height = 7,
  device = cairo_pdf
)

count_data <- cohort_data |>
  dplyr::filter(Age > 100, Sex %in% sexes) |>
  dplyr::mutate(
    death_count = cut(
      Dx,
      breaks = c(-Inf, 0, 9, 99, Inf),
      labels = c("0", "1--9", "10--99", "100+"),
      right = TRUE
    )
  )

count_plot <- ggplot(count_data, aes(x = Year, y = Age, fill = death_count)) +
  geom_tile() +
  facet_grid(sex ~ country) +
  scale_fill_viridis_d(option = "C", na.value = "white", drop = FALSE) +
  labs(x = "Cohort (birth year)", y = "Age", fill = "Deaths") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_directory, "death_heatmap_SI.pdf"),
  count_plot,
  width = 12,
  height = 8,
  device = cairo_pdf
)

selected_countries <- c("Australia", "France", "Japan", "Sweden")
fits <- load_fit_objects(fit_directory)
hyperparameter_draws <- lapply(fits, function(fit_object) {
  if (!fit_object$country %in% selected_countries) return(NULL)
  draws <- rstan::extract(
    fit_object$fit,
    pars = c("tau_a", "tau_b", "tau_gamma", "phi"),
    permuted = TRUE
  )
  data.frame(
    country = fit_object$country,
    sex = fit_object$sex,
    tau_a = draws$tau_a,
    tau_b = draws$tau_b,
    tau_gamma = draws$tau_gamma,
    phi = draws$phi
  ) |>
    tidyr::pivot_longer(
      cols = c(tau_a, tau_b, tau_gamma, phi),
      names_to = "parameter",
      values_to = "value"
    )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    parameter = factor(
      parameter,
      levels = c("tau_a", "tau_b", "tau_gamma", "phi"),
      labels = c(expression(tau[a]), expression(tau[b]), expression(tau[gamma]), expression(phi))
    )
  )

hyperparameter_plot <- ggplot(
  hyperparameter_draws,
  aes(x = value, colour = country)
) +
  geom_density(linewidth = 0.6) +
  facet_grid(sex ~ parameter, scales = "free") +
  labs(x = "Posterior value", y = "Density", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_directory, "posterior_plot.pdf"),
  hyperparameter_plot,
  width = 10,
  height = 5,
  device = cairo_pdf
)

mean_changes <- landmarks |>
  dplyr::group_by(country, sex) |>
  dplyr::arrange(year, .by_group = TRUE) |>
  dplyr::summarise(
    modal_age = mean(diff(mad_mode)),
    deceleration = mean(diff(dec_mode)),
    plateau_onset = mean(diff(onset_mode)),
    .groups = "drop"
  )

utils::write.csv(
  mean_changes,
  file.path(summary_directory, "mean_annual_landmark_changes.csv"),
  row.names = FALSE
)

message("Saved manuscript figures and summary table.")
