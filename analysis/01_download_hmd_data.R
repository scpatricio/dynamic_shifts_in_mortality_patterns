suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(here)
  library(HMDHFDplus)
})

source(here::here("R", "config.R"))
source(here::here("R", "helpers.R"))
hmd_username <- "hmd_username"
hmd_password <- "hmd_password"

if (!nzchar(hmd_username) || !nzchar(hmd_password)) {
  stop(
    "Set HMD_USERNAME and HMD_PASSWORD in ~/.Renviron before downloading data.",
    call. = FALSE
  )
}

raw_directory <- here::here("data", "raw")
derived_directory <- here::here("data", "derived")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_directory, recursive = TRUE, showWarnings = FALSE)

country_lookup <- HMDHFDplus::getHMDcountries() |>
  dplyr::filter(CNTRY %in% country_codes) |>
  dplyr::select(CNTRY, Country)

missing_codes <- setdiff(country_codes, country_lookup$CNTRY)
if (length(missing_codes)) {
  stop("HMD country codes not found: ", paste(missing_codes, collapse = ", "))
}

download_item <- function(item, value_name) {
  purrr::map_dfr(country_codes, function(code) {
    message("Downloading ", item, " for ", code)
    downloaded <- HMDHFDplus::readHMDweb(
      CNTRY = code,
      item = item,
      username = hmd_username,
      password = hmd_password
    )

    downloaded |>
      dplyr::mutate(CNTRY = code) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(sexes),
        names_to = "Sex",
        values_to = value_name
      ) |>
      dplyr::left_join(country_lookup, by = "CNTRY") |>
      dplyr::select(
        dplyr::any_of(c(
          "Year", "Age", "OpenInterval", "CNTRY", "Country", "Sex", value_name
        ))
      )
  })
}

cohort_mx <- download_item("cMx_1x1", "Mx")
cohort_exposure <- download_item("cExposures_1x1", "Ex")

saveRDS(
  list(mx = cohort_mx, exposure = cohort_exposure),
  file.path(raw_directory, "hmd_cohort_download.rds")
)

join_keys <- intersect(
  c("Year", "Age", "OpenInterval", "CNTRY", "Country", "Sex"),
  intersect(names(cohort_mx), names(cohort_exposure))
)

cohort_data <- cohort_exposure |>
  dplyr::inner_join(cohort_mx, by = join_keys) |>
  dplyr::filter(
    Year >= cohort_min,
    Year <= cohort_max,
    Sex %in% sexes,
    is.finite(Ex),
    is.finite(Mx),
    Ex > 0,
    Mx >= 0
  ) |>
  dplyr::mutate(
    Country = standardize_country_name(Country),
    Dx = as.integer(round(Ex * Mx))
  ) |>
  dplyr::arrange(CNTRY, Sex, Year, Age)

saveRDS(cohort_data, file.path(derived_directory, "hmd_cohort_data.rds"))

coverage <- cohort_data |>
  dplyr::group_by(CNTRY, Country, Sex) |>
  dplyr::summarise(
    first_cohort = min(Year),
    last_cohort = max(Year),
    maximum_age = max(Age),
    .groups = "drop"
  )

utils::write.csv(
  coverage,
  file.path(derived_directory, "hmd_cohort_coverage.csv"),
  row.names = FALSE
)

message("Saved analysis-ready data to data/derived/hmd_cohort_data.rds")
