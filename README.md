# Bayesian gamma–Gompertz state-space models for mortality landmark estimation

Silvio C. Patricio

Interdisciplinary Center on Population Dynamics, University of Southern Denmark

Code and reproducibility materials for estimating the modal age at death,
mortality deceleration, plateau onset, and plateau level across birth cohorts.

**Related preprint:** [arXiv:2608.29452](https://doi.org/10.48550/arXiv.2608.29452)

The analysis uses cohort mortality and exposure data from the [Human Mortality Database](https://www.mortality.org/) for 12 low-mortality countries, for females and males. It estimates changes in the modal age at death, mortality deceleration, and plateau onset using a Bayesian gamma–Gompertz state-space model. A spike-and-slab analysis assesses how much evidence the data provide for mortality deceleration or a plateau.

## Repository structure

```text
.
├── analysis/                  # Scripts, run in numerical order
├── R/                         # Shared configuration and functions
├── stan/                      # Main and spike-and-slab models
├── data/
│   ├── raw/                   # Authenticated HMD download (not committed)
│   └── derived/               # Analysis-ready data (not committed)
├── results/
│   ├── fits/                  # Large Stan fit objects (not committed)
│   └── summaries/             # Compact derived results
└── figures/
    ├── diagnostics/           # Population-specific model checks
    └── manuscript/            # Main and supplementary figures
```

## Requirements

R 4.1 or later and the following packages:

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "ggplot2", "patchwork", "here",
  "rstan", "HDInterval", "modeest", "matrixStats", "HMDHFDplus"
))
```

A working C++ toolchain is also required by `rstan`.

## Data access

The script downloads data directly from the HMD (https://www.mortality.org),
which requires a free account. Credentials are read from environment
variables rather than hard-coded:

```bash
export HMD_USERNAME="you@example.com"
export HMD_PASSWORD="your_password"
```

Do not commit or share HMD credentials. Users without an account should
register directly with the Human Mortality Database.

## Run the analysis

Run all commands from the repository root:

```bash
Rscript analysis/01_download_hmd_data.R
Rscript analysis/02_fit_main_model.R
Rscript analysis/03_summarize_landmarks.R
Rscript analysis/04_model_diagnostics.R
Rscript analysis/05_fit_spike_slab.R
Rscript analysis/06_plot_spike_slab.R
Rscript analysis/07_make_figures.R
```

The scripts must be run in numerical order. Stan defaults match the manuscript: four chains, 4,000 warm-up iterations, and 2,000 retained iterations per chain. For a short test run, these settings can be overridden:

```bash
STAN_CHAINS=1 STAN_WARMUP=100 STAN_SAMPLING=100 Rscript analysis/02_fit_main_model.R
```

Set `OVERWRITE=true` to replace existing model fits.

## Main outputs

- `results/summaries/mortality_landmarks.rds`: posterior summaries of the modal age at death, mortality deceleration, plateau onset, equivalent plateau hazard, and survival to the landmarks.
- `results/summaries/spike_slab_summary.rds`: cohort-specific posterior probability of a pure Gompertz tail.
- `figures/manuscript/`: figures used by the manuscript.
- `figures/diagnostics/`: posterior predictive checks and parameter trajectories.

## Data note

HMD data are downloaded by the user under the HMD terms of use. The download script saves a local raw copy and a compact analysis-ready file. Death counts are reconstructed as `round(Ex * Mx)`, following the manuscript analysis.

## License

The code is released under the MIT License. HMD data remain subject to the
Human Mortality Database terms of use and are not redistributed in this
repository or the Zenodo software archive.
