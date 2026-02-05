# Main pipeline for estimating GGM parameters and plateau features

library(foreach)
library(data.table)
library(randomForest)
library(ggplot2)
library(dplyr)
library(splines)
library(rstan)
library(doParallel)
library(bootstrap)
library(reshape2)

# mcprogress is used via mcprogress::pmclapply()
# devtools::install_github("myles-lewis/mcprogress")

source("scripts/plateau_functions.R")

# -----------------------------------------------------------
# Parallel + Stan configuration
# -----------------------------------------------------------

totalCores <- parallel::detectCores()
cluster    <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Compile Stan model once
stan_file <- file.path("stan", "plateau_model.stan")
sm        <- rstan::stan_model(file = stan_file)

# -----------------------------------------------------------
# Settings
# -----------------------------------------------------------

countries      <- c("DNK", "FRA", "ITA", "SWD")
countries_name <- c("Denmark", "France", "Italy", "Sweden")
names(countries_name) <- countries

age_0        <- 60
dir_0        <- getwd()
n_cores_outer <- max(1L, parallel::detectCores() - 2L)

result <- NULL

# -----------------------------------------------------------
# Main loop over countries and sexes
# -----------------------------------------------------------

for (country in countries) {
  country_name <- countries_name[[country]]
  
  # ------------ load cohort data -----------
  mx_file <- file.path(dir_0, "cohort", paste0("mx_", country, ".txt"))
  Ex_file <- file.path(dir_0, "cohort", paste0("Ex_", country, ".txt"))
  
  mx <- read.csv(mx_file, sep = "", skip = 1)
  Ex <- read.csv(Ex_file, sep = "", skip = 1)
  
  Ex <- melt(
    setDT(Ex),
    id.vars      = c("Year", "Age"),
    variable.name = "Sex",
    value.name    = "Ex"
  )
  mx <- melt(
    setDT(mx),
    id.vars      = c("Year", "Age"),
    variable.name = "Sex",
    value.name    = "mx"
  )
  
  df <- merge(as.data.frame(mx), as.data.frame(Ex),
              by = c("Year", "Age", "Sex"))
  
  df$Age[df$Age == "110+"] <- 110
  df$Age <- as.numeric(df$Age)
  
  df <- df |>
    dplyr::filter(mx != ".", Ex != ".") |>
    dplyr::mutate(
      Dx = as.integer(as.numeric(Ex) * as.numeric(mx)),
      Ex = as.numeric(Ex),
      mx = as.numeric(mx)
    )
  
  years_keep <- df |>
    dplyr::group_by(Year, Sex) |>
    dplyr::summarise(
      keep = min(Age) <= age_0 & max(Age) >= 105,
      .groups = "drop"
    )
  
  df <- dplyr::full_join(df, years_keep, by = dplyr::join_by(Year, Sex)) |>
    dplyr::filter(keep)
  
  sexes <- c("Male", "Female")
  years <- df |>
    dplyr::filter(Year >= 1850) |>
    dplyr::select(Year) |>
    dplyr::distinct() |>
    dplyr::pull()
  
  for (sex in sexes) {
    message("Processing: ", country_name, " - ", sex)
    
    res_list <- mcprogress::pmclapply(
      years,
      function(year) {
        # Filter data for this cohort/sex
        df_aux <- df |>
          dplyr::filter(
            Year == year,
            Sex  == sex,
            Dx   > 0,
            Ex   > 0,
            Age  >= age_0
          ) |>
          dplyr::mutate(
            log_mx = log(Dx / Ex)
          ) |>
          dplyr::arrange(Age)
        
        if (nrow(df_aux) == 0L) {
          return(NULL)
        }
        
        # --------------------------------------------------
        # Gamma-Gompertz-Makeham via Stan
        # --------------------------------------------------
        n_chain <- 2000L
        
        init_fun <- function() {
          list(
            a      = runif(1, 0.001, 0.1),
            b      = runif(1, 0.001, 0.1),
            c      = runif(1, 0.001, 0.1),
            sigma2 = runif(1, 0.1, 2)
          )
        }
        
        N <- nrow(df_aux)
        data_list <- list(
          N  = N,
          Dx = as.integer(df_aux$Dx),
          Ex = df_aux$Ex,
          t  = df_aux$Age - age_0 + runif(N, 0, 1)
        )
        
        fit <- sampling(
          sm,
          data    = data_list,
          chains  = 1,
          iter    = n_chain,
          warmup  = n_chain - 100,
          seed    = as.integer(year + ifelse(sex == "Male", 0, 1)),
          refresh = 0,
          init    = init_fun
        )
        
        chain_p <- rstan::extract(fit, pars = c("a", "b", "sigma2", "c"))
        
        # --------------------------------------------------
        # Mortality deceleration (AMD) via random forest
        # --------------------------------------------------
        
        df_aux_mort_dec <- df_aux |>
          dplyr::arrange(Age) |>
          dplyr::mutate(
            deriv_1 = mx - dplyr::lag(mx),
            deriv_2 = deriv_1 - dplyr::lag(deriv_1),
            Age_center = Age,
            Dx_roll3   = (Dx + dplyr::lag(Dx) + dplyr::lead(Dx)) / 3
          ) |>
          stats::na.omit() |>
          dplyr::transmute(
            Age    = Age,
            Dx     = Dx,
            deriv_2 = deriv_1
          )
        
        # bootstrap over ids
        B <- 1000L
        
        boot_amd_candidates <- replicate(
          n = B,
          expr = {
            idx <- sample(seq_len(nrow(df_aux_mort_dec)),
                          nrow(df_aux_mort_dec),
                          replace = TRUE)
            max_func(
              id_1            = idx,
              n_tree          = 100,
              data.frame.full = df_aux_mort_dec
            )
          }
        )
        
        boot_amd_candidates <- unlist(boot_amd_candidates)
        mort_dec_bootstrap  <- bootstrap(boot_amd_candidates, 100, mean)
        
        chain <- data.frame(
          country  = country_name,
          year     = year,
          sex      = sex,
          a        = chain_p$a,
          b        = chain_p$b,
          c        = chain_p$c,
          sigma2   = chain_p$sigma2,
          sum_dx   = sum(df_aux$Dx),
          mort_dec = mort_dec_bootstrap$thetastar
        )
        
        # --------------------------------------------------
        # Plateau onset and level for each posterior draw
        # --------------------------------------------------
        
        chain_posterior <- split(chain, seq_len(nrow(chain)))
        res_plateau     <- lapply(
          chain_posterior,
          function(row_post) plateau_estim(row_post, age0 = age_0)
        )
        
        plateau_df <- as.data.frame(do.call(rbind, res_plateau))
        
        chain$onset   <- plateau_df$onset
        chain$level   <- plateau_df$level
        chain$lastage <- plateau_df$lastage
        
        chain
      },
      mc.cores = n_cores_outer,
      title    = cat("Progress -", country_name, sex, "\n")
    )
    
    chain_full <- do.call(rbind, res_list)
    result     <- rbind(result, chain_full)
  }
}

# Optional: save results to disk
dir.create("results", showWarnings = FALSE)
saveRDS(result, file = file.path("results", "plateau_results.rds"))

# Stop cluster if you want to clean up
stopCluster(cluster)
