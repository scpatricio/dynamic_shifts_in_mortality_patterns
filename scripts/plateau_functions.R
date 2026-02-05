# Core functions for the plateau chapter (GG-Makeham + plateau detection)

library(randomForest)
library(modeest)

# -----------------------------------------------------------
# Force of mortality (Gamma-Gompertz-Makeham)
# theta = (a, b, c, sigma2)
# -----------------------------------------------------------

mu <- function(t, theta) {
  a      <- theta[1]
  b      <- theta[2]
  c      <- theta[3]
  sigma2 <- theta[4]

  a * exp(b * t) / (1 + sigma2 * (a / b) * (exp(b * t) - 1)) + c
}

# For completeness, a general GGM hazard with all four parameters
mu_ggm <- function(t, a = 1, b = 1, c = 1, sigma2 = 0.5) {
  a * exp(b * t) / (1 + sigma2 * (a / b) * (exp(b * t) - 1)) + c
}

# -----------------------------------------------------------
# Random lifetime generators
# -----------------------------------------------------------

rgompertz <- function(n, a = 1, b = 1) {
  u <- runif(n)
  log(1 - b * log(1 - u) / a) / b
}

rgammamak <- function(n, a = 1, b = 1, c = 1, sigma2 = 0.5) {
  z   <- rgamma(n, shape = 1 / sigma2, scale = sigma2)
  out <- apply(
    cbind(
      rgompertz(n, a * z, b),
      rexp(n, c)
    ),
    1,
    min
  )
  out
}

# -----------------------------------------------------------
# Mortality deceleration via random forest on 2nd derivative
# -----------------------------------------------------------

max_func <- function(id_1, n_tree = 100, data.frame.full) {
  df_func <- data.frame.full[id_1, ]
  df_func <- df_func |>
    na.omit() |>
    dplyr::arrange(Age)
  
  # jitter ages
  df_func$Age <- df_func$Age + runif(nrow(df_func), -0.5, 0.5)
  
  grid <- seq(min(df_func$Age), max(df_func$Age), by = 0.1)
  
  rf_fit <- randomForest(
    deriv_2   ~ Age,
    data       = df_func,
    importance = FALSE,
    replace    = TRUE,
    ntree      = n_tree,
    na.action  = na.omit,
    weights    = df_func$Dx
  )
  
  y_pred <- predict(rf_fit, data.frame(Age = grid))
  grid[which.max(y_pred)]
}

# -----------------------------------------------------------
# Plateau onset and level via regression trees on simulated mx
# posterior_sample must have a, b, c, sigma2, sum_dx
# age0: lower age used in the analysis (e.g. 60)
# -----------------------------------------------------------

plateau_estim <- function(posterior_sample, age0) {
  theta <- c(
    posterior_sample$a,
    posterior_sample$b,
    posterior_sample$c,
    posterior_sample$sigma2
  )
  
  n_samp <- posterior_sample$sum_dx
  
  ind_lifetime <- rgammamak(
    n      = n_samp,
    a      = theta[1],
    b      = theta[2],
    c      = theta[3],
    sigma2 = theta[4]
  )
  
  t_vals <- as.numeric(names(table(floor(ind_lifetime))))
  Dx     <- as.vector(table(floor(ind_lifetime)))
  Ex     <- sum(Dx) - cumsum(Dx) + Dx / 2
  
  df_new <- data.frame(
    Age    = t_vals,
    Dx     = Dx,
    Ex     = Ex,
    mx     = Dx / Ex,
    log_mx = log(Dx / Ex)
  )
  
  n_tree <- 100L
  rf_fit <- randomForest(
    log_mx     ~ Age,
    data        = df_new,
    importance  = FALSE,
    replace     = TRUE,
    keep.forest = TRUE,
    ntree       = n_tree,
    na.action   = na.omit,
    weights     = log(df_new$Dx + 2),
    nodesize    = 1
  )
  
  cut   <- numeric(n_tree)
  level <- numeric(n_tree)
  
  for (j in seq_len(n_tree)) {
    tree_j <- getTree(rfobj = rf_fit, k = j, labelVar = FALSE)
    cut[j]   <- max(tree_j[, 4])
    level[j] <- tree_j[which.max(tree_j[, 4]), 6]
  }
  
  c(
    onset   = modeest::mlv(cut) + age0,
    level   = modeest::mlv(exp(level)),
    lastage = max(df_new$Age) + age0
  )
}
