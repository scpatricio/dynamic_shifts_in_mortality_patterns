data {
  int<lower = 1> Y;
  int<lower = 1> N;
  int<lower = 1, upper = Y> yr[N];
  int<lower = 0> Dx[N];
  vector<lower = 0>[N] Ex;
  vector[N] t;
}

transformed data {
  int n_y[Y] = rep_array(0, Y);
  for (n in 1:N) n_y[yr[n]] += 1;

  int max_ny = max(n_y);
  int index_by_cohort[Y, max_ny] = rep_array(0, Y, max_ny);
  {
    int fill[Y] = rep_array(0, Y);
    for (n in 1:N) {
      int y = yr[n];
      fill[y] += 1;
      index_by_cohort[y, fill[y]] = n;
    }
  }
}

parameters {
  vector<upper = log(1)>[Y] log_a;
  vector<lower = log(0.05), upper = log(0.20)>[Y] log_b;
  vector<upper = log(1)>[Y] log_gamma_slab;

  real<lower = 0> tau_a;
  real<lower = 0> tau_b;
  real<lower = 0> tau_gamma;
  real<lower = 0, upper = 1> w;
  real<lower = 0> phi;
}

transformed parameters {
  vector<lower = 0>[Y] a = exp(log_a);
  vector<lower = 0>[Y] b = exp(log_b);
  vector<lower = 0>[Y] gamma_slab = exp(log_gamma_slab);
  vector[Y] log_lik_spike = rep_vector(0, Y);
  vector[Y] log_lik_slab = rep_vector(0, Y);

  for (y in 1:Y) {
    for (j in 1:n_y[y]) {
      int n = index_by_cohort[y, j];
      real growth = exp(b[y] * t[n]);
      real hazard_spike = a[y] * growth;
      real hazard_slab = a[y] * growth /
        (1 + gamma_slab[y] * a[y] * (growth - 1) / b[y]);

      log_lik_spike[y] += neg_binomial_2_lpmf(
        Dx[n] | Ex[n] * hazard_spike, phi
      );
      log_lik_slab[y] += neg_binomial_2_lpmf(
        Dx[n] | Ex[n] * hazard_slab, phi
      );
    }
  }
}

model {
  tau_a ~ normal(0, 0.5);
  tau_b ~ normal(0, 0.5);
  tau_gamma ~ normal(0, 0.5);
  phi ~ normal(0, 1);
  w ~ beta(1, 1);

  for (y in 2:Y) {
    log_a[y] ~ normal(log_a[y - 1], tau_a);
    log_b[y] ~ normal(log_b[y - 1], tau_b);
    log_gamma_slab[y] ~ normal(log_gamma_slab[y - 1], tau_gamma);
  }

  for (y in 1:Y) {
    target += log_mix(w, log_lik_spike[y], log_lik_slab[y]);
  }
}

generated quantities {
  vector[Y] p_no_plateau;
  int Dx_rep[N];

  for (y in 1:Y) {
    real lp_spike = log(w) + log_lik_spike[y];
    real lp_slab = log1m(w) + log_lik_slab[y];
    real p_spike = exp(lp_spike - log_sum_exp(lp_spike, lp_slab));
    int state_spike = bernoulli_rng(p_spike);

    p_no_plateau[y] = p_spike;
    for (j in 1:n_y[y]) {
      int n = index_by_cohort[y, j];
      real growth = exp(b[y] * t[n]);
      real hazard;
      if (state_spike == 1) {
        hazard = a[y] * growth;
      } else {
        hazard = a[y] * growth /
          (1 + gamma_slab[y] * a[y] * (growth - 1) / b[y]);
      }
      Dx_rep[n] = neg_binomial_2_rng(Ex[n] * hazard, phi);
    }
  }
}
