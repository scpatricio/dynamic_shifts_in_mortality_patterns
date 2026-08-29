data {
  int<lower = 1> Y;
  int<lower = 1> N;
  int<lower = 1, upper = Y> yr[N];
  int<lower = 0> Dx[N];
  vector<lower = 0>[N] Ex;
  vector[N] t;
}

parameters {
  vector<upper = log(1)>[Y] log_a;
  vector<lower = log(0.05), upper = log(0.20)>[Y] log_b;
  vector<upper = log(1)>[Y] log_gamma;

  real<lower = 0> tau_a;
  real<lower = 0> tau_b;
  real<lower = 0> tau_gamma;
  real<lower = 0> phi;
}

transformed parameters {
  vector<lower = 0>[Y] a = exp(log_a);
  vector<lower = 0>[Y] b = exp(log_b);
  vector<lower = 0>[Y] gamma = exp(log_gamma);
  vector<lower = 0>[N] hazard;
  vector<lower = 0>[N] expected_deaths;

  for (n in 1:N) {
    int y = yr[n];
    real growth = exp(b[y] * t[n]);
    hazard[n] = a[y] * growth /
      (1 + gamma[y] * a[y] * (growth - 1) / b[y]);
    expected_deaths[n] = Ex[n] * hazard[n];
  }
}

model {
  tau_a ~ normal(0, 0.5);
  tau_b ~ normal(0, 0.5);
  tau_gamma ~ normal(0, 0.5);
  phi ~ normal(0, 1);

  for (y in 2:Y) {
    log_a[y] ~ normal(log_a[y - 1], tau_a);
    log_b[y] ~ normal(log_b[y - 1], tau_b);
    log_gamma[y] ~ normal(log_gamma[y - 1], tau_gamma);
  }

  Dx ~ neg_binomial_2(expected_deaths, phi);
}

generated quantities {
  int Dx_rep[N];
  vector[Y] asymptotic_hazard = b ./ gamma;

  for (n in 1:N) {
    Dx_rep[n] = neg_binomial_2_rng(expected_deaths[n], phi);
  }
}
