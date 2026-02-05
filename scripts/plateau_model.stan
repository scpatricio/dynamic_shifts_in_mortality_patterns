// Gamma-Gompertz-Makeham with Poisson likelihood

data {
  int<lower = 1> N;
  int<lower = 0> Dx[N];
  vector<lower = 0>[N] Ex;
  vector[N] t;
}

parameters {
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = 0> sigma2;
  real<lower = 0> c;
}

transformed parameters {
  vector[N] et;
  vector[N] hazard;
  vector[N] mu;

  et      = exp(b * t);
  hazard  = a .* et ./ (1 + sigma2 * a .* (et - 1) / b) + c; // force of mortality
  mu      = Ex .* hazard;                                    // Poisson mean
}

model {
  // priors
  a      ~ inv_gamma(2, 1);
  b      ~ inv_gamma(2, 1);
  c      ~ inv_gamma(2, 1);
  sigma2 ~ gamma(0.5, 0.5);

  // likelihood
  Dx ~ poisson(mu);
}

generated quantities {
  int Dx_rep[N];

  for (n in 1:N) {
    Dx_rep[n] = poisson_rng(mu[n]);
  }
}
