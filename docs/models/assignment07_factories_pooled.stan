data {
  int<lower=0> N;  // number of data points
  vector[N] y;     // machine quality control data
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  // priors
  mu ~ normal(100, 10);
  sigma ~ inv_chi_square(5);

  // likelihood
  y ~ normal(mu, sigma);
}

generated quantities {
  real ypred;
  vector[N] log_lik;

  ypred = normal_rng(mu, sigma);

  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | mu, sigma);

}
