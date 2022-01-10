data {
  int<lower=0> N;  // number of data points
  int<lower=0> A;  // constant used in model of measurement error
  vector<lower=0>[N] x;  // concentration values
  vector<lower=0>[N] y;  // observed color intensity
}

parameters {
  real<lower=0> beta[4];
  real<lower=0,upper=1> alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[N] g;
  vector<lower=0>[N] tau;

  for (i in 1:N) {
    g[i] = beta[1] + (beta[2] / (1 + (x[i] / beta[3]) ^ (-1 * beta[4])));
    tau[i] = ((g[i] / A) ^ (2.0 * alpha)) * (sigma ^ 2.0);
  }
}

model {
  // Priors
  alpha ~ beta(2, 2);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 5);
  // Likelihood
  y ~ normal(g, tau);
}

generated quantities {
  vector[N] ypred;
  vector[N] log_lik;

  for (i in 1:N) {
    ypred[i] = normal_rng(g[i], tau[i]);
    log_lik[i] = normal_lpdf(y[i] | g[i], tau[i]);
  }
}
