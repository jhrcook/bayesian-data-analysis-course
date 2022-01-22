data {
  int<lower=0> N;           // number of data points
  int<lower=0> A;           // constant used in model of measurement error
  vector<lower=0>[N] x;     // concentration values
  vector<lower=0>[N] y;     // observed color intensity
  int<lower=0> M;           // number of new x values
  vector<lower=0>[M] xnew;  // new x values
}

parameters {
  vector<lower=0>[4] beta;
  real<lower=0,upper=1> alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[N] g;
  vector<lower=0>[N] tau;

  for (i in 1:N) {
    g[i] = beta[1] + beta[2] / (1 + (x[i] / beta[3]) ^ (-beta[4]));
    tau[i] = ((g[i] / A) ^ (2.0 * alpha)) * (sigma ^ 2.0);
  }
}

model {
  // Priors
  alpha ~ beta(1, 1);
  beta[1] ~ normal(10, 2.5);
  beta[2] ~ normal(100, 5);
  beta[3] ~ normal(0, 1);
  beta[4] ~ normal(0, 2.5);
  sigma ~ normal(0, 2.5);

  // Likelihood
  for (i in 1:N) {
    y[i] ~ normal(g[i], tau[i]);
  }
}

generated quantities {
  vector[N] ypred;
  vector[N] log_lik;

  vector[M] g_hat;
  vector[M] tau_hat;
  vector[M] ynew;

  for (i in 1:N) {
    ypred[i] = normal_rng(g[i], tau[i]);
    log_lik[i] = normal_lpdf(y[i] | g[i], tau[i]);
  }

  for (i in 1:M) {
    g_hat[i] = beta[1] + beta[2] / (1 + (xnew[i] / beta[3]) ^ (-beta[4]));
    tau_hat[i] = ((g_hat[i] / A) ^ (2.0 * alpha)) * (sigma ^ 2.0);
    ynew[i] = normal_rng(g_hat[i], tau_hat[i]);
  }
}
