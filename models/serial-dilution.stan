data {
  int<lower=0> N;
  int<lower=0> A;
  vector<lower=0>[N] x;
  vector<lower=0>[N] y;
}

parameters {
  real<lower=0> beta[4];
  real<lower=0,upper=1> alpha;
  real<lower=0> sigma;
}

model {
  vector[N] g;
  alpha ~ beta(1, 1);
  beta ~ normal(0, 5);
  sigma ~ normal(0, 5);

  for (i in 1:N) {
    g[i] = beta[1] + (beta[2] / (1 + (x[i] / beta[3]) ^ (-beta[4])));
  }

  for (i in 1:N) {
    y[i] ~ normal(g[i], (g[i] / A) ^ (2 * alpha) * (sigma ^ 2));
  }
}
