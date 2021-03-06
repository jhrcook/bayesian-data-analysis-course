data {
  int<lower=0> N;  // number of data points per machine
  int<lower=0> J;  // number of machines
  vector[J] y[N];  // quality control data points
}

parameters {
  vector[J] mu;
  real<lower=0> sigma;
  real alpha;
  real<lower=0> tau;
}

model {
  // hyper-priors
  alpha ~ normal(100, 10);
  tau ~ normal(0, 10);

  // priors
  mu ~ normal(alpha, tau);
  sigma ~ inv_chi_square(5);

  // likelihood
  for (j in 1:J){
    y[,j] ~ normal(mu[j], sigma);
  }
}

generated quantities {
  // Compute the predictive distribution for the sixth machine.
  real y6pred;  // Leave for compatibility with earlier assignments.
  vector[J] ypred;
  real mu7pred;
  real y7pred;
  vector[J] log_lik[N];

  y6pred = normal_rng(mu[6], sigma);
  for (j in 1:J) {
    ypred[j] = normal_rng(mu[j], sigma);
  }

  mu7pred = normal_rng(alpha, tau);
  y7pred = normal_rng(mu7pred, sigma);

  for (j in 1:J) {
    for (n in 1:N) {
      log_lik[n,j] = normal_lpdf(y[n,j] | mu[j], sigma);
    }
  }
}
