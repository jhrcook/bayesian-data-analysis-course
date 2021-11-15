data {
  int<lower=0> N;  // number of data points per machine
  int<lower=0> J;  // number of machines
  vector[J] y[N];  // quality control data points
}

parameters {
  vector[J] mu;
  vector<lower=0>[J] sigma;
}

model {
  // priors
  for (j in 1:J) {
    mu[j] ~ normal(100, 10);
    sigma[j] ~ inv_chi_square(5);
  }

  // likelihood
  for (j in 1:J){
    y[,j] ~ normal(mu[j], sigma[j]);
  }
}

generated quantities {
  // Compute the predictive distribution for the sixth machine.
  real y6pred;
  vector[J] log_lik[N];

  y6pred = normal_rng(mu[6], sigma[6]);

  for (j in 1:J) {
    for (n in 1:N) {
      log_lik[n,j] = normal_lpdf(y[n,j] | mu[j], sigma[j]);
    }
  }
}
