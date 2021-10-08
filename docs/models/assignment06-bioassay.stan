
data {
  int<lower=0> N;    // number of data points
  vector[N] x;       // dose
  int<lower=0> n[N]; // number of animals
  int<lower=0> y[N]; // number of deaths

  vector[2] mu;                // prior on mean of theta
  matrix<lower=0>[2, 2] sigma; // prior on covariance matrix of theta
}

parameters {
  vector[2] mdl_params;
}

transformed parameters {
  vector[N] theta;
  theta = mdl_params[1] + mdl_params[2] * x;
}

model {
  mdl_params ~ multi_normal(mu, sigma);
  y ~ binomial_logit(n, theta);
}
