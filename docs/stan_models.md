# (PART) Models {-}

# Stan models

Below are the Stan models built as a part of this course.
The original files are available in the GitHub repo in the "models" directory.


## Model: `8-schools.stan`

```stan
// 8 schools model from 'rstan' documentation.
// https://mc-stan.org/rstan/articles/rstan.html

data {
  int<lower=0> J;          // number of schools
  real y[J];               // estimated treatment effects
  real<lower=0> sigma[J];  // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}
```

## Model: `assignment06-bioassay.stan`

```stan

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
```

## Model: `assignment07_factories_hierarchical.stan`

```stan
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
```

## Model: `assignment07_factories_pooled.stan`

```stan
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
```

## Model: `assignment07_factories_separate.stan`

```stan
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
```

## Model: `assignment07-drownings.stan`

```stan
data {
  int<lower=0> N;  // number of data points
  vector[N] x;     // observation year
  vector[N] y;     // observation number of drowned
  real xpred;      // prediction year
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;  // fix: 'upper' should be 'lower'
}
transformed parameters {
  vector[N] mu = alpha + beta*x;
}
model {
  alpha ~ normal(135, 50);   // prior on `alpha`
  beta ~ normal(0, 26);      // prior on `beta`
  y ~ normal(mu, sigma);     // fix: missing semicolor
}
generated quantities {
  real ypred = normal_rng(alpha + beta*xpred, sigma);  // fix: use `xpred`
}
```

## Model: `serial-dilution.stan`

```stan
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
```
