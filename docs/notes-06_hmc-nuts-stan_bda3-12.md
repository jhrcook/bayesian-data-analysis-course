# Section 6. HMC, NUTS, and Stan

2021-10-04


```r
knitr::opts_chunk$set(echo = TRUE, dpi = 300, comment = "#>")
```

## Resources

- BDA3 chapter 12 and [reading instructions](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/BDA3_ch12_reading-instructions.pdf)
- lectures:
  - ['6.1 HMC, NUTS, dynamic HMC and HMC specific convergence diagnostics'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=1744f6a0-84d3-4218-8a86-aae600ba7e84)
  - ['6.2 probabilistic programming and Stan'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=e60ba1a9-f752-4b0a-88c6-aae600caa61a)
- [slides](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/slides_ch12.pdf)
- [Assignment 6](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-06.pdf)

## Notes

### Reading instructions

#### Outline of the chapter 12 {-}

- 12.1 Efficient Gibbs samplers (not part of the course)
- 12.2 Efficient Metropolis jump rules (not part of the course)
- 12.3 Further extensions to Gibbs and Metropolis (not part of the course)
- 12.4 Hamiltonian Monte Carlo (used in Stan)
- 12.5 Hamiltonian dynamics for a simple hierarchical model (read through)
- 12.6 Stan: developing a computing environment (read through)

#### Hamiltonian Monte Carlo {-}

- review of static HMC (the number of steps in dynamic simulation are not adaptively selected) is @Neal2012-mu
- Stan uses a variant of dynamic Hamiltonian Monte Carlo (using adaptive number of steps in the dynamic simulation), which has been further developed since BDA3 was published
- The first dynamic HMC variant was by @Hoffman2011-nx
- The No-U-Turn Sampler (NUTS) is often associated with Stan, but the current dynamic HMC variant implemented in Stan has some further developments described (mostly) by @Betancourt2017-bp
  - Instead of reading all above, you can also watch a video: [Scalable Bayesian Inference with Hamiltonian Monte Carlo](https://wwwyoutube.com/watch?v=jUSZboSq1zg) by Betancourt

#### Divergences and BFMI {-}

- divergences and Bayesian Fraction of Missing Information (BFMI) are HMC specific convergence diagnostics developed by Betancourt after BDA3 was published
  - Divergence diagnostic checks whether the discretized dynamic simulation has problems due to fast varying density.
    - See more in a [case study](http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html)
  - BFMI checks whether momentum resampling in HMC is sufficiently efficient [@Betancourt2016-je]
  - [Brief Guide to Stan’s Warnings](https://mc-stan.org/misc/warnings.html) provides summary of available convergence diagnostics in Stan and how to interpret them.

### Chapter 12. Computationally efficient Markov chain simulation

(skipping sections 12.1-12.3)

#### 12.4 Hamiltonian Monte Carlo {-}

- random walk of Gibbs sampler and Metropolis algorithm is inherently inefficient
- Hamiltonian Monte Carlo (HMC) uses "momentum" to suppress the local random walk behavior of the Metropolis algorithm
  - such that is moves more rapidly through the target distribution
  - for each component $\theta_j$ in the target space, there is a corresponding *momentum* variable $\phi_j$
  - the posterior density $p(\theta|y)$ is augmented by an *independent* distribution $p(\phi|y)$ on the momentum: $p(\theta, \phi | y) = p(\phi) p(\theta | y)$
  - to compute the momentum, HMC requires the gradient of the log-posterior density
    - in practice, this is computed analytically

##### The momentum distribution $p(\phi)$ {-}

- common to use a multivariate normal distribution with mean 0 and a diagonal *mass matrix* $M$
  - $\phi_j \sim \text{N}(0, M_{jj})$ for each $j = 1, \dots, d$
  - ideally, the mass matrix should scale with the inverse covariance matrix of the posterior distribution $(\text{var}(\theta|y))^{-1}$

##### The three steps of an HMC iteration {-}

1. update $\phi$ by sampling from its posterior (same as its prior): $\phi \sim \text{N}(0, M)$
2. simultaneously update $\theta$ and $\phi$ using a leapfrog algorithm to simulate physical Hamiltonian dynamics where the position and momentum evolve continuously; for $L$ leapfrog steps with scaling factor $\epsilon$:
  a. use the gradient of the log-posterior density of $\theta$ to make a half-step of $\phi$: $\phi \leftarrow \phi + \frac{1}{2} \epsilon \frac{d \log p(\theta|y)}{d \theta}$
  b. use the momentum $\phi$ to update position $\theta$: $\theta \leftarrow \theta + \epsilon M^{-1} \phi$
  c. use the gradient of $\theta$ for the second half-step of $\phi$: $\phi \leftarrow \phi + \frac{1}{2} \epsilon \frac{d \log p(\theta|y)}{d \theta}$
3. calculate the accept/reject ratio $r$ \@ref(eq:accept-reject-r)
4. set $\theta^t = \theta^*$ with probability $\min(r, 1)$, else reject the proposed $\theta$ and set $\theta^t = \theta^{t-1}$

\begin{equation}
  r
    = \frac{p(\theta^*, \phi^* | y)}{p(\theta^{t-1}, \phi^{t-1} | y)}
    = \frac{p(\theta^*|y) p(\phi^*)}{p(\theta^{t-1}|y) p(\phi^{t-1})}
  (\#eq:accept-reject-r)
\end{equation}

> Pause here and watch video on HMC by Betancourt: [Scalable Bayesian Inference with Hamiltonian Monte Carlo](https://youtu.be/jUSZboSq1zg).

#### 12.5 Hamiltonian Monte Carlo for a hierarchical model {-}

(Walks through the process of deciding on model and HMC parameters and tuning $\epsilon$ and $L$ for HMC.)

#### 12.6 Stan: developing a computing environment {-}

(*Very* briefly describes Stan.)

### Lecture notes

#### 6.1 HMC, NUTS, dynamic HMC and HMC specific convergence diagnostics {-}

*Definitely worth looking at the visualizations in this blog post: [Markov Chains: Why Walk When You Can Flow?](http://elevanth.org/blog/2017/11/28/build-a-better-markov-chain/)*

*Interactively play with HMC and NUTS: [MCMC demo](https://chi-feng.github.io/mcmc-demo/)*

- Hamiltonian Monte Carlo
  - uses gradient of log density for more efficient sampling of the posterior
  - parameters: step size $\epsilon$, number of steps in each chain $L$
    - if step size is too large, then the chain can get further and further away from the high density regions leading to "exploding error"
      - can experiment with this in the simulation linked above
  - No U-Turn Sampling (NUTS) and dynamic HMC
    - adaptively selects the number of steps to improve robustness and sampling efficiency
    - "dynamic" HMC refers to dynamic trajectory length
    - to maintain "reversibility" condition for the Markov chain, must simulate in two directions
- dynamic HMC in Stan
  - use a growing tree to increase simulation trajectory until no-U-turn stopping criterion
    - there is a *max tree depth* parameter to control the size of this
    - pick a draw along the trajectory with probability adjusted to account for the error in discretized dynamic simulation and higher probability for parts further away from the start
      - therefore, don't always end up at the end of the trajectory, but usually somewhere near the end
  - Stan can also adjust the mass matrix to "reshape" the posterior to make more circular and reduce correlations between parameters to make sampling faster and more efficient
    - occurs during the initial adaptation in the warm-up
- max tree depth diagnostic
  - indicates inefficiency in sampling leading to higher autocorrelation and lower ESS
  - possibly step size is too small
  - different parameterizations matter
- divergences
  - HMC specific
  - indicates that there are unexpectedly fast changes in log-density
  - comes from exploding error where the chain leaves high-density regions of the posterior and gets lost in other directions
  - occurs in funnel geometries common in hierarchical models because in the neck of the funnel, the high-density region is so thin, the trajectory can easily leave and enter a very low-density region
- problematic distributions
  - *Nonlinear dependencies*
    - simple mass matrix scaling doesn’t help
  - *Funnels*
    - optimal step size depends on location
    - can get stuck in the neck of the funnel causing divergences
  - *Multimodal*
    - difficult to move from one mode to another
    - can be seen if multiple chains end up in different locations
  - *Long-tailed with non-finite variance and mean*
    - efficiency of exploration is reduced
    - central limit theorem doesn’t hold for mean and variance

#### 6.2 probabilistic programming and Stan {-}

example: Binomial model

```stan
data {
  int <lower=0> N; // number of experiments
  int<lower=0,upper=N> y; // number of successes
}

parameters {
  real <lower=0,upper=1> theta ; // parameter of the binomial
}

model {
  theta ~ beta(1,1); //prior
  y ~ binomial (N, theta ); // observation model
}
```

example: running Stan from R

```r
library (rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel ::detectCores())

d_bin <- list(N = 10, y = 7)
fit_bin <- stan(file = "binom.stan", data = d_bin)
```

example: running Stan from Python

```python
import pystan
import stan_utility

data = {"N": 10, "y": 8}
model = stan_utility.compile_model('binom.stan')
fit = model.sampling(data=data)
```

example: Difference between proportions

```stan
data {
  int<lower=0> N1;
  int<lower=0> y1;
  int<lower=0> N2;
  int<lower=0> y2;
}
parameters {
  real <lower=0,upper=1> theta1;
  real <lower=0,upper=1> theta2;
}
model {
  theta1 ~ beta(1,1);
  theta2 ~ beta(1,1);
  y1 ~ binomial(N1,theta1);
  y2 ~ binomial(N2,theta2);
}
generated quantities {
  real oddsratio;
  oddsratio = (theta2 / (1 - theta2)) / (theta1 / (1 - theta1))
}
```

some HMC-specific diagnostics with 'rstan'

```r
check_treedepth(fit_bin2)
check_energy(fit_bin2)
check_div(fit_bin2)
```

example of scaling data in the Stan language

```stan
data {
  int<lower=0> N; // number of data points
  vector[N] x;
  vector[N] y;
  real xpred; // input location for prediction
}

transformed_data {
  vector[N] x_std;
  vector[N] y_std;
  real xpred_std;
  x_std = (x - mean(x)) / sd(x);
  y_std = (y - mean(y)) / sd(y);
  xpred_std = (xpred - mean(x)) / sd(x);
}
```

- other useful R packages worth looking into:
  - 'rstanarm' and 'brms': for quickly building models with the R formula language
  - 'shinystan': interactive diagnostics
  - 'bayesplot': visualization and model checking (see model checking in Ch 6)
  - 'loo': cross-validation model assessment, comparison and averaging (see Ch 7)
  - 'projpred': projection predictive variable selection

---


```r
sessionInfo()
```

```
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] bookdown_0.24     clisymbols_1.2.0  digest_0.6.27     R6_2.5.0         
#>  [5] jsonlite_1.7.2    magrittr_2.0.1    evaluate_0.14     stringi_1.7.3    
#>  [9] rlang_0.4.11      renv_0.14.0       jquerylib_0.1.4   bslib_0.2.5.1    
#> [13] rmarkdown_2.10    tools_4.1.2       stringr_1.4.0     glue_1.4.2       
#> [17] xfun_0.25         yaml_2.2.1        compiler_4.1.2    htmltools_0.5.1.1
#> [21] knitr_1.33        sass_0.4.0
```
