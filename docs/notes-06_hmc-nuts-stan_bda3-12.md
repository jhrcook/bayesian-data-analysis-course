# Section 6. HMC, NUTS, and Stan

2021-10-04



## Resources

- BDA3 chapter 12 and [reading instructions](../reading-instructions/BDA3_ch12_reading-instructions.pdf)
- lectures:
  - ['6.1 HMC, NUTS, dynamic HMC and HMC specific convergence diagnostics'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=1744f6a0-84d3-4218-8a86-aae600ba7e84)
  - ['6.2 probabilistic programming and Stan'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=e60ba1a9-f752-4b0a-88c6-aae600caa61a)
- [slides](../slides/slides_ch12.pdf)
- [Assignment 6](assignments/assignment-06.pdf)

## Notes

### Reading instructions

#### Outline of the chapter 12

- 12.1 Efficient Gibbs samplers (not part of the course)
- 12.2 Efficient Metropolis jump rules (not part of the course)
- 12.3 Further extensions to Gibbs and Metropolis (not part of the course)
- 12.4 Hamiltonian Monte Carlo (used in Stan)
- 12.5 Hamiltonian dynamics for a simple hierarchical model (read through)
- 12.6 Stan: developing a computing environment (read through)

#### Hamiltonian Monte Carlo
- review of static HMC (the number of steps in dynamic simulation are not adaptively selected) is Radford Neal (2011)[^1].
- Stan uses a variant of dynamic Hamiltonian Monte Carlo (using adaptive number of steps in the dynamic simulation), which has been further developed since BDA3 was published
- The first dynamic HMC variant was by Hoffman and Gelman (2014)[^2].
- The No-U-Turn Sampler (NUTS) is often associated with Stan, but the current dynamic HMC variant implemented in Stan has some further developments described (mostly) by Betancourt (2018)[^3]
  - Instead of reading all above, you can also watch a video: [Scalable Bayesian Inference with Hamiltonian Monte Carlo](https://wwwyoutube.com/watch?v=jUSZboSq1zg) by Betancourt

[^1]: Radford Neal (2011). MCMC using Hamiltonian dynamics. In Brooks et al (ed), *Handbook of Markov Chain Monte Carlo*, Chapman & Hall / CRC Press. Preprint https://arxiv.org/pdf/1206.1901.pdf.
[^2]: Matthew D. Hoffman, Andrew Gelman (2014). The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo. *JMLR*, 15:1593–1623 http://jmlr.org/papers/v15/hoffman14a.html.
[^3]: Michael Betancourt (2018). A Conceptual Introduction to Hamiltonian Monte Carlo. arXiv preprint arXiv:1701.02434 https://arxiv.org/abs/1701.02434.

#### Divergences and BFMI

- divergences and Bayesian Fraction of Missing Information (BFMI) are HMC specific convergence diagnostics developed by Betancourt after BDA3 was published
  - Divergence diagnostic checks whether the discretized dynamic simulation has problems due to fast varying density.
    - See more in a [case study](http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html)
  - BFMI checks whether momentum resampling in HMC is sufficiently efficient[^4]
  - [Brief Guide to Stan’s Warnings](https://mc-stan.org/misc/warnings.html) provides summary of available convergence diagnostics in Stan and how to interpret them.

[^4]: Betancourt, Michael. 2016. “Diagnosing Suboptimal Cotangent Disintegrations in Hamiltonian Monte Carlo.” arXiv [stat.ME]. arXiv. http://arxiv.org/abs/1604.00695.

### Chapter 12. Computationally efficient Markov chain simulation

(skipping sections 12.1-12.3)

#### 12.4 Hamiltonian Monte Carlo

- random walk of Gibbs sampler and Metropolis algorithm is inherently inefficient
- Hamiltonian Monte Carlo (HMC) uses "momentum" to suppress the local random walk behavior of the Metropolis algorithm
  - such that is moves more rapidly through the target distribution
  - for each component $\theta_j$ in the target space, there is a corresponding *momentum* variable $\phi_j$
  - the posterior density $p(\theta|y)$ is augmented by an *independent* distribution $p(\phi|y)$ on the momentum: $p(\theta, \phi | y) = p(\phi) p(\theta | y)$
  - to compute the momentum, HMC requires the gradient of the log-posterior density
    - in practice, this is computed analytically

##### The momentum distribution $p(\phi)$

- common to use a multivariate normal distribution with mean 0 and a diagonal *mass matrix* $M$
  - $\phi_j \sim \text{N}(0, M_{jj})$ for each $j = 1, \dots, d$
  - ideally, the mass matrix should scale with the inverse covariance matrix of the posterior distribution $(\text{var}(\theta|y))^{-1}$

##### The three steps of an HMC iteration

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

#### 12.5 Hamiltonian Monte Carlo for a hierarchical model

(Walks through the process of deciding on model and HMC parameters and tuning $\epsilon$ and $L$ for HMC.)

#### 12.6 Stan: developing a computing environment

(*Very* briefly describes Stan.)

### Lecture notes

#### 6.1 HMC, NUTS, dynamic HMC and HMC specific convergence diagnostics

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

#### 6.2 probabilistic programming and Stan

example: Binomial model

```
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

```
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

```
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
