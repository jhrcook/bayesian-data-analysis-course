# Section 11. Normal approximation & Frequency properties

2021-11-19



## Resources

- reading:
  - BDA3 ch 4. *Asymptotics and connections to non-Bayesian approaches*
  - [reading instructions](../reading-instructions/BDA3_ch04_reading-instructions.pdf)
- lectures:
  - [Lecture 11.1. 'Normal approximation (Laplace approximation)'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=e22fedc7-9fd3-4d1e-8318-ab1000ca45a4)
  - [Lecture 11.2. 'Large sample theory and counter examples'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=a8e38a95-a944-4f3d-bf95-ab1000dbdf73)
- [slides](../slides/slides_ch04.pdf)

## Notes

### Reading instructions

- chapter outline:
  - 4.1 Normal approximation (Laplace’s method)
  - 4.2 Large-sample theory
  - 4.3 Counterexamples
  - 4.4 Frequency evaluation (not part of the course, but interesting)
  - 4.5 Other statistical methods (not part of the course, but interesting)

### Chapter 4. Asymptotics and connections to non-Bayesian approaches

- *asymptotic theory*: as sample size increased, the influence of the prior on the posteior decreases
  - used as the justification of non-informative priors in many cases

#### 4.1 Normal approximations of the posterior distribution {-}

- if the posterior distribution $p(\theta|y)$ is unimodal and symmetric, it is useful to approximate it as a normal
  - therefore, the log of the posterior is a quadratic function of $\theta$
- "For a finite sample size $n$, the normal approximation is typically more accurate for conditional and marginal distributions of components of $\theta$ than for the full joint distribution." (pg. 85)
- common to use the normal approximations to quickly debug or sanity-check a model's code

#### 4.2 Large-sample theory {-}

*asymptotic normality of the posterior distribution*: with more data from the same underlying process, the posterior distribution of the parameter vector approaches multivariate normality even if the true distribution of the data is not within the parametric family under consideration (pg. 87)
  - particularly for independent samples from the data-generating process
- summary at the limit of large $n$:
  - the posterior mode $\hat\theta$ approaches the true $\theta_0$
  - the likelihood dominates the prior distribution

#### 4.3 Counterexamples to the theorems {-}

- there are many instances where large amounts of data do not allow for the normal approximation:
- **underidentified models and nonidentified parameters**
  - "the model is *underidentified* given data $y$ if the likelihood $p(\theta|y)$ is equal for a range of $\theta$" (pg. 89)
  - there is no single point $\theta_0$ to which the posterior distribution can converge given infinite data
  - a parameter can be nonidentified if there is no supply of information about it
    - results in its posterior being identical to its prior
- **number of parameters increasing with sample size**
  - in complicated problems, the number of parameters can scale with the amount of data
  - e.g. Gaussian processes or hierarchical models
- **aliasing**
  - the same likelihood function repeats at a discrete set of points
    - a special case of underidentified parameters
  - e.g. a mixture model with two mixed distributions with the same parameters
- **unbounded likelihoods**
  - if the likelihood is unbounded, there might not be any posterior mode with the parameter space
    - invalidates bot the consistency results and the normal approximation
- **improper posterior distributions**
  - an improper posterior integrates to infinity, not to 1 as is required by the normal approximation theory
  - an improper posterior can only occur with an improper prior
- **prior distributions that exclude the point of convergence**
- **convergence to the edge of parameter space**
  - if $\theta_0$ is at the edge of the parameter space, the distribution cannot be symmetric
- **tails of the distribution**
  - the normal approximation can be true for almost all of the mass of the posterior but not be true at the tails

### Lecture notes

#### Lecture 11.1. 'Normal approximation (Laplace approximation)' {-}

(no additional notes)

#### Lecture 11.2. 'Large sample theory and counter examples' {-}

- *large sample theory*:
  - *consistency*: if the true distribution is included in the parametric family then the posterior converges to a point $\theta_0$ when $n \rightarrow \infty$
    - "included in the parametric family": $f(y) = p(y|\theta_0)$ for some $\theta_0$
    - the point does not have uncertainty
    - same result as MLE
  - if the true distribution is not included in the parameteric family, then there is no true $\theta_0$, so replace it with the $\theta_0$ that minimized the KL divergence
