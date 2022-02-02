# Section 3. Multidimensional Posterior

2021-09-02



## Resources

- BDA3 chapter 3 and [reading instructions](../reading-instructions/BDA3_ch03_reading-instructions.pdf)
- lectures:
  - ['Lecture 3. Multiparameter models, joint, marginal and conditional distribution, normal model, bioassay example, grid sampling and grid evaluation'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=ab958b4b-e2c4-4534-8305-aad100ba191f)
- [slides](../slides/slides_ch3.pdf)
- [Assignment 3](assignments/assignment-03.pdf)

## Notes

### Reading instructions

- the trace of a square matrix $tr(A)$ is the sum of the diagonals
  - the following property is used in derivation of 3.11: $tr(ABC) = tr(CAB) = tr(BCA)$

### Chapter 3. Introduction to multiparameter models

#### Averaging over 'nuisance parameters'

- suppose the unknown variable $\theta$ is a vector of length two: $\theta= (\theta_1, \theta_2)$
  - may only care about one of the variables, but the other is still required for a good model
  - example model: $y | \mu, \sigma^2 \sim N(\mu, \sigma^2)$
    - here, $\theta$ would be the unknown values $\mu (=\theta_1)$ and $\sigma (=\theta_2)$, but we really only care about $\mu$
  - we want $p(\theta_1|y)$
    - derive it from the *joint posterior density*: $p(\theta_1, \theta_2) \propto p(y|\theta_1, \theta2) p(\theta_1, \theta_2)$
    - by averaging over $\theta_2$: $p(\theta_1|y) = \int p(\theta_1, \theta_2| y) d\theta_2$
      - "integrate over the uncertainty in $\theta_2$"

#### Summary of elementary modeling and computation

- the following is an outline of a simple Bayesian analysis
  - it will change when we get to more complex models whose posteriors are estimated by more complex sampling processes
1. write the likelihood: $p(y|\theta)$
2. write the posterior density: $p(\theta|y) \propto p(\theta) p(y|\theta)$$
3. estimate the parameters $\theta$ (e.g. using MLE)
4. draw simulations $\theta^1, \dots, \theta^S$ for the posterior distribution (using the results of 3 as a starting point); use the samples to compute any other functions of $\theta$ that are of interest
5. if any predictive quantities $\tilde{y}$ are of interest, simulate $\tilde{y}^1, \dots, \tilde{y}^S$ from $p(\tilde{y} | \theta^s)$

### Lecture notes

(No extra notes were taken — some comments added directly to slides.)
