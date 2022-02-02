# Section 20. Notes on 'Ch 22. Finite mixture models'

2022-01-18



> These are just notes on a single chapter of *BDA3* that were not part of the course.

## Chapter 22. Finite mixture models

- "when measurements of a random variable are taken under two different conditions" (pg. 519)
- where the data contains multiple subpopulations where each has a different, relatively simple model
- basic mixture modeling principle is to introduce *unobserved indicators* $z$ to specify the mixture component for an observation
  - can think of a mixture indicator as missing data

### 22.1 Setting up and interpreting mixture models

#### Finite mixtures

- want to model the distribution of $y = (y_1, \dots, y_n)$ or $y|x$  as a mixture of $H$ components
  - for each component $h \in (1, \dots, H)$, the distribution $f_h(y_i | \theta_h)$ depends on a parameter vector $\theta_h$
  - $\lambda_h$ denotes the proportion of the population in component $h$
    - $\sum_{h=1}^{H} \lambda_h = 1$
  - common to assume all mixture components have the same parametric form
  - thus, the sampling distribution of $y$ is:

$$
p(y_i | \theta, \lambda) = \lambda_1 f(y_i | \theta_1) + \dots + \lambda_H f(y_i | \theta_H)
$$

- can think of the mixture distribution probabilities $\lambda$ as priors over the parameters $\theta_h$
  - or as a description of the variation in $\theta$ in a population
  - akin to a hierarchical model
- introduce the indicator variables $z_{ih}$ where $z_{ih} = 1$ if the $i$th data point is drawn from component $h$ and 0 otherwise
  - the $lambda$ values are used to determine $z$
    - can think of $lambda$ as a hyperprior over $z$
  - joint distribution of the observed data $y$ and the unobserved indicators $z$ conditions on the model parameters:
    - only one $z_{ih}$ can be 1 for each $i$

$$
\begin{aligned}
p(y, z | \theta, \lambda) &= p(z | \lambda) p(y | z, \theta) \\
 &= \prod_{i=1}^n \prod_{h=1}^H (\lambda_h f(y_i | \theta_h))^{z_{i,h}}
\end{aligned}
$$

#### Continuous mixtures

- generalize the finite mixture to allow probability of an observation belongs to a class
- hierarchical models are a form a continuous mixture model
  - each observed value $y_i$ is modeled as coming from a mixture of models defined by the probability of values for $\theta$
- in the book, the focus is on finite mixtures and "minor modifications" are generally required to form a continuous distribution

#### Identifiabilitiy of the mxixture likelihood

- **all finite mixture models are nonidentifiable: the distribution is unchanged if the group labels are permuted**
- in many cases, purposeful, informative priors can solve the issue

#### Priors distributions

- the priors for a finite mixture model's parameters $\theta$ and $\lambda$ are usually the product of the two independent priors on each variable
  - because the vector of mixture indicators $z_i = (z_{i,1}, \dots, z_{i,H})$ is multinomial with parameter $\lambda$, a common prior for $\lambda$ is the Dirichlet
    - $\lambda \sim \text{Dirichlet}(\alpha_1, \dots, \alpha_H)$
  - $\theta = (\theta_1, \dots, \theta_H)$ is the vector of parameters for each component's sub-model
    - some can be shared across components (i.e. equal variance for a group of normal distributions)

#### Number of mixture components

- can model $H$ as unknown but is computationally expensive
- usually can just build models with different $H$ and compare their goodness of fit
  - compare the posterior predictive distributions with a "suitably chosen" test quantity

#### Mixtures as true models or approximating distributions

- two classes of thought[^1]:
  1. *theoretical:* a mixture model is "a realistic characterization of the true data-generating mechanism" (pg. 522)
  2. *pragmatic:* "trying to infer latent subpopulations is an intrinsically ill-defined statistical problem, but finite mixture models are nonetheless useful" (pg. 523)

[^1]: I came up with the names of the two schools of thought as descriptive titles; they do not appear in the book.

### 22.4 Unspecifed number of mixture components

- can assign a Poisson distribution as a on $H$ (the number of groups/components in the mixture model)
  - computationally intensive
  - more common to just fit the model with different $H$ and compare with some statistic and a penalty for model complexity
    - WAIC is theoretically justified, but ignores the uncertainty over $H$
    - LOO-CV may be even better
