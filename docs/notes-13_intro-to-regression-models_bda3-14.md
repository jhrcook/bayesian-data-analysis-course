# Section 13. Notes on 'Ch 14. Introduction to regression models'

2021-12-06



> These are just notes on a single chapter of *BDA3* that were not part of the course.

## Chapter 14. Introduction to regression models

### 14.1 Conditional modeling

- question: how does one quantity $y$ vary as a function of another quantity or vector of quantities $x$?
  - conditional distribution of $y$ given $x$ parameterized as $p(y|\theta,x)$
- key statistical modeling issues:
  1. defining $y$ and $x$ so that $y$ is reasonably linear as a function of the columns of $X$
    - may need to transform $x$
  2. set priors on the model parameters

### 14.2 Bayesian analysis of classical regression

- simplest case: *ordinary linear regression*
  - observation errors are independent and have equal variance

$$
y | \beta, \sigma, X \sim \text{N}(X \beta, \sigma^2 I)
$$

#### Posterior predictive distribution for new data

- posterior predictive distribution has two sources of uncertainty:
  1. the inherent variability in the model represented by $\sigma$ in $y$
  2. posterior uncertainty in $\beta$ and $\sigma$
- draw a random sample $\tilde{y}$ from the posterior predictive distribution:
  - draw $(\beta, \sigma)$ from their posteriors
  - draw $\tilde{y} \sim \text{N}(\tilde{X} \beta, \sigma^2 I)$

## 14.4 Goals of regression analysis

- at least three goals:
  1. understand the behavior of $y$ given $x$
  2. predict $y$ given $x$
  3. causal inference; predict how $y$ would change if $x$ were changed

## 14.5 Assembling the matrix of explanitory variables

### Identifiability and collinearity

- "the parameters in a classical regression cannot be uniquely estimated if there are more parameters than data points or, more generally, if the columns of the matrix $X$ of explanatory variables are not linearly independent" (pg 365)

### Nonlinear relations

- may need to transform variables
- can include more than one transformation in the model as separate covariates
- GLMs and non-linear models are discussed in later chapters

### Indicator variables

- include a categorical variable in a regression using a indicator variable
  - separate effect for each category
  - or model as related with a hierarchical model

### Interactions

- "If the response to a unit change in $x_i$ depends in what value another predictor $x_j$ has been fixed at, then it is necessary to include *interaction* terms in the model" (pg 367)
  - $(x_i - \bar{x_i})(x_j - \bar{x_j})$

## 14.6 Regularization and dimension reduction

- see lecture notes on regularization for more updated recommendations
- "Bayesian regularization":
  - location and scale of the prior
  - analytic form of the prior (e.g. normal vs. Laplacian vs. Cauchy)
  - how the posterior inference is summarized
