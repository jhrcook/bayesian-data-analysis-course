# Section 18. Notes on 'Ch 20. Basis function models'

2022-01-12



> These are just notes on a single chapter of *BDA3* that were not part of the course.

## Chapter 20. Basis function models

- chapter 19 focused on nonlinear models with $\text{E}(y|X,\beta) = \mu(X_i,\phi)$ where $\mu$ is a parametric nonlinear function of unknowns $\phi$
- in this and following chapters, consider models where $\mu$ is also unknown

### 20.1 Splines and weighted sums of basis functions {-}

- replace $X_i \beta$ with $\mu(X_i)$ where $\mu(\cdot)$ is some class of nonlinear functions
  - different options for modeling $\mu$ including with basis function expansions or Gaussian processes (next chapter)
- basis function approach: $\mu(x) = \sum_{h=1}^{H} \beta_h b_h(x)$
  - $b_h$: set of basis functions
  - $\beta_h$: vector of basis coefficients
- common choices for basis functions are:
  - **Gaussian radial basis functions**: multiple centers of the basis functions with a width parameter controlling a set of Gaussian functions
  - **B-spline**: a piecewise continuous function based on a set of knots
    - knots locations control the flexibility of the basis
    - knots cn be placed uniformly or non-uniformly (e.g. based on the density of the data)
    - can use a "free knot approach" with a prior on the number and location of knots, but is computationally demanding
    - instead can use priors on the coefficients $\beta$ to shrink values to near 0

## 20.2 Basis selection and shrinkage coefficients

- common to not know which basis functions are really needed
- can use a variable selection approach to allow the model to estimate the "importance" of each basis function
  - can then either select the best model from the posterior or average over all possible models by weighting each basis by its importance
- possible for some bias based on initial choice of the basis functions
  - implied prior information on the smoothness and shape of the model
  - can include multiple types of basis functions in the initial collection

### Shrinkage priors {-}

- allowing basis function coefficients to be zero with positive probability represents a challenge for sampling from the posterior
  - with many basis functions, it is computationally infeasible to visit all possible states
- may be better to use shrinkage priors instead
  - there are various options discussed in the book and there are likely others recommended now

## 20.3 Non-normal models and regression surfaces

### Other error distributions {-}

- may want to model data that is not a continuous output variable $y$ or does not have Gaussian residuals
  - can modify the residual densities with different prior distributions to accommodate outliers
  - can use the basis function and its coefficients $\eta_i = w_i \beta$ as the linear component in a GLM

### Multivariate regression surfaces {-}

- careful with *curse of dimensionality*
- one option is to assume additive of the covariates so can model as the sum of univariate regression functions
  - this does not always make sense and a different approach using a tensor product is described at the end of the section
- can use informative priors to help restrict the search space
  - is a form of including prior information, so still proper Bayesian results and measures of uncertainty
  - e.g. if we know *a priori* that the mean response variable is non-decreasing, restrict the coefficients $\beta$ to be non-negative
