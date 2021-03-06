# Section 9. Model comparison and selection

2021-11-09



## Resources

- BDA3 chapter 7 and [reading instructions](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/BDA3_ch07_reading-instructions.pdf)
- lectures:
  - ['9.1 PSIS-LOO and K-fold-CV'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=50b2e73f-af0a-4715-b627-ab0200ca7bbd)
  - ['9.2 Model comparison and selection'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=b0299d53-9454-4e33-9086-ab0200db14ee)
  - ['9.3 Variable selection with projection predictive variable selection'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=4b6eeb48-ae64-4860-a8c3-ab0200e40ad8)
- slides:
  - [chapter 7b](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/slides_ch7b.pdf)
- [Assignment 8](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-08.pdf)

## Notes

> See notes for course section 8 for notes on ch. 7 of *BDA3*.

### Lecture notes

#### 9.1 PSIS-LOO and K-fold CV {-}

##### PSIS-LOO CV {-}

- given some data $x$ and observed $y$
- sample from posterior: $\theta^{(s)} \sim p(\theta | x, y)$
- compute the predictive density for new $\tilde{x}$ and $\tilde{y}$: $p(\tilde{y} | \tilde{x}, x, y) \approx \frac{1}{S} \sum_{s=1}^S p(\tilde{y} | \tilde{x} \theta^{(s)})$
- LOO-CV but weighting the draws using importance sampling:
  - importance ratio: $r_i^{(s)} = p(\theta^{(s)} | x_{-i}, y_{-i}) / p(\theta^{(s)} | x,y)$
    - the ratio of the posterior without $x_i$ and $y_i$ to the full posterior distribution
    - get higher $r_i$ (more importance) for draws $s$ where the probabilities are higher without the data points $i$ than with them
  - $r_i^{(s)} = p(\theta^{(s)} | x_{-i}, y_{-i}) / p(\theta^{(s)} | x,y) \propto 1 / p(y_i|x_i, \theta^{(s)})$
    - because $p(\theta^{(s)} | x,y) = \prod_i^n p(y_i|x_i, \theta^{(s)}) p(\theta^{(s)})$ and the difference for the ratio is that the numerator has one less data point ($-i$) than the denominator
    - therefore, everything cancels out to leave $1 / p(y_i|x_i, \theta^{(s)})$ (multiplied against a constant)
  - easier to work with logarithms: $\log(1/p(y_i|x_i, \theta^{(s)})) = -\log(\text{likelihood}[i])$
- importance-weighted predictive distribution: $p(y_i|x_i, x_{-i}, y_{-i}) \approx \sum_{s=1}^S (w_i^{(s)} p(y_i|x_i, \theta^{(s)}))$
  - where weight $w_i \leftarrow \text{PSIS}(r_i)$
    - Pareto smoothing to help stabilize the importance weights
- can have problem where very large importance weights make it that there are very few effective draws for this calculation
  - most draws will have little effect when there are a few dominating importance ratios
  - can look at the distribution of the weights to see how reliable the importance sampling is
  - compute the effective sample size $n_\text{eff}$ from the importance ratios
- Pareto smoothing
  - can approximate the extreme tail of the importance ratios using a generalized Pareto distribution
  - Pareto dist. has two parameters
    - *scale*
    - *shape* ($\hat{k}$): how many finite moments the weight distribution has
      - if the shape is less than 0.5, then the variance is finite
      - else, the variance may be infinite and CLT may not hold - means the approximation is not reliable
      - in practice, $\hat{k} < 0.7$ is okay
- plot the Pareto $\hat{k}$ and $n_\text{eff}$ for each data point (will usually be inversely correlated)

##### K-fold CV {-}

- varieties
  - can approximate LOO by removing some $n$ number of data points
  - leave-one-*group*-out for hierarchical models
  - leave-one-*block*-out for time series
    - remove all data from a certain time point or range

#### 9.2 Model comparison and selection {-}

- CV for model assessment
  - CV is good for model assessment when application specific utility/cost functions are used
  - also useful for PPC
    - model mis-specification diagnostics using Pareto-$\hat{k}$ and $p_\text{LOO}$ (effective number of parameters)
    - checking model calibration and LOO predictive error
- sometimes do not need CV
  - sometimes comparing predictive distribution to observed is enough to compare models

- can compare LLPD between models to see which data points are best modeled by each
- can compare models using ELPD
  - standard error of the difference in ELPD `se_diff` only works for normally-distributed error
    - not always a good assumption
    - only when the models are properly specific and the number of observations is large (>100 data points)

- What is one is not clearly better than others?
  - recommend "continuous expansion including all models"
    - make a model that contains both of the others as special cases
      - instead of a pooled model vs. separate model, use a hierarchical model
      - use a negative binomial instead of a poisson
    - sparse priors like regularized horseshoe prior instead of variable selection
  - model averaging with Bayesian Model Averaging (BMA) or Bayesian stacking
    - when continuous expansion is not an option
  - if nested
    - (e.g. use the model with 2 variables instead of the 10 variable model that contains those 2 variables)
    - use the simpler model assuming some cost for extra variables
    - or use the more complex model if you want to take into account all the uncertainties

- CV and model selection
  - can use CV for model selection if
    1. small number of models
    2. the difference between models is clear
  - do *not* use CV to choose from a large set of models
    - selection process leads to overfitting
    - more likely to chose the model that just happens to best fit the data
  - *overfitting in selection process is not unique for CV*

- **Take-home messages:**
  - it is good to think about predictions of observables, because observables are the only ones we can observe
  - CV can simulate predicting and observing new data
  - CV is good if you do not trust your model
  - different variants of CV are useful in different scenarios
  - CV has high variance, and **if** you trust your model you can beat CV in accuracy

#### 9.3 Variable selection with projection predictive variable selection {-}

- Rich model vs feature selection?
  - If we care only about the predictive performance:
    - include all available prior information
    - integrate over all uncertainties
    - no need for feature selection
  - Variable selection can be useful if:
    - need to reduce measurement or computation cost in the future
    - improve explainability
  - Two options for variable selection:
    1. find a minimal subset of features that yield a good predictive model
    2. identify all features that have predictive information
- shrinkage priors alone do not solve the problem
  - common strategy to fit a model with shrinkage prior and select variable my their marginal posteriors
    - (marginal posteriors = regression coefficients)
  - issues:
    1. marginal posteriors are difficult to interpret if features are correlated
      - their posteriors will be correlated
      - may each have a posterior around 0, but their *joint distribution* is not
    2. how to do post-selection inference correctly?
      - add average cumulative effect
      - set to zero and re-fit
      - etc...
- different approach: focus on predictive performance
  - two-step process:
    1. construct the best predictive model we can $\rightarrow$ *reference model*
    2. variable selection and post-selection inference $\rightarrow$ *projection*
  - instead of focusing on marginal posteriors, find the subset of features with closest predictive performance to the reference model
  - (watch lecture for description with the example)
