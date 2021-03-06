# Section 2. Basics of Bayesian inferences

2021-08-21



## Resources

- BDA3 chapter 2 and [reading instructions](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/BDA3_ch02_reading-instructions.pdf)
- lectures:
  - ['2.1 Basics of Bayesian inference, observation model, likelihood, posterior and binomial model, predictive distribution and benefit of integration'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=9c271082-5a8c-4b66-b6c2-aacc00fc683f)
  - ['2.2 Priors and prior information, and one parameter normal model'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=70655a8a-0eb4-4ddd-9f52-aacc00fc67a2)
  - ['Extra explanations about likelihood, normalization term, density, and conditioning on model M'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=158d119d-8673-4120-8669-ac3900c13304) (optional)
  - ['Summary 2.1. Observation model, likelihood, posterior and binomial model'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=7a297f7d-bb7b-4dd0-9913-a9f500ec822d) (optional)
  - ['Summary 2.2. Predictive distribution and benefit of integration'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=75b9f18f-e379-4557-a5fa-a9f500f11b40) (optional)
  - ['Summary 2.3. Priors and prior information'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=099659a5-f707-473d-8b03-a9f500f39eb5) (optional)
- [slides](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/slides_ch2.pdf) and [extra slides](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/slides_ch2_extra1.pdf)
- [Assignment 2](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-02.pdf)

## Notes

### Chapter instructions

- recommendations about weakly informative priors has changed a bit
  - updated recommendations: [Prior Choice Recommendations](https://github.com/ stan-dev/stan/wiki/Prior-Choice-Recommendations)
  - "5 levels of priors":
    1. **Flat prior** (not usually recommended)
    2. **Super-vague** but proper prior: $N(0, 10^6)$ (not usually recommended)
    3. **Weakly informative prior**: very weak; $N(0, 10)$
    4. **Generic weakly informative prior**: $N(0, 1)$
    5. **Specific informative prior**: $N(0.4, 0.2)$ or whatever; can sometimes be expressed as a scaling followed by a generic prior: $\theta = 0.4 + 0.2z; \text{ } z \sim N(0, 1)$
  - "flat and super-vague priors are not usually recommended"
  - even a seemingly weakly informative prior could be informative
    - e.g. a prior of $N(0, 1)$ could put weight on values too *large* if a large effect size would only be on the scale of 0.1
  - def. **weakly informative**: "if there's a reasonably large amount of data, the likelihood will dominate, and the prior will not be important"
  - section on [General Principles](https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#general-principles); some stand-outs copied here:
    - "Computational goal in Stan: reducing instability which can typically arise from bad geometry in the posterior"
    - "Weakly informative prior should contain enough information to regularize: the idea is that the prior rules out unreasonable parameter values but is not so strong as to rule out values that might make sense"
    - "When using informative priors, be explicit about every choice; write a sentence about each parameter in the model."

### Chapter 2. Single-parameter models

> I took notes in the book, so below are just some main points.

#### 2.2 Posterior as compromise between data and prior information {-}

- "general feature of Bayesian inference: the posterior distribution is centered at a point that represents a compromise between the prior information and the data"

#### 2.3 Estimating a probability from binomial data {-}

- a key benefit of Bayesian modeling is the flexibility of summarizing posterior probabilities
  - can be used to answer the key research questions
- commonly used summary statistics
  - centrality: mean, median, mode
  - variation: standard deviation, interquartile range, highest posterior density

#### 2.4 Informative prior distributions {-}

- *hyperparameter*: parameter of a prior distribution
- *conjugacy*: "the property that the posterior distribution follows the same parameter form as the prior distribution"
  - e.g. the beta prior is a conjugate family for the binomial likelihood
  - e.g. the gamma prior is a conjugate family for the Poisson likelihood
  - convenient because the posterior follows a known parametric family
  - formal definition of conjugacy:

$$
p(\theta | y) \in \mathcal{P} \text{ for all } p(\cdot | \theta) \in \mathcal{F} \text{ and } p(\cdot) \in \mathcal{P}
$$

#### 2.5 Normal distribution with known variance {-}

- *precision* (when discussing normal distributions): the inverse of the variance $\frac{1}{\tau^2}$

#### 2.6 Other standard single-parameter models {-}

- Poisson model for count data
  - data $y$ is the number of positive events
  - unknown rate of the events $\theta$
  - conjugate prior is the gamma distribution
  section 2.7 is a good example of a hierarchical Poisson model

#### 2.8 Noninformative prior distributions {-}

> See more information in the notes from the [chapter instructions](#chapter-instructions).

- rationale: let the data speak for themselves; inferences are unaffected by external information/bias
- problems:
  - can cause the posterior to become improper
  - computationally, makes it harder to sample from the posterior

#### 2.9 Weakly informative prior distributions {-}

> See more information in the notes from the [chapter instructions](#chapter-instructions).

- **weakly informative**: the prior is proper, but intentionally weaker than whatever actual prior knowledge is available
- "in general, any problem has some natural constraints that would allow a weakly informative model"
- small amount of real-world knowledge to ensure the posterior makes sense

### Lecture notes

#### 2.1 Basics of Bayesian inference, observation model, likelihood, posterior and binomial model, predictive distribution and benefit of integration {-}

- predictive distribution
  - "integrate over uncertainties"
  - for the example of pulling red or yellow chips out of a bag:
    - want a predictive distribution for some data point $\tilde{y} = 1$:
    - if we know $\theta$ then it is easy: $p(\tilde{y} = 1 | \theta, y, n, M)$
      - where $n$ is number of draws, $y$ is number of success (red chip), $M$ is model
    - we don't know $\theta$, we weight the probability of the new data for a given $\theta$ by the posterior probability that $\theta$ is that value
      - sum (integrate) over all possible values for $\theta$ ("integrate out the uncertainty of $\theta$")
      - $p(\tilde{y}=1|y, n, M) = \int_0^1 p(\tilde{y} = 1 | \theta, y, n, M) p(\theta | y, n, M) d\theta$
    - now the prediction is not conditioned on $\theta$, just on what was observed
- prior predictive: predictions before seeing any data
  - $p(\tilde{y}=1|M) = \int_o^1 p(\tilde{y}=1 | \theta, y, n, M) p(\theta|M)$

#### 2.2 Priors and prior information, and one parameter normal model {-}

- proper prior: $\int p(\theta) = 1$
  - better to use proper priors
- improper prior density does not have a finite integral
  - the posterior can sometimes still be proper, though
  - uniform distributions to infinity are improper
- a weak prior is *not* non-informative
  - could give a lot of a weight to very unlikely (or impossible) values
  - make sure to check prior values against knowable values
- *sufficient statistic*: the quantity $t(y)$ is a sufficient statistic for $\theta$ because the likelihood for $\theta$ depends on the data $y$ only through the value of $t(y)$
  - smaller dimensional data that fully summarizes the full data
  - can define a Gaussian with just the mean and s.d.

#### Extras: likelihood, normalization term, density, and conditioning on model M {-}

#### Predictive distribution and benefit of integration {-}

- predictive dist.
  - effect of integration
    - predictive dist of new $\hat{y}$ (discrete) with model $M$:
      - if we know $\theta$: $p(\hat{y} = 1| y, n, M) = p(\hat{y} = 1 | \theta, y, n, M)$
      - if we don't know $\theta$: $p(\hat{y} = 1 | \theta, y, n, M) p(\theta| y, n, M)$
        - weight by the probability of the value for $\theta$
    - integrate over all possible values of $\theta$: $p(\hat{y} = 1|y, n, M) = \int_0^1 p(\hat{y}=1| \theta, y, n, M) p(\theta|y, n, M)d\theta$
      - "integrate out the uncertainty over $\theta$"
- prior predictive for new data $\hat{y}$: $p(\hat{y} =1|M) = int_0^1 p(\hat{y}=1|\theta,y,n,M)p(\theta|M)$

#### Priors and prior information {-}

- conjugate priors do not result in any computational benefits in HMC or NUTS
  - can be useful to analytically reduce the size of a model, beforehand, though
