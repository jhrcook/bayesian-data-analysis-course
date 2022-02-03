# Section 1. Course introduction and prerequisites

2021-08-17



## Resources

- BDA chapter 1 and [reading instructions](../reading-instructions/BDA3_ch01_reading-instructions.pdf)
- lectures:
  - ['Computational probabilistic modeling'](https://www.youtube.com/watch?v=ukE5aqdoLZI)
  - ['Introduction to uncertainty and modelling'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=d841f429-9c3d-4d24-8228-a9f400efda7b)
- video: ['Introduction to the course contents'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=13fc7889-cfd1-4d99-996c-a9f400f6e5a2)
- [slides](../slides/bayes_intro.pdf)
- [Assignment 1](assignments/assignment1.pdf)

## Notes

### Reading instructions

- model vs. likelihood for $p(y|\theta, M)$
  - *model* when the function is in terms of $y$
    - should be written as $p_y(y|\theta, M)$
    - used to describe uncertainty about $y$ given values of $\theta$ and $M$
  - *likelihood* when the function is in terms of $\theta$
    - should be written as $p_\theta(y|\theta, M)$
    - the posterior distribution describes the probability for different values of $\theta$ given fixed values for $y$
    - "The likelihood function is unnormalized probability distribution describing uncertainty related to $\theta$ (and that’s why Bayes rule has the normalization term to get the posterior distribution)."
- *exchangeability*
  1. independence is stronger condition than exchangeability
  2. independence implies exchangeability
  3. exchangeability does not imply independence
  4. exchangeability is related to what information is available instead of the properties of unknown underlying data generating mechanism

### Lecture notes

#### Introduction to uncertainty and modelling {-}

- two types of uncertainty:
  - *aleatoric*: due to randomness
  - *epistemic*: due to lack of knowledge
- model vs. likelihood:
  - model
    - $p_y(y|\theta, M)$
    - a function of $y$ given fixed values of $\theta$
    - describes *aleatoric* uncertainty
  - likelihood
    - $p_\theta(y|\theta, M)$
    - function of $\theta$ given fixed values of $y$
    - provides information about the *epistemic* uncertainty
    - is *not* a probability distribution
  - Bayes rule combines the likelihood with prior uncertainty to update the posterior uncertainty
  - example with a bag containing red and yellow chips:
    - probability of red = #red / #red + #yellow = $\theta$
    - $p(y = \text{red} | \theta)$: aleatoric uncertainty
      - predicting the probability of pulling a red chip has uncertainty due to randomness even if we new $\theta$ exactly
    - $p(\theta)$: epistemic uncertainty
      - we don't know $\theta$ but could compute it exactly if we knew the contents of the bag

#### Introduction to the course contents {-}

- benefits of Bayesian approach
  1. integrate over uncertainties to focus on interesting parts
  2. use relevant prior information
  3. hierarchical models
  4. model checking and evaluation
