# Section 10. Decision analysis

2021-11-15



## Resources

- reading
  - BDA3 chapter 9
  - [reading instructions](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/BDA3_ch09_reading-instructions.pdf)
- lectures:
  - ['10.1 Decision analysis'](https://aalto.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=82943720-de0f-4195-8639-ab0900ca2085)
- slides:
  - [Lecture 10.1](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/slides_ch9.pdf)
- [Assignment 9](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-09.pdf)

## Notes

### Reading instructions

- outline of chapter 9
  - 9.1 Context and basic steps (most important part)
  - 9.2 Example
  - 9.3 Multistage decision analysis (you may skip this example)
  - 9.4 Hierarchical decision analysis (you may skip this example)
  - 9.5 Personal vs. institutional decision analysis (important)
- the lectures have simpler examples and discuss  some challenges in selecting utilities or costs
- ch 7 discusses how model selection con be considered as a decision problem

### Chapter 9. Decision analysis

- how can inferences be used in decision making?
- examples in this chapter:
  1. section 9.2: simple example with hierarchical model on how incentives affect survey response rates
    - compare expected response rates of various incentive structures to their expected cost
  2. section 9.3: option of performing a diagnostic test before deciding on a treatment for cancer
    - example of "value of information" and balancing risks of the screening test against the information it would provide
  3. section 9.4: decision and utility analysis of the risk of radon exposure
    - cost of measurement and fixing high exposure
    - example of a full integration if inference with decision analysis

#### 9.1 Bayesian decision theory in difference contexts {-}

- use Bayesian inference in two ways when balancing costs and benefits of decision options under uncertainty:
  1. a decision depends on the predicted quantities which depend on the parameters of the model and type of data
  2. use Bayesian inference within a decisions analysis to estimate outcomes conditional on information from previous decisions

##### Bayesian inference and decision trees {-}

- decision analysis involves optimization over decisions and uncertainties
- **Bayesian decision analysis** is defined as the following steps:
  1. Enumerate the space of all possible decisions $d$ and outcomes $x$.
  2. Determine the probability distribution of $x$ for each decision option $d$.
  3. Define a *utility function* $U(x)$ mapping outcomes onto real numbers (values of interest).
  4. Compute the expected utility $\text{E}(U(x)|d)$ as a function of the decision $d$ and choose the decision with the highest expected utility.
- often, we only do the first two steps and the rest is left to the "decision makers"

### Lecture notes

#### 10.1 Decision analysis {-}

(no new notes)
