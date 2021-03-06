# Section 17. Notes on 'Ch 19. Parametric nonlinear models'

2021-12-09



> These are just notes on a single chapter of *BDA3* that were not part of the course.

## Chapter 19. Parametric nonlinear models

- examples:
  - a ratio; $\text{E}(y = \frac{a_1 + b_1 x_1}{a_2 + b_2 x_2}$
  - sum of nonlinear functions: $\text{E}(y) = A_1 e^{-\alpha_1x} + A_2 e^{-\alpha_2x}$
- parameters of a nonlinear model are often harder to interpret
  - often requires custom visualization techniques
- **"Generally each new modeling problem must be tackled afresh."** (pg 471)
  - these models are less systematic than linear modeling

### 19.1 Example: serial dilution assay {-}

- estimate 10 unknown concentrations of an allergen based off of serial dilutions of a known standard

#### The model {-}

- **Notation**:
  - parameters of interest: concentrations of unknown samples $\theta_1, \dots, \theta_10$
  - known concentration of the standard $\theta_0$
  - dilution of measure $i$ as $x_i$ and color intensity (measurement) as $y_i$
- **Curve of expected measurements given the concentration**
  - use the following equation that is standard in the field
  - parameters:
    - $\beta_1$: color intensity at the limit of 0 concentration
    - $\beta_2$: the increase to saturation
    - $\beta_3$: concentration at which the gradient of the curve turns
    - $\beta_4$: rate at which saturation occurs

$$
\text{E}(y | x, \beta) =
  g(x, \beta) =
  \beta_1 + \frac{\beta_2}{1 + (x / \beta_3)^{-\beta_4}}
$$

- **Measurement error**
  - modeled as normally distributed with unequal variances
  - parameters:
    - $\alpha$: models the pattern that variances are higher for larger measurements
      - restricted $[0, 1]$
    - $A$ is a arbitrary constant to scale the data so $\sigma$ can be interpreted as the deviation from "typical" values
    - $\sigma$: deviation of a measure from the "typical"

$$
y_i \sim \text{N}(g(x_i, \beta), (\frac{g(x_i, \beta)}{A})^{2\alpha} \sigma_y^2)
$$

- **Dilution errors**
  - two possible sources:
    1. *initial dilution*: the accuracy of the creation of the initial standard concentration
    2. *serial dilutions*: error in creation of the subsequent dilutions (low enough to ignore for this analysis)
  - use a normal model on the log scale of the initial dilution error
  - parameters:
    - $\theta_0$: known concentration of the standard solution
    - $d_0^\text{init}$: known initial dilution of the standard that is called for
      - without error, the concentration of the initial solution would be $d_0^\text{init} \theta_0$
    - $x_0^\text{init}$: the *actual* (unknown) concentration of the initial dilution

$$
\log(x_0^\text{init}) \sim \text{N}(\log(d_0^\text{init} \cdot \theta_0), (\sigma^\text{init})^2)
$$

- **Dilution errors** (cont)
  - there is no initial dilution for the unknown samples being tested
    - therefore, the unknown initial concentration for sample $j$ is $x^\text{init} = \theta_j$ for $j = 1, \dots, 10$
    - for the dilutions of the unknown samples, set $x_i = d_i \cdot x_{j(i)}^\text{init}$
      - $j(i)$ is the sample $j$ corresponding to measurement $i$
      - $d_i$ is the dilution of measurement $i$ relative to the initial concentration

#### Prior distributions {-}

- priors used as described by book (are likely different than what would be recommended now):
  - $\log(\beta) \sim U(-\infty, \infty)$
  - $\sigma_y \sim U(0, \infty)$
  - $\alpha \sim U(0,1)$
  - $p(\log \theta_j) \propto 1$ for each unknown $j = 1, \dots, 10$
- cannot estimate $\sigma^\text{init}$ because we only have a single standard
  - use a fixed value of 0.02 based on a previous analysis of different plates

### 19.2 Example: population toxicokinetics {-}

- this is a more complex model
- uses a physiological model with parameters that cannot be solely determined using the data
  - requires informative priors based on previous studies
