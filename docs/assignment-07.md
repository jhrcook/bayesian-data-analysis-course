# Assignment 7

2021-10-18

**[Assignment 7](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-07.pdf)**

## Setup


```r
knitr::opts_chunk$set(echo = TRUE, comment = "#>", dpi = 300)

for (f in list.files(here::here("src"), pattern = "R$", full.names = TRUE)) {
  source(f)
}

library(rstan)
library(tidybayes)
library(magrittr)
library(tidyverse)

theme_set(theme_classic() + theme(strip.background = element_blank()))

options(mc.cores = 2)
rstan_options(auto_write = TRUE)

drowning <- aaltobda::drowning
factory <- aaltobda::factory
```

## 1. Linear model: drowning data with Stan

**The provided data `drowning` in the 'aaltobda' package contains the number of people who died from drowning each year in Finland 1980–2019.**
**A statistician is going to fit a linear model with Gaussian residual model to these data using time as the predictor and number of drownings as the target variable.**
**She has two objective questions:**

i) **What is the trend of the number of people drowning per year? (We would plot the histogram of the slope of the linear model.)**
ii) **What is the prediction for the year 2020? (We would plot the histogram of the posterior predictive distribution for the number of people drowning at $\tilde{x} = 2020$.)**

**Corresponding Stan code is provided in Listing 1.**
**However, it is not entirely correct for the problem.**
**First, there are *three mistakes*.**
**Second, there are no priors defined for the parameters.**
**In Stan, this corresponds to using uniform priors.**

**a) Find the three mistakes in the code and fix them.**
**Report the original mistakes and your fixes clearly in your report.**
**Include the full corrected Stan code in your report.**

1. Declaration of `sigma` on line 10 should be `real<lower=0>`.
2. Missing semicolon at the end of line 16.
3. On line 19, the prediction on new data does not use the new data in `xpred`. This has been changed to `real ypred = normal_rng(alpha + beta*xpred, sigma);`.

Below is a copy of the final model.
The full Stan file is at [models/assignment07-drownings.stan](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment07-drownings.stan).

```
data {
  int<lower=0> N;  // number of data points
  vector[N] x;     // observation year
  vector[N] y;     // observation number of drowned
  real xpred;      // prediction year
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;  // fix: 'upper' should be 'lower'
}
transformed parameters {
  vector[N] mu = alpha + beta*x;
}
model {
  y ~ normal(mu, sigma);  // fix: missing semicolor
}
generated quantities {
  real ypred = normal_rng(alpha + beta*xpred, sigma);  // fix: use `xpred`
}
```

**b) Determine a suitable weakly-informative prior $\text{Normal}(0,\sigma_\beta)$ for the slope $\beta$.**
**It is very unlikely that the mean number of drownings changes more than 50 % in one year.**
**The approximate historical mean yearly number of drownings is 138.**
**Hence, set $\sigma_\beta$ so that the following holds for the prior probability for $\beta$: $Pr(−69 < \beta < 69) = 0.99$.**
**Determine suitable value for $\sigma_\beta$ and report the approximate numerical value for it.**


```r
x <- rnorm(1e5, 0, 26)
print(mean(-69 < x & x < 69))
```

```
#> [1] 0.99212
```

```r
plot_single_hist(x, alpha = 0.5, color = "black") + geom_vline(xintercept = c(-69, 69)) + labs(x = "beta")
```

```
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="assignment-07_files/figure-html/unnamed-chunk-1-1.png" width="2100" />

**c) Using the obtained σβ, add the desired prior in the Stan code.**

From some trial and error, it seems that a prior of $\text{Normal}(0, 26)$ should work.
I have added this prior distribution to `beta` in the model at line 17.

```
beta ~ normal(0, 26);   // prior on `beta`
```

**d) In a similar way, add a weakly informative prior for the intercept alpha and explain how you chose the prior.**

To use the year directly as the values for $x$ would lead to a massive value of $\alpha$ because the values for $x$ range from 1980 to 2019.
Thus, it would be advisable to first center the year, meaning at the prior distribution for $\alpha$ can be centered around the average of the number of drownings per year and a standard deviation near that of the actual number of drownings.


```r
head(drowning)
```

```
#>   year drownings
#> 1 1980       149
#> 2 1981       127
#> 3 1982       139
#> 4 1983       141
#> 5 1984       122
#> 6 1985       120
```


```r
print(mean(drowning$drownings))
```

```
#> [1] 134.35
```

```r
print(sd(drowning$drownings))
```

```
#> [1] 28.48441
```

Therefore, I add the prior $\text{Normal}(135, 50)$ to $\alpha$ on line 16.

```
alpha ~ normal(135, 50);   // prior on `alpha`
```


```r
data <- list(
  N = nrow(drowning),
  x = drowning$year - mean(drowning$year),
  y = drowning$drownings,
  xpred = 2020 - mean(drowning$year)
)
drowning_model <- stan(
  here::here("models", "assignment07-drownings.stan"),
  data = data
)
```


```r
variable_post <- spread_draws(drowning_model, alpha, beta) %>%
  pivot_longer(c(alpha, beta), names_to = "variable", values_to = "value")
head(variable_post)
```

```
#> # A tibble: 6 × 5
#>   .chain .iteration .draw variable   value
#>    <int>      <int> <int> <chr>      <dbl>
#> 1      1          1     1 alpha    135.   
#> 2      1          1     1 beta      -1.08 
#> 3      1          2     2 alpha    134.   
#> 4      1          2     2 beta      -0.918
#> 5      1          3     3 alpha    135.   
#> 6      1          3     3 beta      -1.58
```


```r
variable_post %>%
  ggplot(aes(x = .iteration, y = value, color = factor(.chain))) +
  facet_grid(rows = vars(variable), scales = "free_y") +
  geom_path(alpha = 0.5) +
  scale_x_continuous(expand = expansion(c(0, 0))) +
  scale_y_continuous(expand = expansion(c(0.02, 0.02))) +
  labs(x = "iteration", y = "value", color = "chain")
```

<img src="assignment-07_files/figure-html/unnamed-chunk-6-1.png" width="2100" />


```r
variable_post %>%
  ggplot(aes(x = value)) +
  facet_grid(cols = vars(variable), scales = "free_x") +
  geom_histogram(color = "black", alpha = 0.3, bins = 30) +
  scale_x_continuous(expand = expansion(c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02)))
```

<img src="assignment-07_files/figure-html/unnamed-chunk-7-1.png" width="2100" />


```r
spread_draws(drowning_model, ypred) %$%
  plot_single_hist(ypred, alpha = 0.3, color = "black") +
  labs(x = "predicted number of drownings in 2020")
```

```
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="assignment-07_files/figure-html/unnamed-chunk-8-1.png" width="2100" />


```r
red <- "#C34E51"

bayestestR::describe_posterior(drowning_model, ci = 0.89, test = c()) %>%
  as_tibble() %>%
  filter(str_detect(Parameter, "mu")) %>%
  select(Parameter, Median, CI_low, CI_high) %>%
  janitor::clean_names() %>%
  mutate(idx = row_number()) %>%
  left_join(drowning %>% mutate(idx = row_number()), by = "idx") %>%
  ggplot(aes(x = year)) +
  geom_point(aes(y = drownings), data = drowning, color = "#4C71B0") +
  geom_line(aes(y = median), color = red, size = 1.2) +
  geom_smooth(
    aes(y = ci_low),
    method = "loess",
    formula = "y ~ x",
    linetype = 2,
    se = FALSE,
    color = red,
    size = 1
  ) +
  geom_smooth(
    aes(y = ci_high),
    method = "loess",
    formula = "y ~ x",
    linetype = 2,
    se = FALSE,
    color = red,
    size = 1
  ) +
  labs(x = "year", y = "number of drownings (mean ± 89% CI)")
```

<img src="assignment-07_files/figure-html/unnamed-chunk-9-1.png" width="2100" />

## 2. Hierarchical model: factory data with Stan

**The `factory` data in the 'aaltobda' package contains quality control measurements from 6 machines in a factory (units of the measurements are irrelevant here).**
**In the data file, each column contains the measurements for a single machine.**
**Quality control measurements are expensive and time-consuming, so only 5 measurements were done for each machine.**
**In addition to the existing machines, we are interested in the quality of another machine (the seventh machine).**

**For this problem, you’ll use the following Gaussian models:**

- **a separate model, in which each machine has its own model**
- **a pooled model, in which all measurements are combined and there is no distinction between machines**
- **a hierarchical model, which has a hierarchical structure as described in BDA3 Section 11.6**

**As in the model described in the book, use the same measurement standard deviation $\sigma$ for all the groups in the hierarchical model.**
**In the separate model, however, use separate measurement standard deviation $\sigma_j$ for each group $j$.**
**You should use weakly informative priors for all your models.**

**Complete the following questions for each of the three models (separate, pooled, hierarchical).**

**a) Describe the model with mathematical notation.**
**Also describe in words the difference between the three models.**

*Separate model*: The separate model is described below where each machine has its own centrality $\mu$ and dispersion $\sigma$ parameters that do not influence the parameters of the other machines.

$$
y_{ij} \sim N(\mu_j, \sigma_j) \\
\mu_j \sim N(0, 1) \\
\sigma_j \sim \text{Inv-}\chi^2(10)
$$

*Pooled model*: The pool model is described below where there is no distinction between the models but instead a single set of parameters for all of the data.

$$
y_{i} \sim N(\mu, \sigma) \\
\mu \sim N(0, 1) \\
\sigma \sim \text{Inv-}\chi^2(10)
$$

*Hierarchical model*: The hierarchical model is described below where each machine has its own centrality $\mu$ parameter which are linked through a hyper-prior distribution from which they are drawn.
The machines will all share a common dispersion paramete $\sigma$

$$
y_{ij} \sim N(\mu_j, \sigma_j) \\
\mu_j \sim N(\alpha, \tau) \\
\alpha \sim N(0, 1) \\
\tau \sim \text{HalfNormal}(2.5) \\
\sigma \sim \text{Inv-}\chi^2(10)
$$

The separate model is effectively building a different linear model for each machine where as the pooled model treats all the measurements as coming from the same model.
The hierarchical model is treating the machines as having come from a single, shared distribution.

**b) Implement the model in Stan and include the code in the report.**
**Use weakly informative priors for all your models.**


```r
print_model_code <- function(path) {
  for (l in readLines(path)) {
    cat(l, "\n")
  }
}
```

*Separate model*


```r
separate_model_code <- here::here(
  "models", "assignment07_factories_separate.stan"
)
print_model_code(separate_model_code)
```

```
#> data { 
#>   int<lower=0> N;  // number of data points per machine 
#>   int<lower=0> J;  // number of machines 
#>   vector[J] y[N];  // quality control data points 
#> } 
#>  
#> parameters { 
#>   vector[J] mu; 
#>   vector<lower=0>[J] sigma; 
#> } 
#>  
#> model { 
#>   // priors 
#>   for (j in 1:J) { 
#>     mu[j] ~ normal(100, 10); 
#>     sigma[j] ~ inv_chi_square(5); 
#>   } 
#>  
#>   // likelihood 
#>   for (j in 1:J){ 
#>     y[,j] ~ normal(mu[j], sigma[j]); 
#>   } 
#> } 
#>  
#> generated quantities { 
#>   // Compute the predictive distribution for the sixth machine. 
#>   real y6pred; 
#>   vector[J] log_lik[N]; 
#>  
#>   y6pred = normal_rng(mu[6], sigma[6]); 
#>  
#>   for (j in 1:J) { 
#>     for (n in 1:N) { 
#>       log_lik[n,j] = normal_lpdf(y[n,j] | mu[j], sigma[j]); 
#>     } 
#>   } 
#> }
```


```r
separate_model_data <- list(
  y = factory,
  N = nrow(factory),
  J = ncol(factory)
)
separate_model <- rstan::stan(
  separate_model_code,
  data = separate_model_data,
  verbose = FALSE,
  refresh = 0
)
knitr::kable(
  bayestestR::describe_posterior(separate_model, ci = 0.89, test = NULL),
  digits = 3
)
```



|   |Parameter    |  Median|   CI|  CI_low| CI_high|      ESS|  Rhat|
|:--|:------------|-------:|----:|-------:|-------:|--------:|-----:|
|31 |mu[1]        |  85.111| 0.89|  74.306|  96.601| 2713.665| 1.000|
|32 |mu[2]        | 105.263| 0.89|  98.325| 111.841| 3975.737| 0.999|
|33 |mu[3]        |  90.187| 0.89|  82.651|  97.276| 3363.204| 1.000|
|34 |mu[4]        | 110.726| 0.89| 106.159| 115.729| 3507.330| 1.001|
|35 |mu[5]        |  91.315| 0.89|  84.852|  97.457| 3387.003| 1.000|
|36 |mu[6]        |  90.724| 0.89|  81.465| 100.714| 3947.610| 1.000|
|37 |y6pred       |  90.870| 0.89|  62.936| 121.718| 4077.420| 1.001|
|1  |log_lik[1,1] |  -3.847| 0.89|  -4.405|  -3.393| 2050.755| 1.000|
|7  |log_lik[2,1] |  -3.979| 0.89|  -4.431|  -3.576| 4879.696| 1.000|
|13 |log_lik[3,1] |  -3.979| 0.89|  -4.431|  -3.576| 4429.343| 0.999|
|19 |log_lik[4,1] |  -6.266| 0.89|  -8.130|  -4.896| 5557.199| 1.001|
|25 |log_lik[5,1] |  -4.386| 0.89|  -5.030|  -3.759| 5096.882| 1.000|
|2  |log_lik[1,2] |  -4.006| 0.89|  -4.863|  -3.250| 6597.596| 0.999|
|8  |log_lik[2,2] |  -3.348| 0.89|  -3.909|  -2.880| 2587.751| 0.999|
|14 |log_lik[3,2] |  -3.688| 0.89|  -4.358|  -3.117| 2452.719| 1.000|
|20 |log_lik[4,2] |  -3.281| 0.89|  -3.774|  -2.801| 3074.626| 1.001|
|26 |log_lik[5,2] |  -4.972| 0.89|  -6.871|  -3.669| 5487.854| 0.999|
|3  |log_lik[1,3] |  -3.901| 0.89|  -4.723|  -3.337| 3825.434| 1.000|
|9  |log_lik[2,3] |  -3.423| 0.89|  -3.880|  -2.959| 3092.288| 1.000|
|15 |log_lik[3,3] |  -3.394| 0.89|  -3.868|  -2.943| 2587.751| 0.999|
|21 |log_lik[4,3] |  -3.437| 0.89|  -4.010|  -2.947| 3505.274| 1.000|
|27 |log_lik[5,3] |  -5.626| 0.89|  -7.689|  -4.092| 2961.958| 1.001|
|4  |log_lik[1,4] |  -3.277| 0.89|  -3.979|  -2.647| 3786.606| 1.000|
|10 |log_lik[2,4] |  -3.696| 0.89|  -4.675|  -2.832| 4896.409| 1.000|
|16 |log_lik[3,4] |  -3.199| 0.89|  -3.900|  -2.614| 4317.514| 1.000|
|22 |log_lik[4,4] |  -3.775| 0.89|  -5.014|  -2.892| 6923.344| 0.999|
|28 |log_lik[5,4] |  -3.199| 0.89|  -3.900|  -2.614| 2505.018| 1.000|
|5  |log_lik[1,5] |  -4.127| 0.89|  -5.201|  -3.309| 2777.028| 1.001|
|11 |log_lik[2,5] |  -3.410| 0.89|  -3.964|  -2.926| 5870.349| 1.000|
|17 |log_lik[3,5] |  -4.020| 0.89|  -5.144|  -3.216| 5096.882| 1.000|
|23 |log_lik[4,5] |  -4.127| 0.89|  -5.201|  -3.309| 3577.546| 1.000|
|29 |log_lik[5,5] |  -3.190| 0.89|  -3.677|  -2.722| 3005.680| 1.000|
|6  |log_lik[1,6] |  -5.888| 0.89|  -7.758|  -4.540| 5509.757| 0.999|
|12 |log_lik[2,6] |  -3.771| 0.89|  -4.234|  -3.314| 6725.445| 1.000|
|18 |log_lik[3,6] |  -4.147| 0.89|  -4.723|  -3.642| 3786.606| 1.000|
|24 |log_lik[4,6] |  -4.139| 0.89|  -4.700|  -3.565| 2685.671| 1.001|
|30 |log_lik[5,6] |  -3.978| 0.89|  -4.411|  -3.486| 3718.073| 1.000|

*Pooled model*


```r
pooled_model_code <- here::here("models", "assignment07_factories_pooled.stan")
print_model_code(pooled_model_code)
```

```
#> data { 
#>   int<lower=0> N;  // number of data points 
#>   vector[N] y;     // machine quality control data 
#> } 
#>  
#> parameters { 
#>   real mu; 
#>   real<lower=0> sigma; 
#> } 
#>  
#> model { 
#>   // priors 
#>   mu ~ normal(100, 10); 
#>   sigma ~ inv_chi_square(5); 
#>  
#>   // likelihood 
#>   y ~ normal(mu, sigma); 
#> } 
#>  
#> generated quantities { 
#>   real ypred; 
#>   vector[N] log_lik; 
#>  
#>   ypred = normal_rng(mu, sigma); 
#>  
#>   for (i in 1:N) 
#>     log_lik[i] = normal_lpdf(y[i] | mu, sigma); 
#>  
#> }
```


```r
pooled_model_data <- list(
  y = unname(unlist(factory)),
  N = length(unlist(factory))
)
pooled_model <- rstan::stan(
  pooled_model_code,
  data = pooled_model_data,
  verbose = FALSE,
  refresh = 0
)

knitr::kable(
  bayestestR::describe_posterior(pooled_model, ci = 0.89, test = NULL),
  digits = 3
)
```



|   |Parameter   | Median|   CI| CI_low| CI_high|      ESS|  Rhat|
|:--|:-----------|------:|----:|------:|-------:|--------:|-----:|
|31 |mu          | 93.516| 0.89| 88.619|  98.538| 2792.680| 1.000|
|32 |ypred       | 93.338| 0.89| 63.567| 121.323| 4090.508| 1.000|
|1  |log_lik[1]  | -3.978| 0.89| -4.190|  -3.753| 2834.900| 1.000|
|12 |log_lik[2]  | -3.796| 0.89| -3.990|  -3.586| 2552.985| 1.002|
|23 |log_lik[3]  | -3.796| 0.89| -3.990|  -3.586| 2552.985| 1.002|
|25 |log_lik[4]  | -7.500| 0.89| -8.975|  -6.111| 3175.958| 1.002|
|26 |log_lik[5]  | -4.947| 0.89| -5.454|  -4.497| 3189.124| 1.001|
|27 |log_lik[6]  | -4.694| 0.89| -5.126|  -4.305| 2997.643| 1.000|
|28 |log_lik[7]  | -4.190| 0.89| -4.444|  -3.942| 2833.338| 1.001|
|29 |log_lik[8]  | -4.482| 0.89| -4.830|  -4.161| 2951.303| 1.000|
|30 |log_lik[9]  | -3.977| 0.89| -4.195|  -3.776| 2707.695| 1.003|
|2  |log_lik[10] | -3.863| 0.89| -4.075|  -3.662| 2664.152| 1.001|
|3  |log_lik[11] | -3.888| 0.89| -4.075|  -3.680| 2665.731| 1.004|
|4  |log_lik[12] | -3.793| 0.89| -3.992|  -3.586| 2549.054| 1.003|
|5  |log_lik[13] | -3.796| 0.89| -3.990|  -3.586| 2552.985| 1.002|
|6  |log_lik[14] | -3.887| 0.89| -4.102|  -3.681| 2710.840| 1.001|
|7  |log_lik[15] | -4.947| 0.89| -5.454|  -4.497| 3189.124| 1.001|
|8  |log_lik[16] | -4.013| 0.89| -4.226|  -3.798| 2728.733| 1.003|
|9  |log_lik[17] | -4.853| 0.89| -5.367|  -4.429| 3019.529| 1.000|
|10 |log_lik[18] | -4.621| 0.89| -5.016|  -4.248| 2984.177| 1.000|
|11 |log_lik[19] | -3.914| 0.89| -4.119|  -3.717| 2676.073| 1.003|
|13 |log_lik[20] | -4.621| 0.89| -5.016|  -4.248| 2984.177| 1.000|
|14 |log_lik[21] | -4.142| 0.89| -4.399|  -3.910| 2963.757| 1.000|
|15 |log_lik[22] | -3.814| 0.89| -4.013|  -3.618| 2539.292| 1.003|
|16 |log_lik[23] | -3.944| 0.89| -4.156|  -3.746| 2690.057| 1.003|
|17 |log_lik[24] | -4.142| 0.89| -4.399|  -3.910| 2963.757| 1.000|
|18 |log_lik[25] | -3.796| 0.89| -3.990|  -3.586| 2552.985| 1.002|
|19 |log_lik[26] | -5.986| 0.89| -6.897|  -5.180| 3193.539| 1.002|
|20 |log_lik[27] | -3.796| 0.89| -3.990|  -3.586| 2552.985| 1.002|
|21 |log_lik[28] | -3.977| 0.89| -4.195|  -3.776| 2707.695| 1.003|
|22 |log_lik[29] | -4.242| 0.89| -4.517|  -3.984| 3028.028| 1.000|
|24 |log_lik[30] | -3.864| 0.89| -4.049|  -3.654| 2628.695| 1.004|

*Hierarchical model*


```r
hierarchical_model_code <- here::here(
  "models", "assignment07_factories_hierarchical.stan"
)
print_model_code(hierarchical_model_code)
```

```
#> data { 
#>   int<lower=0> N;  // number of data points per machine 
#>   int<lower=0> J;  // number of machines 
#>   vector[J] y[N];  // quality control data points 
#> } 
#>  
#> parameters { 
#>   vector[J] mu; 
#>   real<lower=0> sigma; 
#>   real alpha; 
#>   real<lower=0> tau; 
#> } 
#>  
#> model { 
#>   // hyper-priors 
#>   alpha ~ normal(100, 10); 
#>   tau ~ normal(0, 10); 
#>  
#>   // priors 
#>   mu ~ normal(alpha, tau); 
#>   sigma ~ inv_chi_square(5); 
#>  
#>   // likelihood 
#>   for (j in 1:J){ 
#>     y[,j] ~ normal(mu[j], sigma); 
#>   } 
#> } 
#>  
#> generated quantities { 
#>   // Compute the predictive distribution for the sixth machine. 
#>   real y6pred;  // Leave for compatibility with earlier assignments. 
#>   vector[J] ypred; 
#>   real mu7pred; 
#>   real y7pred; 
#>   vector[J] log_lik[N]; 
#>  
#>   y6pred = normal_rng(mu[6], sigma); 
#>   for (j in 1:J) { 
#>     ypred[j] = normal_rng(mu[j], sigma); 
#>   } 
#>  
#>   mu7pred = normal_rng(alpha, tau); 
#>   y7pred = normal_rng(mu7pred, sigma); 
#>  
#>   for (j in 1:J) { 
#>     for (n in 1:N) { 
#>       log_lik[n,j] = normal_lpdf(y[n,j] | mu[j], sigma); 
#>     } 
#>   } 
#> }
```


```r
hierarchical_model_data <- list(
  y = factory,
  N = nrow(factory),
  J = ncol(factory)
)
hierarchical_model <- rstan::stan(
  hierarchical_model_code,
  data = hierarchical_model_data,
  verbose = FALSE,
  refresh = 0
)
```

```
#> Warning: There were 32 divergent transitions after warmup. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
```

```
#> Warning: Examine the pairs() plot to diagnose sampling problems
```

```r
knitr::kable(
  bayestestR::describe_posterior(hierarchical_model, ci = 0.89, test = NULL),
  digits = 3
)
```



|   |Parameter    |  Median|   CI| CI_low| CI_high|      ESS|  Rhat|
|:--|:------------|-------:|----:|------:|-------:|--------:|-----:|
|32 |mu[1]        |  81.502| 0.89| 71.656|  91.170| 2666.504| 1.002|
|33 |mu[2]        | 102.520| 0.89| 92.735| 111.327| 3058.328| 1.001|
|34 |mu[3]        |  89.846| 0.89| 81.341|  98.753| 3882.504| 1.001|
|35 |mu[4]        | 106.568| 0.89| 96.595| 116.746| 2209.765| 1.001|
|36 |mu[5]        |  91.266| 0.89| 82.593|  99.788| 3592.183| 0.999|
|37 |mu[6]        |  88.523| 0.89| 79.551|  97.316| 4059.488| 0.999|
|1  |alpha        |  94.256| 0.89| 87.420| 102.841| 3417.937| 0.999|
|39 |tau          |  10.574| 0.89|  4.933|  17.387| 1643.252| 1.003|
|40 |y6pred       |  88.177| 0.89| 63.554| 112.444| 4098.415| 1.000|
|42 |ypred[1]     |  81.369| 0.89| 56.341| 105.538| 3544.553| 1.001|
|43 |ypred[2]     | 103.058| 0.89| 77.523| 127.228| 3675.682| 1.000|
|44 |ypred[3]     |  89.805| 0.89| 64.398| 112.847| 4492.295| 1.000|
|45 |ypred[4]     | 106.496| 0.89| 82.866| 131.399| 3414.140| 1.000|
|46 |ypred[5]     |  91.079| 0.89| 66.153| 115.148| 4135.596| 1.000|
|47 |ypred[6]     |  88.358| 0.89| 63.394| 113.503| 4218.295| 0.999|
|38 |mu7pred      |  94.467| 0.89| 74.222| 113.732| 3826.575| 1.001|
|41 |y7pred       |  94.836| 0.89| 65.364| 125.617| 3463.076| 1.000|
|2  |log_lik[1,1] |  -3.649| 0.89| -3.937|  -3.369| 2271.617| 1.001|
|8  |log_lik[2,1] |  -3.879| 0.89| -4.447|  -3.464| 3246.640| 1.001|
|14 |log_lik[3,1] |  -3.879| 0.89| -4.447|  -3.464| 4042.348| 1.000|
|20 |log_lik[4,1] |  -6.710| 0.89| -8.676|  -5.034| 1622.617| 1.001|
|26 |log_lik[5,1] |  -4.113| 0.89| -4.811|  -3.528| 3865.045| 0.999|
|3  |log_lik[1,2] |  -4.112| 0.89| -4.769|  -3.490| 4840.755| 1.000|
|9  |log_lik[2,2] |  -3.711| 0.89| -4.161|  -3.400| 3519.330| 1.002|
|15 |log_lik[3,2] |  -3.916| 0.89| -4.482|  -3.416| 2107.559| 1.002|
|21 |log_lik[4,2] |  -3.642| 0.89| -3.931|  -3.364| 2807.278| 1.000|
|27 |log_lik[5,2] |  -4.194| 0.89| -5.003|  -3.598| 2275.083| 1.001|
|4  |log_lik[1,3] |  -3.918| 0.89| -4.455|  -3.471| 2627.950| 1.001|
|10 |log_lik[2,3] |  -3.655| 0.89| -3.939|  -3.382| 2405.063| 1.001|
|16 |log_lik[3,3] |  -3.643| 0.89| -3.934|  -3.389| 3519.330| 1.002|
|22 |log_lik[4,3] |  -3.658| 0.89| -3.962|  -3.376| 2765.118| 1.001|
|28 |log_lik[5,3] |  -4.873| 0.89| -5.868|  -3.895| 2634.266| 1.001|
|5  |log_lik[1,4] |  -3.647| 0.89| -3.939|  -3.367| 1987.050| 1.001|
|11 |log_lik[2,4] |  -3.977| 0.89| -4.646|  -3.435| 3317.433| 1.000|
|17 |log_lik[3,4] |  -3.821| 0.89| -4.410|  -3.408| 3859.162| 1.000|
|23 |log_lik[4,4] |  -3.691| 0.89| -4.031|  -3.401| 4217.883| 1.000|
|29 |log_lik[5,4] |  -3.821| 0.89| -4.410|  -3.408| 1814.465| 1.002|
|6  |log_lik[1,5] |  -3.970| 0.89| -4.519|  -3.481| 2446.320| 1.003|
|12 |log_lik[2,5] |  -3.702| 0.89| -4.060|  -3.402| 2201.263| 1.001|
|18 |log_lik[3,5] |  -3.942| 0.89| -4.500|  -3.494| 3865.045| 0.999|
|24 |log_lik[4,5] |  -3.970| 0.89| -4.519|  -3.481| 3729.297| 0.999|
|30 |log_lik[5,5] |  -3.630| 0.89| -3.906|  -3.377| 2944.595| 1.000|
|7  |log_lik[1,6] |  -6.046| 0.89| -7.683|  -4.598| 4021.581| 1.000|
|13 |log_lik[2,6] |  -3.665| 0.89| -3.944|  -3.379| 4637.809| 1.000|
|19 |log_lik[3,6] |  -4.186| 0.89| -4.998|  -3.612| 1987.050| 1.001|
|25 |log_lik[4,6] |  -3.930| 0.89| -4.470|  -3.453| 2039.388| 1.001|
|31 |log_lik[5,6] |  -3.930| 0.89| -4.524|  -3.496| 3890.526| 1.000|

**c) Using the model (with weakly informative priors) report, comment on and, if applicable, plot histograms for the following distributions:**

i) **the posterior distribution of the mean of the quality measurements of the sixth machine.**
ii) **the predictive distribution for another quality measurement of the sixth machine.**
iii) **the posterior distribution of the mean of the quality measurements of the seventh machine.**


```r
plot_hist_mean_of_sixth <- function(vals) {
  plot_single_hist(vals, bins = 30, color = "black", alpha = 0.3) +
    labs(x = "mean of 6th machine", y = "posterior density")
}

plot_hist_sixth_predictions <- function(vals) {
  plot_single_hist(vals, bins = 30, color = "black", alpha = 0.3) +
    labs(x = "posterior predictions for 6th machine", y = "posterior density")
}

plot_hist_mean_of_seventh <- function(vals) {
  plot_single_hist(vals, bins = 30, color = "black", alpha = 0.3) +
    labs(x = "mean of 7thth machine", y = "posterior density")
}
```


*Separate model*


```r
plot_hist_mean_of_sixth(rstan::extract(separate_model)$mu[, 6])
```

<img src="assignment-07_files/figure-html/unnamed-chunk-18-1.png" width="2100" />


```r
plot_hist_sixth_predictions(rstan::extract(separate_model)$y6pred)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-19-1.png" width="2100" />

It is not possible to estimate the posterior for the mean of some new 7th machine because all machines are treated separately.

*Pooled model*


```r
plot_hist_mean_of_sixth(rstan::extract(pooled_model)$mu)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-20-1.png" width="2100" />


```r
plot_hist_sixth_predictions(rstan::extract(pooled_model)$ypred)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-21-1.png" width="2100" />

The predicted mean for a new machine is the same as the pooled mean $mu$.


```r
plot_hist_mean_of_seventh(rstan::extract(pooled_model)$mu)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-22-1.png" width="2100" />

*Hierarchical model*


```r
plot_hist_mean_of_sixth(rstan::extract(hierarchical_model)$mu[, 6])
```

<img src="assignment-07_files/figure-html/unnamed-chunk-23-1.png" width="2100" />


```r
plot_hist_sixth_predictions(rstan::extract(hierarchical_model)$y6pred)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-24-1.png" width="2100" />


```r
plot_hist_mean_of_seventh(rstan::extract(hierarchical_model)$mu7pred)
```

<img src="assignment-07_files/figure-html/unnamed-chunk-25-1.png" width="2100" />

**d) Report the posterior expectation for $\mu_1$ with a 90% credible interval but using a $\text{Normal}(0,10)$ prior for the $\mu$ parameter(s) and a $\text{Gamma}(1,1)$ prior for the $\sigma$ parameter(s).**
**For the hierarchical model, use the $\text{Normal}(0, 10)$ and $\text{Gamma}(1, 1)$ as hyper-priors.**

(I'm going to skip this one, but come back to it if it is needed for future assignments.)

---


```r
sessionInfo()
```

```
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#>  [1] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.7         
#>  [4] purrr_0.3.4          readr_2.0.1          tidyr_1.1.3         
#>  [7] tibble_3.1.3         tidyverse_1.3.1      magrittr_2.0.1      
#> [10] tidybayes_3.0.1      rstan_2.21.2         ggplot2_3.3.5       
#> [13] StanHeaders_2.21.0-7
#> 
#> loaded via a namespace (and not attached):
#>  [1] nlme_3.1-153         matrixStats_0.61.0   fs_1.5.0            
#>  [4] lubridate_1.7.10     insight_0.14.4       httr_1.4.2          
#>  [7] rprojroot_2.0.2      tensorA_0.36.2       tools_4.1.2         
#> [10] backports_1.2.1      bslib_0.2.5.1        utf8_1.2.2          
#> [13] R6_2.5.0             mgcv_1.8-38          DBI_1.1.1           
#> [16] colorspace_2.0-2     ggdist_3.0.0         withr_2.4.2         
#> [19] tidyselect_1.1.1     gridExtra_2.3        prettyunits_1.1.1   
#> [22] processx_3.5.2       curl_4.3.2           compiler_4.1.2      
#> [25] cli_3.0.1            rvest_1.0.1          arrayhelpers_1.1-0  
#> [28] xml2_1.3.2           bayestestR_0.11.0    labeling_0.4.2      
#> [31] bookdown_0.24        posterior_1.1.0      sass_0.4.0          
#> [34] scales_1.1.1         checkmate_2.0.0      aaltobda_0.3.1      
#> [37] callr_3.7.0          digest_0.6.27        rmarkdown_2.10      
#> [40] pkgconfig_2.0.3      htmltools_0.5.1.1    highr_0.9           
#> [43] dbplyr_2.1.1         rlang_0.4.11         readxl_1.3.1        
#> [46] rstudioapi_0.13      jquerylib_0.1.4      farver_2.1.0        
#> [49] generics_0.1.0       svUnit_1.0.6         jsonlite_1.7.2      
#> [52] distributional_0.2.2 inline_0.3.19        loo_2.4.1           
#> [55] Matrix_1.3-4         Rcpp_1.0.7           munsell_0.5.0       
#> [58] fansi_0.5.0          abind_1.4-5          lifecycle_1.0.0     
#> [61] stringi_1.7.3        yaml_2.2.1           snakecase_0.11.0    
#> [64] pkgbuild_1.2.0       grid_4.1.2           parallel_4.1.2      
#> [67] crayon_1.4.1         lattice_0.20-45      splines_4.1.2       
#> [70] haven_2.4.3          hms_1.1.0            knitr_1.33          
#> [73] ps_1.6.0             pillar_1.6.2         codetools_0.2-18    
#> [76] clisymbols_1.2.0     stats4_4.1.2         reprex_2.0.1        
#> [79] glue_1.4.2           evaluate_0.14        V8_3.4.2            
#> [82] renv_0.14.0          RcppParallel_5.1.4   modelr_0.1.8        
#> [85] vctrs_0.3.8          tzdb_0.1.2           cellranger_1.1.0    
#> [88] gtable_0.3.0         datawizard_0.2.1     assertthat_0.2.1    
#> [91] xfun_0.25            janitor_2.1.0        broom_0.7.9         
#> [94] coda_0.19-4          ellipsis_0.3.2       here_1.0.1
```
