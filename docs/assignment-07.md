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

```stan
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
#> [1] 0.99239
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
#> 1      1          1     1 alpha    132.   
#> 2      1          1     1 beta      -0.424
#> 3      1          2     2 alpha    134.   
#> 4      1          2     2 beta      -1.51 
#> 5      1          3     3 alpha    133.   
#> 6      1          3     3 beta      -0.758
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
|31 |mu[1]        |  85.174| 0.89|  73.289|  95.922| 3962.294| 1.000|
|32 |mu[2]        | 105.190| 0.89|  98.463| 112.044| 3659.706| 1.000|
|33 |mu[3]        |  90.004| 0.89|  82.904|  97.480| 4128.795| 1.001|
|34 |mu[4]        | 110.647| 0.89| 105.726| 115.085| 3922.168| 1.000|
|35 |mu[5]        |  91.384| 0.89|  85.014|  97.826| 3421.298| 1.000|
|36 |mu[6]        |  90.850| 0.89|  81.474| 101.167| 4594.646| 1.000|
|37 |y6pred       |  89.881| 0.89|  62.157| 119.343| 3864.304| 1.000|
|1  |log_lik[1,1] |  -3.859| 0.89|  -4.382|  -3.349| 2880.706| 1.000|
|7  |log_lik[2,1] |  -3.983| 0.89|  -4.450|  -3.539| 4302.312| 1.000|
|13 |log_lik[3,1] |  -3.983| 0.89|  -4.450|  -3.539| 4666.343| 1.001|
|19 |log_lik[4,1] |  -6.258| 0.89|  -8.164|  -4.781| 5303.575| 1.000|
|25 |log_lik[5,1] |  -4.396| 0.89|  -5.095|  -3.780| 6571.813| 1.000|
|2  |log_lik[1,2] |  -4.008| 0.89|  -4.822|  -3.235| 7073.528| 1.000|
|8  |log_lik[2,2] |  -3.346| 0.89|  -3.885|  -2.863| 3151.064| 0.999|
|14 |log_lik[3,2] |  -3.684| 0.89|  -4.307|  -3.051| 2558.008| 1.000|
|20 |log_lik[4,2] |  -3.287| 0.89|  -3.779|  -2.840| 2982.817| 1.002|
|26 |log_lik[5,2] |  -4.984| 0.89|  -6.822|  -3.608| 5852.231| 1.000|
|3  |log_lik[1,3] |  -3.914| 0.89|  -4.775|  -3.330| 3720.636| 1.000|
|9  |log_lik[2,3] |  -3.418| 0.89|  -3.901|  -2.967| 3378.874| 0.999|
|15 |log_lik[3,3] |  -3.391| 0.89|  -3.886|  -2.948| 3151.064| 0.999|
|21 |log_lik[4,3] |  -3.426| 0.89|  -3.971|  -2.946| 3604.859| 1.000|
|27 |log_lik[5,3] |  -5.647| 0.89|  -7.856|  -4.104| 2865.449| 1.003|
|4  |log_lik[1,4] |  -3.255| 0.89|  -3.976|  -2.692| 3799.897| 1.000|
|10 |log_lik[2,4] |  -3.714| 0.89|  -4.681|  -2.849| 4830.580| 1.000|
|16 |log_lik[3,4] |  -3.201| 0.89|  -3.867|  -2.610| 4642.050| 1.000|
|22 |log_lik[4,4] |  -3.762| 0.89|  -5.035|  -2.964| 7079.724| 0.999|
|28 |log_lik[5,4] |  -3.201| 0.89|  -3.867|  -2.610| 2390.203| 1.000|
|5  |log_lik[1,5] |  -4.139| 0.89|  -5.126|  -3.290| 2913.758| 1.003|
|11 |log_lik[2,5] |  -3.431| 0.89|  -3.969|  -2.918| 5957.479| 0.999|
|17 |log_lik[3,5] |  -4.015| 0.89|  -5.186|  -3.243| 6571.813| 1.000|
|23 |log_lik[4,5] |  -4.139| 0.89|  -5.126|  -3.290| 3954.586| 1.000|
|29 |log_lik[5,5] |  -3.208| 0.89|  -3.721|  -2.745| 3639.411| 1.000|
|6  |log_lik[1,6] |  -5.915| 0.89|  -7.770|  -4.497| 4970.547| 1.000|
|12 |log_lik[2,6] |  -3.770| 0.89|  -4.215|  -3.311| 6311.004| 0.999|
|18 |log_lik[3,6] |  -4.136| 0.89|  -4.714|  -3.627| 3799.897| 1.000|
|24 |log_lik[4,6] |  -4.131| 0.89|  -4.729|  -3.588| 2680.640| 1.000|
|30 |log_lik[5,6] |  -3.964| 0.89|  -4.424|  -3.517| 4191.982| 1.000|

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
|31 |mu          | 93.641| 0.89| 88.251|  98.204| 2149.604| 1.002|
|32 |ypred       | 93.435| 0.89| 64.898| 121.674| 3804.014| 0.999|
|1  |log_lik[1]  | -3.982| 0.89| -4.211|  -3.773| 2506.400| 1.000|
|12 |log_lik[2]  | -3.802| 0.89| -4.022|  -3.608| 2839.754| 1.001|
|23 |log_lik[3]  | -3.802| 0.89| -4.022|  -3.608| 2839.754| 1.001|
|25 |log_lik[4]  | -7.497| 0.89| -9.029|  -6.120| 3273.206| 1.002|
|26 |log_lik[5]  | -4.950| 0.89| -5.474|  -4.504| 3097.572| 1.002|
|27 |log_lik[6]  | -4.682| 0.89| -5.108|  -4.300| 2559.576| 1.000|
|28 |log_lik[7]  | -4.185| 0.89| -4.436|  -3.938| 2540.364| 1.002|
|29 |log_lik[8]  | -4.470| 0.89| -4.805|  -4.143| 2549.949| 1.001|
|30 |log_lik[9]  | -3.976| 0.89| -4.193|  -3.773| 2842.081| 1.002|
|2  |log_lik[10] | -3.866| 0.89| -4.082|  -3.663| 2597.352| 1.001|
|3  |log_lik[11] | -3.886| 0.89| -4.098|  -3.695| 2904.081| 1.002|
|4  |log_lik[12] | -3.798| 0.89| -4.012|  -3.598| 2865.422| 1.001|
|5  |log_lik[13] | -3.802| 0.89| -4.022|  -3.608| 2839.754| 1.001|
|6  |log_lik[14] | -3.889| 0.89| -4.111|  -3.690| 2567.709| 1.000|
|7  |log_lik[15] | -4.950| 0.89| -5.474|  -4.504| 3097.572| 1.002|
|8  |log_lik[16] | -4.012| 0.89| -4.238|  -3.808| 2750.360| 1.002|
|9  |log_lik[17] | -4.840| 0.89| -5.323|  -4.407| 2568.657| 1.000|
|10 |log_lik[18] | -4.607| 0.89| -4.994|  -4.237| 2554.527| 1.000|
|11 |log_lik[19] | -3.913| 0.89| -4.127|  -3.718| 2888.267| 1.002|
|13 |log_lik[20] | -4.607| 0.89| -4.994|  -4.237| 2554.527| 1.000|
|14 |log_lik[21] | -4.151| 0.89| -4.407|  -3.921| 2591.254| 1.001|
|15 |log_lik[22] | -3.813| 0.89| -4.024|  -3.622| 2909.574| 1.002|
|16 |log_lik[23] | -3.944| 0.89| -4.169|  -3.755| 2867.237| 1.002|
|17 |log_lik[24] | -4.151| 0.89| -4.407|  -3.921| 2591.254| 1.001|
|18 |log_lik[25] | -3.802| 0.89| -4.022|  -3.608| 2839.754| 1.001|
|19 |log_lik[26] | -5.980| 0.89| -6.888|  -5.166| 3229.785| 1.002|
|20 |log_lik[27] | -3.802| 0.89| -4.022|  -3.608| 2839.754| 1.001|
|21 |log_lik[28] | -3.976| 0.89| -4.193|  -3.773| 2842.081| 1.002|
|22 |log_lik[29] | -4.251| 0.89| -4.503|  -3.972| 2774.348| 1.001|
|24 |log_lik[30] | -3.862| 0.89| -4.077|  -3.673| 2914.081| 1.002|

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
#> Warning: There were 24 divergent transitions after warmup. See
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
|32 |mu[1]        |  81.345| 0.89| 70.957|  91.416| 1570.827| 1.002|
|33 |mu[2]        | 102.578| 0.89| 92.830| 111.344| 2783.436| 1.000|
|34 |mu[3]        |  89.607| 0.89| 81.637|  98.634| 3589.968| 0.999|
|35 |mu[4]        | 106.470| 0.89| 96.713| 116.634| 2088.096| 1.000|
|36 |mu[5]        |  91.070| 0.89| 82.058|  99.991| 4352.602| 0.999|
|37 |mu[6]        |  88.586| 0.89| 79.606|  97.331| 3656.506| 1.000|
|1  |alpha        |  94.153| 0.89| 86.746| 101.832| 3187.614| 1.000|
|39 |tau          |  10.641| 0.89|  4.939|  17.460| 1400.133| 1.002|
|40 |y6pred       |  88.233| 0.89| 64.426| 113.590| 3835.410| 1.000|
|42 |ypred[1]     |  81.315| 0.89| 55.074| 105.141| 3397.156| 1.001|
|43 |ypred[2]     | 102.308| 0.89| 79.549| 129.232| 3614.848| 0.999|
|44 |ypred[3]     |  89.656| 0.89| 66.374| 115.144| 4141.933| 1.000|
|45 |ypred[4]     | 106.861| 0.89| 81.677| 131.402| 4063.645| 0.999|
|46 |ypred[5]     |  90.948| 0.89| 66.010| 114.804| 4091.596| 1.000|
|47 |ypred[6]     |  88.700| 0.89| 65.257| 113.780| 4314.855| 1.000|
|38 |mu7pred      |  93.868| 0.89| 74.147| 114.332| 3978.326| 1.000|
|41 |y7pred       |  93.879| 0.89| 64.091| 125.733| 4121.636| 1.000|
|2  |log_lik[1,1] |  -3.655| 0.89| -3.955|  -3.361| 2097.381| 1.002|
|8  |log_lik[2,1] |  -3.885| 0.89| -4.544|  -3.485| 3163.668| 1.000|
|14 |log_lik[3,1] |  -3.885| 0.89| -4.544|  -3.485| 3860.247| 1.000|
|20 |log_lik[4,1] |  -6.718| 0.89| -8.641|  -4.968| 2147.503| 1.002|
|26 |log_lik[5,1] |  -4.106| 0.89| -4.821|  -3.466| 4097.038| 0.999|
|3  |log_lik[1,2] |  -4.108| 0.89| -4.760|  -3.530| 4922.625| 1.000|
|9  |log_lik[2,2] |  -3.709| 0.89| -4.129|  -3.388| 3629.553| 1.000|
|15 |log_lik[3,2] |  -3.920| 0.89| -4.487|  -3.474| 2250.320| 1.000|
|21 |log_lik[4,2] |  -3.634| 0.89| -3.924|  -3.359| 2683.134| 1.003|
|27 |log_lik[5,2] |  -4.194| 0.89| -4.986|  -3.582| 2067.028| 1.000|
|4  |log_lik[1,3] |  -3.921| 0.89| -4.413|  -3.465| 2831.905| 1.000|
|10 |log_lik[2,3] |  -3.654| 0.89| -3.947|  -3.381| 3278.739| 1.001|
|16 |log_lik[3,3] |  -3.641| 0.89| -3.927|  -3.377| 3629.553| 1.000|
|22 |log_lik[4,3] |  -3.650| 0.89| -3.957|  -3.370| 2789.647| 1.000|
|28 |log_lik[5,3] |  -4.857| 0.89| -5.887|  -3.945| 2492.188| 1.004|
|5  |log_lik[1,4] |  -3.647| 0.89| -3.946|  -3.377| 1805.693| 1.000|
|11 |log_lik[2,4] |  -3.984| 0.89| -4.605|  -3.445| 3738.544| 1.000|
|17 |log_lik[3,4] |  -3.819| 0.89| -4.351|  -3.396| 4141.907| 1.000|
|23 |log_lik[4,4] |  -3.689| 0.89| -4.026|  -3.403| 4046.580| 1.001|
|29 |log_lik[5,4] |  -3.819| 0.89| -4.351|  -3.396| 2092.830| 1.000|
|6  |log_lik[1,5] |  -3.954| 0.89| -4.501|  -3.476| 2316.620| 1.003|
|12 |log_lik[2,5] |  -3.711| 0.89| -4.056|  -3.389| 3332.573| 1.001|
|18 |log_lik[3,5] |  -3.949| 0.89| -4.551|  -3.518| 4097.038| 0.999|
|24 |log_lik[4,5] |  -3.954| 0.89| -4.501|  -3.476| 3433.696| 1.000|
|30 |log_lik[5,5] |  -3.632| 0.89| -3.917|  -3.370| 1655.024| 1.002|
|7  |log_lik[1,6] |  -6.042| 0.89| -7.729|  -4.636| 3281.137| 1.001|
|13 |log_lik[2,6] |  -3.664| 0.89| -3.961|  -3.381| 4986.803| 1.000|
|19 |log_lik[3,6] |  -4.188| 0.89| -4.998|  -3.629| 1805.693| 1.000|
|25 |log_lik[4,6] |  -3.923| 0.89| -4.495|  -3.499| 2097.078| 1.000|
|31 |log_lik[5,6] |  -3.924| 0.89| -4.517|  -3.486| 4313.828| 1.000|

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
