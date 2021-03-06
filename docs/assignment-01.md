# Assignment 1

2021-08-19

**[Assignment 1](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/course-material/assignment-01.pdf)**

## Setup


```r
knitr::opts_chunk$set(echo = TRUE, comment = "#>", dpi = 300)

library(glue)
```

## Exercise 1

**(Basic probability theory notation and terms). This can be trivial or you may need to refresh your memory on these concepts. Note that some terms may be different names for the same concept. Explain each of the following terms with one sentence:**

a) probability: how likely some assertion is to be true
b) probability mass: how likely a discrete random variable is to be some value
c) probability density: how likely a continuous random variable is to be some value
d) probability mass function (pmf): a function describing how likely a discrete random variable is to be any possible value
e) probability density function (pdf): a function describing how likely a continuous random variable is to be any possible value
f) probability distribution: a function describing how likely a random variable is to be of some value
g) discrete probability distribution: a probability distribution over discrete values
h) continuous probability distribution: a probability distribution over continuous values
i) cumulative distribution function (cdf): the cumulative sum of probabilities from a probability distribution
j) likelihood: how probable some event is under a given hypothesis

## Exercise 3

**(Bayes' theorem) A group of researchers has designed a new inexpensive and painless test for detecting lung cancer.**
**The test is intended to be an initial screening test for the population in general.**
**A positive result (presence of lung cancer) from the test would be followed up immediately with medication, surgery or more extensive and expensive test.**
**The researchers know from their studies the following facts:**

- Test gives a positive result in 98% of the time when the test subject has lung cancer.
- Test gives a negative result in 96 % of the time when the test subject does not have lung cancer.
- In general population approximately 1 person in 1000 has lung cancer.

**The researchers are happy with these preliminary results (about 97% success rate), and wish to get the test to market as soon as possible.**
**How would you advise them?**
**Base your answer on Bayes??? rule computations.**

**Hint: Relatively high false negative (cancer doesn???t get detected) or high false positive (unnecessarily administer medication) rates are typically bad and undesirable in tests.**

**Hint: Here are some probability values that can help you figure out if you copied the right conditional probabilities from the question.**

- $P(\text{Test gives positive} | \text{Subject does not have lung cancer}) = 0.04$
- $P(\text{Test gives positive and Subject has lung cancer}) = 0.00098$ this is also referred to as the joint probability of test being positive and the subject having lung cancer.

We are interested in the false positive and false negative rates.

The false positive rate, a positive test when the patient does not have cancer, is calculated below using Bayes' rule:

$$
\Pr(\text{no cancer} | \text{positive test}) = \frac{\Pr(\text{no cancer}) \Pr(\text{positive test} | \text{no cancer})}{\Pr(\text{positive test})}
$$

Each of the components:

$$
\Pr(\text{no cancer}) = 999/1000\\
\Pr(\text{positive test} | \text{no cancer}) = 4\% = 4/100\\
\begin{aligned}
\Pr(\text{positive test}) &= \Pr(\text{positive test AND cancer}) + \Pr(\text{positive test AND no cancer}) \\
&= \frac{98}{100} \frac{1}{1000} + \frac{4}{100} \frac{999}{1000} = \frac{4094}{100000}
\end{aligned}
$$

thus

$$
\begin{aligned}
\Pr(\text{no cancer} | \text{positive test}) &= \frac{\Pr(\text{no cancer}) \Pr(\text{positive test} | \text{no cancer})}{\Pr(\text{positive test})} \\
&= \frac{\frac{999}{1000} \frac{4}{100}}{\frac{4094}{100000}} \\
&= 0.976 \\
&= 97.6 \%
\end{aligned}
$$

The false positive rate, a negative test when the patient does have cancer, is calculated below using Bayes' rule:

$$
\Pr(\text{cancer} | \text{negative test}) = \frac{\Pr(\text{cancer}) \Pr(\text{negative test} | \text{cancer})}{\Pr(\text{negative test})}
$$

Each of the components:

$$
\Pr(\text{cancer}) = 1/1000\\
\Pr(\text{negative test} | \text{cancer}) = 1-0.98 = 2/100 \\
\begin{aligned}
\Pr(\text{negative test}) &= \Pr(\text{negative test AND cancer}) + \Pr(\text{negative test AND no cancer}) \\
&= \frac{2}{100} \frac{1}{1000} + \frac{96}{100} \frac{999}{1000} = \frac{95906}{100000}
\end{aligned}
$$

Thus,

$$
\begin{aligned}
\Pr(\text{cancer} | \text{negative test}) &= \frac{\Pr(\text{cancer}) \Pr(\text{negative test} | \text{cancer})}{\Pr(\text{negative test})} \\
&= \frac{\frac{1}{1000} \frac{2}{100}}{\frac{95906}{100000}} \\
&= 0.0000209 \\
&= 0.00209 \%
\end{aligned}
$$

The false negative rate is quite low, primarily because the cancer is relatively rare to begin with.
But, for the same reason, the false positive rate is very high.
Therefore, it could be advised that positive tests are confirmed with an independently-conducted second test or another testing procedure (e.g. CT scan).
Taking the test twice will only help if the false positive is caused by randomness and is not due to some other factor of the patient (e.g. the test is picking up some metabolite this patient produces due to their specific diet).

## Exercise 4

**We have three boxes, A, B, and C. There are**

- **2 red balls and 5 white balls in the box A,**
- **4 red balls and 1 white ball in the box B, and**
- **1 red ball and 3 white balls in the box C.**

**Consider a random experiment in which one of the boxes is randomly selected and from that box, one ball is randomly picked up.**
**After observing the color of the ball it is replaced in the box it came from.**
**Suppose also that on average box A is selected 40% of the time and box B 10% of the time (i.e. P (A) = 0.4).**

a) **What is the probability of picking a red ball?**
b) **If a red ball was picked, from which box it most probably came from?**

**Implement two functions in R that computes the probabilities.**

### 4.a

$$
\Pr(red) = \Sigma_{b}^{boxes} \Pr(\text{red} | \text{box}_b) \Pr(\text{box}_b)
$$


```r
# Contents of the boxes.
boxes <- matrix(
  c(2, 4, 1, 5, 1, 3),
  ncol = 2,
  dimnames = list(c("A", "B", "C"), c("red", "white"))
)
# Probability of selecting each box.
box_probs <- matrix(
  c(0.4, 0.1, 0.5),
  ncol = 1,
  dimnames = list(c("A", "B", "C"))
)

# Calculate the probability of red per box.
p_red_per_box <- function(boxes) {
  return(boxes[, "red"] / apply(boxes, 1, sum))
}

# Calculate the probability of pulling red.
p_red <- function(boxes, box_probs) {
  red_probs <- p_red_per_box(boxes)
  return(sum(red_probs * box_probs))
}

prob_red <- p_red(boxes, box_probs)
print(glue("probability of red: {prob_red}"))
```

```
#> probability of red: 0.319285714285714
```


```r
# Calculate the prob. that each box was used given a red ball was selected.
p_box <- function(boxes, box_probs) {
  prob_red <- p_red(boxes, box_probs)
  red_probs <- p_red_per_box(boxes)
  box_probs_red <- c()
  for (box in rownames(box_probs)) {
    p <- unlist(box_probs[box, ]) * unlist(red_probs[box]) / prob_red
    box_probs_red <- c(box_probs_red, p)
  }
  names(box_probs_red) <- rownames(box_probs)
  return(box_probs_red)
}

probs_of_boxes <- p_box(boxes = boxes, box_probs = box_probs)
stopifnot(sum(probs_of_boxes) == 1)
probs_of_boxes
```

```
#>         A         B         C 
#> 0.3579418 0.2505593 0.3914989
```

## Exercise 5

**Assume that on average fraternal twins (two fertilized eggs and then could be of different sex) occur once in 150 births and identical twins (single egg divides into two separate embryos, so both have the same sex) once in 400 births (Note! This is not the true values, see Exercise 1.6, page 28, in BDA3).**
**American male singer-actor Elvis Presley (1935 ??? 1977) had a twin brother who died in birth.**
**Assume that an equal number of boys and girls are born on average.**
**What is the probability that Elvis was an identical twin?**
**Show the steps how you derived the equations to compute that probability.**

State the problem as a conditional probability and use Bayes' rule to decompose into simpler terms.

$$
\begin{aligned}
\Pr(\text{Elvis was an identical twin}) &= \\
\Pr(\text{identical} | \text{twin AND brother}) &= \frac{\Pr(\text{identical})\Pr(\text{twin AND bro} |
\text{identical})}{\Pr(\text{twin AND bro})} \\
\end{aligned}
$$

Each component separately.

$$
\Pr(\text{identical}) = \frac{1}{400} \\
\Pr(\text{twin AND bro} | \text{identical}) = \frac{1}{2} \\
$$
$$
\begin{aligned}
\Pr(\text{twin AND bro}) &= \Pr(\text{identical AND bro}) + \Pr(\text{fraternal AND bro}) \\
&= \frac{1}{400} \frac{1}{2} + \frac{1}{150} \frac{1}{4} \\
&= \frac{7}{2400} \\
\end{aligned}
$$

Thus,

$$
\Pr(\text{identical} | \text{twin AND brother}) = \frac{\frac{1}{400} \frac{1}{2}}{\frac{7}{2400}} = 0.429
$$

**Implement this as a function in R that computes the probability.**


```r
p_identical_twin <- function(fraternal_prob, identical_prob) {
  (identical_prob * 0.5) / (identical_prob * 0.5 + fraternal_prob * 0.25)
}

# Tests from provided examples.
stopifnot(round(p_identical_twin(1 / 125, 1 / 300), 7) == 0.4545455)
stopifnot(round(p_identical_twin(1 / 100, 1 / 500), 7) == 0.2857143)

p_identical_twin(fraternal_prob = 1 / 150, identical_prob = 1 / 400)
```

```
#> [1] 0.4285714
```

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
#> [1] glue_1.4.2
#> 
#> loaded via a namespace (and not attached):
#>  [1] bookdown_0.24     clisymbols_1.2.0  crayon_1.4.1      digest_0.6.27    
#>  [5] R6_2.5.0          jsonlite_1.7.2    magrittr_2.0.1    evaluate_0.14    
#>  [9] stringi_1.7.3     rlang_0.4.11      renv_0.14.0       jquerylib_0.1.4  
#> [13] bslib_0.2.5.1     rmarkdown_2.10    tools_4.1.2       stringr_1.4.0    
#> [17] xfun_0.25         yaml_2.2.1        compiler_4.1.2    htmltools_0.5.1.1
#> [21] knitr_1.33        sass_0.4.0
```
