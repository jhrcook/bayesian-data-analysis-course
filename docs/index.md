---
bibliography: ["citations.bib"]
biblio-style: "apalike"
link-citations: true
---

# Bayesian Data Analysis course {-}

This book serves as a collection of my notes and exercises completed for the [*Bayesian Data Analysis*](https://avehtari.github.io/BDA_course_Aalto/) course taught by Aki Vehtari.
It is normally taught as CS-E5710 at Aalto University, but the lectures and assignments have been made freely available online and the course is based around the text book [*Bayesian Data Analysis* (3rd ed.)](http://www.stat.columbia.edu/~gelman/book/) by Gelman *et al* [@bda3].

The source code for my course work is available on GitHub: [jhrcook/bayesian-data-analysis-course](https://github.com/jhrcook/bayesian-data-analysis-course).
If you see any problems (big or small), please let me know by opening an [Issue](https://github.com/jhrcook/bayesian-data-analysis-course/issues).
All R Markdown files are self-contained an can be executed on their own (i.e. do not need to be run in a series or specific order).

## Resources

- [Course website](https://avehtari.github.io/BDA_course_Aalto/)
- [2021 Schedule](https://avehtari.github.io/BDA_course_Aalto/Aalto2021.html#Schedule_2021)
- [GitHub repo](https://github.com/avehtari/BDA_course_Aalto) ([my fork](https://github.com/jhrcook/BDA_course_Aalto))
- free PDF of [*Bayesian Data Analysis* (3e)](https://users.aalto.fi/~ave/BDA3.pdf) (a.k.a "BDA3") ([exercise solutions](https://github.com/jhrcook/bayesian-data-analysis-course/blob/master/reading/bda3-exercise-solutions.pdf))
- [Chapter Notes](https://avehtari.github.io/BDA_course_Aalto/chapter_notes/BDA_notes.pdf)
- [Video lectures](https://aalto.cloud.panopto.eu/Panopto/Pages/Sessions/List.aspx#folderID=%22f0ec3a25-9e23-4935-873b-a9f401646812%22) or individually lists [here](https://avehtari.github.io/BDA_course_Aalto/#videos)
- [Lecture slides](https://github.com/avehtari/BDA_course_Aalto/tree/master/slides)

## How to study

> The following are recommendations from the course creators on how to take the course.

The recommended way to go through the material is:

1. Read the reading instructions for a chapter in the [chapter notes](https://avehtari.github.io/BDA_course_Aalto/chapter_notes/BDA_notes.pdf).
2. Read the chapter in BDA3 and check that you find the terms listed in the reading instructions.
3. Watch the corresponding [video lecture](https://aalto.cloud.panopto.eu/Panopto/Pages/Sessions/List.aspx#folderID=%22f0ec3a25-9e23-4935-873b-a9f401646812%22) to get explanations for most important parts.
4. Read corresponding additional information in the chapter notes.
5. Run the corresponding demos in [R demos](https://github.com/avehtari/BDA_R_demos) or [Python demos](https://github.com/avehtari/BDA_py_demos).
6. Read the exercise instructions and make the corresponding assignments. Demo codes in R demos and Python demos have a lot of useful examples for handling data and plotting figures. If you have problems, visit TA sessions or ask in course slack channel.
7. If you want to learn more, make also self study exercises listed below.

## Course sections

Below are the main sections of the course where each section should take about a week to complete.

1. Course Introduction
2. Basics of Bayesian Inference
3. Multidimensional Posterior
4. Monte Carlo
5. Markov chain Monte Carlo
6. HMC, NUTS, and Stan
7. Hierarchical models and exchangeability
8. Model checking & Cross-validation
9. Model comparison and selection
10. Decision analysis
11. Normal approximation & Frequency properties
12. Extended topics

## Additional notes

I took additional notes on other chapters in the book.
I appended these chapters after the primary course sections.

1. Ch. 14 Introduction to regression models
1. Ch. 15 Hierarchical linear models
1. Ch. 19 Parametric nonlinear models
1. Ch. 20 Basis function models
1. Ch. 21 Gaussian process models
1. Ch. 22 Finite mixture models
1. Ch. 23 Dirichlet process models

I plan to also read through and take notes on chapters 16 and 18 on "Generalized linear models" and "Models for missing data" in the future.

## Assignments and exercises

I completed the assignments that accompany the course and also did some of additional exercises from the book.
These are available in their own sections of the book: [Assignments](#assignments-intro) and [Exercises](#exercises-intro).

## Stan models

Below is a list of the models built during the course.
The original code is available in the GitHub repo ["models"](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/) directory and they are also copied in the [Models](#stan-models-1) section of this website.

- [Drug bioassay model](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment06-bioassay.stan) for Assignment 6
- [8 school SAT model](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/8-schools.stan)
- [Drownings](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment07-drownings.stan) for Assignment 7
- Factory machine measurements for Assignments 7 & 8:
  - [pooled](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment07_factories_pooled.stan)
  - [separate](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment07_factories_separate.stan)
  - [hierarchical](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/assignment07_factories_hierarchical.stan) (also used in Assignment 9)
- [serial dilution assay](https://github.com/jhrcook/bayesian-data-analysis-course/tree/master/models/serial-dilution.stan) for chapter 19 exercises
