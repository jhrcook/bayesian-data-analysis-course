# Chapter 5 Exercises

2021-10-31



> Complete question 5.1.

## Question 1

**Exchangeability with known model parameters: For each of the following three examples, answer: (i) Are observations $y_1$ and $y_2$ exchangeable? (ii) Are observations $y_1$ and $y_2$ independent? (iii) Can we act *as if* the two observations are independent?**

**a) A box has one black ball and one white ball.**
**We pick a ball $y_1$ at random, put it back, and pick another ball $y_2$ at random.**

The observations are exchangeable, independent, and we can act as if they are independent.

**b) A box has one black ball and one white ball.**
**We pick a ball $y_1$ at random, we do not put it back, then we pick ball $y_2$.**

The observations are exchangeable because we don't have information about which ball is most likely to be picked first, but they are not independent because with $y_1$, we know the result for $y_2$.
I don't think we can treat the observations as independent

**c) A box has a million black balls and a million white balls.**
**We pick a ball $y_1$ at random, we do not put it back, then we pick ball $y_2$ at random.**

The observations are exchangeable for the same reason as in the answer to (b).
The observations are not independent because we will know that, after $y_1$, the other color is slightly more likely to be picked.
Since there are so many balls, we can likely treat the observations are independent.
