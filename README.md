# IPDMA_power_binary
Riley RD, Hattle M, Collins GS, RWhittle R, Ensor J. Calculating the power to examine treatment-covariate interactions when planning an 
individual participant data meta-analysis of randomised trials with a binary outcome.

ABSTRACT
Before embarking on an individual participant data meta-analysis (IPDMA) project, researchers and funders need assurance it is worth their 
time and cost. This should include consideration of how many studies are promising their IPD and, given the characteristics of 
these studies, the power of an IPDMA including them. Here, we show how to estimate the power of a planned IPDMA of randomised 
trials to examine treatment-covariate interactions at the participant-level (i.e. treatment effect modifiers). We focus on a binary 
outcome with binary or continuous covariates, and propose a three-step approach, which assumes the true interaction size is common 
to all trials. In step one, the user must specify a minimally important interaction size and, for each trial separately (e.g. as 
obtained from trial publications), the following aggregate data: the number of participants and events in control and treatment 
groups, the mean and standard deviation for each continuous covariate, and the proportion of participants in each category for
each binary covariate. This allows the variance of the interaction estimate to be calculated for each trial, using an analytic 
solution for Fisher’s information matrix from a logistic regression model. Step two calculates the variance of the summary interaction 
estimate from the planned IPDMA (equal to the inverse of the sum of the inverse trial variances from step one), and step three 
calculates the corresponding power based on a two-sided Wald test. Stata and R code are provided, and two examples given for 
illustration. Extension to allow for between-study heterogeneity is also considered.

(1) IPDAMA_SS_binary.R : binary covariate
(2) IPDAMA_SS_continuous.R : continuous covariate
