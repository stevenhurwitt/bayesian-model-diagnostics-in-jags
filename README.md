# bayesian-model-diagnostics-in-jags

This R code consists of a bayesian model run in jags 
for a meta-analysis on 22 clinical trials 
of the effectiveness of beta blockers compared to a placebo.

Model is defined in jags with prior distributions on parameters.
MCMC is performed to simulate draws from the posterior distribution.
Diagnostic statistics such as AIC, DIC, WAIC & LOO-CV are calculated
for the resulting model.
