# Variational Bayesian Inference for Group Testing Data — Simulation Code

This repository contains the R scripts for the simulation study associated with the paper “Variational Bayesian Inference for Group Testing Data: Comparing a New Hybrid Variational Inference Approach to Traditional Methods” (Zhao et al.).

The simulation study evaluates and compares the performance of Variational Bayesian (VB) and Markov chain Monte Carlo (MCMC) algorithms across multiple group-testing protocols under varying prevalence levels and sample sizes.

## Overview

This repository includes R code used to:

- Implement VB and MCMC algorithms for several group-testing protocols:

   - Master Pooling Testing (MPT) with known sensitivity and specificity

   - Individual Testing (IT) with known sensitivity and specificity

   - Dorfman Testing (DT) with unknown sensitivity and specificity

- Summarize key performance metrics, including bias, standard deviation, coverage probabilities, and computation time

- Generate the boxplot of computation time across different simulation scenarios
