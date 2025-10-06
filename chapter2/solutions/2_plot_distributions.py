# Introduction to probability, statistics and hypothesis testing
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 2: Predictive probability



# Exercise 1: 
# Represent the PMF of the Bernoulli, Binomial and Poisson distribution
# Represent the PDF of the Gaussian distribution



# Import libraries ............................................................

import numpy as np
import matplotlib.pyplot as plt
from math import comb, exp, factorial, erf, sqrt
from scipy import stats



# Bernoulli distribution ......................................................

print("\nBernoulli distribution")

# Simulate a flip of coin
# x represents the face observed: heads (success, 1) or tails (failure, 0)
x = np.array([0, 1])

# Parameters
p = 1/2

# Compute Bernoulli probabilities
pmf_values = stats.bernoulli(p).pmf(x)
print(f"\nRandom variable x: {x}, format: {type(x)}, shape: {x.shape}")
print(f"Probability P(x): {pmf_values}, format: {type(pmf_values)}, shape: {pmf_values.shape}")

# Plot distribution
plt.bar(x, pmf_values, color = 'skyblue', edgecolor = 'black', width = 0.4)
plt.title(f"Bernoulli distribution (p = {p})")
plt.xlabel("x")
plt.ylabel("Bern(x; p)")
plt.show()



# Binomial distribution .......................................................

print("\nBinomial distribution")

# Simulate 10 rolls of dice
# x represents the possible number of sixes: 0, 1, 2, ..., 10
x = np.arange(11)

# Parameters
n, p = 10, 1/6

# Compute Binomial probability
pmf_values = stats.binom(n, p).pmf(x)
print(f"\nRandom variable x: {x}, format: {type(x)}, shape: {x.shape}")
print(f"Probability(x): {pmf_values}, format: {type(pmf_values)}, shape: {pmf_values.shape}")

# Plot distribution
plt.bar(x, pmf_values, color = 'skyblue', edgecolor = 'black', width = 0.6)
plt.title(f"Binomial distribution (n = {n}, p = {p:.5f})")
plt.xlabel("x")
plt.ylabel("B(x; n, p)")
plt.show()



# Poisson distribution ........................................................

print("\nPoisson distribution")


# Simulate a series of Poisson observations
# x represents the number of observations
x = np.arange(0, 15)

# Parameters
lmbda = 3

# Compute Poisson probability
pmf_values = stats.poisson(lmbda).pmf(x)
print(f"\nRandom variable x: {x}, format: {type(x)}, shape: {x.shape}")
print(f"Probability P(x): {pmf_values}, format: {type(pmf_values)}, shape: {pmf_values.shape}")

# Plot distribution
plt.bar(x, pmf_values, color = 'skyblue', edgecolor = 'black', width = 0.6)
plt.title(rf"Poisson distribution ($\lambda = {lmbda}$)")
plt.xlabel("x")
plt.ylabel(r"$P(x; \lambda)$")
plt.show()



# Gaussian distribution ........................................................

print("\nGaussian distribution")

# Parameters for gaussian
mu, sigma = 0, 1

# Prepare grid of points for plot
x = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 100)

# Compute Gaussian probability
pdf_values = stats.norm(mu, sigma).pdf(x)
print(f"\nRandom variable x: {x.shape}")
print(f"Probability density f(x): {pdf_values.shape}")

# Plot distribution
plt.plot(x, pdf_values, lw = 2)
plt.title(rf"Gaussian distribution ($\mu = {mu}, \sigma = {sigma}$)")
plt.xlabel("x")
plt.ylabel(r"$f\,(x; \mu, \sigma)$")
plt.show()


