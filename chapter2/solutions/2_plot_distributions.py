# Introduction to probability, statistics and hypothesis testing
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 2: Predictive probability



# Exercise 1: 
# Represent the PMF of the Bernoulli, Binomial and Poisson distribution
# Represent the PDF of the Gaussian distribution
# Compute mean and variance and add it to the plots


# Import libraries ............................................................

import numpy as np
import matplotlib.pyplot as plt
from math import comb, exp, factorial, erf, sqrt
from scipy import stats



# Binomial distribution .......................................................

# Simulate 10 rolls of dice and plot probability distribution
x = np.arange(11)
sixes = stats.binom(10, 1/6)
plt.plot(x, sixes.pmf(x), "ro", ms = 8)
plt.vlines(x, 0, sixes.pmf(x), colors = "r", lw = 4)
plt.xlabel("x"); plt.ylabel("B(x)")
plt.title("Binomial distribution")
plt.show()

# Compute mean and variance
mean = sixes.mean()
var = sixes.var()
print("mean = ", round(mean, 2), "\nvariance = ", round(var, 2))



# Poisson distribution ........................................................

# Plot probability distribution
impacts = stats.poisson(4) # e.g. an average of 4 meteorite impacts per year.
x = np.arange(16)
plt.plot(x, impacts.pmf(x), "ro", ms = 8)
plt.vlines(x, 0, impacts.pmf(x), colors = "r", lw = 4)
plt.xlabel("x"); plt.ylabel("P(x)")
plt.title("Poisson distribution")
plt.show()

# Compute mean and variance
mean = impacts.mean()
var = impacts.var()
print("mean = ", round(mean, 2), "\nvariance = ", round(var, 2))



# Gaussian distribution .......................................................

# Plot gaussian distrubtion
x = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 1000)
gaussian_distribution = stats.norm.pdf(x, mu, sigma)

# Plot the Gaussian distribution
plt.plot(x, gaussian_distribution, label = f"mu = {mu}, sigma = {sigma}")
plt.xlabel('x'); plt.ylabel('Probability Density')
plt.legend(); plt.grid(True)
plt.title('Gaussian Distribution')
plt.show()
