# Introduction to probability, statistics and hypothesis testing
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 2: Predictive probability



# Exercise 1: 
# Implement a manual calculation of the Bernoulli, Binomial, Poisson and Gaussian probabilty
# Check calculation with numpy / scipy implementation



# Import libraries ............................................................

import numpy as np
import matplotlib.pyplot as plt
from math import comb, exp, factorial, erf, sqrt
from scipy import stats



# Bernoulli distribution ......................................................

def bernoulli_probability(x, p):
    """
    Calculate the Bernoulli probability of success (x=1) or failure (x=0)
    with a probability of success p.
    """
    if x == 1:  # success
        return p
    elif x == 0:  # failure
        return 1 - p
    else:
        raise ValueError("x must be 0 (failure) or 1 (success)")

print("\nBernoulli distribution")

# Probability of success (1) in a coin flip
x = 1; p = 1/2
prob_success = bernoulli_probability(x, p)
prob_success_scipy = stats.bernoulli.pmf(x, p)
print(f"\nProbability of heads (success) (manual): {prob_success}")
print(f"Probability of heads (success) (scipy): {prob_success_scipy}")

# Probability of failure (0) in a coin flip
x = 0; p = 1/2
prob_failure = bernoulli_probability(x, p)
prob_failure_scipy = stats.bernoulli.pmf(x, p)
print(f"\nProbability of tails (failure) (manual): {prob_failure}")
print(f"Probability of tails (failure) (scipy): {prob_failure_scipy}")



# Binomial distribution .......................................................

def binomial_probability(x, n, p):
    """
    Calculate the binomial probability of getting k successes in n trials
    with a probability of success p.
    """
    return comb(n, x) * (p ** x) * ((1 - p) ** (n - x))

print("\nBinomial distribution")

# Probability of 5 heads in 10 flips of a coin
x = 5; n = 10; p = 1/2
probability1 = binomial_probability(x, n, p)
probability2 = stats.binom.pmf(x, n, p)
print(f"\nProbability {x} heads in {n} coin tosses: {probability1}")
print(f"Probability {x} heads in {n} coin tosses: {probability2}")

# Probability of 3 times a 6 in 10 rolls of dice
x = 3; n = 10; p = 1/6
probability1 = binomial_probability(x, n, p)
probability2 = stats.binom.pmf(x, n, p)
print(f"\nProbability {x} 6s in {n} dice rolls: {probability1}")
print(f"Probability {x} 6s in {n} dice rolls: {probability2}")

# Probability of passing an (A, B, C) exam answering randomly
x = 5; n = 10; p = 1/3
probability1 = binomial_probability(x, n, p)
probability2 = stats.binom.pmf(x, n, p)
print(f"\nProbability {x} 5 questions out of {n}: {probability1}")
print(f"Probability {x} 5 questions out of {n}: {probability2}")



# Poisson distribution ........................................................

def poisson_probability(x, lmbda):
    """
    Calculate the probability of observing k events in a Poisson distribution
    with rate parameter lmbda.
    """
    return (lmbda ** x) * exp(-lmbda) / factorial(x)

print("\nPoisson distribution")

# Probability of 3 cancer patients with average 5
x = 3; lmbda = 5
probability1 = poisson_probability(x, lmbda)
probability2 = stats.poisson.pmf(x, lmbda)
print(f"\nProbability {x} cancer patients with average {lmbda}: {probability1}")
print(f"Probability {x} cancer patients with average {lmbda}: {probability2}")

# Probability of 5 or less patients, with same average
x = 5; lmbda = 5
probability1 = stats.poisson.pmf(0, lmbda) + stats.poisson.pmf(1, lmbda) + stats.poisson.pmf(2, lmbda) + stats.poisson.pmf(3, lmbda) + stats.poisson.pmf(4, lmbda) + stats.poisson.pmf(5, lmbda)
probability2 = stats.poisson.cdf(x, lmbda)
print(f"\nProbability of observing {x} or less cancer patients with average {lmbda}: {probability1}")
print(f"Probability of observing {x} or less cancer patients with average {lmbda}: {probability2}")

# Probability of more than 5 patients, with same average
x = 5; lmbda = 5
probability1 = 1 - probability1
probability2 = 1 - probability2
print(f"\nProbability of observing more than {x} cancer patients with average {lmbda}: {probability1}")
print(f"Probability of observing more than {x} cancer patients with average {lmbda}: {probability2}")



# Gaussian distribution .......................................................

def gaussian_density(x, mu, sigma):
    """
    Calculate the probability of a Gaussian distribution falling between x1 and x2.
    """
    return erf((x - mu) / (sqrt(2) * sigma)) / 2

print("\nGaussian distribution")

# Calculating probability
x1 = 1; x2 = 2; mu = 1; sigma = 2
probability1 = gaussian_density(x2, mu, sigma) - gaussian_density(x1, mu, sigma)
probability2 = stats.norm.cdf(x2, mu, sigma) - stats.norm.cdf(x1, mu, sigma)
print(f"\nProbability of finding an event between {x1} and {x2}: {probability1}")
print(f"Probability of finding an event between {x1} and {x2}: {probability2}")

# Calculating probability
x1 = -2; x2 = 2; mu = 1; sigma = 2
probability1 = gaussian_density(x2, mu, sigma) - gaussian_density(x1, mu, sigma)
probability2 = stats.norm.cdf(x2, mu, sigma) - stats.norm.cdf(x1, mu, sigma)
print(f"\nProbability of finding an event between {x1} and {x2}: {probability1}")
print(f"Probability of finding an event between {x1} and {x2}: {probability2}")


