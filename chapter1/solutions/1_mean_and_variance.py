# Introduction to probability, statistics and hypothesis testing
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 1: Descriptive statistics



# Exercise 1: 
# Simulate random data following a gaussian distribution
# Implement a manual computation of mean, median, variance and std
# Check calculation with numpy / scipy implementation



# Import libraries ............................................................

import numpy as np
import matplotlib.pyplot as plt



# Simulate data ...............................................................

# Random seed
np.random.seed(123)

# Simulate data
data = np.random.normal(loc = 5, scale = 2, size = 10)
print(f"\nRandom data:\n{data}")
print(f"Data format: {type(data)}")
print(f"Data shape: {data.shape}")



# Compute mean value ..........................................................

# Manual calculation
def manual_mean(data):
    
    return sum(data) / len(data)

# Check implementation
mean_manual = manual_mean(data)
mean_np = np.mean(data)
print("\nMean (manual):", mean_manual)
print("Mean (numpy):", mean_np)



# Compute variance ............................................................

# Manual calculation
def manual_variance(data):

    mean_val = manual_mean(data)
    squared_diffs = [(x - mean_val) ** 2 for x in data]

    return sum(squared_diffs) / len(data) # population variance
    # return sum(squared_diffs) / (len(data) - 1)  # sample variance

# Check implementation
var_manual = manual_variance(data)
var_np = np.var(data, ddof = 0) # population variance
# var_np = np.var(data, ddof = 1) # sample variance
print("\nVariance (manual):", var_manual)
print("Variance (numpy):", var_np)



# Compute median ..............................................................

# Manual calculation
def manual_median(data):

    # Sort data and find middle point
    sorted_data = sorted(data)
    n = len(sorted_data)
    middle = n // 2
    
    # Check if even / odd data
    if n % 2 == 0:
        return (sorted_data[middle - 1] + sorted_data[middle]) / 2
    else:
        return sorted_data[middle]

# Check implementation
median_manual = manual_median(data)
median_np = np.median(data)
print("\nMedian (manual):", median_manual)
print("Median (numpy):", median_np)



# Compute std .................................................................

# Manual calculation
def manual_std(data):
    var = manual_variance(data)
    return var ** 0.5

# Check implementation
std_manual = manual_std(data)
std_np = np.std(data, ddof=0)
print("\nStandard deviation (manual):", std_manual)
print("Standard deviation (numpy):", std_np)


