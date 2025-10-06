# Introduction to probability, statistics and hypothesis testing
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 1: Descriptive statistics



# Exercise 2: 
# Simulate random samples following a gaussian distribution
# Compute for each the mean and variance
# Represent data as histogram, box and violin plots



# Import libraries ............................................................

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



# Simulate data ...............................................................

# Random seed
np.random.seed(42)

# Generate gaussian distributed samples
sample1 = np.random.normal(loc = 5, scale = 1, size = 500)
sample2 = np.random.normal(loc = 2, scale = 2, size = 500)
sample3 = np.random.normal(loc = 4, scale = 0.5, size = 500)

# Collect samples in a named list
samples = [sample1, sample2, sample3]
labels = ['Sample 1', 'Sample 2', 'Sample 3']

# Compute mean and variance
for i, s in enumerate(samples):

    mean = np.mean(s)
    var = np.var(s)
    print(f"{labels[i]} Mean: {mean:.3f}, Variance: {var:.3f}")



# Box plot ....................................................................

plt.figure(figsize = (8, 5))
sns.boxplot(data = samples)
plt.title("Box plot of gaussian samples")
plt.xticks(ticks = range(3), labels = labels)
plt.ylabel("Value")
plt.grid(True, alpha = 0.3)
plt.show()



# Violin plot ....................................................................

plt.figure(figsize = (8, 5))
sns.violinplot(data = samples)
plt.title("Violin plot of gaussian samples")
plt.xticks(ticks = range(3), labels = labels)
plt.ylabel("Value")
plt.grid(True, alpha = 0.3)
plt.show()


# Histogram ....................................................................

# Compute common bin size
all_data = np.concatenate(samples)
bins = np.linspace(all_data.min(), all_data.max(), 30)

# Plot histogram
plt.figure(figsize = (8, 5))

for s, label in zip(samples, labels):
    plt.hist(s, bins = bins, edgecolor = 'black', linewidth = 0.8, label = label, alpha = 0.5)

plt.title("Histogram of gaussian samples")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, alpha = 0.3)
plt.show()


