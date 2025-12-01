# Introduction to Computational Statistics
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 1: Descriptive statistics



# Exercise 2: 

# 1. Read expression data from EM-seq experiment
# Compute methylation score s = nC / (nC + nT)
# Average score per replicate

# 2. Implement a manual computation of mean, median, variance and std
# Check calculation with numpy / scipy implementation

# 3. Plot data as a box plot, violin and histogram
# Plot statistics summary on top of the histogram



# Import libraries ............................................................

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



# Load EM-seq .................................................................

# Read EM-seq data
emseq = pd.read_csv("data/emseq_counts.cleaned.csv")
print("\nEM-seq data\n", emseq.head())
print("\nEM-seq dimensions: ", emseq.shape)

# Compute methylation score
def score(numCs, numTs):
    return numCs / (numCs + numTs + 1e-9)

# Average scores for 3uM deprived
s_3uM_rep1 = score(emseq["3uM_rep1_numCs"], emseq["3uM_rep1_numTs"])
s_3uM_rep2 = score(emseq["3uM_rep2_numCs"], emseq["3uM_rep2_numTs"])
score_3uM = (s_3uM_rep1 + s_3uM_rep2) / 2

# Average scores for 200uM control
s_200uM_rep1 = score(emseq["200uM_rep1_numCs"], emseq["200uM_rep1_numTs"])
s_200uM_rep2 = score(emseq["200uM_rep2_numCs"], emseq["200uM_rep2_numTs"])
score_200uM = (s_200uM_rep1 + s_200uM_rep2) / 2



# Compute mean value ..........................................................

# Manual calculation
def mean_manual(x):
    return sum(x) / len(x)

# Check implementation
print("\nMean value")
print("3uM (manual):", mean_manual(score_3uM.values))
print("3uM (Numpy):", np.mean(score_3uM.values))
print("200uM (manual):", mean_manual(score_200uM.values))
print("200uM (Numpy):", np.mean(score_200uM.values))


# Compute variance ............................................................

# Manual calculation
def var_manual(x):
    m = mean_manual(x)
    return sum((xi - m)**2 for xi in x) / len(x)

# Check implementation
print("\nVariance")
print("3uM (manual):", var_manual(score_3uM.values))
print("3uM (Numpy):", np.var(score_3uM.values))
print("200uM (manual):", var_manual(score_200uM.values))
print("200uM (Numpy):", np.var(score_200uM.values))



# Compute median ..............................................................

# Manual calculation
def median_manual(x):
    s = sorted(x)
    n = len(s)
    mid = n // 2
    if n % 2 == 0:
        return (s[mid-1] + s[mid]) / 2
    return s[mid]

# Check implementation
print("\nMedian")
print("3uM (manual):", median_manual(score_3uM.values))
print("3uM (Numpy):", np.median(score_3uM.values))
print("200uM (manual):", median_manual(score_200uM.values))
print("200uM (Numpy):", np.median(score_200uM.values))



# Compute std .................................................................

# Manual calculation
def std_manual(x):
    return var_manual(x)**0.5

# Check implementation
print("\nStandard deviation")
print("3uM (manual):", std_manual(score_3uM.values))
print("3uM (Numpy):", np.std(score_3uM.values))
print("200uM (manual):", std_manual(score_200uM.values))
print("200uM (Numpy):", np.std(score_200uM.values))



# Plot data and summary statistics ............................................

# Prepare for plot
samples = [score_3uM.values, score_200uM.values]
labels = ["3uM Deprived", "200uM Control"]

# Box plot
plt.figure(figsize = (8, 5))
sns.boxplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Score")
plt.grid(True, alpha = 0.3)
# plt.savefig("emseq_box.png", dpi = 300, bbox_inches = "tight")
# plt.show()

# Violin plot
plt.figure(figsize = (8, 5))
sns.violinplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Score")
plt.grid(True, alpha = 0.3)
# plt.savefig("emseq_violin.png", dpi = 300, bbox_inches = "tight")
# plt.show()

# Histogram

# Extract mean per group
mean_3uM = mean_manual(samples[0])
mean_200uM = mean_manual(samples[1])

# Prepare binsize
all_data = np.concatenate(samples)
bins = np.linspace(all_data.min(), all_data.max(), 30)

# Plot histogram
plt.figure(figsize = (8, 5))
plt.hist(samples[0], bins = bins, label = labels[0], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.hist(samples[1], bins = bins, label = labels[1], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.axvline(mean_3uM, color = "blue", linestyle = "--", label = f"3uM mean = {mean_3uM:.3f}")
plt.axvline(mean_200uM, color = "darkorange", linestyle = "--", label = f"200uM mean = {mean_200uM:.3f}")
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, alpha = 0.3)
# plt.savefig("emseq_histogram.png", dpi = 300, bbox_inches = "tight")
# plt.show()
