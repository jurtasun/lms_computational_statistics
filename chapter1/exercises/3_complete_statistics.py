# Introduction to Computational Statistics
# Jesús Urtasun Elizari: MRC LMS 2026
# Chapter 1: Descriptive statistics



# Exercise 3: 
# Read expression data from RNA-seq experiment
# Read expression data from EM-seq experiment
# Average counts per replicate
# Plot data as a box plot, violin and histogram
# Implement a manual computation of mean, median, variance and std
# Check calculation with numpy / scipy implementation



# Import libraries ............................................................

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



# Load data ...................................................................

emseq = pd.read_csv("data/emseq_counts.cleaned.csv")
print("\nEM-seq data\n", emseq.head())
print("\nEM-seq dimensions: ", emseq.shape)

rnaseq = pd.read_csv("data/rnaseq_counts.cleaned.csv")
print("\nRNA-seq data\n", rnaseq.head())
print("\nRNA-seq dimensions: ", rnaseq.shape)



# Process data ................................................................

# EM-seq: compute methylation score
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


#  RNA-seq expression: extract deprived vs control

rnaseq_3uM = rnaseq[[c for c in rnaseq.columns if "3uM" in c]]
rnaseq_200uM = rnaseq[[c for c in rnaseq.columns if "200uM" in c]]

exp_3uM = rnaseq_3uM.values.flatten()
exp_200uM = rnaseq_200uM.values.flatten()


# Plot EM-seq data ............................................................

# Prepare for plot
samples = [score_3uM.values, score_200uM.values]
labels = ["3uM Deprived", "200uM Control"]

# Box plot
plt.figure(figsize = (8, 5))
sns.boxplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Score")
plt.grid(True, alpha = 0.3)
plt.savefig("emseq_box.png", dpi=300, bbox_inches="tight")
plt.show()

# Violin plot
plt.figure(figsize = (8, 5))
sns.violinplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Score")
plt.grid(True, alpha = 0.3)
plt.savefig("emseq_violin.png", dpi=300, bbox_inches="tight")
plt.show()

# Histogram
all_data = np.concatenate(samples)
bins = np.linspace(all_data.min(), all_data.max(), 30)
plt.figure(figsize = (8, 5))
plt.hist(samples[0], bins = bins, label = labels[0], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.hist(samples[1], bins = bins, label = labels[1], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, alpha = 0.3)
plt.savefig("emseq_histogram.png", dpi=300, bbox_inches="tight")
plt.show()



# Plot RNA-seq data (raw values) ........................................................

# Prepare for plot
samples = [exp_3uM, exp_200uM]
labels = ["3uM Deprived", "200uM Control"]

# Log transform for proper visualization
log_3uM = np.log1p(exp_3uM)
log_200uM = np.log1p(exp_200uM)
samples = [log_3uM, log_200uM]
labels = ["3uM Deprived (log)", "200uM Control (log)"]

# Box plot
plt.figure(figsize = (8, 5))
sns.boxplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Expression")
plt.grid(True, alpha = 0.3)
plt.savefig("rnaseq_box.png", dpi=300, bbox_inches="tight")
plt.show()

# Violin plot
plt.figure(figsize = (8, 5))
sns.violinplot(data = samples)
plt.xticks(range(2), labels)
plt.ylabel("Expression")
plt.grid(True, alpha = 0.3)
plt.savefig("rnaseq_violin.png", dpi=300, bbox_inches="tight")
plt.show()

# Histogram
all_data = np.concatenate(samples)
bins = np.linspace(all_data.min(), all_data.max(), 30)
plt.figure(figsize = (8, 5))
plt.hist(samples[0], bins = bins, label = labels[0], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.hist(samples[1], bins = bins, label = labels[1], alpha = 0.5, edgecolor = "black", linewidth = 0.8)
plt.xlabel("Expression")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, alpha = 0.3)
plt.savefig("rnaseq_histogram.png", dpi=300, bbox_inches="tight")
plt.show()



# Compute summary statistics ..................................................

# Manual calculation
def manual_mean(x):
    return sum(x) / len(x)

# Manual calculation
def manual_variance(x):
    m = manual_mean(x)
    return sum((xi - m)**2 for xi in x) / len(x)

# Manual calculation
def manual_median(x):
    s = sorted(x)
    n = len(s)
    mid = n // 2
    if n % 2 == 0:
        return (s[mid-1] + s[mid]) / 2
    return s[mid]

# Manual calculation
def manual_std(x):
    return manual_variance(x)**0.5



# Check implementation ........................................................

def print_stats(x, label):

    # Get current comparison
    print(f"\n{label} summmary stats\n")

    # Check mean
    print("Manual mean:", manual_mean(x))
    print("Numpy mean:", np.mean(x))

    # Check variance
    print("Manual variance:", manual_variance(x))
    print("Numpy variance:", np.var(x))

    # Check median
    print("Manual median:", manual_median(x))
    print("Numpy median:", np.median(x))

    # Check std
    print("Manual std:", manual_std(x))
    print("Numpy std:", np.std(x))

print_stats(score_3uM.values, "EM-seq score 3uM (deprived)")
print_stats(score_200uM.values, "EM-seq score 200uM (control)")
print_stats(exp_3uM, "RNA-seq expression 3uM (deprived)")
print_stats(exp_200uM, "RNA-seq expression 200uM (control)")
