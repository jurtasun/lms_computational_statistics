# MRC LMS Introduction to probability, statistics and hypothesis testing
# Chapter 3: Parameter estimation

# Exercise 1



# Import libraries ............................................................

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binom, poisson, norm, nbinom

# Set default plotting style
plt.style.use('seaborn-v0_8-whitegrid')

# Function plotting nice histogram
def plot_hist(data, bins = 30, label = "Observed", color = "skyblue"):
    plt.hist(data, bins = bins, density = True, alpha = 0.6, color = color,
             edgecolor = 'black', linewidth = 0.8, label = label)



# 1. BINOMIAL: Allelic counts .................................................

# Read data
binom_df = pd.read_csv("00_data/allelic_counts.csv")
alt = binom_df["alt_reads"]
tot = binom_df["ref_reads"] + binom_df["alt_reads"]

# Compute meand and variance
mean_obs = alt.mean()
var_obs = alt.var()

# Use overall average n, p for approximate fit
n = int(np.round(tot.mean()))
p = np.clip(mean_obs / n, 0, 1)

x = np.arange(0, n + 1)
plt.figure(figsize=(6,4))
plot_hist(alt, bins=50)
plt.plot(x, binom.pmf(x, n, p), 'r-', lw=2,
         label=f'Approx. Binomial(n={n}, p={p:.2f})')
plt.title(f"Binomial-like SNP Allelic Counts (mixture)\nMean={mean_obs:.2f}, Var={var_obs:.2f}")
plt.xlabel("Alt allele count")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()



# 2. POISSON: ChIP-seq counts .................................................

# Read data
poiss_df = pd.read_csv("00_data/chipseq_counts.bed", 
                        sep = "\t", header = None, names = ["chrom", "start", "end", "read_count"])
counts = poiss_df["read_count"]

# Compute meand and variance
mean_obs = counts.mean()
var_obs = counts.var()

x = np.arange(0, counts.max() + 1)
plt.figure(figsize=(6,4))
plot_hist(counts, bins=50)
plt.plot(x, poisson.pmf(x, mean_obs), 'r-', lw=2,
         label=f'Approx. Poisson(λ={mean_obs:.2f})')
plt.title(f"Poisson-like Read Counts (mixture)\nMean={mean_obs:.2f}, Var={var_obs:.2f}")
plt.xlabel("Read count")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()



# 3. GAUSSIAN: Log expression ................................................

# Read data
gauss_df = pd.read_csv("00_data/log_expression.csv")
expr = gauss_df.drop(columns="gene_id").values.flatten()

# Compute meand and variance
mean_obs = expr.mean()
var_obs = expr.var()

x = np.linspace(expr.min(), expr.max(), 200)
plt.figure(figsize=(6,4))
plot_hist(expr, bins=50)
plt.plot(x, norm.pdf(x, mean_obs, np.sqrt(var_obs)), 'r-', lw=2,
         label=f'N({mean_obs:.2f}, {var_obs:.2f})')
plt.title(f"Gaussian – Log Expression\nMean={mean_obs:.2f}, Var={var_obs:.2f}")
plt.xlabel("Expression (log)")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()



# 4. NEGATIVE BINOMIAL: RNA-seq counts ........................................

# Read data
nb_df = pd.read_csv("00_data/rna_counts.csv")
counts = nb_df.drop(columns="gene_id").values.flatten()

# Compute mean and variance
mean_obs = counts.mean()
var_obs = counts.var()

# Method-of-moments parameter estimates
r_est = mean_obs**2 / (var_obs - mean_obs) if var_obs > mean_obs else 10
p_est = r_est / (r_est + mean_obs)

x = np.arange(0, np.percentile(counts, 99).astype(int))  # cut extreme tail
plt.figure(figsize=(6,4))
plot_hist(counts, bins=100)
plt.plot(x, nbinom.pmf(x, r_est, p_est), 'r-', lw=2,
         label=f'NegBin(r={r_est:.1f}, p={p_est:.2f})')
plt.title(f"Negative Binomial – RNA-seq Counts\nMean={mean_obs:.2f}, Var={var_obs:.2f}")
plt.xlabel("Read count")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()

print("\n\u2705 All distributions processed and plotted (with black edges).")