# Import libraries
import numpy as np
import pandas as pd

# Random seed
np.random.seed(42)



# 1. Binomial: Allelic counts at SNPs .............................................................

n_snps = 10_000  # Increased from 10
total_reads = np.random.randint(20, 200, size=n_snps)     # coverage per SNP
allele_freq = np.random.uniform(0.05, 0.95, size=n_snps)  # allele frequency
alt_reads = np.random.binomial(total_reads, allele_freq)

binomial_df = pd.DataFrame({
    "chrom": np.random.choice(["chr1", "chr2", "chr3"], size=n_snps),
    "position": np.arange(1_000_000, 1_000_000 + n_snps * 10, 10),
    "ref_reads": total_reads - alt_reads,
    "alt_reads": alt_reads
})
binomial_df.to_csv("allelic_counts.csv", index=False)
print("Saved: allelic_counts.csv (Binomial, n=10,000 SNPs)")



# 2. Poisson: Read counts per region (e.g., ChIP-seq peaks) .......................................

n_regions = 5_000  # Increased from 10
lambda_values = np.random.uniform(5, 30, size=n_regions)
region_counts = np.random.poisson(lam=lambda_values)

poisson_df = pd.DataFrame({
    "chrom": np.random.choice(["chr1", "chr2", "chr3"], size=n_regions),
    "start": np.arange(2_000_000, 2_000_000 + n_regions * 200, 200),
    "end": np.arange(2_000_100, 2_000_100 + n_regions * 200, 200),
    "read_count": region_counts
})
poisson_df.to_csv("chipseq_counts.bed", sep="\t", header=False, index=False)
print("Saved: chipseq_counts.bed (Poisson, n=5,000 regions)")



# 3. Gaussian: Normalized log expression levels ...................................................

n_genes = 5_000   # Increased from 10
n_samples = 10
samples = [f"Sample_{i}" for i in range(1, n_samples + 1)]
expression = np.random.normal(loc=0, scale=1, size=(n_genes, n_samples))

gaussian_df = pd.DataFrame(expression, columns=samples)
gaussian_df.insert(0, "gene_id", [f"Gene_{i}" for i in range(1, n_genes + 1)])
gaussian_df.to_csv("log_expression.csv", index=False)
print("Saved: log_expression.csv (Gaussian, 5,000 genes × 10 samples)")



# 4. Negative Binomial: Raw RNA-seq counts ........................................................

# Parameters: mean (mu), dispersion (r)
n_genes = 5_000
n_samples = 10
mu = np.random.uniform(10, 500, size=n_genes)
r = 10  # dispersion (higher = less overdispersion)
p = r / (r + mu)
nb_counts = np.random.negative_binomial(r, p[:, None], size=(n_genes, n_samples))

neg_binom_df = pd.DataFrame(nb_counts, columns=samples)
neg_binom_df.insert(0, "gene_id", [f"Gene_{i}" for i in range(1, n_genes + 1)])
neg_binom_df.to_csv("rna_counts.csv", index=False)
print("Saved: rna_counts.csv (Negative Binomial, 5,000 genes × 10 samples)")

print("\n\u2705 All large simulated genomic datasets saved successfully.")
