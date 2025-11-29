#!/usr/bin/env python3
import pandas as pd

# -------------------------------------------------------------------
#                  LOAD INPUT FILES
# -------------------------------------------------------------------

emseq_file = "emseq_counts.csv"
rnaseq_file = "rnaseq_counts.csv"

# ============================================================
#                        EM-seq
# ============================================================

print("\nLoading EM-seq...")
emseq = pd.read_csv(emseq_file)

# ---- 1. Create Tile_ID = chr_start_end ----
emseq["Tile_ID"] = (
    emseq["chr"].astype(str)
    + "_"
    + emseq["start"].astype(str)
    + "_"
    + emseq["end"].astype(str)
)

# Move Tile_ID to first column
emseq = emseq[["Tile_ID"] + [c for c in emseq.columns if c != "Tile_ID"]]

# ---- 2. Remove coverage columns ----
coverage_cols = [c for c in emseq.columns if "coverage" in c]
emseq = emseq.drop(columns=coverage_cols)

# ---- 3. Rename columns with _rep1 and _rep2 consistently ----
rename_map = {
    # 3uM treatment
    "Met_1_numCs": "3uM_rep1_numCs",
    "Met_1_numTs": "3uM_rep1_numTs",
    "Met_2_numCs": "3uM_rep2_numCs",
    "Met_2_numTs": "3uM_rep2_numTs",

    # 6uM treatment
    "Met_3_numCs": "6uM_rep1_numCs",
    "Met_3_numTs": "6uM_rep1_numTs",
    "Met_4_numCs": "6uM_rep2_numCs",
    "Met_4_numTs": "6uM_rep2_numTs",

    # 200uM treatment
    "Met_5_numCs": "200uM_rep1_numCs",
    "Met_5_numTs": "200uM_rep1_numTs",
    "Met_6_numCs": "200uM_rep2_numCs",
    "Met_6_numTs": "200uM_rep2_numTs",
}

emseq = emseq.rename(columns=rename_map)

# ---- 4. Remove duplicated rows (exact duplicates) ----
numeric_cols = [c for c in emseq.columns if c not in ["Tile_ID", "chr", "start", "end", "strand"]]

before = len(emseq)
emseq = emseq.drop_duplicates(subset=numeric_cols, keep="first")
after = len(emseq)

print(f"Removed {before - after} duplicated EM-seq windows")

# ---- Save cleaned EM-seq ----
emseq.to_csv("emseq_counts.cleaned.csv", index=False)
print("✔ EM-seq cleaned → emseq_counts.cleaned.csv")


# ============================================================
#                         RNA-seq
# ============================================================

print("\nLoading RNA-seq...")
rnaseq = pd.read_csv(rnaseq_file)

# ---- 1. Rename first column to Ensemble_ID ----
rnaseq = rnaseq.rename(columns={rnaseq.columns[0]: "Ensemble_ID"})

# ---- 2. Rename sample intensities ----
new_names = {}

for c in rnaseq.columns:
    if "_3_" in c:
        new_names[c] = c.replace("_3_", "_3uM_")
    elif "_6_" in c:
        new_names[c] = c.replace("_6_", "_6uM_")
    elif "_200_" in c:
        new_names[c] = c.replace("_200_", "_200uM_")

rnaseq = rnaseq.rename(columns=new_names)

# ---- Save cleaned RNA-seq ----
rnaseq.to_csv("rnaseq_counts.cleaned.csv", index=False)
print("✔ RNA-seq cleaned → rnaseq_counts.cleaned.csv\n")

print("All preprocessing completed.\n")
