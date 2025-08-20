#!/usr/bin/env python3

import pandas as pd
import glob
import os

# Collect all quant files
files = glob.glob("*.Oarfish.byReads.quant")

# Dictionary to hold each sample's counts
dfs = {}

for file in files:
    # Extract sample name (remove extension)
    sample_name = os.path.basename(file).replace(".Oarfish.byReads.quant", "")

    # Read file
    df = pd.read_csv(file, sep="\t")

    # Store num_reads with tname as index
    dfs[sample_name] = df.set_index("tname")["num_reads"]

# Combine into a single dataframe (counts matrix)
counts_matrix = pd.DataFrame(dfs)

# Compute CPM (normalize each column to sum to 1e6)
cpm_matrix = counts_matrix.div(counts_matrix.sum(axis=0), axis=1) * 1e6

# Save outputs
counts_matrix.to_csv("counts_matrix.csv")
cpm_matrix.to_csv("cpm_matrix.csv")

print(
    "Counts matrix and CPM matrix have been written to 'counts_matrix.csv' and 'cpm_matrix.csv'"
)
