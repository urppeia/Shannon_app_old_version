import pandas as pd
import numpy as np
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy, ks_2samp, f_oneway
from scipy.special import rel_entr
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def calculate_divergences_and_ks(df, smoothing=0.00001):
    """Calculates divergences (including normalized JSD) and the KS statistic for methylation data."""

    # Add smoothing
    df["methylated"] = df["methylated"] + smoothing
    df["unmethylated"] = df["unmethylated"] + smoothing

    # Calculate probabilities
    total = df["methylated"] + df["unmethylated"]
    p = df["methylated"] / total
    q = df["unmethylated"] / total

    # Divergence calculations
    js_divergence = jensenshannon(p, q)
    normalized_js_divergence = js_divergence / np.log(2)
    kl_divergence = entropy(p, q)  # KL from methylated to unmethylated
    gjs_divergence = np.sqrt(
        0.5 * (entropy(p, q, base=2) + entropy(q, p, base=2))
    )  # Use base=2 for bits
    
    # Symmetric geometric js divergence
    m = (p + q) / 2
    sgjs_divergence = 0.5 * (jensenshannon(p, m, base=2) + jensenshannon(q, m, base=2))  

    # Kolmogorov-Smirnov Test
    sorted_methylated = np.sort(df["methylated"])
    sorted_unmethylated = np.sort(df["unmethylated"])

    ecdf_methylated = np.arange(1, len(sorted_methylated) + 1) / len(sorted_methylated)
    ecdf_unmethylated = np.arange(1, len(sorted_unmethylated) + 1) / len(sorted_unmethylated)

    ks_statistic, _ = ks_2samp(sorted_methylated, sorted_unmethylated)

    return normalized_js_divergence, kl_divergence, gjs_divergence, sgjs_divergence, ks_statistic


def read_bismark_file(filename):
    """Reads Bismark methylation data from a file."""
    column_names = ["chr", "start", "end", "coverage", "methylated", "unmethylated"]
    df = pd.read_csv(filename, sep='\t', header=None, names=column_names, compression='gzip')

    # Convert columns to correct data types
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    df['coverage'] = pd.to_numeric(df['coverage'])
    df['methylated'] = pd.to_numeric(df['methylated'])
    df['unmethylated'] = pd.to_numeric(df['unmethylated'])

    return df

def clean_data(results_df):
    """Cleans the results DataFrame."""
    results_df['Sample'] = results_df['Sample'].astype(str).str.replace('.bedgraph.gz', '', regex=False)
    return results_df




# Data processing pipeline
data_directory = "/shares/grossniklaus.botinst.uzh/eharputluoglu/test_run_datasets/"
results = []

# Iterate through files
for file_path in glob.glob(os.path.join(data_directory, "**/*.bismark.cov.gz"), recursive=True):  
    filename = os.path.basename(file_path)
    parts = filename.split("_")
    context_type = parts[0]
    sample_name = "_".join(parts[1:3])  # Combine SRX and se parts
    chromosome = parts[-1].split(".")[0][3:]  # Extract chromosome (e.g., '1' from 'chr1')

    df = read_bismark_file(file_path)
    normalized_js, kl, gjs, sgjs, ks = calculate_divergences_and_ks(df)


# Store results
    results.append({
        "Sample": sample_name,
        "Chromosome": chromosome,
        "Context": context_type,
        "Normalized JS Divergence": normalized_js,
        "KL Divergence": kl,
        "GJS Divergence": gjs,
        "SGJS Divergence": sgjs,
        "KS Statistic": ks,
    })

# Create a DataFrame from results
results_df = pd.DataFrame(results)

# Clean data
results_df = clean_data(results_df)



