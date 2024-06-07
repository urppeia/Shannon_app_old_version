import pandas as pd
import numpy as np
import glob
import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy, ks_2samp, f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


######################### function definitions ######################

def calculate_divergences_and_ks(df, smoothing=1e-10):
    """Calculates divergences (JSD, KL, GJS, SGJS) and KS test for methylation data.

    Args:
        df (pandas.DataFrame): DataFrame with 'methylated' and 'unmethylated' columns.
        smoothing (float, optional): Smoothing factor to add to counts. Defaults to 1e-10.

    Returns:
        tuple: Normalized divergences (JSD, KL, GJS, SGJS) and KS test statistic and p-value.
    """
    
    # Add smoothing 
    df["methylated"] += smoothing
    df["unmethylated"] += smoothing
    
    # Calculate probabilities
    total = df["methylated"] + df["unmethylated"]
    p = df["methylated"] / total
    q = df["unmethylated"] / total
    
    # Divergence calculations (with normalizations)
    normalized_js_divergence = float(jensenshannon(p, q, base=2))  # JSD

    normalized_kl_divergence = float(entropy(p, q, base=2) - entropy(p, base=2))  # KLD

    # Geometric Jensen-Shannon Divergence
    normalized_gjs_divergence = float(np.sqrt(
        0.5 * (np.sqrt(entropy(p, q, base=2)) + np.sqrt(entropy(q, p, base=2)))
    ) / np.sqrt(np.log(2)))
    
    # Symmetric Geometric Jensen-Shannon Divergence
    m = (p + q) / 2
    normalized_sgjs_divergence = float( 0.5 * (
        np.sqrt(jensenshannon(p, m, base=2)) + np.sqrt(jensenshannon(q, m, base=2))
    ) / np.sqrt(np.log(2)))

    # Kolmogorov-Smirnov Test
    ks_statistic, ks_pvalue = ks_2samp(df["methylated"], df["unmethylated"])  # KST

    return normalized_js_divergence, normalized_kl_divergence, normalized_gjs_divergence, normalized_sgjs_divergence, ks_statistic, ks_pvalue

def read_bismark_file(filename):
    column_names = ["chr", "start", "end", "coverage", "methylated", "unmethylated"]
    df = pd.read_csv(filename, sep='\t', header=None, names=column_names, compression='gzip')

    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    df['coverage'] = pd.to_numeric(df['coverage'])
    df['methylated'] = pd.to_numeric(df['methylated'])
    df['unmethylated'] = pd.to_numeric(df['unmethylated'])

    return df

def clean_data(results_df):
    """Cleans the results DataFrame by removing file extensions from sample names."""
    results_df['Sample'] = results_df['Sample'].astype(str).str.replace('.bedgraph.gz', '', regex=False)
    return results_df



################# RUN CODE HERE #################


data_directory = "/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets_gz"
results = []

# Iterate through files
for file_path in glob.glob(os.path.join(data_directory, "**/*.bismark.cov.gz"), recursive=True):
    try:
        filename = os.path.basename(file_path)
        parts = filename.split("_")
        context_type = parts[0]
        sample_name = "_".join(parts[1:3])
        chromosome = parts[-1].split(".")[0][3:]

        df = read_bismark_file(file_path)
        divergences = calculate_divergences_and_ks(df) 

        results.append(
            {
                "Sample": sample_name,
                "Chromosome": chromosome,
                "Context": context_type,
                "JS Divergence": divergences[0],  
                "KL Divergence": divergences[1],
                "GJS Divergence": divergences[2],
                "SGJS Divergence": divergences[3],
                "KS Statistic": divergences[4],
                "KS P-value": divergences[5]
            }
        )

    except Exception as e:
        logging.error(f"Error processing file {filename}: {e}")

results_df = pd.DataFrame(results)

results_df = clean_data(results_df)

# save to csv
results_df.to_csv("/shares/grossniklaus.botinst.uzh/eharputluoglu/calculation_output/methylation_analysis_results.csv", index=False)
