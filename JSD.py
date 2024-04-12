
import os
import logging
import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf
from tqdm import tqdm

# Setup basic logging
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)

################## Core functions.#############################################

def divergence(sample, chrom=None, data_columns=None, outfile=None, chunksize=None):
    """Compute within-group divergence using Jensen-Shannon divergence.

    Args:
        sample (dict): Contains 'url' for data file path and 'label' for sample labels.
        chrom (str, optional): Chromosome number to limit analysis. Processes all if None.
        data_columns (list[str], optional): Data columns to use for analysis.
        outfile (str, optional): Output file path. If None, doesn't write to file.
        chunksize (int, optional): Sites per chunk for processing. Uses default size if None.
    """
    # Get regions and their data
    regions_pct, regions = gpf.get_regions(sample['url'], chrom, chunksize)
    regions_data = gpf.get_data(sample['url'], sample['label'], data_columns, regions)

    progress_bar = tqdm(total=len(regions_pct))

    for progress, data in zip(regions_pct, regions_data):
        if data.empty:
            progress_bar.set_description('Skipped empty region')
        else:
            div = est.js_divergence(data)
            if not div.empty:
                # Check and write header if needed
                header = not os.path.isfile(outfile) or os.stat(outfile).st_size == 0
                div.round({'JSD_bit_': 3, 'HMIX_bit_': 3}).to_csv(
                    outfile, header=header, sep='\t', mode='a', index=True)
            else:
                progress_bar.set_description('Skipped low-quality region')
        progress_bar.update()

    progress_bar.close()
    logging.info("The File is created!")


################## Estimators.###################################################


def shannon_entropy(countmatrix, axis=1):
    """Calculate Shannon entropy of feature frequency profile in a count matrix."""
    logger.debug("Calculating Shannon entropy...")
    count_sum = countmatrix.sum(axis=axis, keepdims=True)
    prob = countmatrix / count_sum
    entropy = np.sum(np.where(prob > 0, -prob * np.log(prob), 0), axis=axis)
    logger.debug("Shannon entropy calculation complete.")
    return entropy

def js_divergence(indata, weights=None):
    """Calculate Jensen-Shannon divergence among groups in a DataFrame."""
    LOG2E = np.log2(math.e)
    logger.debug("Calculating JS divergence...")

    count_per_unit = indata.groupby(level=0).sum()
    samplesize = count_per_unit.notnull().sum(axis=1)
    data_filtered = indata[(count_per_unit >= 3).any(axis=1) & (samplesize >= 2)]

    if data_filtered.empty:
        return pd.DataFrame()  # Returning empty DataFrame if filters remove all data

    data_unit_sum = count_per_unit.loc[data_filtered.index]
    data_feature_sum = data_filtered.groupby(level=1).sum().astype(np.int32)

    counts = data_filtered.values.reshape(data_unit_sum.shape[0], data_unit_sum.shape[1], -1)
    mix_entropy = shannon_entropy(data_feature_sum.values)
    avg_entropy = np.average(shannon_entropy(counts, axis=2), weights=data_unit_sum.fillna(0), axis=1)

    div = pd.DataFrame({
        'JSD_bit_': LOG2E * (mix_entropy - avg_entropy),
        'sample size': samplesize.loc[data_filtered.index],
        'HMIX_bit_': LOG2E * mix_entropy
    }, index=data_filtered.index)

    logger.debug("JS divergence calculation complete.")
    return div
