
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









###################################################### ARGS  ++++++++++++++########################



import argparse
import io
import os
import sys

import pandas as pd

from shannonlib.core import divergence


def run_divergence(args):

    if os.path.isfile(args.output) and not os.stat(args.output).st_size == 0:
        msg = "-- Stopped!\n-- Output file exists and is not empty."
        sys.exit(msg)

    meta = args.metadata.read()

    try:
        sample = pd.read_csv(io.StringIO(meta), comment='#', header=0)
        _ = sample['url']
        _ = sample['label']
    except KeyError:
        try:
            sample = pd.read_csv(io.StringIO(
                meta), comment='#', header=0, sep='\t')
            _ = sample['url']
            _ = sample['label']
        except KeyError:
            msg = ('-- Stopped!\n'
                   '-- Could not parse metadata.\n'
                   '-- 1. Ensure that fields are tab- or comma-separated\n'
                   '-- 2. Ensure that columns "url" and "label" are present')
            sys.exit(msg)

    # GPF data columns
    try:
        assert(len(args.dcols) == len(args.dnames))
    except AssertionError:
        msg = ('-- Stopped!\n'
               '-- Length of --dcols and --dnames must match')
        sys.exit(msg)

    dtypes = [float if args.prob else int] * len(args.dcols)
    dcols = [col - 1 for col in args.dcols]
    gpf_data = [list(zip(dcols, args.dnames, dtypes))]

    # if groupby is not None:
    #     subsample_family = sample.groupby(groupby)
    #     for key, subsample in subsample_family:
    #         filename = groupname(by=groupby, name=key, fname=args.output)
    #         divergence(subsample, chrom=args.sequence, data_columns=gpf_data,
    #                    outfile=filename, chunksize=args.chunksize)
    # else:
    print('processing sequence {} ...'.format(args.sequence))
    divergence(sample, chrom=args.sequence, data_columns=gpf_data,
               outfile=args.output, chunksize=args.chunk)

    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # main parser
    parser.prog = 'shannon'
    parser.description = 'Command-line interface to %(prog)s.'

    # cluster
    parser_cluster = subparsers.add_parser('cluster')

    # segment
    parser_segment = subparsers.add_parser('seg')

    # divergence
    parser_div = subparsers.add_parser(
        'div', formatter_class=argparse.RawTextHelpFormatter)

    parser_div.set_defaults(func=run_divergence)
    parser_div.help = 'JS Divergence for genome position files (GPFs).'
    parser_div.description = parser_div.help
    parser_div_required = parser_div.add_argument_group('required arguments')

    parser_div.add_argument(
        '--prob', action='store_true',
        help='indicate that data are probabilites (default: counts)')

    parser_div.add_argument(
        '--chunk', metavar='SIZE', default=1e4, type=int,
        help=('set size of data to process in-memory (default: %(default)d)\n'
              '- in terms of expected number of genome positions\n'
              '- higher numbers lead to more memory-hungry, faster computations'))

    parser_div_required.add_argument(
        '-m', '--metadata', metavar='FILE', type=argparse.FileType('r'),
        required=True, help=('metadata for GPFs\n'
                             '- comment lines (#) are ignored\n'
                             '- values must be comma- or tab-separated\n'
                             '- first non-comment line must be header\n'
                             '- "url" and "label" columns are required\n'
                             '- if stdin is metadata use "--metadata -"'))

    parser_div_required.add_argument(
        '-o', '--output', metavar='FILE', required=True, help='output filepath')

    parser_div_required.add_argument(
        '-s', '--sequence', metavar='ID', required=True, type=str,
        help='query sequence (chromosome/scaffold) in GPF')

    parser_div_required.add_argument(
        '-c', '--dcols', metavar='COLN', nargs='+', required=True, type=int,
        help='column numbers (1-based) in GPFs that hold the data')

    parser_div_required.add_argument(
        '-n', '--dnames', metavar='NAME', nargs='+', required=True, type=str,
        help='names of data columns following the order in --dcols')

    # parser.add_argument('-g', '--groupby', metavar='STR', nargs='+', type=str,
    #                     help='''
    #                     The factor according to which the selected set is
    #                     partitioned; has to match column names in the input
    #                     metadata. One or more factors can be given, which
    #                     produces a file for each combination of factors.
    #                     ''')

    args = parser.parse_args()
    args.func(args)























































