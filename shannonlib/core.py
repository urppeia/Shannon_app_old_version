# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf

from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')



def divergence(sample, chrom=None, data_columns=None, outfile=None, chunksize=None):
    """Computes within-group divergence for population.
    """

    regions_pct, regions = gpf.get_regions(
        sample['url'], chrom=chrom, exp_numsites=chunksize)

    regions_data = gpf.get_data(sample['url'], labels=sample['label'],
                                data_columns=data_columns, regions=regions)

    progress_bar = tqdm(total=len(regions_pct))

    for progress, data in zip(regions_pct, regions_data):

        if data.empty:
            progress_bar.set_description('Skipped empty region')
            progress_bar.update()
            continue

        div = est.js_divergence(data)

        if div.empty:
            progress_bar.set_description('Skipped low-quality region')
            progress_bar.update()
            continue

        if not os.path.isfile(outfile) or os.stat(outfile).st_size == 0:
            header = True
        else:
            header = False

        (div
         .round({'JSD_bit_': 3, 'HMIX_bit_': 3})
         .to_csv(outfile, header=header, sep='\t', index=True, mode='a'))

        progress_bar.set_description('Progress')
        progress_bar.update()

    progress_bar.close()
    logging.info("The File is created!")

    return None
