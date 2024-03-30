# -*- coding:utf-8 -*-
# core.py

"""Core functions.
"""

import os

import shannonlib.estimators as est
import shannonlib.gpf_utils as gpf

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

    for progress, data in zip(regions_pct, regions_data):

        if data.empty:
            print('...{:>5} % (skipped empty region)'.format(progress))
            continue

        # data is pandas DF
        #logging.info(f"Column name of the data: {data.columns}, info: {data.info()}")
        
        div = est.js_divergence(data)

        if div.empty:
            print('...{:>5} % (skipped low-quality region)'.format(progress))
            continue

        if not os.path.isfile(outfile):
            header = True
        elif os.stat(outfile).st_size == 0:
            header = True
        else:
            header = False

        (div
         .round({'JSD_bit_': 3, 'HMIX_bit_': 3})
         .to_csv(outfile, header=header, sep='\t', index=True, mode='a'))

        print('...{:>5} %'.format(progress))

    logging.info("The File is created!")

    return None
