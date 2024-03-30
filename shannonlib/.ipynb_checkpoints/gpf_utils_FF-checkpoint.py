# -*- coding:utf-8 -*-
# gpf_utils.py

"""
Utility functions for handling genome position files.
"""

import logging
import math
import subprocess
import numpy as np
import pandas as pd

logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')


class GPFUtilsError(Exception):
    """Base exception class for GPF utility errors."""
    pass


def get_regions(files, chrom=None, exp_numsites=1e3):
    """Get stepsize and list of regions for tabix-indexed files."""

    sup_position = supremum_position(files, chrom)
    if sup_position is None:
        logging.info("Skipping because chromosome is missing.")
        return False

    sup_numsites = supremum_numsites(files, chrom)
    if sup_numsites in (None, 0):
        logging.info("Skipping because there are no entries.")
        return False

    stepsize = math.ceil(sup_position / sup_numsites * exp_numsites)
    stepsize = min(stepsize, sup_position)

    pos_start = list(range(0, sup_position, stepsize + 1))
    pos_end = list(range(stepsize, sup_position, stepsize + 1)) + [sup_position]
    progress = [round(100 * pos / sup_position, 1) for pos in pos_end]

    regions = zip([chrom] * len(pos_start), pos_start, pos_end)
    return progress, regions


def get_data(files, keys, data_columns, regions, join='outer', preset='bed'):
    """Combines tabix-indexed genome position files."""

    validate_input(files, data_columns)
    index, index_col, levels = get_preset(preset)
    columns = [index + cols for cols in data_columns]

    dframes = []
    for file_, region in zip(files, regions):
        for i, col in enumerate(columns):
            dframe = read_tabix_output(file_, region, index_col, col, keys, levels, join)
            dframes.append(dframe)

    merged_dframe = pd.concat(dframes, axis=1, keys=keys, names=levels, join=join)
    return merged_dframe


def supremum_numsites(files, chrom):
    """Return the least upper bound for the number of covered sites."""
    return process_tabix_files(files, chrom, command=["wc", "-l"])


def supremum_position(files, chrom):
    """Return the least upper bound for the chrom end coordinate."""
    return process_tabix_files(files, chrom, command=["tail", "-1"], cut_field=3)


def process_tabix_files(files, chrom, command, cut_field=None):
    """General function to process tabix file output for supremum calculations."""
    preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"
    results = []

    for f in files:
        full_path = f"{preF}{f.replace('myProcessed', '').lstrip('/')}"
        tabix = subprocess.Popen(["tabix", full_path, chrom], stdout=subprocess.PIPE)
        process = subprocess.Popen(command, stdin=tabix.stdout, stdout=subprocess.PIPE)
        tabix.stdout.close()

        output, _ = process.communicate()
        if cut_field is not None:
            output = output.decode('utf-8').strip().split()[cut_field-1]

        try:
            result = int(output)
            results.append(result)
        except ValueError:
            continue

    return max(results, default=None)


def validate_input(files, data_columns):
    if len(data_columns) not in (1, len(files)):
        raise GPFUtilsError('Mismatch in the number of files and data columns provided.')


def get_preset(preset):
    if preset == 'bed':
        index = [(0, '#chrom', str), (1, 'start', np.int64), (2, 'end', np.int64)]
    # Implement other presets as needed
    return index, [i[0] for i in index], ['sample', 'event']


def read_tabix_output(file_, region, index_col, columns, keys, levels, join):
    """Process individual tabix file output for a given region."""
    query = '{0}:{1}-{2}'.format(*region)
    preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"
    full_path = f"{preF}{file_.replace('myProcessed', '').lstrip('/')}"

    with subprocess.Popen(["tabix", full_path, query], stdout=subprocess.PIPE, universal_newlines=True) as process:
        return pd.read_table(process.stdout, header=None, index_col=index_col, comment='#',
                             usecols=[f[0] for f in columns], names=[f[1] for f in columns], 
                             dtype={f[1]: f[2] for f in columns})