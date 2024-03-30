# -*- coding:utf-8 -*-
# gpf_utils.py

"""Utility functions for handling genome position files.
"""

import logging
import math
import subprocess

import numpy as np
import pandas as pd
logger = logging.getLogger(__name__)

class InputMismatchError(Exception):
    pass


class MissingInputError(Exception):
    pass


logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')


def get_regions(files, chrom=None, exp_numsites=1e3):
    """Get stepsize and list of regions for tabix-indexed files.
    """

    sup_position = supremum_position(files, chrom)

    #print("sup_position", sup_position)

    if sup_position == None:
        logging.info("Skipping because chromosome is missing.")
        return False

    sup_numsites = supremum_numsites(files, chrom)
    print("sup_numsites", sup_numsites)

    if sup_numsites == None or sup_numsites == 0:
        logging.info("Skipping because there are no entries.")
        return False

    step = math.ceil(sup_position / sup_numsites * exp_numsites)

    if step < sup_position:
        stepsize = step
    else:
        stepsize = sup_position

    pos_start = list(
        range(0, sup_position, stepsize + 1))
    pos_end = list(
        range(stepsize, sup_position, stepsize + 1)) + [sup_position]

    progress = [round(100 * pos / sup_position, 1) for pos in pos_end]

    regions = zip([chrom] * len(pos_start), pos_start, pos_end)

    #print("REGIONS", list(regions))
    #print("PROGRESS", progress)

    return progress, regions


def get_data(files, keys=None, data_columns=None,
             regions=None, join='outer', preset='bed'):
    """Combines tabix-indexed genome position files.

    data_columns : list of tuple list
    """

    if data_columns is None:
        raise MissingInputError(
            'The list of data_colums must have at least one entry!')
    elif len(data_columns) == len(files):
        pass
    elif len(data_columns) == 1:
        data_columns = data_columns * len(files)
    else:
        raise InputMismatchError(
            'Either supply a single entry in data_columns or'
            'the number of entries must match the number of files!')

    if preset == 'bed':
        index = [(0, '#chrom', str),
                 (1, 'start', np.int64),
                 (2, 'end', np.int64)]
        index_col = [i[0] for i in index]
    if preset == 'gff':
        # TODO
        pass
    if preset == 'vcf':
        # TODO
        pass
    if preset == 'sam':
        # TODO
        pass

    levels = ['sample', 'event']

    columns = [index + cols for cols in data_columns]

    for region in regions:
        query = '{0}:{1}-{2}'.format(*region)

        logger.info(f'QUERY for REGION {region}: {query}')


        """
        tabix = enumerate(
            subprocess.Popen(
                ['tabix', file_, query],
                stdout=subprocess.PIPE,
                universal_newlines=True)
            for file_ in files)
        """

        preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"

        
        tabix = enumerate(
            subprocess.Popen(
                ['tabix', f"{preF}{file_.replace('myProcessed', '').lstrip('/')}", query],
                stdout=subprocess.PIPE,
                universal_newlines=True)
            for file_ in files)

        
        dframes = (
            pd.read_table(
                tbx.stdout,
                header=None,
                index_col=index_col,
                comment='#',
                usecols=[f[0] for f in columns[i]],
                names=[f[1] for f in columns[i]],
                dtype={f[1]: f[2] for f in columns[i]})
            for i, tbx in tabix)

        merged_dframe = pd.concat(
            dframes, axis=1, keys=keys, names=levels, join=join)

        print("#######################################")
        print(merged_dframe.head())
        print("#######################################")

        yield merged_dframe


def supremum_numsites(tabixfiles, chrom):
    '''Return the least upper bound for the number of covered sites.'''

    preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"

        
    sites = list()

    #print('\033[94m' + "THIS IS FILE"+ tabixfiles)
    #print('\x1b[0m')

    for f in tabixfiles[0:3]:
        f_modified = f.replace("myProcessed", "").lstrip('/')
        full_path = preF + f_modified

        #print('\033[94m' + "THIS IS PATH", full_path + '\x1b[0m')
       # print('\033[94m' + "THIS IS FILE"+ full_path)
        
        tabix = subprocess.Popen(["tabix", full_path, chrom], stdout=subprocess.PIPE)
        wcl = subprocess.Popen(
            ["wc", "-l"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        tabix.stdout.close()  # Allow tabix to receive a SIGPIPE if wcl exits.
        try:
            site_count = int(wcl.communicate()[0])
            sites.append(site_count)
        except ValueError:
            continue

    try:
        out = np.max(sites)
    except ValueError:
        out = None

    return out


def supremum_position(tabixfiles, chrom):
    """Return the least upper bound for the chrom end coordinate."""

    end_coordinate = []

    preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"

    for f in tabixfiles:
        f_modified = f.replace("myProcessed", "").lstrip('/')
        full_path = preF + f_modified

        #print('\033[93m' + "THIS IS PATH", full_path + '\x1b[0m')
        
        # Execute tabix and pipe its output to tail
        tabix = subprocess.Popen(["tabix", full_path, chrom], stdout=subprocess.PIPE)
        tail = subprocess.Popen(["tail", "-1"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        cut = subprocess.Popen(["cut", "-f3"], stdin=tail.stdout, stdout=subprocess.PIPE)
        
        # Close the stdout of tabix and tail to signal EOF to the subprocesses
        tabix.stdout.close()
        tail.stdout.close()
        
        # Wait for the cut process to finish and capture its output
        cut_output, _ = cut.communicate()
        cut_output = cut_output.decode('utf-8').strip()
        
        if cut_output:  # Check if the output is not empty
            try:
                base_position = int(cut_output)  
                end_coordinate.append(base_position)  
            except ValueError as e:
                print(f"Failed to convert base position to int for file: {f}, Error: {e}")
        else:
            print(f"No data for specified chromosome in file: {f}")
            
    
    print("End Coordinates:", end_coordinate)
    try:
        out = np.max(end_coordinate)
    except ValueError:
        out = None

    print(out)

    return out




def supremum_position_____(tabixfiles, chrom):
    """Return the least upper bound for the chrom end coordinate.
    """

    end_coordinate = list()

    for f in tabixfiles:
        preF = "/shares/grossniklaus.botinst.uzh/dkt/scienceCloud/okartal_marcws_processed/"
        f_modified = f.replace("myProcessed", "").lstrip('/')
        f = preF + f_modified
        
        tabix = subprocess.Popen(["tabix", f, chrom], stdout=subprocess.PIPE)
        tabix.wait()

        tabix_output = tabix.communicate()[0].decode('utf-8').strip()
        print("Tabix Command Output:", tabix_output)

        tail = subprocess.Popen(["tail", "-1"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        tail.wait()

        tail_output = tail.communicate()[0].decode('utf-8').strip()
        print("Tail Command Output:", tail_output)

        
        cut = subprocess.Popen(["cut", "-f3"], stdin=tail.stdout, stdout=subprocess.PIPE)
        cut.wait()

        cut_output = cut.communicate()[0].decode('utf-8').strip()  
        print("Cut Command Output:", cut_output) 
        
        try:
            base_position = int(cut_output)  
            end_coordinate.append(base_position)  
        except ValueError:
            print("Failed to convert base position to int for file:", f)  
            continue  

        
        # Allow first process to receive a SIGPIPE if process 2 exits.
        tabix.stdout.close()
        try:
            base_position = int(cut.communicate()[0])
            end_coordinate.append(base_position)
        except ValueError:
            continue

    print("End Coordinates:", end_coordinate)
    try:
        out = np.max(end_coordinate)
    except ValueError:
        out = None

    return out
