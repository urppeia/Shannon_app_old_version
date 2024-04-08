# -*- coding:utf-8 -*-
# # estimators.py

"""Estimators of information-theoretic quantities.
"""

import numexpr as ne
import numpy as np
import pandas as pd
import math

#import shannonlib.constants as constant

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s",
                    level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')


def shannon_entropy(countmatrix, axis=1, method='plug-in'):
    """Shannon entropy (in nat) of the feature frequency profile."""
    
    #logger.debug("Calculating Shannon entropy...")
    
    if method == 'plug-in':
        expression = ("sum(where(prob > 0, -prob * log(prob), 0), axis={})"
                      .format(axis))
        count_distribution = countmatrix.sum(axis)[..., np.newaxis]
        prob = countmatrix / count_distribution
        entropy = ne.evaluate(expression)
        
        #logger.debug("Shannon entropy calculation complete.")
        return entropy


def js_divergence(indata, weights=None): # indata is DF with multiindex: 

    LOG2E = np.log2(math.e)
    
    #logger.debug("Calculating count per unit...")
    count_per_unit = indata.groupby(level=0, axis=1).sum()
    samplesize = count_per_unit.notnull().sum(axis=1)
    
    #logger.debug("Applying QC filters...")
    min_samplesize = 2
    min_count = 3
    count_filter = (count_per_unit >= min_count).any(axis=1)
    samplesize_filter = (samplesize >= min_samplesize)
    combined_filter = (count_filter & samplesize_filter)
    data = indata[combined_filter]

    if data.empty:
        return data
    else:
       # logger.debug("Calculating JS divergence...")
        data_unit = count_per_unit[combined_filter]
        data_feature = data.groupby(level=1, axis=1).sum().astype(np.int32)
        
        counts = data.values.reshape(
            data_unit.shape[0],
            data_unit.shape[1],
            data_feature.shape[1])

        mix_entropy = shannon_entropy(data_feature.values)
        avg_entropy = np.average(
            shannon_entropy(counts, axis=2),
            weights=data_unit.fillna(0),
            axis=1)

        div = data_feature
        div.insert(0, 'JSD_bit_', LOG2E * (mix_entropy - avg_entropy))
        div.insert(1, 'sample size', samplesize[combined_filter])
        div.insert(2, 'HMIX_bit_', LOG2E * mix_entropy)

        #logger.debug("JS divergence calculation complete.")
        
        return div
