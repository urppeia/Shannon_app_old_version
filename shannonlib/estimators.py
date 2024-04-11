# -*- coding:utf-8 -*-
# estimators.py


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
    """Shannon entropy (in nat) of the feature frequency profile.
    """

    if method == 'plug-in':
        expression = ("sum(where(prob > 0, -prob * log(prob), 0), axis={})"
                      .format(axis))
        count_distribution = countmatrix.sum(axis)[..., np.newaxis]
        prob = countmatrix / count_distribution
        return ne.evaluate(expression)


def js_divergence(indata, weights=None): # indata is DF with multiindex: 
    """
    """

    LOG2E = np.log2(math.e)

    # if imputation:
    #     impute_value = impute(indata, method='pseudocount')
    #     indata.fillna(impute_value, inplace=True)


    count_per_unit = indata.groupby(axis=1, level=['sampling_unit']).sum()
    
    #count_per_unit = indata.sum(axis=1, index='sampling_unit')
    #count_per_unit = indata.groupby(level=0, axis=1).sum()
    #count_per_unit = indata.T.groupby(level=0).sum()
    #logging.info(f"count_per_unit: {count_per_unit.head()}")
    samplesize = count_per_unit.notnull().sum(axis=1)
    #logging.info(f"samplesize for count_per_unit: {samplesize.head()}")

    # QC filter
    min_samplesize = 2
    min_count = 3
    count_filter = (count_per_unit >= min_count).any(axis=1)
    samplesize_filter = (samplesize >= min_samplesize)
    combined_filter = (count_filter & samplesize_filter)
    data = indata[combined_filter]


    if data.empty:
        return data
    else:
        data_unit = count_per_unit[combined_filter]
        #data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
        #data_feature = data.groupby(level=1, axis=1).sum().astype(np.int32)
        data_feature = data.groupby(axis=1, level=['feature']).sum().astype(np.int32)
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
        div.insert(0, 'JSD_bit_', LOG2E *
                   (mix_entropy - avg_entropy))
        div.insert(1, 'sample size', samplesize[combined_filter])
        div.insert(2, 'HMIX_bit_', LOG2E * mix_entropy)
        # div['members'] = (data_unit.apply(lambda x: ','.join(x.dropna().index), axis=1))

        return div




def js_divergence_(indata, weights=None):
    """
    """

    LOG2E = np.log2(math.e)
    # if imputation:
    #     impute_value = impute(indata, method='pseudocount')
    #     indata.fillna(impute_value, inplace=True)

    count_per_unit = indata.groupby(axis=1, level=['sampling_unit']).sum()
    samplesize = count_per_unit.notnull().sum(axis=1)

    # QC filter
    min_samplesize = 2
    min_count = 3
    count_filter = (count_per_unit >= min_count).any(axis=1)
    samplesize_filter = (samplesize >= min_samplesize)
    combined_filter = (count_filter & samplesize_filter)
    data = indata[combined_filter]

    if data.empty:
        return data
    else:
        data_unit = count_per_unit[combined_filter]
       # data_feature = (data.sum(axis=1, level='feature').astype(np.int32))
        data_feature = data.groupby(axis=1, level=['feature']).sum().astype(np.int32)
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
        div.insert(0, 'JSD_bit_', LOG2E *
                   (mix_entropy - avg_entropy))
        div.insert(1, 'sample size', samplesize[combined_filter])
        div.insert(2, 'HMIX_bit_', LOG2E * mix_entropy)
        # div['members'] = (data_unit.apply(lambda x: ','.join(x.dropna().index), axis=1))

        return div






















