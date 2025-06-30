#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: zmk214
"""

import os
import numpy as np
import pandas as pd


def extract_filenames(dir_path, keep_word=None, drop_word=None):

    filename_ls = os.listdir(dir_path)

    if keep_word:
        filename_ls = [name for name in filename_ls if all(substring in name for substring in keep_word)]
    elif keep_word is not None:
        raise ValueError("Invalid parameter type for keep_word. Expected list of strings.")

    if drop_word:
        filename_ls = [name for name in filename_ls if all(substring not in name for substring in drop_word)]
    elif drop_word is not None:
        raise ValueError("Invalid parameter type for drop_word. Expected list of strings.")

    if '.DS_Store' in filename_ls:
        filename_ls.remove('.DS_Store')

    filename_ls.sort()

    return(filename_ls)

def extract_ranges(df):

    """
    Extracts chromosome, start, end, and strand from DataFrame index.

    Args:
    df (pandas.DataFrame): DataFrame with index containing strings in the format "chrom:start-end;strand".

    Returns:
    pandas.DataFrame: DataFrame with 'chrom', 'chromStart', 'chromEnd', and 'strand' columns.
    """
    # Split index and expand into separate columns
    split_index = df.index.str.split(':|-|;', expand=True)

    # Untuple the split_index
    chrom, chromStart, chromEnd, strand = zip(*split_index.values)

    # Create a DataFrame with column names
    result_df = pd.DataFrame({
        'chrom': chrom,
        'chromStart': chromStart,
        'chromEnd': chromEnd,
        'strand': strand
    })

    return result_df

def extract_profiles(df):
    """
    Extract profile vectors from a DataFrame.

    Parameters:
    df (pandas.DataFrame): DataFrame containing profile data.

    Returns:
    numpy.ndarray: NumPy array representing profile vectors.
    """
    # Convert DataFrame to NumPy array
    numpy_array = df.values

    # Convert NaN values to 0
    # numpy_array[np.isnan(numpy_array)] = 0

    return numpy_array
