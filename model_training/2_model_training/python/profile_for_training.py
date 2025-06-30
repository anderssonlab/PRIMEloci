import os
import numpy as np
import pandas as pd
from python.extraction import extract_profiles


def profile_for_training(profile_dir, profile_file_ls, label):
    """
    Loads profile data from CSV files, concatenates them into a DataFrame,
    extracts features and labels, and returns them in numpy array format.

    Parameters:
    profile_dir (str): Directory containing the profile files.
    profile_file_ls (list): List of profile filenames to process.
    label (int): Label to assign to the profiles 
    (e.g., 1 for positive, 0 for negative).

    Returns:
    tuple: A tuple containing two numpy arrays:
           - profile_np: Numpy array of profile features.
           - label_np: Numpy array of labels corresponding to the profiles.
    """
    data_in = pd.DataFrame()

    for filename in profile_file_ls:
        file_path = os.path.join(profile_dir, filename)
        read_in = pd.read_csv(file_path, header=0, index_col="rownames")
        print(f"{filename} contains {len(read_in)} rows")
        data_in = pd.concat([data_in, read_in])

    data_in.reset_index(inplace=True, drop=True)

    profile_np = extract_profiles(data_in)
    label_np = np.array([label] * len(profile_np))

    return profile_np, label_np
