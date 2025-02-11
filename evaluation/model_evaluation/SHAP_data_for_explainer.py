import numpy as np
import pandas as pd
import os

from python.extraction import extract_filenames, extract_profiles
from python.helpfn_for_prediction import set_index_if_exists


def prepare_shap_data(profile_dir, 
                      file_filter,
                      output_path,
                      sample_size=None,
                      random_seed=42):
    """
    Extracts and processes test profiles based on a filename filter.
    If sample_size is specified, randomly selects that many samples; otherwise, returns all samples.
    Saves the processed data for later use.

    Parameters:
    - profile_dir (str): Path to the directory containing profile files.
    - file_filter (str): Substring to filter files (e.g., "_pos_" or "_neg_").
    - output_path (str): Path to save the processed test data.
    - sample_size (int or None): Number of samples to extract (if None, take all).
    - random_seed (int): Seed for random sampling (default=42 for reproducibility).

    Returns:
    - test_X (numpy array): Processed and sampled test data.
    """

    # Step 1: Extract filenames in the directory
    profile_file_ls = extract_filenames(profile_dir)

    # Step 2: Filter files based on the provided substring
    filtered_files = [f for f in profile_file_ls if file_filter in f]

    if not filtered_files:
        print(f"No files matching '{file_filter}' found in {profile_dir}. Returning an empty array.")
        return np.array([])

    def load_profile(file_path):
        """Dynamically loads CSV or Parquet profile files."""
        if file_path.endswith(".csv"):
            df = pd.read_csv(file_path, header=0, index_col=None)
        elif file_path.endswith(".parquet"):
            df = pd.read_parquet(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path}")
        return set_index_if_exists(df)

    # Step 3: Load and process profiles
    profiles = [extract_profiles(load_profile(os.path.join(profile_dir, f))) for f in filtered_files]

    # Step 4: Concatenate profiles if available
    if profiles:
        concatenated_profiles = np.concatenate(profiles, axis=0)
        print(f"Processed {len(filtered_files)} files. Final shape: {concatenated_profiles.shape}")
    else:
        concatenated_profiles = np.array([])
        print("No valid profiles were processed.")

    # Step 5: Handle Sampling
    total_samples = concatenated_profiles.shape[0]

    if sample_size is None or sample_size >= total_samples:
        test_X = concatenated_profiles  # Take all samples
        print(f"Returning all {total_samples} samples.")
    else:
        np.random.seed(random_seed)  # Ensure reproducibility
        selected_indices = np.random.choice(total_samples, sample_size, replace=False)
        test_X = concatenated_profiles[selected_indices]  # Randomly select samples
        print(f"Randomly selected {sample_size} samples.")

    # Step 6: Save Processed Data
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure directory exists
    np.savez_compressed(output_path, test_X=test_X)  # Save as compressed .npz file
    print(f"Saved processed data to {output_path}")

    return test_X


profile_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres/K562_endres_profiles/profiles_subtnorm/"
output_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP"

# Save positive samples
test_X_pos = prepare_shap_data(profile_dir,
                               file_filter="K562_N",
                               output_path=os.path.join(output_dir, "K562_N_endres.npz"),
                               sample_size=None)
