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
    Ensures that the final result retains the correct indices that match the profile.

    If sample_size is specified, randomly selects that many samples; otherwise, returns all samples.
    Saves the processed data and related indices for later use.

    Parameters:
    - profile_dir (str): Path to the directory containing profile files.
    - file_filter (str): Substring to filter files (e.g., "_pos_" or "_neg_").
    - output_path (str): Path to save the processed test data.
    - sample_size (int or None): Number of samples to extract (if None, take all).
    - random_seed (int): Seed for random sampling (default=42 for reproducibility).

    Returns:
    - test_X (numpy array): Processed and sampled test data.
    - final_indices (list): Indices that correspond to the selected profiles.
    """

    # Step 1: Extract filenames in the directory
    profile_file_ls = extract_filenames(profile_dir)

    # Step 2: Filter files based on the provided substring
    filtered_files = [f for f in profile_file_ls if file_filter in f]

    if not filtered_files:
        print(f"No files matching '{file_filter}' found in {profile_dir}. Returning an empty array.")
        return np.array([]), []

    def load_profile(file_path):
        """Dynamically loads CSV or Parquet profile files and retains indices."""
        if file_path.endswith(".csv"):
            df = pd.read_csv(file_path, header=0, index_col=None)
        elif file_path.endswith(".parquet"):
            df = pd.read_parquet(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path}")
        return set_index_if_exists(df)

    # Step 3: Load and process profiles while tracking indices
    profiles = []
    index_list = []

    for f in filtered_files:
        file_path = os.path.join(profile_dir, f)
        pre_profile = load_profile(file_path)
        extracted_profile = extract_profiles(pre_profile)

        if len(extracted_profile) > 0:
            profiles.append(extracted_profile)
            index_list.append(pre_profile.index)  # Keep track of indices

    # Step 4: Ensure we have valid profiles
    if not profiles:
        print("No valid profiles were processed.")
        return np.array([]), []

    # Step 5: Flatten profiles and corresponding indices
    concatenated_profiles = np.concatenate(profiles, axis=0)
    matched_indices = np.concatenate(index_list, axis=0)

    print(f"Processed {len(filtered_files)} files. Final shape: {concatenated_profiles.shape}")

    # Step 6: Handle Sampling
    total_samples = concatenated_profiles.shape[0]

    if sample_size is None or sample_size >= total_samples:
        test_X = concatenated_profiles  # Take all samples
        final_indices = matched_indices  # Keep all indices
        print(f"Returning all {total_samples} samples.")
    else:
        np.random.seed(random_seed)  # Ensure reproducibility
        selected_indices = np.random.choice(total_samples, sample_size, replace=False)
        test_X = concatenated_profiles[selected_indices]  # Select profiles
        final_indices = matched_indices[selected_indices]  # Select corresponding indices
        print(f"Randomly selected {sample_size} samples.")

    # Step 7: Save Processed Data and Indices
    os.makedirs(os.path.dirname(output_path), exist_ok=True)  # Ensure directory exists

    # Save the processed data as compressed .npz file
    np.savez_compressed(output_path+".npz", test_X=test_X)

    # Save the corresponding indices as a text file
    index_output_path = output_path + ".index"
    np.savetxt(index_output_path, final_indices, fmt="%s", delimiter="\n")

    print(f"Saved processed data to {output_path}")
    print(f"Saved indices to {index_output_path}")

    return test_X, final_indices


profile_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres/K562_endres_profiles/profiles_subtnorm/"
output_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP"

# Save positive samples
test_X_pos = prepare_shap_data(profile_dir,
                               file_filter="K562_C",
                               output_path=os.path.join(output_dir, "K562_C_endres_v2"),
                               sample_size=None)
