import os
import pickle
import pandas as pd
from lightgbm import LGBMClassifier

import argparse

#import pyPRIMEloci

from python.extraction import extract_filenames
from python.extraction import extract_ranges
from python.extraction import extract_profiles


# Function to load environment variables from .env file
def load_env_vars(script_dir, profile_main_dir, subdir_name, model_path):
    
    # Set the working directory
    os.chdir(script_dir)

    # Import model
    model = pickle.load(open(model_path, 'rb'))

    # Ensure the output directory exists
    if not os.path.exists(profile_main_dir+"/predictions"):
        os.makedirs(profile_main_dir+"/predictions")
        print(f"Directory '{profile_main_dir}/predictions' created successfully.")
    else:
        print(f"Directory '{profile_main_dir}/predictions' already exists.")

    # Ensure the output sub-directory exists
    if not os.path.exists(profile_main_dir+"/predictions/"+subdir_name):
        os.makedirs(profile_main_dir+"/predictions/"+subdir_name)
        print(f"Directory '{profile_main_dir}/predictions/{subdir_name}' created successfully.")
    else:
        print(f"Directory '{profile_main_dir}/predictions/{subdir_name}' already exists.") 

    return script_dir, profile_main_dir, subdir_name, model


def wrapup_model_prediction(script_dir, profile_dir, profile_filename, metadata_dir, metadata_filename, model, output_dir, name_prefix, threshold):

    # Example usage of extract_filenames

    filenames_without_extensions = os.path.splitext(profile_filename)[0]
    print(filenames_without_extensions)

    # Define the input paths
    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
    input_metadata_path = os.path.join(metadata_dir, metadata_filename)

    # Read CSV files
    subtnorm_df = pd.read_csv(input_profiles_subtnorm_path,
                              header=0,
                              index_col=0)
    metadata_df = pd.read_csv(input_metadata_path,
                              header=0,
                              index_col=None)

    if len(subtnorm_df) != len(metadata_df):
        raise ValueError("The metadata and profile subtnorm dataframes do not have the same length")
        sys.exit(1)

    # Extract ranges_df and profiles
    ranges_df = extract_ranges(subtnorm_df)
    subtnorm_np = extract_profiles(subtnorm_df)

    # Predict probabilities
    y_proba = model.predict_proba(subtnorm_np)
    ranges_df['score'] = y_proba[:, 1]

    ranges_df['name'] = filenames_without_extensions

    # Reorder columns
    ranges_df = ranges_df.reindex(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

    # Extract extra column from metadata
    columns_to_keep = ["seqnames", "start", "end", "width", "strand", 'chrom', 'chromStart', 'chromEnd', 'name', 'score']
    columns_to_move = [col for col in metadata_df.columns if col not in columns_to_keep]

    # Move columns to ranges_df
    for col in columns_to_move:
        ranges_df[col] = metadata_df[col]

    # Save results to .bed files
    output_all_results = os.path.join(output_dir, f'{name_prefix}_pred_all_{filenames_without_extensions}.bed')
    ranges_df.to_csv(output_all_results, sep='\t', header=True, index=False)

    # Save selected results if score >= threshold
    output_slt_results = os.path.join(output_dir, f'{name_prefix}_pred_slt{threshold}_{filenames_without_extensions}.bed')
    selected_ranges = ranges_df[ranges_df['score'] >= threshold]
    selected_ranges.to_csv(output_slt_results, sep='\t', header=True, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to predict genomewide TC normalized data.')

    parser.add_argument('-w', '--script_dir', default=".", type=str, 
                        help='Path to the script directory')
    parser.add_argument('-p', '--profile_main_dir', type=str, required=True, 
                        help='Path to the profile directory')
    parser.add_argument('-r', '--profile_sub_dir', type=str, default="tcs", 
                        help='Sub-directory name in the main profile directory')
    parser.add_argument('-m', '--model_path', type=str, required=True, 
                        help='Path to the input model')
    parser.add_argument('-n', '--name_prefix', type=str, required=True, 
                        help='Name added to the output files, indicate model name and library (celltype) name') 
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Threshold value for prediction')

    args = parser.parse_args()

    # Load environment variables
    script_dir = args.script_dir
    profile_main_dir = args.profile_main_dir
    profile_sub_dir = args.profile_sub_dir.split(',')
    model_path = args.model_path

    name_prefix = args.name_prefix

    threshold = args.threshold

    for subdir_name in profile_sub_dir:

        # Load environment variables and model
        script_dir, profile_main_dir, subdir_name, model = load_env_vars(script_dir, profile_main_dir, subdir_name, model_path)

        profile_dir = profile_main_dir+"/profiles_subtnorm/"+subdir_name
        metadata_dir = profile_main_dir+"/metadata/"+subdir_name
        output_dir = profile_main_dir+"/predictions/"+subdir_name

        profile_file_ls = extract_filenames(profile_dir)
        metadata_file_ls = extract_filenames(metadata_dir)

        # Execute data processing  
        for i in range(len(profile_file_ls)):
            wrapup_model_prediction(script_dir,
                                    profile_dir,
                                    profile_file_ls[i],
                                    metadata_dir,
                                    metadata_file_ls[i],
                                    model,
                                    output_dir,
                                    name_prefix,
                                    threshold)
