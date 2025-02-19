import pandas as pd
import os

def filter_profiles(prediction_path, profile_path, output_path, annotation_filter=None, score_threshold=None):
    """
    Filters profile data based on conditions in the prediction file.
    If no annotation or score is specified, it takes all annotations or scores.

    Parameters:
    - prediction_path (str): Path to the BED file containing annotation and score.
    - profile_path (str): Path to the profile file (CSV or Parquet).
    - output_path (str): Path to save the filtered profiles.
    - annotation_filter (str or None): Annotation type to filter (default=None, meaning take all).
    - score_threshold (float or None): Score threshold for filtering (default=None, meaning take all).

    Returns:
    - Saves the filtered profile data to a CSV file.
    """

    # Load prediction file (assuming it's tab-separated with a header)
    pred_df = pd.read_csv(prediction_path, sep="\t", header=0)

    # Load profile file (detect CSV or Parquet)
    if profile_path.endswith(".csv"):
        profile_df = pd.read_csv(profile_path, header=0)
    elif profile_path.endswith(".parquet"):
        profile_df = pd.read_parquet(profile_path)
    else:
        raise ValueError(f"Unsupported file format: {profile_path}")

    # Ensure row counts match before merging
    if pred_df.shape[0] != profile_df.shape[0]:
        raise ValueError("Row counts do not match between prediction and profile files!")

    # Merge by column (row-wise alignment)
    merged_df = pd.concat([pred_df, profile_df], axis=1)
    print(f"Merged dataset with {merged_df.shape[0]} rows and {merged_df.shape[1]} columns.")

    # Apply filtering conditions
    if annotation_filter:
        merged_df = merged_df[merged_df["annotation"] == annotation_filter]
        print(f"Filtered by annotation: '{annotation_filter}', remaining rows: {merged_df.shape[0]}")
    
    if score_threshold is not None:
        merged_df = merged_df[merged_df["score"] >= score_threshold]
        print(f"Filtered by score >= {score_threshold}, remaining rows: {merged_df.shape[0]}")

    # Extract only profile data (columns starting with "Pos")
    profile_columns = [col for col in merged_df.columns if col.startswith("Pos")]
    filtered_profiles_df = merged_df[profile_columns]

    # Ensure output directory exists before saving
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Save filtered profiles
    filtered_profiles_df.to_csv(output_path, index=False)
    print(f"Filtered profiles saved to: {output_path}")



# Define paths
prediction_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/GM12878_wt10M_profiles_te_withAnno/predictions/PRIMEloci_pred_all_profiles_subtnorm_pos_te_GM12878_reseq_C_1mio_1.bed"
profile_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/GM12878_wt10M_profiles_te/profiles_subtnorm/profiles_subtnorm_pos_te_GM12878_reseq_C_1mio_1.csv"  # Supports Parquet too
output_path = "/Users/natsudanav/Documents/filtered_promoter_0.75_profiles_subtnorm_pos_te_GM12878_reseq_C_1mio_1.csv"

# Run filtering
filter_profiles(prediction_path,
                profile_path,
                output_path,
                annotation_filter="promoter",
                score_threshold=0.75)
