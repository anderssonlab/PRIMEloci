import os
import pandas as pd
import numpy as np
import re


def extract_exclude_patterns(directory, extract_pattern):
    """
    Extracts patterns from filenames in a directory based on the given regex.

    Parameters:
        directory (str): The directory to scan for files.
        extract_pattern (re.Pattern): The regex pattern to extract parts from filenames.

    Returns:
        list: A list of extracted patterns.
    """
    exclude_patterns = []
    for filename in os.listdir(directory):
        match = extract_pattern.search(filename)
        if match:
            exclude_patterns.append(match.group(1))
    return exclude_patterns


def set_score_to_nan_where_sum_count_zero(df, set_nan):
    """
    Sets 'score' to NaN where 'sum_count' is 0 in the dataframe if set_nan is True.

    Parameters:
        df (pd.DataFrame): The dataframe to modify.
        set_nan (bool): Whether to set 'score' to NaN.

    Returns:
        pd.DataFrame: The modified dataframe.
    """
    if set_nan:
        print("Filtering: Setting 'score' to NaN where 'sum_count' = 0.")
        df['score'] = df.apply(
            lambda row: np.nan if row['sum_count'] == 0 else row['score'],
            axis=1
        )
    return df


def process_prediction_data(directory, pattern, outname, noCAGEtoNaN=True, save_csv=True):
    """
    Processes BED files, applies filtering, and splits data into 'score' and 'sum_count' dataframes.
    Saves results to files. Row names are formatted as chrom:chromStart+1-chromEnd;strand (R-format).

    Parameters:
        directory (str): The directory containing the files.
        pattern (re.Pattern): The regex pattern to extract sample names from filenames.
        outname (str): Output filename prefix for saving results.
        noCAGEtoNaN (bool): Whether to set 'score' to NaN where 'sum_count' = 0.
        save_csv (bool): Whether to save the resulting dataframes to TSV files.
    """
    # Extract sample names from the directory
    samples = extract_exclude_patterns(directory, pattern)

    # Initialize dataframes for scores and sum_counts
    df_scores = pd.DataFrame()
    df_sum_counts = pd.DataFrame()

    # Loop through each sample
    for sample in samples:
        # Construct the expected filename pattern
        filename_pattern = f"FANTOM5-cellfacet-on-PRIMEloci_pred_all_{sample}_combined.bed"

        # Look for the file in the directory that matches the pattern
        for filename in os.listdir(directory):
            if filename == filename_pattern:
                # Construct full file path
                file_path = os.path.join(directory, filename)

                # Read the BED file into a dataframe
                df = pd.read_csv(file_path, sep='\t', header=0)

                # Create row names in the format chrom:chromStart+1-chromEnd;strand
                df['row_name'] = df.apply(
                    lambda row: f"{row['chrom']}:{row['chromStart'] + 1}-{row['chromEnd']};{row['strand']}" 
                    if 'strand' in df.columns else f"{row['chrom']}:{row['chromStart'] + 1}-{row['chromEnd']}",
                    axis=1
                )
                df.set_index('row_name', inplace=True)

                # Apply noCAGEtoNaN filtering to scores
                if noCAGEtoNaN:
                    print("Filtering: Setting 'score' to NaN where 'sum_count' = 0.")
                    df['score'] = df.apply(
                        lambda row: np.nan if row['sum_count'] == 0 else row['score'],
                        axis=1
                    )

                # Extract and add score column to `df_scores`
                score_col = df['score'].rename(sample)
                df_scores = pd.concat([df_scores, score_col], axis=1)

                # Extract and add sum_count column to `df_sum_counts`
                sum_count_col = df['sum_count'].rename(sample)
                df_sum_counts = pd.concat([df_sum_counts, sum_count_col], axis=1)

                break

    # Save to files if specified
    if save_csv:
        # Adjust score file name based on noCAGEtoNaN
        score_filename = outname + ('_score_setnoCAGEtoNaN.tsv' if noCAGEtoNaN else '_score.tsv')
        df_scores.to_csv(score_filename, sep='\t', index=True)  # Save with row names

        # Save sum_count and position files
        df_sum_counts.to_csv(outname + '_sumcount.tsv', sep='\t', index=True)  # Save with row names
        df_position = df.reset_index()[['chrom', 'chromStart', 'chromEnd']]  # Extract original BED positions
        df_position.to_csv(outname + '_position.bed', sep='\t', index=False, header=False)

    # Print the resulting dataframes (for verification)
    print("Scores DataFrame:")
    print(df_scores.head())
    print("\nSum Counts DataFrame:")
    print(df_sum_counts.head())

if __name__ == "__main__":
    directory = "/Users/natsudanav/Documents/testfile"
    pattern = re.compile(r"FANTOM5-cellfacet-on-PRIMEloci_pred_all_(.*?)_combined\.bed")
    outname = "FANTOM5-cellfacet"

    # Load and filter data with options
    process_prediction_data(directory,
                            pattern,
                            outname,
                            noCAGEtoNaN=True,
                            save_csv=True)