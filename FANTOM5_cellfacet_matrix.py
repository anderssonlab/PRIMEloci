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


def set_score_to_nan_where_sum_count_zero(df):
    """
    Sets 'score' to NaN where 'sum_count' is 0 in the dataframe.

    Parameters:
        df (pd.DataFrame): The dataframe to modify.

    Returns:
        pd.DataFrame: The modified dataframe.
    """
    print("Filtering: Setting 'score' to NaN where 'sum_count' = 0.")
    df['score'] = df.apply(
        lambda row: np.nan if row['sum_count'] == 0 else row['score'],
        axis=1
    )
    return df


def parse_row_name_to_bed(row_names):
    """
    Parses row names in the format chrom:chromStart-chromEnd;strand to extract BED fields.

    Parameters:
        row_names (pd.Index): Row names in the specified format.

    Returns:
        pd.DataFrame: DataFrame with columns [chrom, chromStart, chromEnd, strand].
    """
    parsed = row_names.str.extract(r"(?P<chrom>[^:]+):(?P<chromStart>\d+)-(?P<chromEnd>\d+)(?:;(?P<strand>.+))?")
    parsed['chromStart'] = parsed['chromStart'].astype(int) - 1  # Convert to 0-based start
    parsed['chromEnd'] = parsed['chromEnd'].astype(int)
    return parsed


def process_prediction_data(directory, pattern, outname, noCAGEtoNaN=True, save_tsv=True):
    """
    Processes BED files, applies filtering, and splits data into 'score' and 'sum_count' dataframes.
    Saves results to files. Row names are formatted as chrom:chromStart+1-chromEnd;strand (R-format).

    Parameters:
        directory (str): The directory containing the files.
        pattern (re.Pattern): The regex pattern to extract sample names from filenames.
        outname (str): Output filename prefix for saving results.
        noCAGEtoNaN (bool): Whether to set 'score' to NaN where 'sum_count' = 0.
        save_tsv (bool): Whether to save the resulting dataframes to TSV files.
    """
    # Extract sample names from the directory
    samples = extract_exclude_patterns(directory, pattern)

    # Initialize dataframes for scores and sum_counts
    score_frames = []
    sum_count_frames = []

    # Process each sample
    for sample in samples:
        print(sample)
        filename_pattern = f"FANTOM5-cellfacet-on-PRIMEloci_pred_all_{sample}_combined.bed"
        file_path = os.path.join(directory, filename_pattern)

        if not os.path.isfile(file_path):
            print(f"File not found for sample: {sample}. Skipping.")
            continue

        # Read the BED file into a dataframe
        df = pd.read_csv(file_path, sep='\t', header=0)

        # Create row names in the format chrom:chromStart+1-chromEnd;strand
        df['row_name'] = df.apply(
            lambda row: f"{row['chrom']}:{row['chromStart'] + 1}-{row['chromEnd']};{row['strand']}"
            if 'strand' in df.columns else f"{row['chrom']}:{row['chromStart'] + 1}-{row['chromEnd']}",
            axis=1
        )
        df.set_index('row_name', inplace=True)

        # Apply noCAGEtoNaN filtering to scores if specified
        if noCAGEtoNaN:
            df = set_score_to_nan_where_sum_count_zero(df)

        # Extract and collect score and sum_count columns
        score_frames.append(df['score'].rename(sample))
        sum_count_frames.append(df['sum_count'].rename(sample))

    # Combine all scores and sum_counts into respective dataframes
    df_scores = pd.concat(score_frames, axis=1)
    df_sum_counts = pd.concat(sum_count_frames, axis=1)

    # Check if row names in both DataFrames match
    if not df_scores.index.equals(df_sum_counts.index):
        raise ValueError("Row names in scores and sum counts do not match. Ensure consistent data alignment.")

    # Generate position BED file from row names
    position_df = parse_row_name_to_bed(df_scores.index)
    
    # Save to files if specified
    if save_tsv:
        score_filename = outname + ('_score_setnoCAGEtoNaN.tsv' if noCAGEtoNaN else '_score.tsv')
        df_scores.to_csv(score_filename, sep='\t', index=True)

        sum_count_filename = outname + '_sumcount.tsv'
        df_sum_counts.to_csv(sum_count_filename, sep='\t', index=True)

        position_filename = outname + '_position.bed'
        position_df.to_csv(position_filename, sep='\t', index=False, header=False)

    # Print summary of results
    print("Scores DataFrame:")
    print(df_scores.head())
    print("\nSum Counts DataFrame:")
    print(df_sum_counts.head())
    print("\nPosition BED DataFrame:")
    print(position_df.head())


if __name__ == "__main__":
    # Define parameters
    directory = "/home/zmk214/data/projects/nucleiCAGEproject/8.Genomewide_prediction/FANTOM5_rmSingletons_cellfacet/prediction"
    pattern = re.compile(r"FANTOM5-cellfacet-on-PRIMEloci_pred_all_(.*?)_combined\.bed")
    outname = "FANTOM5-cellfacet"

    # Process data
    process_prediction_data(directory, pattern, outname, noCAGEtoNaN=True)

