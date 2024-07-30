import os
import pandas as pd
import numpy as np


def set_score_to_nan_where_sum_count_zero(df):
    """
    Sets 'score' to NaN where 'sum_count' is 0 in the dataframe.
    """
    print("Filtering: Setting 'score' to NaN where 'sum_count' = 0.")
    df['score'] = df.apply(
        lambda row: np.nan if row['sum_count'] == 0 else row['score'],
        axis=1
    )
    return df


def load_and_filter_prediction_data(directory,
                                    models,
                                    samples,
                                    pattern,
                                    y_true):
    """
    Loads data files for specified models and samples, applies filtering,
    adds a uniform label, and outputs dataframes.

    Parameters:
        directory (str): The directory containing the files.
        models (list): List of models to process.
        samples (list): List of sample names in the desired order.
        pattern (str): The pattern to match filenames.
        label (int): A uniform label to apply across all rows.

    Returns:
        dict: Dataframes for scores and sum_counts for each model.
    """
    # Initialize dictionaries to store dataframes for scores and sum_counts
    # for each model
    df_scores = {model: pd.DataFrame() for model in models}
    df_sum_counts = {model: pd.DataFrame() for model in models}

    # Loop through each model and sample name
    for model in models:
        for sample in samples:
            # Construct the expected filename pattern
            filename_pattern = pattern.format(model=model, sample=sample)

            # Look for the file in the directory that matches the pattern
            for filename in os.listdir(directory):
                if filename_pattern in filename:
                    # Construct full file path
                    file_path = os.path.join(directory, filename)

                    # Read the CSV file into a dataframe
                    df = pd.read_csv(file_path, sep='\t', header=0)

                    # Filter scores based on sum_count
                    df = set_score_to_nan_where_sum_count_zero(df)

                    # Append the score and sum_count to
                    # the respective dataframes
                    df_scores[model][sample] = df['score']
                    df_sum_counts[model][sample] = df['sum_count']
                    break

    # Add uniform label and reset index to include chrom, chromStart,
    # and chromEnd for each model's dataframe
    for model in models:
        if not df_scores[model].empty:
            df_scores[model] = df[['chrom',
                                   'chromStart',
                                   'chromEnd']].join(df_scores[model])
            df_sum_counts[model] = df[['chrom',
                                       'chromStart',
                                       'chromEnd']].join(df_sum_counts[model])

            # Add uniform label column
            df_scores[model]['y_true'] = y_true
            df_sum_counts[model]['y_true'] = y_true

            # Save the resulting dataframes to CSV files (optional)
            df_scores[model].to_csv(f'scores_{model}.csv', index=False)
            df_sum_counts[model].to_csv(f'sum_counts_{model}.csv', index=False)

            # Print the resulting dataframes (for verification)
            print(f"Scores {model} DataFrame:")
            print(df_scores[model].head())
            print(f"\nSum Counts {model} DataFrame:")
            print(df_sum_counts[model].head())

    return df_scores, df_sum_counts
