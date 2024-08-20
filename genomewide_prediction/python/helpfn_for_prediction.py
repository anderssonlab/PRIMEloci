def set_index_if_exists(df, index_name='rownames'):
    """
    Sets a specified column as the index of the DataFrame if it exists,
    and removes the index name to prevent it from appearing as an additional
    column in the DataFrame.

    Parameters:
    df (pd.DataFrame): The DataFrame in which the index needs to be set.
    index_name (str): The name of the column to set as the index.

    Returns:
    pd.DataFrame: The DataFrame with the specified column set as the index.
    """
    if index_name in df.columns:
        df = df.set_index(index_name)
        df.index.name = None  # Remove the name of the index to avoid an extra column named 'rownames'
    return df


def adjust_genomic_positions_for_bed(ranges_df):
    """
    Adjusts genomic positions for BED file format.

    Parameters:
    ranges_df (pd.DataFrame): DataFrame with columns 'chrom', 'chromStart', 'chromEnd', 'strand'.

    Returns:
    pd.DataFrame: Adjusted DataFrame ready for BED file output.
    """
    ranges_df['chromStart'] = ranges_df['chromStart'] - 1  # Convert to 0-based start
    return ranges_df


def move_metadata_columns_to_ranges(ranges_df, metadata_df, columns_to_keep):
    """
    Moves columns from the metadata DataFrame to the ranges DataFrame,
    excluding the columns specified in `columns_to_keep`.

    Parameters:
    ranges_df (pd.DataFrame): The DataFrame containing genomic ranges.
    metadata_df (pd.DataFrame): The DataFrame containing additional metadata.
    columns_to_keep (list): List of columns that should not be moved from metadata_df.

    Returns:
    pd.DataFrame: The updated ranges DataFrame with additional metadata columns.
    """
    columns_to_move = [col for col in metadata_df.columns if col not in columns_to_keep]
    for col in columns_to_move:
        ranges_df[col] = metadata_df[col]
    return ranges_df
