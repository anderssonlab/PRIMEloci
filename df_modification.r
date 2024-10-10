# This file are the manipulation of the input DataFrame
#
# Parameters: The input DataFrame to be manipulated.
# Returns: The manipulated DataFrame.
#
# These functions takes a DataFrame as input, performs some manipulation,
# such as adding or removing columns, filtering rows, or computing new
# columns based on existing ones, and returns the manipulated DataFrame.

#' Calculate Midpoint ("Thick") Positions in a DataFrame
#'
#' This function calculates the midpoint positions, referred to as "thick" 
#' positions, for each row in a DataFrame containing "start" and "end" columns. 
#' These columns must be numeric and represent the boundaries of genomic 
#' intervals.
#'
#' @param df A data frame containing genomic intervals with "start" and "end" 
#'   columns.
#' @return A data frame with an additional "thick" column containing the 
#'   calculated midpoint positions.
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @examples
#' df <- data.frame(start = c(100, 200, 300), end = c(150, 250, 350))
#' thick_df <- calculate_thick_df(df)
#' @export
calculate_thick_df <- function(df) {
  # Check if the start and end columns exist and are numeric
  assertthat::assert_that(exists("start", where = df) && exists("end", where = df), # nolint: line_length_linter.
                          msg = "The 'start' and 'end' columns must exist.")
  assertthat::assert_that(is.numeric(df$start) && is.numeric(df$end),
                          msg = "The 'start' and 'end' columns must be numeric.") # nolint: line_length_linter.

  # Calculate the thick position
  thick_df <- df %>%
    dplyr::mutate(thick = .data$start + round((.data$end - .data$start) / 2))

  return(thick_df)
}



#' Extend Intervals from the Center
#'
#' Extends the genomic intervals from their center points, known as "thick"
#' positions, by a specified distance. It can optionally maintain the original
#' interval lengths and group intervals with exactly the same positions.
#'
#' @param df A data frame containing "start" and "end" columns representing
#'   genomic intervals.
#' @param dis A numeric value specifying the distance to extend the intervals
#'   from the center.
#' @param keep_same_length A logical value indicating whether to maintain the
#'   same length for all intervals after extension. Defaults to TRUE.
#' @param group_exact_positions A logical value indicating whether to group
#'   intervals that have the exact same positions. Defaults to TRUE.
#' @return A data frame with updated "start", "end", and "length" columns, and
#'   possibly an updated "thick" column.
#' @importFrom dplyr mutate case_when
#' @importFrom magrittr %>%
#' @examples
#' df <- data.frame(start = c(10, 20, 30), end = c(15, 25, 35))
#' extended_df <- extend_from_center_thick_df(df, dis = 5)
#' extended_df_same_length <- extend_from_center_thick_df(df, dis = 5, keep_same_length = TRUE) # nolint: line_length_linter.
#' @export
extend_from_center_thick_df <- function(df,
                                        dis,
                                        keep_same_length = TRUE,
                                        group_exact_positions = TRUE) {
  # Check if the start and end columns exist and are numeric
  assertthat::assert_that(exists("start", where = df) && exists("end", where = df), # nolint: line_length_linter.
                          msg = "The 'start' and 'end' columns must exist.")
  assertthat::assert_that(is.numeric(df$start) && is.numeric(df$end),
                          msg = "The 'start' and 'end' columns must be numeric.") # nolint: line_length_linter.

  # Check if the thick column exists
  if ("thick" %in% colnames(df)) {
    # Check if the thick column is numeric
    if (!is.numeric(df$thick)) {
      df <- calculate_thick_df(df)
    }
  } else {
    df <- calculate_thick_df(df)
  }

  if (keep_same_length) {
    ext_df <- df %>%
      dplyr::mutate(start = .data$thick - dis,
                    end = .data$thick + dis,
                    start = dplyr::case_when(start < 1 ~ 1,
                                             TRUE ~ start),
                    end = dplyr::case_when(start == 1 ~ (dis * 2 + 1),
                                           TRUE ~ end),
                    length = end - start + 1)
    # Assert that all intervals have the same length
    assertthat::assert_that(length(unique(ext_df$length)) == 1,
                            msg = "The length of the intervals is not the same.") # nolint: line_length_linter.
    # Recalculate the thick position to match the adjusted end position
    ext_df <- calculate_thick_df(ext_df)
  } else {
    ext_df <- df %>%
      dplyr::mutate(start = .data$thick - dis,
                    end = .data$thick + dis,
                    start = dplyr::case_when(.data$start < 1 ~ 1,
                                             TRUE ~ start),
                    length = end - start + 1)
    print("Length was recalculated to match the new start/end")
  }

  if (group_exact_positions) {
    ext_df <- group_exact_positions_df(ext_df)
    print("Grouped by exact positions")
  }

  return(ext_df)
}



#' Filter Genomic Intervals by Specified Chromosomes
#'
#' This function filters genomic intervals in a DataFrame to retain only those 
#' intervals located on specified chromosomes. The DataFrame must include a 
#' "seqnames" column indicating the chromosome names.
#'
#' @param df A data frame containing genomic intervals, with a "seqnames" 
#'   column specifying the chromosome names.
#' @param given_chromosomes A character vector of chromosome names to retain 
#'   in the DataFrame.
#' @return A data frame containing only the genomic intervals located on the 
#'   specified chromosomes.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @examples
#' df <- data.frame(seqnames = c("chr1", "chr2", "chrX", "chrY"),
#'                  start = c(100, 200, 300, 400),
#'                  end = c(150, 250, 350, 450))
#' filtered_df <- keep_given_chromosomes_df(df, c("chr1", "chr2"))
#' @export
keep_given_chromosomes_df <- function(df, given_chromosomes) {
  assertthat::assert_that(exists("seqnames", where = df),
                          msg = "The 'seqnames' column must exist in the DataFrame.") # nolint: line_length_linter.
  return(df %>% dplyr::filter(.data$seqnames %in% given_chromosomes))
}



#' Group Genomic Intervals by Exact Positions
#'
#' This function groups genomic intervals in a DataFrame that have identical
#' "seqnames", "start", and "end" positions. Optionally, the "thick" column can
#' also be considered if present. For each group of intervals with the same
#' positions, the function concatenates unique values from other columns, using
#' a semicolon as the separator.
#'
#' @param df A data frame containing genomic intervals with "seqnames",
#'   "start", and "end" columns.
#' @return A data frame with intervals grouped by identical positions, with
#'   concatenated unique values from other columns.
#' @importFrom dplyr group_by summarise_all across
#' @importFrom magrittr %>%
#' @examples
#' df <- data.frame(seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'                  start = c(100, 100, 200, 200),
#'                  end = c(150, 150, 250, 250),
#'                  value = c("A", "B", "C", "D"))
#' grouped_df <- group_exact_positions_df(df)
#' @export
group_exact_positions_df <- function(df) {
  assertthat::assert_that(exists("seqnames", where = df),
                          msg = "The 'seqnames' column must exist in the DataFrame.")
  assertthat::assert_that(exists("start", where = df) && exists("end", where = df),
                          msg = "The 'start' and 'end' columns must exist in the DataFrame.")
  assertthat::assert_that(is.numeric(df$start) && is.numeric(df$end),
                          msg = "The 'start' and 'end' columns must be numeric.")

  # Group by identical positions
  group_cols <- c("seqnames", "start", "end")
  if ("thick" %in% colnames(df)) {
    assertthat::assert_that(is.numeric(df$thick), msg = "The 'thick' column must be numeric.")
    group_cols <- c(group_cols, "thick")
  }

  group_df <- df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise_all(~ paste(unique(.), collapse = ";"))

  return(group_df)
}