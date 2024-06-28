# This file are the manipulation of the input DataFrame
#
# Parameters: The input DataFrame to be manipulated.
# Returns: The manipulated DataFrame.
#
# These functions takes a DataFrame as input, performs some manipulation,
# such as adding or removing columns, filtering rows, or computing new
# columns based on existing ones, and returns the manipulated DataFrame.


suppressPackageStartupMessages({
  library(assertthat)
  library(magrittr)
  library(dplyr)
  library(rlang)
})


#' Calculate the "thick" position in a dataframe
#'
#' This function calculates the "thick" position in a dataframe (`df`)
#' based on the midpoint of the "start" and "end" positions.
#' The dataframe must contain columns named "start" and "end",
#' and these columns must be numeric.
#'
#' @param df Dataframe containing genomic intervals.
#' @return Dataframe with an additional column named "thick"
#' containing the calculated "thick" positions.
#'
#' @examples
#' df <- data.frame(start = c(100, 200, 300), end = c(150, 250, 350))
#' thick_df <- calculate_thick_df(df)
#'
#' #' @export
calculate_thick_df <- function(df) {
  # check if the start and end columns exist and are numeric
  assert_that(exists("start", where = df) && exists("end", where = df),
              msg = "The 'start' and 'end' columns must exist.")
  assert_that(is.numeric(df$start) && is.numeric(df$end),
              msg = "The 'start' and 'end' columns must be numeric.")

  # calculate the thick position
  # Add the missing import for the `%>%` operator
  thick_df <- df %>%
    dplyr::mutate(thick = .data$start + round((.data$end - .data$start) / 2))

  return(thick_df)
}

#' Extend genomic intervals from a central "thick" position in a dataframe
#'
#' This function extends the "start" and "end" positions of genomic intervals
#' in a dataframe (`df`) based on a central "thick" position.
#' The dataframe must contain columns named "start" and "end",
#' and these columns must be numeric.
#' If a column named "thick" does not exist or is not numeric,
#' it calculates the "thick" column using the calculate_thick_df function.
#' Then, it extends the intervals symmetrically around the "thick" position
#' by a specified distance (`dis`).
#' If the resulting "start" positions become negative, they are set to 0.
#'
#' @param df Dataframe containing genomic intervals.
#' @param dis Distance by which to extend the intervals
#' from the "thick" position.
#' @return Extended dataframe with modified start and end positions.
#'
#' @examples
#' df <- data.frame(start = c(100, 200, 300), end = c(150, 250, 350))
#' extended_df <- extend_from_center_thick_df(df, dis = 20)
#'
#' @export
extend_from_center_thick_df <- function(df, dis) {
  # check if the start and end columns exist and are numeric
  assert_that(exists("start", where = df) && exists("end", where = df),
              msg = "The 'start' and 'end' columns must exist.")
  assert_that(is.numeric(df$start) && is.numeric(df$end),
              msg = "The 'start' and 'end' columns must be numeric.")

  # check if the thick column exists
  if ("thick" %in% colnames(df)) {
    # check if the thick column is numeric
    if (!is.numeric(df$thick)) {
      df <- calculate_thick_df(df)
    }
  } else {
    df <- calculate_thick_df(df)
  }

  # extend the start and end positions
  ext_df <- df %>%
    dplyr::mutate(start = .data$thick - dis,
                  end = .data$thick + dis,
                  start = dplyr::case_when(.data$start < 0 ~ 0,
                                           TRUE ~ .data$start))

  return(ext_df)
}



#' Filter genomic intervals by given chromosomes in a dataframe
#'
#' This function filters genomic intervals in a dataframe (`df`)
#' to retain only those intervals located on chromosomes
#' specified by the user.
#' The dataframe need to contain a column named "seqnames"
#' that specifies the chromosome names.
#' It filters the dataframe to keep only the intervals
#' located on the chromosomes specified in the `given_chromosomes` vector.
#'
#' @param df Dataframe containing genomic intervals.
#' @param given_chromosomes A vector specifying
#' the chromosome names to keep.
#' @return Dataframe containing only the genomic intervals
#' located on the specified chromosomes.
#'
#' @examples
#' df <- data.frame(seqnames = c("chr1", "chr2", "chrX", "chrY"),
#'                  start = c(100, 200, 300, 400),
#'                  end = c(150, 250, 350, 450))
#' given_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5",
#'                           "chr6", "chr7", "chr8", "chr9", "chr10")
#' filtered_df <- keep_given_chromosomes_df(df, given_chromosomes)
#'
#' @export
keep_given_chromosomes_df <- function(df, given_chromosomes) {
  assert_that(exists("seqnames", where = df),
              msg = "The 'seqnames' column must exist in the DataFrame.")
  return(df %>% dplyr::filter(.data$seqnames %in% given_chromosomes))
}



#' Group genomic intervals with exactly the same positions in a dataframe
#'
#' This function groups genomic intervals in a dataframe (`df`)
#' that have exactly the same positions based on
#' the "seqnames" and "thick" columns.
#' The dataframe must contain a column named "seqnames"
#' that specifies the chromosome names and a column named "thick"
#' that specifies the thick positions.
#' If no thick column exists or is not numeric,
#' it calculates the "thick" column using the calculate_thick_df function.
#' It groups the intervals with the same positions
#' and concatenates the unique values of all other columns
#' using a semicolon as the separator.
#'
#' @param df Dataframe containing genomic intervals.
#' @return Dataframe with grouped intervals
#' that have exactly the same positions.
#'
#' @examples
#' df <- data.frame(seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'                  thick = c(100, 100, 200, 200),
#'                  value = c("A", "B", "C", "D"))
#' grouped_df <- group_exact_positions_df(df)
#'
#' @export
group_exact_positions_df <- function(df) {
  assert_that(exists("seqnames", where = df),
              msg = "The 'seqnames' column must exist in the DataFrame.")

  # check if the thick column exists
  if ("thick" %in% colnames(df)) {
    # check if the thick column is numeric
    if (!is.numeric(df$thick)) {
      df <- calculate_thick_df(df)
    }
  } else {
    df <- calculate_thick_df(df)
  }

  # group the 'exactly' same positions
  group_df <- df %>%
    dplyr::group_by(.data$seqnames, .data$thick) %>%
    dplyr::summarise_all(function(x) paste(unique(x), collapse = ";"))

  return(group_df)

}
