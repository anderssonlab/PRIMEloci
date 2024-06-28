suppressPackageStartupMessages({
  library(tidyverse)
})

#' Measure and Report Execution Time of a Function
#'
#' The `report_time_execution` function measures 
#' and reports the execution time of a given function.
#'
#' @param fun A function call to be executed and timed.
#' @return The output of the executed function.
#' @export
report_time_execution <- function(fun) {
  start_time <- Sys.time()          # Capture the start time
  output <- fun                     # Execute the function and store the output
  print(Sys.time() - start_time)    # Calculate and print the elapsed time
  return(output)                    # Return the function's output
}



#' Convert SummarizedExperiment to GRanges with Assay Data
#'
#' This function converts a `SummarizedExperiment` object to a `GRanges` object
#' and adds a specified assay data column to the resulting `GRanges` object.
#'
#' @param rse A `SummarizedExperiment` object
#' containing the genomic ranges and assay data.
#' @param assay A character string specifying the assay name to extract from
#' the `SummarizedExperiment` object. Default is "counts".
#' @param coln_assay An integer specifying the column index of the assay to use.
#' Default is 1.
#' @param colname A character string specifying the name of the column
#' to be added to the `GRanges` object. Default is "score".
#'
#' @return A `GRanges` object with the specified assay data
#' added as a metadata column.
#'
#' @importFrom SummarizedExperiment rowRanges assay
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom tidyverse %>%
#' @export
cast_rse_to_granges <- function(rse,
                                assay = "counts",
                                coln_assay = 1,
                                colname = "score") {
  # Extract the row ranges from the SummarizedExperiment object
  # Ensure it's a GRanges object
  gr <- SummarizedExperiment::rowRanges(rse) %>% GenomicRanges::GRanges() # nolint: pipe_operator_linter

  # Extract the assay data
  assay_data <- SummarizedExperiment::assay(rse, assay)

  # Assign the assay data to the specified column name in the GRanges object
  GenomicRanges::mcols(gr)[[colname]] <- assay_data[, coln_assay]

  return(gr)
}