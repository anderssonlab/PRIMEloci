#' Get Highest Non-Overlapping Ranges for a Chromosome
#'
#' This function takes a GenomicRanges::GRanges object
#' and returns the highest non-overlapping ranges
#' based on the scores in the metadata column.
#'
#' @param gr A GenomicRanges::GRanges object with a 'score' metadata column.
#' @return A GenomicRanges::GRanges object with
#' the highest non-overlapping ranges.
#' @import GenomicRanges
#' @import IRanges
#' @export
get_highest_non_overlap_chr <- function(gr) {

  gr <- GenomicRanges::sort(gr, by = ~score, decreasing = TRUE)
  result <- GenomicRanges::GRanges()

  while (length(gr) > 0) {
    highest <- gr[1]
    result <- c(result, highest)
    gr <- gr[!IRanges::overlapsAny(gr, highest)]
  }

  return(result)

}

#' Get Highest Non-Overlapping Ranges for All Chromosomes in Parallel
#'
#' This function takes a `GenomicRanges::GRanges` object and returns the highest 
#' non-overlapping ranges for all chromosomes, processed in parallel.
#'
#' The function splits the input `GRanges` object by chromosome and then 
#' processes each chromosome in parallel to find the highest non-overlapping 
#' ranges based on the 'score' metadata column. The results from each 
#' chromosome are then combined and sorted.
#'
#' @param gr A `GenomicRanges::GRanges` object with a 'score' metadata column.
#' @param num_cores The number of cores to use for parallel processing. 
#' Default is `parallel::detectCores() - 1`.
#'
#' @return A `GenomicRanges::GRanges` object with the highest non-overlapping 
#' ranges for all chromosomes.
#'
#' @import GenomicRanges
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
get_highest_non_overlap <- function(gr,
                                    num_cores = parallel::detectCores() - 1) {

  gr_list <- GenomicRanges::split(gr, GenomicRanges::seqnames(gr))

  doParallel::registerDoParallel(cores = num_cores)

  results <- foreach::foreach(i = seq_along(gr_list), .combine = c) %dopar% {
    get_highest_non_overlap_chr(gr_list[[i]])
  }

  doParallel::stopImplicitCluster()

  results <- GenomicRanges::sort(results, by = ~start)

  return(results)

}


#' Process a BED file and save the highest non-overlapping ranges
#'
#' This function loads a BED file, converts it to a `GRanges` object, and selects the highest
#' non-overlapping ranges using the `get_highest_non_overlap` function. The results are saved
#' as both RDS and BED files.
#'
#' @param input_bed Character. The path to the input BED file.
#' @param output_dir Character. The directory where the output files will be saved.
#'
#' @return A `GRanges` object containing the highest non-overlapping ranges.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom GenomicRanges sort
#' @export
wrapup_filter_bed_to_reduce <- function(input_bed, output_dir = NULL) {
  # Load the bed file and create GRanges object
  bed_file <- load_bed_file(input_bed)
  gr <- create_granges_from_bed(bed_file)

  # Inform the user that processing is starting
  writeLines("Starting to get highest non-overlapping ranges...")
  writeLines("It may take a while depending on the size of the input file...")
  start_time <- Sys.time()
  print(start_time)

  # Get the highest non-overlapping ranges
  selected_gr <- get_highest_non_overlap(gr)
  selected_gr <- GenomicRanges::sort(selected_gr)
  end_time <- Sys.time()
  print(end_time - start_time)

  # If output_dir is specified, save the GRanges object to the directory
  if (!is.null(output_dir) && output_dir != FALSE) {
    input_basename <- tools::file_path_sans_ext(basename(input_bed))
    save_granges_to_bed(selected_gr, output_dir, input_basename, bed_file)
  }

  writeLines("Done!")

  # Return the selected GRanges object
  return(selected_gr)
}