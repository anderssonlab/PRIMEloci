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
#' @importFrom stringr str_replace_all
#' @importFrom dplyr %>%
#' @export
wrapup_filter_bed_to_reduce <- function(input_bed,
                                        output_dir = NULL,
                                        threshold) {
  # Load the bed file and create GRanges object
  bed_file <- load_bed_file(input_bed)
  gr <- create_granges_from_bed(bed_file)

  filtered_gr <- gr[gr$score >= threshold]

  # Inform the user that processing is starting
  writeLines("Starting to get highest non-overlapping ranges...")
  writeLines("It may take a while depending on the size of the input file...")
  start_time <- Sys.time()
  print(start_time)

  # Get the highest non-overlapping ranges
  selected_gr <- get_highest_non_overlap(filtered_gr)
  selected_gr <- GenomicRanges::sort(selected_gr)
  end_time <- Sys.time()
  print(end_time - start_time)

  # If output_dir is specified, save the GRanges object to the directory
  if (!is.null(output_dir) && output_dir != FALSE) {
    input_basename <- tools::file_path_sans_ext(basename(input_bed))
    input_basename <- input_basename %>%
      stringr::str_replace_all("all", as.character(threshold)) %>%  # Replace "all" with threshold # nolint: line_length_linter.
      stringr::str_replace_all("[^[:alnum:]]", "_")                # Replace non-alphanumeric characters with "_" # nolint: line_length_linter.

    save_granges_to_bed(selected_gr, output_dir, input_basename, bed_file)
  }

  writeLines("Done!")

  # Return the selected GRanges object
  return(selected_gr)
}


#' Process a BED file, filter by score, extract core regions, extend, and save the reduced, non-overlapping ranges
#'
#' This function processes a BED file by first converting it to a `GRanges` object, filtering it based on a
#' specified score threshold, and extracting core regions using the `extract_core` function. The core regions 
#' are extended left and right to a fixed width of 401 base pairs, and overlapping or adjacent regions are 
#' reduced using the `reduce()` function. The function saves two outputs: 
#' 1) the selected non-overlapping core regions with the suffix "_reduce_core"
#' 2) the extended and reduced ranges with the suffix "_extended_reduced". 
#' The results are saved as both RDS and BED files in the specified output directory.
#'
#' @param input_bed Character. The path to the input BED file.
#' @param output_dir Character. The directory where the output files will be saved.
#' @param threshold Numeric. The minimum score threshold to filter ranges.
#' @param ext_core Numeric. The extension to apply to core regions (default is 75 base pairs).
#'
#' @return A `GRanges` object containing the reduced, extended ranges.
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom GenomicRanges sort reduce resize
#' @importFrom stringr str_replace_all
#' @importFrom dplyr %>%
#' @export
wrapup_filter_bed_to_reduce_core <- function(input_bed,
                                             output_dir = NULL,
                                             threshold,
                                             width = 401,
                                             ext_core = 75) {
  # Load the bed file and create GRanges object
  bed_file <- load_bed_file(input_bed)
  gr <- create_granges_from_bed(bed_file)

  # Filter by score
  filtered_gr <- gr[gr$score >= threshold]

  # Extract core regions
  #core_gr <- extract_core(filtered_gr, ext_core = ext_core)
  core_gr <- GenomicRanges::resize(filtered_gr, width = (ext_core*2+1), fix = "center")

  # Inform the user that processing is starting
  writeLines("Starting to get highest non-overlapping ranges -- core ...")
  writeLines("It may take a while depending on the size of the input file...")
  start_time <- Sys.time()
  print(start_time)

  # Get the highest non-overlapping ranges
  selected_gr <- get_highest_non_overlap(core_gr)
  selected_gr <- GenomicRanges::sort(selected_gr)
  end_time <- Sys.time()
  print(end_time - start_time)

  # Save the selected GRanges object with name _reduce_core
  if (!is.null(output_dir) && output_dir != FALSE) {
    input_basename <- tools::file_path_sans_ext(basename(input_bed)) %>%
      stringr::str_replace_all("all", as.character(threshold)) %>%
      stringr::str_replace_all("[^[:alnum:]]", "_") %>%
      paste0("_reduce_core")

    save_granges_to_bed(selected_gr, output_dir, input_basename, bed_file)
  }

  # Extend the core regions left/right to a fixed width of 401 bp
  extended_gr <- GenomicRanges::resize(core_gr, width = width, fix = "center")
  extended_gr
  # Use the reduce() function on the extended GRanges to merge overlapping or adjacent ranges
  reduced_gr <- GenomicRanges::reduce(extended_gr)
  reduced_gr <- GenomicRanges::sort(reduced_gr)
  reduced_gr
  # Save the reduced, extended GRanges
  if (!is.null(output_dir) && output_dir != FALSE) {
    extended_basename <- tools::file_path_sans_ext(basename(input_bed)) %>%
      stringr::str_replace_all("all", as.character(threshold)) %>%
      stringr::str_replace_all("[^[:alnum:]]", "_") %>%
      paste0("_extended_reduced")
    bed_file_path <- file.path(output_dir, paste0(extended_basename, ".bed"))
    rtracklayer::export.bed(reduced_gr, bed_file_path)
    rds_file_path <- file.path(output_dir, paste0(extended_basename, ".rds"))
    saveRDS(reduced_gr, rds_file_path)
  }

  writeLines("Done!")

  # Return the reduced GRanges object
  return(reduced_gr)
}