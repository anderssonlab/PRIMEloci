#!/usr/bin/env Rscript

writeLines("\n### Running _filter_bed_to_reduce_gr.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(assertthat)
  library(data.table)
  library(argparse)
  library(foreach)
  library(doParallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(PRIMEloci)
})

# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input_bed", type = "character",
                    help = "Input BED file")
parser$add_argument("-o", "--output_dir", type = "character",
                    default = NULL,
                    help = "Output directory, default will be the same as the input file") # nolint: line_length_linter.
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")

args <- parser$parse_args()


# Execute the main function
wrapup_filter_bed_to_reduce_core(args$input_bed,
                                 args$output_dir,
                                 args$score_threshold,
                                 width = 401,
                                 ext_core = 75)

writeLines("Done!")



input_bed = "/Users/natsudanav/Desktop/PRIMEloci/example/results/K562-on-PRIMEloci1.0-model_pred_all_profiles_subtnorm_tcs_K562_C.bed"
output_dir = "/Users/natsudanav/Desktop/PRIMEloci/example/results"
threshold = 0.7
width = 401
ext_core = 75
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