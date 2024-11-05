#!/usr/bin/env Rscript

writeLines("\n### Running NEW POSTPROCESS ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(GenomicRanges)
  library(parallel)
  library(tools)
  library(stringr)
  library(PRIMEloci)
})

# Create the argparse object
parser <- ArgumentParser(description = "Process GRanges with selective core merging based on score and width thresholds.") # nolint: line_length_linter.
parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-d", "--score_diff", type = "double", default = 10,
                    help = "Score difference threshold for merging.")
parser$add_argument("-w", "--max_width", type = "integer", default = 500,
                    help = "Maximum width for merged regions.")
parser$add_argument("-c", "--use_max_cores", action = "store_true",
                    default = TRUE,
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")

# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- args$score_threshold
score_diff <- args$score_diff
max_width <- args$max_width
bed_file <- args$input_file
use_max_cores <- args$use_max_cores
output_dir <- args$output_dir

# Set the number of cores
num_cores <- if (use_max_cores) { parallel::detectCores() } else { 1 }

# Load and prepare data
bed <- load_bed_file(bed_file)
gr <- create_granges_from_bed(bed)

# Filter GRanges by score threshold
filtered_gr <- gr[gr$score >= score_threshold]

# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))

# Function to selectively merge cores based on score difference and max width
selective_merge_cores <- function(overlap_set, score_diff, max_width) {
  # Rank the cores by score (highest first)
  overlap_set <- overlap_set[order(-overlap_set$score)]
  merged_cores <- GenomicRanges::GRanges()  # Initialize an empty GRanges object
  
  while (length(overlap_set) > 0) {
    # Take the highest-ranked core
    x <- overlap_set[1]
    score_x <- x$score
    
    # Identify cores to merge based on score difference
    merge_candidates <- overlap_set[overlap_set$score >= score_x - score_diff]
    
    # Merge only if the resulting width is within the max_width
    merged_region <- GenomicRanges::reduce(merge_candidates)
    
    if (width(merged_region) <= max_width) {
      # Assign the score of x to the merged region
      mcols(merged_region)$score <- score_x
      
      # Add the merged core to the final list
      merged_cores <- c(merged_cores, merged_region)
    }
    
    # Remove merged candidates from the overlap set
    overlap_set <- overlap_set[!(seqnames(overlap_set) %in% seqnames(merge_candidates) &
                                 start(overlap_set) %in% start(merge_candidates) &
                                 end(overlap_set) %in% end(merge_candidates))]
  }
  
  return(merged_cores)
}

# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list, function(chr) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]
    core_gr <- extract_core(chr_gr)
    
    overlaps <- GenomicRanges::findOverlaps(core_gr)
    collapsed_ranges <- GenomicRanges::reduce(core_gr[unique(S4Vectors::queryHits(overlaps))])
    
    # Apply selective merging on overlapping sets of cores
    final_merged_gr <- selective_merge_cores(collapsed_ranges, score_diff, max_width)
    
    # Calculate and add metadata for the final merged cores
    thick_vals <- numeric(length(final_merged_gr))
    max_scores <- numeric(length(final_merged_gr))
    all_scores <- character(length(final_merged_gr))
    
    for (i in seq_along(final_merged_gr)) {
      metadata <- get_metadata(final_merged_gr[i], filtered_gr)
      thick_vals[i] <- metadata$thick
      max_scores[i] <- metadata$max_score
      all_scores[i] <- metadata$all_scores
    }
    
    # Add metadata to the final merged cores
    final_merged_gr$thick <- thick_vals
    final_merged_gr$max_score <- max_scores
    final_merged_gr$all_scores <- all_scores
    
    return(final_merged_gr)
    
  }, error = function(e) {
    message(paste("Error processing chromosome", chr, ":", e$message))
    return(NULL)
  })
}, mc.cores = num_cores)

# Combine all results
collapsed_gr <- do.call(c, collapsed_gr_list)
collapsed_gr <- sort(collapsed_gr)

# If output_dir is specified, save the GRanges object to the directory
if (!is.null(output_dir) && output_dir != FALSE) {
  input_basename <- tools::file_path_sans_ext(basename(bed_file))
  input_basename <- input_basename %>%
    stringr::str_replace_all("all", as.character(score_threshold)) %>%
    stringr::str_replace_all("[^[:alnum:]]", "_")

  save_ovlcorereduced_to_bed(collapsed_gr, output_dir, input_basename)
}

writeLines("Done!")
