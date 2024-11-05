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

library(GenomicRanges)
library(IRanges)
library(S4Vectors)

# Create the argparse object
parser <- ArgumentParser(description = "Process GRanges with selective core merging based on score and width thresholds.") # nolint: line_length_linter.
parser$add_argument("-i", "--input_file", required = TRUE,
                    help = "Path to the input BED file.")
parser$add_argument("-t", "--score_threshold", type = "double", default = 0.7,
                    help = "Score threshold for filtering GRanges.")
parser$add_argument("-d", "--score_diff", type = "double", default = 0.1,
                    help = "Score difference threshold for merging.")
parser$add_argument("-w", "--max_width", type = "integer", default = 1000,
                    help = "Maximum width for merged regions.")
parser$add_argument("-m", "--use_max_cores", action = "store_true",
                    default = TRUE,
                    help = "Flag to use the maximum number of available cores.")
parser$add_argument("-o", "--output_dir", default = NULL,
                    help = "Directory to save the output files.")

# Parse arguments
args <- parser$parse_args()

# Define parameters based on parsed arguments
score_threshold <- as.numeric(args$score_threshold)
score_diff <- as.numeric(args$score_diff)
max_width <- as.numeric(args$max_width)
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
print(head(filtered_gr))
# Get a list of unique chromosomes
chr_list <- unique(as.character(seqnames(filtered_gr)))
print(chr_list)

selective_merge_cores <- function(core_gr, score_diff, max_width) {
  # Rank the cores by score (highest first)
  core_gr <- core_gr[order(-core_gr$score)]
  
  # Initialize an empty GRanges object to store the final merged cores
  merged_cores <- GenomicRanges::GRanges()
  
  # Initialize containers for storing metadata
  thick_vals <- IRanges::IRanges()
  max_scores <- numeric(0)
  
  while (length(core_gr) > 0) {
    # Take the highest-ranked core
    x <- core_gr[1]
    
    # Calculate thick as the midpoint of the range of x
    thick_x <- IRanges::IRanges(start = start(x) + floor(width(x) / 2),
                                width = 1)  # Create an IRanges object for thick
    
    score_x <- x$score  # Use the existing score
    
    # Find overlapping cores with x
    overlaps <- GenomicRanges::findOverlaps(x, core_gr)
    overlap_set <- core_gr[subjectHits(overlaps)]
    
    # Filter overlap set based on score difference
    merge_candidates <- overlap_set[overlap_set$score >= score_x - score_diff]
    
    # Attempt to merge candidates if there is more than one
    if (length(merge_candidates) > 1) {
      merged_region <- GenomicRanges::reduce(merge_candidates)
      
      # Append metadata
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      
      # Add the merged core to the final result
      merged_cores <- c(merged_cores, merged_region)
      
      # Remove the merged cores from core_gr
      core_gr <- core_gr[!(seqnames(core_gr) %in% seqnames(merge_candidates) &
                           start(core_gr) %in% start(merge_candidates) &
                           end(core_gr) %in% end(merge_candidates))]
    } else {
      # If no valid merge candidates, just add the highest core x
      thick_vals <- c(thick_vals, thick_x)
      max_scores <- c(max_scores, score_x)
      
      merged_cores <- c(merged_cores, x)
      core_gr <- core_gr[-1]  # Remove only the highest core x
    }
  }
  
  # Add metadata to the final merged cores
  mcols(merged_cores)$thick <- thick_vals
  mcols(merged_cores)$max_score <- max_scores
  
  return(merged_cores)
}




# Run in parallel across chromosomes using the specified number of cores
collapsed_gr_list <- mclapply(chr_list, function(chr) {
  tryCatch({
    # Subset GRanges by chromosome
    chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]
    core_gr <- GenomicRanges::resize(chr_gr, width = 151, fix = "center")
    print(head(core_gr))
    # Apply selective merging on overlapping sets of cores
    final_merged_gr <- selective_merge_cores(core_gr, score_diff, max_width)
    print(head(final_merged_gr))

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