#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(GenomicRanges)
  library(parallel)
  library(argparse)
  library(PRIME)
}))

# Create argument parser
parser <- ArgumentParser(description = "Sliding window on genomic ranges")
parser$add_argument("-i", "--infile", required = TRUE,
                    help = "Path to the input .rds file (tc_grl)")

parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "sld_tc_grl.rds",
                    help = "Output file name for sld_tc_grl object")

parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")
parser$add_argument("-s", "--sld_by", type = "integer", default = 20,
                    help = "Slide by parameter")
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")


# Parse arguments
args <- parser$parse_args()

# Load the input data
tc_grl <- readRDS(args$infile)

# Output
output_dir <- args$output_dir
PRIME::plc_create_output_dir(output_dir)
outfile_sld_tc_grl <- args$outfile

primeloci_tmp <- PRIME::plc_setup_tmp_dir(output_dir)

# parameters
sld_by <- as.integer(args$sld_by)
num_cores <- args$num_cores
ext_dis <- as.integer(args$ext_dis)

assertthat::assert_that(
  is.numeric(sld_by),
  sld_by %% 1 == 0,
  sld_by > 0,
  msg = "`sld_by` must be a positive integer."
)

assertthat::assert_that(
  is.numeric(ext_dis),
  ext_dis %% 1 == 0,
  ext_dis > 0,
  msg = "`ext_dis` must be a positive integer."
)

if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "❌ `num_cores` must be a positive integer or NULL."
  )
}

if (is.null(num_cores)) {
  num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
}
if (num_cores == 1) {
  processing_method <- "callr"
  plc_message("⚠️ num_workers was set to 1. Using callr backend: tasks will run sequentially (despite using multiple R sessions).") # nolint: line_length_linter.
} else {
  processing_method <- PRIME::plc_detect_parallel_plan()
}


plc_message("🚀 Running PRIMEloci -3: sliding windows covering reduced TC regions") # nolint: line_length_linter.

if (inherits(tc_grl, "GenomicRanges::GRanges")) {
  start_time <- Sys.time()
  tc_sliding_window_grl <- PRIME::plc_tc_sliding_window(tc_grl,
                                                        sld_by = sld_by,
                                                        ext_dis = ext_dis,
                                                        num_cores = num_cores,
                                                        processing_method = processing_method) # nolint: line_length_linter.
  plc_message(sprintf("⏱️ Time taken: %.2f minutes",
                      as.numeric(difftime(Sys.time(),
                                          start_time,
                                          units = "mins"))))
} else if (inherits(tc_grl, "GenomicRanges::GRangesList") ||
             inherits(tc_grl, "CompressedGRangesList")) {
  tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
    start_time <- Sys.time()
    gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
    plc_message(sprintf("🔹 Processing: %s", gr_name))
    result <- plc_tc_sliding_window(tc_grl[[i]],
                                    sld_by = sld_by,
                                    ext_dis = ext_dis,
                                    num_cores = num_cores,
                                    processing_method = processing_method)
    plc_message(sprintf("⏱️ Time taken: %.2f minutes",
                        as.numeric(difftime(Sys.time(),
                                            start_time,
                                            units = "mins"))))
    result
  })
  if (!is.null(tc_sliding_window_grl) &&
        length(tc_sliding_window_grl) > 0) {
    tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
    names(tc_sliding_window_grl) <- names(tc_grl)
  } else {
    plc_error("❌ Processed TC object list is empty. Ensure tc_grl contains valid data.") # nolint: line_length_linter.
  }
} else {
  plc_error("❌ tc_grl must be either a GRanges, GRangesList, or CompressedGRangesList object.") # nolint: line_length_linter.
}

plc_message(sprintf("Saving TC objects to PRIMEloci_tmp .."))
saveRDS(tc_sliding_window_grl,
        file.path(primeloci_tmp, outfile_sld_tc_grl))

plc_message("✅ DONE :: Sliding window TC object is saved to PRIMEloci_tmp")
