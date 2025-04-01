#!/usr/bin/env Rscript

writeLines("\n### Running _3_get_sld_window_from_tc.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(parallel)
  library(argparse)
  library(PRIME)
})

# Create argument parser
parser <- ArgumentParser(description = "Sliding window on genomic ranges")
parser$add_argument("-i", "--infile", required = TRUE,
                    help = "Path to the input .rds file (tc_grl)")

parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "sld_tc_grl.rds",
                    help = "Output file name for sld_tc_grl object")
parser$add_argument("-l", "--log", default = NULL,
                    help = "Log file name e.g. PRIMEloci-3.log")

parser$add_argument("-s", "--sld_by", type = "integer", default = 20,
                    help = "Slide by parameter")
parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")


# Parse arguments
args <- parser$parse_args()

# Load the input data
tc_grl <- readRDS(args$infile)

# Output
output_dir <- args$output_dir
create_output_dir(output_dir)
outfile_sld_tc_grl <- args$outfile

primeloci_tmp <- setup_tmp_dir(output_dir)

log <- if (args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

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

if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "`num_cores` must be a positive integer or NULL."
  )
}

assertthat::assert_that(
  is.numeric(ext_dis),
  ext_dis %% 1 == 0,
  ext_dis > 0,
  msg = "`ext_dis` must be a positive integer."
)


plc_log("\n\n\n 🚀 Running PRIMEloci -3: sliding windows covering reduced TC regions", log_target) # nolint: line_length_linter.

if (inherits(tc_grl, "GenomicRanges::GRanges")) {
  start_time <- Sys.time()
  tc_sliding_window_grl <- PRIME::tc_sliding_window(tc_grl,
                                                    sld_by = sld_by,
                                                    ext_dis = ext_dis,
                                                    log_file = log_target,
                                                    num_cores = num_cores)
  plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                  as.numeric(difftime(Sys.time(),
                                      start_time,
                                      units = "mins"))),
          log_target)
} else if (inherits(tc_grl, "GenomicRanges::GRangesList") ||
             inherits(tc_grl, "CompressedGRangesList")) {
  tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
    start_time <- Sys.time()
    gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
    plc_log(sprintf("🔹 Processing: %s", gr_name), log_target)
    result <- tc_sliding_window(tc_grl[[i]],
                                sld_by = sld_by,
                                ext_dis = ext_dis,
                                log_file = log_target,
                                num_cores = num_cores)
    plc_log(sprintf("⏱️ Time taken: %.2f minutes",
                    as.numeric(difftime(Sys.time(),
                                        start_time,
                                        units = "mins"))),
            log_target)
    result
  })
  if (!is.null(tc_sliding_window_grl) &&
        length(tc_sliding_window_grl) > 0) {
    tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
    names(tc_sliding_window_grl) <- names(tc_grl)
  } else {
    msg <- "Processed TC object list is empty. Ensure tc_grl contains valid data." # nolint: line_length_linter.
    plc_log(msg, log_target, level = "❌ ERROR")
    stop(msg)
  }
} else {
  msg <- "tc_grl must be either a GRanges, GRangesList, or CompressedGRangesList object." # nolint: line_length_linter.
  plc_log(msg, log_target, level = "❌ ERROR")
  stop(msg)
}

plc_log(sprintf("Saving TC objects to PRIMEloci_tmp .."),
        log_target)
saveRDS(tc_sliding_window_grl,
        file.path(primeloci_tmp, outfile_sld_tc_grl))

plc_log("✅ DONE :: Sliding window TC object is saved to PRIMEloci_tmp.",
        log_target)
tc_for_profile <- tc_sliding_window_grl
assertthat::assert_that(!is.null(tc_for_profile),
                        msg = "❌ tc_for_profile is NULL at return. Something failed.") # nolint: line_length_linter.
