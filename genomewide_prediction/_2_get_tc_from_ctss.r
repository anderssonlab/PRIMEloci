#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(argparse)
  library(CAGEfightR)
  library(GenomicRanges)
  library(PRIME)
  library(assertthat)
}))

### ARGPARSE
parser <- ArgumentParser()

# Input
parser$add_argument("-i", "--ctss_rse", default = "./ctss_rse.rds",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "tc_grl.rds",
                    help = "Output file name for tc_grl object")

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

args <- parser$parse_args()

# Setting up variables
ext_dis <- as.integer(args$ext_dis)

infile_ctss_rse <- args$ctss_rse
ctss_rse <- readRDS(infile_ctss_rse)

output_dir <- args$output_dir
PRIME:::create_output_dir(output_dir)
outfile_tc_grl <- args$outfile


assertthat::assert_that(
  methods::is(ctss_rse, "RangedSummarizedExperiment"),
  msg = "`❌ ctss_rse` must be a RangedSummarizedExperiment object."
)

plc_message("🚀 Running PRIMEloci -2: get extended the tc object provided and validated") # nolint: line_length_linter.
plc_message(sprintf("🕒 Pipeline started at: %s", Sys.time()))

plc_message("Creating tc object ...")
tc_grl <- PRIME::plc_get_tcs_and_extend_fromthick(ctss_rse,
                                                  ext_dis = ext_dis)
plc_message("Saving tc object ...")
saveRDS(tc_grl, file = file.path(output_dir, outfile_tc_grl))

plc_message("Validating tc object ...")
validate_tc <- PRIME::plc_validate_tc_object(tc_grl,
                                             ctss_rse,
                                             ext_dis = ext_dis)
if (!validate_tc) {
  plc_error("TC object validation failed. Ensure the TC object is valid.")
}
plc_message("✅ DONE :: TC object is created, validated, and ready to use.")
