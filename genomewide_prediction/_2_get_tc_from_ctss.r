#!/usr/bin/env Rscript

writeLines("\n### Running _2_get_tc_from_ctss.r ###")

writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(argparse)
  library(CAGEfightR)
  library(GenomicRanges)
  library(PRIME)
  library(assertthat)
})

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
parser$add_argument("-l", "--log", default = NULL,
                    help = "Log file name e.g. PRIMEloci-2.log")

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

args <- parser$parse_args()

# Setting up variables
ext_dis <- as.integer(args$ext_dis)

infile_ctss_rse <- args$ctss_rse
ctss_rse <- readRDS(infile_ctss_rse)

output_dir <- args$output_dir
create_output_dir(output_dir)
outfile_tc_grl <- args$outfile

log <- if (args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

assertthat::assert_that(
  methods::is(ctss_rse, "RangedSummarizedExperiment"),
  msg = "`ctss_rse` must be a RangedSummarizedExperiment object."
)

plc_log("\n\n\n 🚀 Running PRIMEloci -2: get extended tc and validate the tc object provided", # nolint: line_length_linter.
        log_target)
plc_log(sprintf("🕒 Pipeline started at: %s", Sys.time()), log_target)

plc_log("🔹 Creating tc object\n", log_target)
tc_grl <- PRIME::get_tcs_and_extend_fromthick(ctss_rse,
                                              ext_dis = ext_dis)
plc_log("🔹 Saving tc object\n", log_target)
saveRDS(tc_grl, file = file.path(output_dir, outfile_tc_grl))

plc_log("🔹 Validating tc object\n", log_target)
validate_tc <- PRIME::validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis)
if (!validate_tc) {
  msg <- "\nTC object validation failed. Ensure the TC object is valid."
  plc_log(msg, log_target, level = "❌ ERROR")
  stop(msg)
}
plc_log("✅ DONE :: TC object is created, validated, and ready to use.",
        log_target)
