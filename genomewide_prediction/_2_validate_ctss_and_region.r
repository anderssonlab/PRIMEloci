#!/usr/bin/env Rscript

writeLines("\n### Running _2_validate_ctss_and_region.r ###")

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
parser$add_argument("-t", "--tc_grl", default = "./tc_grl.rds",
                    help = "FULL PATH to input file name for input regions rds object") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-l", "--log", default = "PRIMEloci-2.log",
                    help = "Log file name e.g. PRIMEloci-2.log")

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

args <- parser$parse_args()

# Setting up variables
ext_dis <- as.integer(args$ext_dis)

infile_ctss_rse <- args$ctss_rse
ctss_rse <- readRDS(infile_ctss_rse)

infile_tc_grl <- args$tc_grl
tc_grl <- readRDS(infile_tc_grl)

output_dir <- args$output_dir
create_output_dir(args$output_dir)
outfile_ctss_rse <- args$outfile

log <- if (is.null(args$log) || args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

assertthat::assert_that(
  methods::is(ctss_rse, "RangedSummarizedExperiment"),
  msg = "`ctss_rse` must be a RangedSummarizedExperiment object."
)

plc_log("\n\n\n 🚀 Running PRIMEloci -2: validating the ctss and region object provided", # nolint: line_length_linter.
        log_target)

plc_log("🔹 Validating tc object\n", log_target)
validate_tc <- PRIME::validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis)
if (!validate_tc) {
  msg <- "\nTC object validation failed. Ensure the TC object is valid."
  plc_log(msg, log_target, level = "❌ ERROR")
  stop(msg)
}
plc_log("✅ DONE :: Region object is validated, and ready to use.",
        log_target)
