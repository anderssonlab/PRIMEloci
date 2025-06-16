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
parser$add_argument("--ctss_rse", default = "./ctss_rse.rds",
                    help = "FULL PATH to input file name for ctss_rse rds object") # nolint: line_length_linter.
parser$add_argument("--region", default = "./tc_grl.rds",
                    help = "FULL PATH to input file name for input regions rds object") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")

# Parameters
parser$add_argument("-e", "--ext_dis", default = 200,
                    help = "Extension distance")

args <- parser$parse_args()

# Setting up variables
ext_dis <- as.integer(args$ext_dis)

infile_ctss_rse <- args$ctss_rse
ctss_rse <- readRDS(infile_ctss_rse)

infile_region_gr <- args$region
region_gr <- readRDS(infile_region_gr)

output_dir <- args$output_dir
PRIME::plc_create_output_dir(args$output_dir)
outfile_ctss_rse <- args$outfile

assertthat::assert_that(
  methods::is(ctss_rse, "RangedSummarizedExperiment"),
  msg = "`ctss_rse` must be a RangedSummarizedExperiment object."
)

plc_message("🚀 Validating the region gr object provided")
ext_dis <- as.integer(args$ext_dis)
len_vec <- ext_dis * 2 + 1

# Check if region object is a GRanges
assertthat::assert_that(
  inherits(region_gr, "GRanges"),
  msg = "❌ The object must be a GRanges object."
)

plc_message("🚀 Running PRIMEloci -0: validating the ctss and region object provided") # nolint: line_length_linter.

# Ensure all regions have the correct width
if (all(GenomicRanges::width(region_gr) != len_vec)) {
  msg <- paste("⚠️ All regions in the object (GRanges) must have width",
               len_vec,
               " : extend 401 bp from thick if existed. It was saved as an extended object.") # nolint: line_length_linter.
  region_gr <- PRIME::plc_extend_fromthick(tc_gr = region_gr,
                                           ext_dis = ext_dis)
  saveRDS(region_gr, file = paste0(tools::file_path_sans_ext(infile_region_gr),
                                   "_extended.rds"))

} else {
  msg <- paste("✅ All regions in the object (GRanges) have width", len_vec) # nolint: line_length_linter.
}
plc_message(msg)


validate_tc <- PRIME::plc_validate_tc_object(region_gr,
                                             ctss_rse,
                                             ext_dis = ext_dis)
if (!validate_tc) {
  plc_error("TC object validation failed. Ensure the TC object is valid.")
}

plc_message("✅ DONE :: Region object is validated, and ready to use.")
