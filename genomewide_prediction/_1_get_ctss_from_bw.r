#!/usr/bin/env Rscript


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
parser$add_argument("-i", "--input_dir", default = "./",
                    help = "CAGE BigWig directory")
parser$add_argument("-m", "--design_matrix", default = "design.matrix.tsv",
                    help = "Design matrix of all CAGE data in .tsv format, example can be found in the 'example/resources' directory") # nolint: line_length_linter.

# Output
parser$add_argument("-o", "--output_dir", default = "./",
                    help = "Output directory")
parser$add_argument("-n", "--outfile", default = "ctss_rse.rds",
                    help = "Output file name for ctss_rse object")
parser$add_argument("-l", "--log", default = NULL,
                    help = "Log file name or NULL to log to console")

# Parameters
parser$add_argument("-k", "--keep_standard_chr", action = "store_true",
                    help = "Keep only standard chromosomes")

args <- parser$parse_args()


# Setting up variables
input_dir <- args$input_dir
design_matrix_file <- args$design_matrix

output_dir <- args$output_dir
create_output_dir(args$output_dir)
outfile_ctss_rse <- args$outfile

log <- if (is.null(args$log) || args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)


plc_log("\n\n\n 🚀 Running PRIMEloci -1: Getting CTSS from bw files",
        log_target)
# Read in CAGE BigWig, based on design matrix provided,
# no filtering/subsetBySupport steps are performed
plc_log("# Reading design matrix and importing BigWig files...", log_target)
design_matrix <- read.table(design_matrix_file,
                            header = TRUE,
                            sep = "\t",
                            row.names = "Name")
plc_log("🔹 Running get_ctss_from_bw() ...", log_target)
ctss_log_output <- capture.output({
  ctss_rse <- get_ctss_from_bw(input_dir, design_matrix)
})
plc_log(paste(ctss_log_output, collapse = "\n"), log_target)

# Keep only standard chromosomes if specified
if (args$keep_standard_chr) {
  plc_log("\n# Keeping only standard chromosomes..", log_target)
  ctss_rse <- GenomeInfoDb::keepStandardChromosomes(ctss_rse,
                                                    pruning.mode = "coarse")
}

# Save the ctss_rse object
file_path <- file.path(output_dir, outfile_ctss_rse)
plc_log(sprintf("📝 Saving CTSS RSE to: %s", file_path), log_target)
saveRDS(ctss_rse, file_path)

plc_log(sprintf("✅ DONE :: CTSS RSE saved to: %s",
                file.path(output_dir, outfile_ctss_rse)),
        log_target)
