# Setting up variables
infile_ctss_rse <- "/Users/natsudanav/Desktop/PRIMEloci/example/results/ctss_rse.rds"
infile_tc_grl <- "/Users/natsudanav/Desktop/PRIMEloci/example/results/tc_grl.rds"

outdir_dir <- "/Users/natsudanav/Desktop/PRIMEloci/example/results"
outdir_dir_name <- "testtest"
file_format <- "parquet"
save_count_profiles <- TRUE

ext_dis <- 200



# Read in .RDS files
writeLines("\nReading in data..")
ctss_rse <- readRDS(infile_ctss_rse)
tc_grl <- readRDS(infile_tc_grl)


# Create output directory
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name)



PRIMEloci_profile(ctss_rse,
                  tc_grl,
                  outdir_dir,
                  outdir_dir_name,
                  ext_dis,
                  save_count_profiles = save_count_profiles, # nolint: line_length_linter.
                  file_type = file_format)


lapply(outdir_subdir_name, create_profiles)

writeLines("\n### Finished get_tc_profiles.r ###\n")

library(GenomicRanges)
library(S4Vectors)
library(parallel)
library(SummarizedExperiment)
library(arrow)
library(PRIMEloci)
library(PRIME)

# Split the GRanges object into individual ranges
regions_list <- split(tc_grl[[1]], seqnames(tc_grl[[1]]))

# Check the result
print(length(regions_list))  # This should be equal to the number of ranges in regions_gr
print(regions_list[1])



tc_sliding_window(gr,
                  slide_by = 20,
                  expand_by = 200,
                  use_max_cores = TRUE)
