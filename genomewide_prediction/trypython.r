writeLines("\n# Importing R libraries..")
suppressPackageStartupMessages({
  library(os)
  library(assertthat)
  library(data.table)
  library(argparse)
  library(foreach)
  library(doParallel)

  library(GenomicRanges)
  library(rtracklayer)
  library(PRIMEloci)
  library(reticulate)
})


CTSS_RSE_NAME="ctss_rse.rds"
TC_GRL_NAME="tc_grl.rds"

# NEED TO DEFINE THESE
OUTPUT_DIR="../example/results"
PROFILE_MAIN_DIR="example_tc_profiles"
PYTHON_SCRIPT_DIR="."
MODEL_PATH="../example/resources/model_Meena_v3_gm12878_fullModel_C.sav"
PREFIX_OUT_NAME="K562-on-GM12878-C-model"


# CAN BE USE THIS SETTING AS A DEFAULT
PROFILE_SUB_DIR="tcs"
PROFILE_FILE_TYPE="parquet"
# add -s if you want to save count profiles
THRESHOLD=0.2
PREDICTION_DIR=file.path(OUTPUT_DIR, PROFILE_MAIN_DIR, "predictions", PROFILE_SUB_DIR)
PARTIAL_NAME="pred_slt.*\\.bed"

outdir_dir <- OUTPUT_DIR
outdir_dir_name <- PROFILE_MAIN_DIR
file_format <- PROFILE_FILE_TYPE
save_count_profiles <- FALSE

ext_dis <- 200
outdir_main_name <- c("metadata",
                      "profiles",
                      "profiles_subtnorm",
                      "predictions")
outdir_subdir_name <- unlist(strsplit(PROFILE_SUB_DIR, ","))

ctss_rse <- readRDS(file.path(outdir_dir, CTSS_RSE_NAME))
tc_grl <- readRDS(file.path(outdir_dir, TC_GRL_NAME))


# Create output directory
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)


# Define a function to create profiles for a single subdir_name
create_profiles <- function(subdir_name) {
  writeLines(paste0("\nCreating profiles for ", subdir_name, ".."))
  report_time_execution(wrapup_make_profiles(ctss_rse,
                                             tc_grl,
                                             outdir_dir,
                                             outdir_dir_name,
                                             subdir_name,
                                             ext_dis,
                                             save_count_profiles = save_count_profiles, # nolint: line_length_linter.
                                             file_type = file_format))
}

lapply(outdir_subdir_name, create_profiles)

writeLines("\n### Finished get_tc_profiles.r ###\n")

run_python_script <- function(script_path, script_dir, profile_main_dir, profile_sub_dir, model_path, name_prefix, threshold = 0.5, file_format = 'parquet') {
  library(reticulate)
  
  # Ensure the correct Python environment is used
  # use_python("/path/to/python", required = TRUE)  # Adjust this path as needed
  
  # Check if the script exists
  if (file.exists(script_path)) {
    # Create a Python command with the required arguments
    command <- sprintf(
      "import sys; sys.argv = ['%s', '-w', '%s', '-p', '%s', '-r', '%s', '-m', '%s', '-n', '%s', '-t', '%f', '-f', '%s']; exec(open('%s').read())",
      script_path, script_dir, profile_main_dir, profile_sub_dir, model_path, name_prefix, threshold, file_format, script_path
    )
    
    # Run the Python script with arguments
    py_run_string(command)
  } else {
    stop("The specified script does not exist.")
  }
}

run_python_script(
  script_path = "_predict_profile_probabilities.py",
  script_dir = PYTHON_SCRIPT_DIR,
  profile_main_dir = file.path(OUTPUT_DIR, PROFILE_MAIN_DIR),
  profile_sub_dir = PROFILE_SUB_DIR,
  model_path = MODEL_PATH,
  name_prefix = PREFIX_OUT_NAME,
  threshold = THRESHOLD,
  file_format = PROFILE_FILE_TYPE
)



process_all_files <- function(prediction_dir, partial_name) {
  # List all files matching the partial name in the specified directory
  files <- list.files(path = prediction_dir, pattern = partial_name, full.names = TRUE)

  if (length(files) == 0) {
    stop("No files found matching the pattern in the specified directory.")
  }

  gr_list <- GRangesList()  # Initialize an empty GRangesList

  for (file in files) {
    writeLines(paste("Processing", file, "..."))
    
    # Call wrapup_filter_bed_to_reduce and store the returned GRanges in the list
    gr <- wrapup_filter_bed_to_reduce(file, prediction_dir)
    
    # Add the GRanges object to the GRangesList
    gr_list <- c(gr_list, GRangesList(gr))
  }

  # Return the combined GRangesList
  return(gr_list)
}


# Call the main function with the parsed argument
gr_list <- process_all_files(PREDICTION_DIR, PARTIAL_NAME)

writeLines("Done!")