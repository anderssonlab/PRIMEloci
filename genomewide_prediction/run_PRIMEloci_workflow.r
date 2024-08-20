#' Run the PRIMEloci Workflow
#'
#' This function encapsulates the entire PRIMEloci workflow,
#' including creating profiles,
#' running a Python script for profile prediction,
#' and processing the results into a `GRangesList`.
#'
#' @param ctss_rse A `SummarizedExperiment` object containing CTSS data.
#' @param tc_grl A `GRangesList` object containing TC data.
#' @param config_file Character. Path to the configuration file in YAML format.
#' Default is "config.yaml".
#'
#' @return A `GRangesList` containing the processed results.
#'
#' @importFrom yaml read_yaml
#' @importFrom tools file_path_sans_ext
#' @importFrom assertthat assert_that
#' @importFrom data.table fread fwrite
#' @importFrom argparse ArgumentParser
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom GenomicRanges GRangesList sort
#' @importFrom rtracklayer import
#' @importFrom PRIMEloci prep_profile_dir wrapup_make_profiles wrapup_filter_bed_to_reduce # nolint: line_length_linter.
#' @importFrom reticulate py_run_string
#' @export
#' #' @examples
#' Assuming you have already created the ctss_rse and tc_grl objects
#' ctss_rse <- loadRDS("path_to_ctss_rse.rds")
#' tc_grl <- loadRDS("path_to_tc_grl.rds")
#'
#' # Run the PRIMEloci workflow with a YAML configuration file
#' gr_list <- run_PRIMEloci_workflow(
#'   ctss_rse = ctss_rse,
#'   tc_grl = tc_grl,
#'   config_file = "config_PRIMEloci.yaml"
#' )
#' 
run_PRIMEloci_workflow <- function(ctss_rse,
                                   tc_grl,
                                   config_file = "config_PRIMEloci.yaml") {

  # Load configuration from YAML file
  config <- yaml::read_yaml(config_file)

  # Set up directories and file paths based on config
  prediction_dir <- file.path(config$output_dir,
                              config$profile_main_dir,
                              "predictions",
                              config$profile_sub_dir)
  outdir_main_name <- c("metadata",
                        "profiles",
                        "profiles_subtnorm",
                        "predictions")
  outdir_subdir_name <- config$profile_sub_dir

  # Create output directories
  PRIMEloci::prep_profile_dir(output_dir = config$output_dir,
                              output_dir_name = config$profile_main_dir,
                              output_main_name = outdir_main_name,
                              output_subdir_name = outdir_subdir_name)

  # Create profiles for the specified subdir_name
  writeLines(paste0("\nCreating profiles for ", outdir_subdir_name, ".."))
  PRIMEloci::wrapup_make_profiles(ctss_rse,
                                  tc_grl,
                                  config$output_dir,
                                  config$profile_main_dir,
                                  outdir_subdir_name,
                                  config$ext_dis,
                                  save_count_profiles = config$save_count_profiles, # nolint: line_length_linter.
                                  file_type = config$profile_file_type)

  writeLines("\n### Finished profile creation ###\n")

  # Run Python script for profile prediction
  run_python_script <- function(script_path,
                                script_dir,
                                profile_main_dir,
                                profile_sub_dir,
                                model_path, name_prefix,
                                threshold = 0.5,
                                file_format = "parquet") {
    # Check if the script exists
    if (file.exists(script_path)) {
      # Create a Python command with the required arguments
      command <- sprintf(
        "import sys; sys.argv = ['%s', '-w', '%s', '-p', '%s', '-r', '%s', '-m', '%s', '-n', '%s', '-t', '%f', '-f', '%s']; exec(open('%s').read())", # nolint: line_length_linter.
        script_path, script_dir, profile_main_dir, profile_sub_dir, model_path, name_prefix, threshold, file_format, script_path # nolint: line_length_linter.
      )

      # Run the Python script with arguments
      reticulate::py_run_string(command)
    } else {
      stop("The specified script does not exist.")
    }
  }

  run_python_script(
    script_path = "_predict_profile_probabilities.py",
    script_dir = config$python_script_dir,
    profile_main_dir = file.path(config$output_dir, config$profile_main_dir),
    profile_sub_dir = config$profile_sub_dir,
    model_path = config$model_path,
    name_prefix = config$prefix_out_name,
    threshold = config$threshold,
    file_format = config$profile_file_type
  )

  process_all_files <- function(prediction_dir, partial_name) {
    # List all files matching the partial name in the specified directory
    files <- list.files(path = prediction_dir,
                        pattern = partial_name,
                        full.names = TRUE)
    
    if (length(files) == 0) {
      stop("No files found matching the pattern in the specified directory.")
    }
    
    gr_list <- GenomicRanges::GRangesList()  # Initialize an empty GRangesList
    
    for (file in files) {
      writeLines(paste("Processing", file, "..."))
      
      # Call wrapup_filter_bed_to_reduce and store the returned GRanges in the list
      gr <- PRIMEloci::wrapup_filter_bed_to_reduce(file, prediction_dir)
      
      # Extract the base name of the file without the extension to use as the name
      file_name <- tools::file_path_sans_ext(basename(file))
      
      # Add the GRanges object to the GRangesList with the file name as the name
      gr_list[[file_name]] <- gr
    }
    
    # Return the combined GRangesList with names
    return(gr_list)
  }
  

  # Call the main function with the parsed argument
  gr_list <- process_all_files(prediction_dir, config$partial_name)

  writeLines("Done!")

  return(gr_list)
}
