#predict_PRIMEloci <- function(script_dir, model_path, output_dir, profile_main_dir, prefix_out_name, profile_file_type, calibrate = FALSE) {
#  # Ensure reticulate is loaded and Python is set
#  library(reticulate)
#  use_python("/usr/bin/python3", required = TRUE)
#  py_config()
#
#  # Define the path to the Python script
#  python_script <- file.path(script_dir, "_5_predict_profile_probability.py")
#
#  # Ensure required directories exist
#  if (!dir.exists(output_dir)) {
#    stop("Error: Output directory '", output_dir, "' does not exist. Please check the path.")
#  }
#  if (!file.exists(python_script)) {
#    stop("Error: Python script '", python_script, "' not found in the specified script directory.")
#  }
#
#  message("\n### Running _5_predict_profile_probability.py ###")
#
#  # Construct arguments list
#  args <- c(
#    "-w", script_dir,
#    "-m", model_path,
#    "-p", file.path(output_dir, "PRIMEloci_tmp", profile_main_dir),
#    "-n", prefix_out_name,
#    "-f", profile_file_type
#  )
#
#  # Run the Python script using system2()
#  result <- system2("python3", args = c(python_script, args), stdout = TRUE, stderr = TRUE)
#
#  # Print output from the Python script
#  cat(result, sep = "\n")
#
#  message("DONE :: Prediction completed successfully.")
#
#  ### **Step 5 Part 2: Merge BED Files Using R**
#  message("\n### Merging BED Files in R ###")
#
#  # Define the directory where BED files are stored
#  bed_input_dir <- file.path(output_dir, "PRIMEloci_tmp", profile_main_dir, "predictions")
#
#  # List all BED files
#  bed_files <- list.files(bed_input_dir, pattern = "\\.bed$", full.names = TRUE)
#
#  if (length(bed_files) == 0) {
#    warning("No BED files found in: ", bed_input_dir)
#    return()
#  }
#
#  # Read and merge all BED files
#  merged_bed <- rbindlist(lapply(bed_files, fread), use.names = TRUE, fill = TRUE)
#
#  # Define output BED file path
#  merged_bed_file <- file.path(output_dir, paste0(prefix_out_name, "_merged.bed"))
#
#  # Save the merged BED file
#  fwrite(merged_bed, merged_bed_file, sep = "\t", quote = FALSE, col.names = FALSE)
#
#  message("DONE :: Merged BED file saved to: ", merged_bed_file)
#}

predict_PRIMEloci <- function(script_dir, model_path, output_dir, profile_main_dir, prefix_out_name, profile_file_type, calibrate = FALSE) {
  # Ensure Python is set correctly
  reticulate::use_python("/usr/bin/python3", required = TRUE)
  reticulate::py_config()

  # Define the path to the Python script
  python_script <- file.path(script_dir, "_5_predict_profile_probability.py")

  # Ensure required directories exist
  if (!dir.exists(output_dir)) {
    stop("Error: Output directory '", output_dir, "' does not exist. Please check the path.")
  }
  if (!file.exists(python_script)) {
    stop("Error: Python script '", python_script, "' not found in the specified script directory.")
  }

  message("\n### Running _5_predict_profile_probability.py ###")

  # Construct arguments list
  args <- c(
    "-w", script_dir,
    "-m", model_path,
    "-p", file.path(output_dir, "PRIMEloci_tmp", profile_main_dir),
    "-n", prefix_out_name,
    "-f", profile_file_type
  )

  # Run the Python script using system2()
  result <- system2("python3", args = c(python_script, args), stdout = TRUE, stderr = TRUE)
  cat(result, sep = "\n")
  message("DONE :: Prediction completed successfully.")
}

combine_bed_files <- function(input_dir, output_dir) {
  message("\n### Running combine_bed_files ###")

  # Ensure output directory exists
  if (!base::dir.exists(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Step 1: Get all BED files and extract unique prefixes
  bed_files <- base::list.files(input_dir, pattern = "\\.bed$", full.names = TRUE)

  # Extract sample prefixes before `_chrXXX_`
  prefixes <- base::unique(base::sub("(_chr[^_]+)_.*\\.bed$", "", base::basename(bed_files)))

  # Debugging: Print found prefixes
  message("Detected prefixes:")
  base::print(prefixes)

  # Step 2: Combine files per prefix
  for (prefix in prefixes) {
    # Get all matching files for this prefix
    matching_files <- bed_files[base::grepl(base::paste0("^", prefix, "_chr[^_]+_"), base::basename(bed_files))]

    if (base::length(matching_files) == 0) {
      message("No matching files found for prefix: ", prefix)
      next  # Skip if no matching files
    }

    message("Combining files for prefix: ", prefix)

    # Define output file path
    combined_file <- base::file.path(output_dir, base::paste0(prefix, "_combined.bed"))

    # Read first file to extract the header
    first_file <- matching_files[1]
    header <- base::readLines(first_file, n = 1)
    col_names <- base::strsplit(header, "\t")[[1]]  # Extract column names from first line

    # Read and merge all BED files while ensuring consistent column names
    bed_data <- base::lapply(matching_files, function(file) {
      data <- data.table::fread(file, skip = 1)  # Skip first line (header)
      data.table::setnames(data, col_names)  # Apply column names dynamically
      return(data)
    })

    # Combine all files into one
    merged_bed <- data.table::rbindlist(bed_data, use.names = TRUE, fill = TRUE)

    # Ensure column types are correctly formatted before writing
    if ("chromStart" %in% col_names) merged_bed[, chromStart := base::as.integer(chromStart)]
    if ("chromEnd" %in% col_names) merged_bed[, chromEnd := base::as.integer(chromEnd)]
    if ("sum_count" %in% col_names) merged_bed[, sum_count := base::as.integer(sum_count)]
    if ("score" %in% col_names) merged_bed[, score := base::as.numeric(format(score, scientific = FALSE))]

    # Write the BED file ensuring proper formatting
    data.table::fwrite(merged_bed, combined_file, sep = "\t", quote = FALSE, col.names = FALSE, append = TRUE, 
                       scipen = 999)  # Prevent scientific notation

    message("✅ Combined files into: ", combined_file)
  }

  message("✅ DONE :: All BED files successfully merged.")
}







predict_PRIMEloci(
  script_dir = "/Users/natsudanav/Desktop/PRIMEloci/genomewide_prediction",
  model_path = "/Users/natsudanav/Desktop/PRIMEloci/model/PRIMEloci_GM12878_model_1.0.sav",
  output_dir = "/Users/natsudanav/Desktop/PRIMEloci/example/results",
  profile_main_dir = "PRIMEloci_profiles",
  prefix_out_name = "PRIMEloci",
  profile_file_type = "parquet"
)

combine_bed_files(
  input_dir = "/Users/natsudanav/Desktop/PRIMEloci/example/results/PRIMEloci_tmp/PRIMEloci_profiles/predictions",
  output_dir = "/Users/natsudanav/Desktop/PRIMEloci/example/results"
)
