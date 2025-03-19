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



apply_coreovl_postprocessing <- function(input_dir, output_dir, input_file = NULL, core_width = 151, 
                                         score_threshold = 0.75, score_diff = 0.10, num_cores = NULL) {
  message("\n### Running coreovlwith-d postprocess ###")

  # If input_file is not provided, find all matching files in input_dir
  if (is.null(input_file)) {
    input_file <- list.files(input_dir, pattern = "pred_all.*_combined\\.bed$", full.names = FALSE)

    if (length(input_file) == 0) {
      stop("Error: No matching 'pred_all*_combined.bed' files found in: ", input_dir)
    }

    message("Auto-detected input files: ")
    print(input_file)
  }

  # Convert file names to full paths
  input_file <- file.path(input_dir, input_file)

  # Ensure at least one file exists
  if (length(input_file) == 0 || !all(file.exists(input_file))) {
    stop("Error: Some input BED files do not exist. Please check the paths.")
  }

  # Set number of cores
  if (is.null(num_cores)) {
    num_cores <- min(25, parallel::detectCores() %/% 2)
  }

  options(future.globals.maxSize = 4 * 1024^3)  # Set memory limit (4GB)
  future::plan(future::multisession, workers = num_cores)
  on.exit(future::plan(future::sequential))  # Reset plan after execution

  # Loop over multiple input files
  for (file in input_file) {
    message("\nProcessing file: ", file)

    # Load and prepare data
    bed <- PRIMEloci::load_bed_file(file)
    gr <- PRIMEloci::create_granges_from_bed(bed)  # Convert to GenomicRanges object

    # Ensure `gr` is not NULL or empty before accessing columns
    if (is.null(gr) || length(gr) == 0) {
      message("Skipping empty file: ", file)
      next
    }

    # Ensure 'score' column exists before filtering
    if (!"score" %in% colnames(GenomicRanges::mcols(gr))) {
      stop("Error: 'score' column is missing in the GRanges object.")
    }

    # Filter GRanges by score threshold (FIXED)
    message("Filtering regions by score threshold...")
    filtered_gr <- gr[GenomicRanges::mcols(gr)[, "score"] >= score_threshold]

    # If no regions remain after filtering, skip file
    if (length(filtered_gr) == 0) {
      message("Skipping file: ", file, " (No regions meet the score threshold)")
      next
    }

    # Get a list of unique chromosomes
    chr_list <- unique(as.character(GenomicRanges::seqnames(filtered_gr)))

    # If no chromosomes are found, skip processing
    if (length(chr_list) == 0) {
      message("Skipping file: ", file, " (No valid chromosomes found)")
      next
    }

    message("Processing chromosomes with core merging...")

    # Process each chromosome using future.apply
    collapsed_gr_list <- future.apply::future_lapply(chr_list, function(chr) {
      tryCatch({
        chr_gr <- filtered_gr[GenomicRanges::seqnames(filtered_gr) == chr]
        
        # Ensure chr_gr is not empty
        if (length(chr_gr) == 0) {
          message("Skipping chromosome: ", chr, " (No valid regions)")
          return(NULL)
        }

        core_gr <- GenomicRanges::resize(chr_gr, width = core_width, fix = "center")
        core_gr$thick <- IRanges::start(core_gr) + floor(core_width / 2)

        # Selective merging based on score difference
        final_merged_gr <- PRIMEloci::selective_merge_cores(core_gr, score_diff)

        return(final_merged_gr)
      }, error = function(e) {
        message(paste("Error processing chromosome", chr, ":", e$message))
        return(NULL)
      })
    })

    # Remove NULL elements from collapsed_gr_list
    collapsed_gr_list <- collapsed_gr_list[!sapply(collapsed_gr_list, is.null)]

    # If all chromosomes failed, skip file
    if (length(collapsed_gr_list) == 0) {
      message("Skipping file: ", file, " (No valid regions after processing)")
      next
    }

    # Combine results and sort
    collapsed_gr <- do.call(c, collapsed_gr_list)
    collapsed_gr <- GenomicRanges::sort(collapsed_gr)

    # Save processed data
    if (!is.null(output_dir) && output_dir != FALSE) {
      input_basename <- tools::file_path_sans_ext(basename(file))
      input_basename <- stringr::str_replace_all(input_basename, "all", as.character(score_threshold))
      input_basename <- stringr::str_replace_all(input_basename, "[^[:alnum:]]", "_")

      PRIMEloci::write_granges_to_bed_coreovlwithd(collapsed_gr, output_dir, input_basename, score_diff)
    }
  }

  message("DONE! Post-processing complete for all files.")
}


library(GenomicRanges)
library(PRIME)
library(PRIMEloci)
library(assertthat)
PRIMEloci <- function(ctss_rse,
                      tc_object = NULL,
                      outdir,
                      ext_dis = 200,
                      save_to_log = TRUE,
                      save_tc = TRUE,
                      tc_object_name = "tc_grl.rds",
                      sld_needed = TRUE,
                      sld_object_name = "sld_tc_grl.rds",
                      sld_by = 20,
                      keep_tmp_dir = TRUE,
                      save_count_profiles = FALSE,
                      profile_dir_name = "PRIMEloci_profiles",
                      profile_file_format = "parquet",
                      python_script_dir,
                      model_path,
                      prefix_out_name,
                      run_postprocessing = FALSE,
                      core_width = 151,
                      score_threshold = 0.75,
                      score_diff = 0.10) {

  # Ensure directories exist
  if (dir.exists(outdir)) {
    warning("Warning: Output directory '", outdir, "' already exists. Files may be overwritten.")
  } else {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  
  primeloci_tmp <- file.path(outdir, "PRIMEloci_tmp")
  if (!dir.exists(primeloci_tmp)) dir.create(primeloci_tmp, recursive = TRUE, showWarnings = FALSE)

  log_dir <- file.path(outdir, "PRIMEloci_log")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # Ensure numeric parameters are integers
  assertthat::assert_that(is.numeric(ext_dis), msg = "ext_dis must be numeric.")
  assertthat::assert_that(ext_dis %% 1 == 0, msg = "ext_dis must be an integer.")
  
  assertthat::assert_that(is.numeric(sld_by), msg = "sld_by must be numeric.")
  assertthat::assert_that(sld_by %% 1 == 0, msg = "sld_by must be an integer.")

  # _2_ Running TC Extraction
  message("\n### Running _2_get_tc_from_ctss.r ###")

  ext_dis <- as.integer(ext_dis)

  if (is.null(tc_object)) {
    tc_grl <- get_tcs_and_extend_fromthick(ctss_rse, ext_dis = ext_dis)
    print(class(tc_grl))
    validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis)
    if (save_tc) saveRDS(tc_grl, file = file.path(outdir, tc_object_name))
  } else {
    validate_tc_object(tc_object, ctss_rse, ext_dis = ext_dis)
    tc_grl <- tc_object
  }

  message("DONE :: TC object is validated and ready to use.")

  # _3_ Running Sliding Window
  if (sld_needed) {
    message("\n### Running _3_get_sld_window_from_tc.r ###")
    cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

    if (inherits(tc_grl, "GenomicRanges::GRanges")) {
      start_time <- Sys.time()
      message("Processing single GRanges object")
      tc_sliding_window_grl <- tc_sliding_window(tc_grl, slide_by = sld_by, expand_by = ext_dis, num_cores = NULL)
      cat("Time taken:", difftime(Sys.time(), start_time, units = "mins"), "minutes\n")
    } else if (inherits(tc_grl, "GenomicRanges::GRangesList") || inherits(tc_grl, "CompressedGRangesList")) {
      tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
        start_time <- Sys.time()
        gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i)
        message("Processing:", gr_name)
        result <- tc_sliding_window(tc_grl[[i]], slide_by = sld_by, expand_by = ext_dis, num_cores = NULL)
        cat("Time taken for", gr_name, ":", difftime(Sys.time(), start_time, units = "mins"), "minutes\n")
        return(result)
      })

      if (length(tc_sliding_window_grl) > 0) {
        tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl)
      } else {
        warning("Processed TC object list is empty. Ensure tc_grl contains valid data.")
      }
    } else {
      stop("tc_grl must be either a GenomicRanges::GRanges or GenomicRanges::GRangesList object.")
    }

    message("\nSaving TC objects to PRIMEloci_tmp ..")
    saveRDS(tc_sliding_window_grl, file.path(primeloci_tmp, sld_object_name))
  }

  # _4_ Running Profile Generation
  message("\n### Running _4_get_profile.r ###")
  region_grl <- if (sld_needed) tc_sliding_window_grl else tc_grl
  PRIMEloci_profile_2(ctss_rse, region_grl, outdir, profile_dir_name, ext_dis, save_count_profiles, profile_file_format)

  # **Step 5 - Run Prediction**
  message("\n### Running Step 5: Prediction ###")
  predict_PRIMEloci(python_script_dir, model_path, outdir, profile_dir_name, prefix_out_name, profile_file_format)
  message("DONE :: Step 5 (Prediction) completed successfully.")

  # **Step 5 Part 2 - Merge BED Files**
  message("\n### Running Step 5 Part 2: Merging BED Files ###")
  combine_bed_files(file.path(outdir, "PRIMEloci_tmp", profile_dir_name, "predictions"), outdir)
  message("DONE :: Step 5 Part 2 (Merging) completed successfully.")

  # **Step 6 - Apply Post-Processing**
  if (run_postprocessing) {
    message("\n### Running Step 6: Post-Processing ###")
    apply_coreovl_postprocessing(outdir, outdir, core_width, score_threshold, score_diff)
    message("DONE :: Step 6 (Post-Processing) completed successfully.")
  }

  # Cleanup temporary files
  if (!keep_tmp_dir) {
    unlink(primeloci_tmp, recursive = TRUE, force = TRUE)
    message("Temporary directory removed:", primeloci_tmp)
  }
}

ctss_rse <- readRDS("/Users/natsudanav/Desktop/PRIMEloci/example/results/filtered_ctss_rse.rds")

PRIMEloci(
  ctss_rse = ctss_rse,
  outdir = "/Users/natsudanav/Desktop/PRIMEloci/example/results",
  ext_dis = 200,
  save_to_log = TRUE,
  save_tc = TRUE,
  tc_object_name = "tc_grl.rds",
  sld_needed = TRUE,
  sld_object_name = "sld_tc_grl.rds",
  sld_by = 20,
  keep_tmp_dir = TRUE,
  save_count_profiles = FALSE,
  profile_dir_name = "PRIMEloci_profiles",
  profile_file_format = "parquet",
  python_script_dir = "/Users/natsudanav/Desktop/PRIMEloci/genomewide_prediction",
  model_path = "/Users/natsudanav/Desktop/PRIMEloci/model/PRIMEloci_GM12878_model_1.0.sav",
  prefix_out_name = "K562-on-PRIMEloci",
  run_postprocessing = TRUE,
  core_width = 151,
  score_threshold = 0.75,
  score_diff = 0.10
)
