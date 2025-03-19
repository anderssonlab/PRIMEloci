PRIMEloci <- function(ctss_rse,
                      tc_object = NULL,
                      outdir = outdir,
                      ext_dis = 200,
                      save_to_log = TRUE,
                      save_tc = TRUE,
                      tc_object_name = "tc_grl.rds",
                      sld_needed = TRUE,
                      sld_object_name = "sld_tc_grl.rds",
                      sld_by = 20,
                      keep_tmp_dir = FALSE) {


  # Ensure directories exist, warn if it already exists (may overwrite files)
  if (dir.exists(outdir)) {
    warning("Warning: Output directory '",
            outdir,
            "' already exists. Files may be overwritten.")
  } else {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  primeloci_tmp <- file.path(outdir, "PRIMEloci_tmp")
  if (dir.exists(primeloci_tmp)) {
    warning("Warning: Temporary directory 'PRIMEloci_tmp' already exists. Files may be overwritten.") # nolint: line_length_linter.
  } else {
    dir.create(primeloci_tmp, recursive = TRUE, showWarnings = FALSE)
  }
  log_dir <- file.path(outdir, "PRIMEloci_log")
  if (dir.exists(log_dir)) {
    warning("Warning: Log directory 'PRIMEloci_log' already exists. Files may be overwritten.") # nolint: line_length_linter.
  } else {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }


  numeric_vars <- c("ext_dis", "sld_by")
  for (var in numeric_vars) {
    if (!is.numeric(get(var, envir = parent.frame()))) {
      stop(var, " must be a numeric value.")
    }
    assign(var,
           as.integer(get(var, envir = parent.frame())),
           envir = parent.frame())
    if (get(var, envir = parent.frame()) != as.numeric(get(var, envir = parent.frame()))) { # nolint: line_length_linter.
      warning(var,
              " must be an integer. Automatically converting to ",
              get(var, envir = parent.frame()))
    }
  }


  # _2_
  message("\n### Running _2_get_tc_from_ctss.r ###")

  ext_dis <- as.integer(ext_dis)

  if (is.null(tc_object)) {
    # If tc_object is not provided, generate it
    tc_grl <- get_tcs_and_extend_fromthick(ctss_rse, ext_dis = ext_dis)
    validate_tc_object(tc_grl, ctss_rse, ext_dis = ext_dis)
    if (save_tc) {
      saveRDS(tc_grl, file = file.path(outdir, tc_object_name))
      message(paste("TC object saved to", file.path(outdir, tc_object_name)))
    }
  } else {
    # Validate the provided TC object and assign it to tc_grl
    validate_tc_object(tc_object, ctss_rse, ext_dis = ext_dis)
    tc_grl <- tc_object
  }

  message("DONE :: TC object is validated and ready to use.")


  # _3_
  if (sld_needed) {

    message("\n### Running _3_get_sld_window_from_tc.r ###")
    cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    overall_start_time <- Sys.time()

    if (inherits(tc_grl, "GenomicRanges::GRanges")) {
      # If tc_grl is a single GRanges object
      start_time <- Sys.time()
      message("Processing single GRanges object")
      tc_sliding_window_grl <- tc_sliding_window(tc_grl,
                                                 slide_by = sld_by,
                                                 expand_by = ext_dis,
                                                 num_cores = NULL)
      cat("Time taken for GRanges:", difftime(Sys.time(),
                                              start_time,
                                              units = "mins"),
          "minutes\n")

    } else if (inherits(tc_grl, "GenomicRanges::GRangesList")) {

      tc_sliding_window_grl <- lapply(seq_along(tc_grl), function(i) {
        start_time <- Sys.time()
        gr_name <- if (!is.null(names(tc_grl))) names(tc_grl)[i] else paste0("Sample_", i) # nolint: line_length_linter.
        message("Processing:", gr_name)
        result <- tc_sliding_window(tc_grl[[i]],
                                    slide_by = sld_by,
                                    expand_by = ext_dis,
                                    num_cores = NULL)
        cat("Time taken for", gr_name, ":", difftime(Sys.time(),
                                                     start_time,
                                                     units = "mins"),
            "minutes\n")
        return(result)
      })

      if (length(tc_sliding_window_grl) > 0) {
        tc_sliding_window_grl <- GenomicRanges::GRangesList(tc_sliding_window_grl) # nolint: line_length_linter.
      } else {
        warning("Processed TC object list is empty. Ensure tc_grl contains valid data.") # nolint: line_length_linter.
      }

    } else {
      stop("tc_grl must be either a GenomicRanges::GRanges or GenomicRanges::GRangesList object.") # nolint: line_length_linter.
    }

    # Save the processed data
    message("\nSaving TC objects to PRIMEloci_tmp ..")
    saveRDS(tc_sliding_window_grl, file.path(primeloci_tmp, sld_object_name))

    cat("Overall time taken:", difftime(Sys.time(),
                                        overall_start_time,
                                        units = "mins"),
        "minutes\n")
    cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  }

  # _4_
  message("\n### Running _4_get_profile.r ###")
  # if sld = TRUE, run sld_tc_grl
  PRIMEloci_profile_2(ctss_rse,
                      tc_grl,
                      outdir_dir,
                      outdir_dir_name,
                      ext_dis,
                      save_count_profiles = save_count_profiles, # nolint: line_length_linter.
                      file_type = file_format
                      )
  #else tc_grl






  # cleans up temporary files unless the user explicitly wants to keep them
  if (!keep_tmp_dir) {
    unlink(primeloci_tmp, recursive = TRUE, force = TRUE)
    message("Temporary directory removed:", primeloci_tmp)
  }


}





# make sure profile hangle both gr and grl
