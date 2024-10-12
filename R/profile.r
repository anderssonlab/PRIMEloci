#' Process profiles for a single chromosome and save results
#'
#' This function processes the profiles
#' for a specific chromosome (GRanges object),
#' computes the count profiles, combines the results,
#' and saves them in the specified format (CSV or Parquet).
#' The filename is provided by the main function,
#' and this function appends the chromosome name to the filename.
#'
#' @param current_region_gr GRanges object for the current chromosome.
#' @param chr_name The name of the chromosome being processed.
#' @param filtered_ctss_gr A GRanges object containing CTSS data
#' filtered for the current chromosome.
#' @param file_path The base file path (including filename)
#' where the output will be saved.
#' The chromosome name will be appended to this base file path.
#' @param ext_dis Integer value for extending the distance
#' used in profile computations.
#' @param save_count_profiles Logical flag indicating
#' whether count profiles should be saved.
#' @param file_type The format in which to save the files ("csv" or "parquet").
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom arrow write_parquet
#'
PRIMEloci_profile_chr <- function(current_region_gr, # nolint: object_name_linter
                                  filtered_ctss_gr,
                                  chr_name,
                                  file_path,
                                  base_file_name,
                                  ext_dis,
                                  save_count_profiles,
                                  file_type) {

  current_region_gr <- convert_strand_to_nostrand_gr(current_region_gr)
  current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

  # Make count profiles
  count_profiles <- PRIME::heatmapData(current_region_gr, filtered_ctss_gr)
  rm(current_region_gr, filtered_ctss_gr)

  len_vec <- ext_dis * 2 + 1

  # Combine the count profiles
  combined_count_profiles <- combine_plus_minus_profiles(count_profiles,
                                                         len_vec)
  rm(count_profiles)

  combined_subtnorm_profiles <- strands_norm_subtraction_all(combined_count_profiles, # nolint: line_length_linter.
                                                             ext_dis,
                                                             len_vec)

  # Create metadata
  combined_count_metadata <- create_granges_from_rownames(rownames(combined_count_profiles)) # nolint: line_length_linter.
  sum_count <- data.frame(rowSums(combined_count_profiles))
  colnames(sum_count) <- "sum_count"
  combined_count_metadata$sum_count <- sum_count

  # Add rownames as a column and save data

  ## Save profiles
  combined_subtnorm_profiles$rownames <- rownames(combined_subtnorm_profiles)
  save_to_file(combined_subtnorm_profiles, "_profiles_subtnorm",
               file_type, chr_name,
               file.path(file_path, "profiles_subtnorm", base_file_name))
  rm(combined_subtnorm_profiles)

  ## Save count profiles if requested
  if (save_count_profiles) {
    combined_count_profiles$rownames <- rownames(combined_count_profiles)
    save_to_file(combined_count_profiles, "_profiles",
                 file_type, chr_name,
                 file.path(file_path, "profiles", base_file_name))
  }
  rm(combined_count_profiles)

  ## Save metadata per chromosome
  combined_count_metadata$rownames <- rownames(combined_count_metadata)
  save_to_file(combined_count_metadata, "_metadata",
               file_type, chr_name,
               file.path(file_path, "metadata", base_file_name))
  rm(combined_count_metadata)

  gc()
}

#PRIMEloci_profile <- function(ctss_rse, # nolint: object_name_linter.
#                              regions_gr,
#                              output_dir,
#                              output_dir_name,
#                              ext_dis,
#                              addtn_to_filename = "",
#                              save_count_profiles = FALSE,
#                              file_type = "parquet") {
#
#  prep_profile_dir(output_dir = output_dir,
#                   output_dir_name = output_dir_name)
#
#  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {
#
#    print(paste("Processing:", SummarizedExperiment::colnames(ctss_rse)[i]))
#    current_datetime <- Sys.time()
#    print(current_datetime)
#
#    # Temporarily turn off parallelization for debugging
#    # Detect number of cores (commented out)
#    # num_cores <- parallel::detectCores() - 1 # Use one less core to avoid overloading the system
#
#    if (inherits(regions_gr, "GRangesList")) {
#      current_region_gr <- regions_gr[[i]]
#    } else if (inherits(regions_gr, "GRanges")) {
#      current_region_gr <- regions_gr
#    } else {
#      stop("regions_gr is neither GRanges nor GRangesList")
#    }
#
#    # Extract the GRanges for the current column of the ctss_rse
#    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)
#
#    regions_list <- split(current_region_gr,
#                          GenomicRanges::seqnames(current_region_gr))
#
#    # Generate the base file path to be passed
#    # to the chromosome-specific function
#    file_path <- file.path(output_dir, output_dir_name)
#    base_file_path <- paste0(SummarizedExperiment::colnames(ctss_rse)[i],
#                             addtn_to_filename)
#
#    # Run sequentially for debugging (using lapply instead of mclapply)
#    lapply(names(regions_list), function(chr_name) {
#
#      # Filter ctss_gr for the current chromosome
#      filtered_ctss_gr <- ctss_gr[GenomicRanges::seqnames(ctss_gr) == chr_name]
#
#      # Call the chromosome-specific function
#      PRIMEloci_profile_chr(
#        current_region_gr = regions_list[[chr_name]],
#        chr_name = chr_name,
#        filtered_ctss_gr = filtered_ctss_gr,
#        file_path = file_path,
#        base_file_name = base_file_path,
#        ext_dis = ext_dis,
#        save_count_profiles = save_count_profiles,
#        file_type = file_type
#      )
#    })
#
#    current_datetime <- Sys.time()
#    print(current_datetime)
#
#  }
#}


#' Main function to process profiles for each chromosome and save results
#'
#' This function processes the profiles for all chromosomes using a GRanges
#' object, computes the count profiles, and calls a chromosome-specific function
#' to handle individual chromosomes and save the results. The filenames are
#' generated and passed to the chromosome-specific function.
#'
#' If `regions_gr` is a `GRangesList`, the function expects paired profiles
#' between `ctss_rse` and `regions_gr` such that for each element `i`, the
#' `ctss_rse` column and the corresponding element of `regions_gr` will be
#' processed together. It is assumed that the order of the elements in
#' `regions_gr` corresponds to the order in `ctss_rse`.
#'
#' @param ctss_rse A RangedSummarizedExperiment object containing CTSS counts.
#' @param regions_gr A GRanges or GRangesList object containing genomic regions.
#' @param output_dir The directory where output files will be saved.
#' @param output_dir_name The directory name for the output.
#' @param ext_dis Integer value for extending the distance
#' used in profile computations.
#' @param addtn_to_filename A string to add to the output filename.
#' @param save_count_profiles Logical flag indicating whether count profiles
#' should be saved.
#' @param file_type The format in which to save the files ("csv" or "parquet").
#' @param num_cores The number of cores to use for parallel processing.
#' If not provided, it defaults to using all available cores minus one.
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment colnames
#' @importFrom arrow write_parquet
#' @importFrom parallel mclapply detectCores
#' @export
PRIMEloci_profile <- function(ctss_rse, # nolint: object_name_linter.
                              regions_gr,
                              output_dir,
                              output_dir_name,
                              ext_dis,
                              addtn_to_filename = "",
                              save_count_profiles = FALSE,
                              file_type = "parquet",
                              num_cores = NULL) {

  # Prepare the output directory
  prep_profile_dir(output_dir = output_dir,
                   output_dir_name = output_dir_name)

  # Set number of cores to use for parallel processing
  if (is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 1
  }

  # Process each column in ctss_rse
  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {
    print(paste("Processing:", SummarizedExperiment::colnames(ctss_rse)[i]))
    current_datetime <- Sys.time()
    print(current_datetime)

    # Handle GRangesList and GRanges objects
    if (inherits(regions_gr, "GRangesList")) {
      current_region_gr <- regions_gr[[i]]
    } else if (inherits(regions_gr, "GRanges")) {
      current_region_gr <- regions_gr
    } else {
      stop("regions_gr is neither GRanges nor GRangesList")
    }

    # Extract the GRanges for the current column of ctss_rse
    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)

    # Split regions by chromosome
    regions_list <- split(current_region_gr,
                          GenomicRanges::seqnames(current_region_gr))

    # Generate the base file path
    file_path <- file.path(output_dir, output_dir_name)
    base_file_path <- paste0(SummarizedExperiment::colnames(ctss_rse)[i],
                             addtn_to_filename)

    # Parallel execution for each chromosome using mclapply
    parallel::mclapply(names(regions_list), function(chr_name) {
      # Filter ctss_gr for the current chromosome
      filtered_ctss_gr <- ctss_gr[GenomicRanges::seqnames(ctss_gr) == chr_name]

      # Call the chromosome-specific function
      PRIMEloci_profile_chr(
        current_region_gr = regions_list[[chr_name]],
        chr_name = chr_name,
        filtered_ctss_gr = filtered_ctss_gr,
        file_path = file_path,
        base_file_name = base_file_path,
        ext_dis = ext_dis,
        save_count_profiles = save_count_profiles,
        file_type = file_type
      )
    }, mc.cores = num_cores)

    current_datetime <- Sys.time()
    print(current_datetime)
  }
}
