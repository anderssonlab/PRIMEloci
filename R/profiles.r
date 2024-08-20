#' Convert Row Name Strands to No-strand
#'
#' This function converts the strand information
#' in row names from "+" or "-" to "*".
#'
#' @param strand_str A character vector of row names
#' containing strand information.
#' @return A character vector with strand information converted to "*".
#'
convert_rowname_to_nostrand <- function(strand_str) {
  strand_str <- gsub("[+-]$", "*", strand_str)
  return(strand_str)
}



#' Modify Profile Row Names
#'
#' This function modifies the row names of profile data by converting strand
#' information to no-strand and ensuring consistency between forward and
#' reverse strands.
#'
#' @param profiles A data frame of profile data where row names represent
#' genomic coordinates.
#' @param count_profiles A list containing count profiles for both forward (`+`)
#' and reverse (`-`) strands.
#' @return A data frame with modified row names
#' if the forward and reverse strand lists
#' are identical after strand conversion.
#' If the lists are not identical,
#' a message is printed, and the row names are not modified.
#'
modify_profile_rownames <- function(profiles, count_profiles) {
  plus_converted <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`+`)) # nolint: line_length_linter.
  minus_converted <- convert_rowname_to_nostrand(rownames(count_profiles$`*`$`-`)) # nolint: line_length_linter.
  if (identical(sort(plus_converted), sort(minus_converted))) {
    message("Both lists are the same after strand conversion.")
    rownames(profiles) <- plus_converted
  } else {
    message("The lists are not the same after strand conversion.")
  }
  return(profiles)
}



#' Combine Plus and Minus Strand Profiles
#'
#' This function combines forward (plus) and reverse (minus) strand count
#' profiles and modifies their row names.
#'
#' @param count_profiles A list containing count profiles for both forward (`+`)
#' and reverse (`-`) strands.
#' @param len_vec An integer specifying the length of the profiles.
#' @return A data frame containing combined profiles with modified row names.
#'
combine_plus_minus_profiles <- function(count_profiles, len_vec) {
  combined_profiles <- cbind(data.frame(count_profiles$`*`$`+`),
                             data.frame(count_profiles$`*`$`-`))
  colnames(combined_profiles) <- c(paste("Plus_", 1:len_vec, sep = ""),
                                   paste("Minus_", 1:len_vec, sep = ""))
  combined_profiles <- modify_profile_rownames(combined_profiles,
                                               count_profiles)
  return(combined_profiles)
}



#' Extract Components from Row Names
#'
#' This function extracts the chromosome, start, end, and strand components 
#' from row names formatted as "chr:start-end;strand".
#'
#' @param row_names_cpn A character vector of row names with each name formatted
#' as "chr:start-end;strand".
#' @return A character vector containing the chromosome, start, end, and strand
#' extracted from the row name.
#'
extract_rowname_components <- function(row_names_cpn) {
  # Split by ':' to separate chromosome and rest
  parts <- strsplit(row_names_cpn, ":")[[1]]
  chr <- parts[1]
  # Split the rest by '-' and ';' to get start, end, and strand
  rest <- strsplit(parts[2], "[-;]")[[1]]
  start <- as.numeric(rest[1])
  end <- as.numeric(rest[2])
  strand <- rest[3]
  return(c(chr, start, end, strand))
}



#' Create GRanges Object from Row Names
#'
#' This function creates a `GRanges` object from row names by extracting
#' the chromosome, start, end, and strand components.
#'
#' @param row_names_str A character vector of row names with each name formatted
#' as "chr:start-end;strand".
#' @return A `GRanges` object created from the extracted components of
#' the rownames.
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
create_granges_from_rownames <- function(row_names_str) {
  # Apply the function to each row name
  components <- t(sapply(row_names_str, extract_rowname_components))
  # Create GRanges object
  gr <- GenomicRanges::GRanges(seqnames = components[, 1],
                               ranges = IRanges::IRanges(start = as.numeric(components[, 2]), # nolint: line_length_linter.
                                                         end = as.numeric(components[, 3])),  # nolint: line_length_linter.
                               strand = components[, 4])
  return(gr)
}



#' Normalized Strand Subtraction
#'
#' This function performs a normalized subtraction between forward (plus) and 
#' reverse (minus) strand signals.
#'
#' @param vec A numeric vector containing the combined forward and reverse 
#' strand signals.
#' @param len_vec An integer specifying the length of the forward strand 
#' signal within `vec`.
#' @return A numeric vector containing the normalized subtraction result of the
#' forward and reverse strand signals.
#'
strands_norm_subtraction <- function(vec, len_vec) {
  p <- as.numeric(vec[1:len_vec])
  m <- as.numeric(vec[(len_vec + 1):(len_vec * 2)])
  # Normalized strand subtraction
  return((p - m) / max(abs(c(p, m))))
}



#' Apply Normalized Strand Subtraction to All Windows
#'
#' This function applies the normalized subtraction between forward (plus) and
#' reverse (minus) strand signals to all windows.
#'
#' @param windows A data frame where each row represents a window containing
#' combined forward and reverse strand signals.
#' @param ext_dis An integer specifying the extension distance,
#' used for adjusting column names.
#' @param len_vec An integer specifying the length of the forward strand signal
#' within each window.
#' @return A data frame containing the normalized subtraction results 
#' for all windows.
#'
strands_norm_subtraction_all <- function(windows, ext_dis, len_vec) {
  #' Apply normalized forward and reverse subtraction to all windows
  #'
  p_min_m_norm_df <- as.data.frame(t(apply(windows, 1,
                                           strands_norm_subtraction, len_vec)))
  pos <- seq(1, ncol(p_min_m_norm_df))
  colnames(p_min_m_norm_df) <- paste0("Pos", pos - ext_dis - 1)
  return(p_min_m_norm_df)
}



#' Wrap up profile for input of the model
#'
#' This function processes multiple columns of
#' data in \code{ctss_rse} and generates
#' count profiles and subtracted normalized profiles,
#' saving them along with metadata
#' as CSV or Parquet files in specified directories.
#'
#' @param ctss_rse SummarizedExperiment object containing count data.
#' @param regions_gr GRanges or GRangesList object specifying genomic regions.
#' @param output_dir Directory where output profiles and metadata will be saved.
#' @param output_dir_name Output directory name within \code{dir_results}.
#' @param output_subdir_name Subdirectory name within \code{outdir_dir_name}
#' where output files will be stored.
#' @param ext_dis Numeric value specifying the extension distance.
#' @param addtn_to_filename Additional text to append to the output file names.
#' @param save_count_profiles Logical value indicating
#' whether to save count profiles.
#' @param file_type Character string indicating the file type for output 
#' ("csv" or "parquet"). Default is "parquet".
#'
#' @details
#' This function iterates through each column of \code{ctss_rse},
#' processes the data to generate count and subtracted normalized profiles,
#' and saves them as CSV or Parquet files.
#' It also creates corresponding metadata files based on
#' row names of the profiles.
#'
#' @importFrom PRIME heatmapData
#' @importFrom arrow write_parquet
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' wrapup_make_profiles(ctss_rse, regions_gr, dir_results, outdir_dir_name,
#'                      outdir_subdir_name, ext_dis, file_type = "parquet")
#' }
#'
#' @export
wrapup_make_profiles <- function(ctss_rse,
                                 regions_gr,
                                 output_dir,
                                 output_dir_name,
                                 output_subdir_name,
                                 ext_dis,
                                 addtn_to_filename = "",
                                 save_count_profiles = FALSE,
                                 file_type = "parquet") {
  for (i in seq_along(SummarizedExperiment::colnames(ctss_rse))) {

    print(SummarizedExperiment::colnames(ctss_rse)[i])

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(formatted_datetime)

    if (inherits(regions_gr, "GRangesList")) {
      current_region_gr <- regions_gr[[i]]
    } else if (inherits(regions_gr, "GRanges")) {
      current_region_gr <- regions_gr
    } else {
      stop("regions_gr is neither GRanges nor GRangesList")
    }

    ctss_gr <- cast_rse_to_granges(ctss_rse, assay = "counts", coln_assay = i)

    current_region_gr <- convert_strand_to_nostrand_gr(current_region_gr)
    current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

    print("Start making profile")
    count_profiles <- PRIME::heatmapData(current_region_gr, ctss_gr)
    print("Finish making profile")
    rm(current_region_gr, ctss_gr)

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(formatted_datetime)

    len_vec <- ext_dis * 2 + 1

    combined_count_profiles <- combine_plus_minus_profiles(count_profiles, len_vec)
    rm(count_profiles)

    combined_subtnorm_profiles <- strands_norm_subtraction_all(combined_count_profiles, ext_dis, len_vec)

    combined_count_metadata <- create_granges_from_rownames(rownames(combined_count_profiles))
    sum_count <- data.frame(rowSums(combined_count_profiles))
    colnames(sum_count) <- "sum_count"
    combined_count_metadata$sum_count <- sum_count

    # Add rownames as a column before saving
    combined_count_metadata$rownames <- rownames(combined_count_metadata)

    if (file_type == "csv") {
      # Save as CSV without rownames
      write.csv(as.data.frame(combined_count_metadata),
                file = file.path(output_dir,
                                 output_dir_name,
                                 "metadata",
                                 output_subdir_name,
                                 paste0("metadata_count_", output_subdir_name,
                                        addtn_to_filename, "_",
                                        SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                row.names = FALSE)
    } else if (file_type == "parquet") {
      # Save as Parquet
      arrow::write_parquet(as.data.frame(combined_count_metadata),
                           file.path(output_dir,
                                     output_dir_name,
                                     "metadata",
                                     output_subdir_name,
                                     paste0("metadata_count_", output_subdir_name,
                                            addtn_to_filename, "_",
                                            SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
    }
    rm(combined_count_metadata)

    if (save_count_profiles) {
      combined_count_profiles$rownames <- rownames(combined_count_profiles)

      if (file_type == "csv") {
        # Save as CSV without rownames
        write.csv(as.data.frame(combined_count_profiles),
                  file = file.path(output_dir,
                                   output_dir_name,
                                   "profiles",
                                   output_subdir_name,
                                   paste0("profiles_count_", output_subdir_name,
                                          addtn_to_filename, "_",
                                          SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                  row.names = FALSE)
      } else if (file_type == "parquet") {
        # Save as Parquet
        arrow::write_parquet(as.data.frame(combined_count_profiles),
                             file.path(output_dir,
                                       output_dir_name,
                                       "profiles",
                                       output_subdir_name,
                                       paste0("profiles_count_", output_subdir_name,
                                              addtn_to_filename, "_",
                                              SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
      }
    }
    rm(combined_count_profiles)

    combined_subtnorm_profiles$rownames <- rownames(combined_subtnorm_profiles)

    if (file_type == "csv") {
      # Save as CSV without rownames
      write.csv(as.data.frame(combined_subtnorm_profiles),
                file = file.path(output_dir,
                                 output_dir_name,
                                 "profiles_subtnorm",
                                 output_subdir_name,
                                 paste0("profiles_subtnorm_", output_subdir_name,
                                        addtn_to_filename, "_",
                                        SummarizedExperiment::colnames(ctss_rse)[i], ".csv")),
                row.names = FALSE)
    } else if (file_type == "parquet") {
      # Save as Parquet
      arrow::write_parquet(as.data.frame(combined_subtnorm_profiles),
                           file.path(output_dir,
                                     output_dir_name,
                                     "profiles_subtnorm",
                                     output_subdir_name,
                                     paste0("profiles_subtnorm_", output_subdir_name,
                                            addtn_to_filename, "_",
                                            SummarizedExperiment::colnames(ctss_rse)[i], ".parquet")))
    }
    rm(combined_subtnorm_profiles)

    gc()

    current_datetime <- Sys.time()
    formatted_datetime <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
    print(current_datetime)
  }
}
