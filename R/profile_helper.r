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
#' 
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



#' Save data to file with an appropriate file extension based on file type
#'
#' This function saves the given data either
#' in CSV or Parquet format based on the specified `file_type`.
#' The function dynamically determines the correct file extension
#' (`.csv` or `.parquet`) and appends the chromosome name and suffix
#' to the provided file path to create the full file name.
#'
#' @param data A data frame or object that can be converted to
#' a data frame to be saved.
#' @param suffix A string to be appended to the file path
#' (e.g., "_metadata", "_profiles").
#' @param file_type The format in which to save the file.
#' Should be either "csv" or "parquet".
#' @param chr_name A string representing the chromosome name
#' that will be appended to the file name.
#' @param file_path The base file path where the file will be saved,
#' excluding the chromosome name and file extension.
#' @importFrom arrow write_parquet
#' @importFrom utils write.csv
#'
save_to_file <- function(data, suffix, file_type, chr_name, file_path) {
  # Determine the full file path with the correct extension based on file_type
  extension <- ifelse(file_type == "csv", ".csv", ".parquet")
  full_file_path <- paste0(file_path, "_", chr_name, suffix, extension)

  # Save based on the file type
  if (file_type == "csv") {
    # Save as CSV
    utils::write.csv(as.data.frame(data), full_file_path, row.names = FALSE)
  } else if (file_type == "parquet") {
    # Save as Parquet
    arrow::write_parquet(as.data.frame(data), full_file_path)
  } else {
    stop("Unsupported file type. Use 'csv' or 'parquet'.")
  }
}
