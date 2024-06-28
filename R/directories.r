#' Prepare Directory Structure for Profile Output
#'
#' This function creates a specified directory structure
#' for storing profile output data.
#' It generates a main output directory
#' containing several subdirectories for different types of profile data.
#' If the specified main output directory already exists,
#' the function will not recreate it or its subdirectories.
#'
#' @param output_dir A character string specifying the path
#' to the output directory. Defaults to the current directory (".").
#' @param output_dir_name A character string specifying
#' the name of the main output directory. Defaults to "profile_output".
#' @param output_main_name A character vector specifying
#' the names of the main subdirectories to be created
#' within the main output directory. Defaults to
#' c("metadata", "metadata_subtnorm", "profiles", "profiles_subtnorm").
#' @param output_subdir_name A character vector specifying
#' the names of the subdirectories to be created within each main subdirectory.
#' Defaults to c("pos", "neg").
#'
#' @return The function does not return any value.
#' It creates the specified directory structure on the file system.
#'
#' @export
prep_profile_dir <- function(output_dir = ".",
                             output_dir_name = "profile_output",
                             output_main_name = c("metadata",
                                                  "metadata_subtnorm",
                                                  "profiles",
                                                  "profiles_subtnorm"),
                             output_subdir_name = c("pos", "neg")) {

  # Create the output dir and main dir
  new_path <- file.path(output_dir, output_dir_name)

  if (!file.exists(new_path)) {
    dir.create(new_path)

    # Create main directories and their subdirectories
    lapply(output_main_name, function(main) {
      main_path <- file.path(new_path, main)
      dir.create(main_path)

      if (length(output_subdir_name) > 0) {
        lapply(output_subdir_name, function(subdir) {
          subdir_path <- file.path(main_path, subdir)
          dir.create(subdir_path)
        })
      }
    })

    cat("New folder created:", new_path, "\n")
  } else {
    cat("Folder already exists:", new_path, "\n")
  }
}