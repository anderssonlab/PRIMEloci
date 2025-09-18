#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(argparse)
  library(CAGEfightR)
  library(parallel)
  library(GenomicRanges)
  library(PRIME)
  library(PRIMEloci)
  library(future.apply)
  library(SummarizedExperiment)
}))

# ---- Argument Parsing ----
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir",
                    default = "./PRIMEloci_profiles",
                    help = "Path to the input directory containing the main PRIMEloci profiles (default: ./PRIMEloci_profiles)") # nolint: line_length_linter.
#parser$add_argument("-o", "--output_dir", default = "./",
#                    help = "Path to the output directory")
parser$add_argument("--profile_dir_name", default = "PRIMEloci_profiles",
                    help = "Name of the profile main directory (default: PRIMEloci_profiles)") # nolint: line_length_linter.

parser$add_argument("--python_path", default = "~/.virtualenvs/prime-env",
                    help = "Path to Python executable. If not provided, the system's default Python will be used.") # nolint: line_length_linter.

parser$add_argument("-m", "--model_path",
                    default = file.path(system.file("model", package = "PRIME"),
                                        "PRIMEloci_GM12878_model_1.0.sav"),
                    help = "Model full path")
parser$add_argument("--name_prefix", default = "PRIMEloci",
                    help = "Prefix for prediction output")

parser$add_argument("-p", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores to use for parallel processing")

args <- parser$parse_args()

# ---- Setup Directories ----
# Setup
profile_main_dir <- args$input_dir
profiles_subtnorm_dir <- file.path(profile_main_dir, "profiles_subtnorm")

if (!dir.exists(profiles_subtnorm_dir)) {
  plc_error("❌ Directory 'profiles_subtnorm' does not exist in the specified profile_main_dir.") # nolint: line_length_linter.
}

profile_files <- list.files(profiles_subtnorm_dir,
                            pattern = "\\.(npz|parquet|csv)$")
assertthat::assert_that(length(profile_files) > 0,
                        msg = paste("❌ No profile files found in:",
                                    profiles_subtnorm_dir))

model_path <- args$model_path
predict_script_path <- file.path(system.file("python", package = "PRIMEloci"),
                                 "main.py")

assertthat::assert_that(file.exists(predict_script_path),
                        msg = paste("❌ Prediction script not found at:",
                                    predict_script_path))
assertthat::assert_that(file.exists(model_path),
                        msg = paste("❌ Model file not found at:",
                                    model_path))

name_prefix <- args$name_prefix

num_cores <- args$num_cores

# Assertions
if (!is.null(num_cores)) {
  assertthat::assert_that(
    is.numeric(num_cores),
    num_cores %% 1 == 0,
    num_cores > 0,
    msg = "`num_cores` must be a positive integer or NULL."
  )
}

if (is.null(num_cores)) {
  num_cores <- max(1, min(25, parallel::detectCores() %/% 2))
}
if (num_cores == 1) {
  processing_method <- "callr"
  plc_message("⚠️ num_workers was set to 1. Using callr backend: tasks will run sequentially (despite using multiple R sessions).") # nolint: line_length_linter.
} else {
  processing_method <- PRIMEloci:::plc_detect_parallel_plan()
}

# Python config
plc_message("🚀 Setting up Python environment")

if (is.null(args$python_path)) {
  py <- reticulate::import("sys")
  python_path <- py$executable
} else {
  python_path <- args$python_path
}
py_conf <- PRIMEloci:::plc_configure_python(python_path = python_path)


plc_message("🚀 Running PRIMEloci -5: Prediction using PRIMEloci model")

prediction_cmd <- c(
  python_path, predict_script_path,
  "--script_dir", dirname(predict_script_path),
  "--profile_main_dir", profile_main_dir,
  "--combined_outdir", dirname(profile_main_dir),
  "--model_path", model_path,
  "--log_file", file.path(profile_main_dir, "PRIMEloci_prediction.log"),
  "--name_prefix", name_prefix
)

if (!is.null(num_cores)) {
  assertthat::assert_that(is.numeric(num_cores), num_cores > 0)
  prediction_cmd <- c(prediction_cmd, "--num_core", as.character(num_cores))
}

plc_message(paste("🔧 Python command:",
                  paste(shQuote(prediction_cmd), collapse = " ")))

plc_message("🔹 Running Python prediction script...")
result <- tryCatch(
  {
    output <- system2(python_path,
                      args = prediction_cmd[-1],
                      stdout = TRUE,
                      stderr = TRUE)
    attr(output, "status") <- 0
    output
  },
  error = function(e) {
    msg <- paste("❌ ERROR during prediction execution:", e$message)
    plc_message(msg)
    attr(msg, "status") <- 1
    msg
  }
)

if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
  plc_error("❌ Prediction script failed. Check PRIMEloci.log for details.")
} else {
  plc_message("✅ DONE :: Prediction script executed successfully.")
}
