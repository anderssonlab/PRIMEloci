#!/usr/bin/env Rscript

writeLines("\n### Running PRIMEloci Step 5: Prediction from profiles ###")

suppressPackageStartupMessages({
  library(argparse)
  library(assertthat)
  library(PRIME)
})

# ---- Argument Parsing ----
parser <- ArgumentParser()

parser$add_argument("-o", "--output_dir", required = TRUE,
                    help = "Path to the main output directory (where PRIMEloci_tmp lives)")
parser$add_argument("--profile_dir_name", default = "PRIMEloci_profiles",
                    help = "Name of the profile main directory (default: PRIMEloci_profiles)")
parser$add_argument("--python_path", default = "~/.virtualenvs/prime-env/bin/python",
                    help = "Full path to the Python executable")
parser$add_argument("-m", "--model_name", default = "PRIMEloci_GM12878_model_1.0.sav",
                    help = "Model file name (inside the PRIME R package 'model' folder)")
parser$add_argument("--name_prefix", default = "PRIMEloci",
                    help = "Prefix for prediction output")
parser$add_argument("-l", "--log", default = NULL,
                    help = "Log file name or NULL to log to console")
parser$add_argument("-n", "--num_cores", type = "integer", default = NULL,
                    help = "Number of cores for prediction. Default is NULL")

args <- parser$parse_args()

# ---- Setup Directories and Logging ----
output_dir <- create_output_dir(args$output_dir)

primeloci_tmp <- setup_tmp_dir(output_dir)

log <- if (is.null(args$log) || args$log == "NULL") NULL else args$log
log_target <- setup_log_target(log, output_dir)

# ---- Derived Paths ----
profile_main_dir <- file.path(primeloci_tmp, args$profile_dir_name)
profiles_subtnorm_dir <- file.path(profile_main_dir, "profiles_subtnorm")
python_path <- path.expand(args$python_path)
model_path <- file.path(system.file("model", package = "PRIME"), args$model_name)

name_prefix <- args$name_prefix
num_cores <- args$num_cores

# ---- Input Checks ----
profile_files <- list.files(profiles_subtnorm_dir, pattern = "\\.(npz|parquet|csv)$")
assert_that(length(profile_files) > 0,
            msg = paste("❌ No profile files found in:", profiles_subtnorm_dir))

model_path <- file.path(system.file("model", package = "PRIME"), model_name)
predict_script_path <- file.path(system.file("python", package = "PRIME"), "main.py")

assert_that(file.exists(predict_script_path),
            msg = paste("❌ Prediction script not found at:", predict_script_path))
assert_that(file.exists(model_path),
            msg = paste("❌ Model file not found at:", model_path))

# ---- Python Configuration ----
py_conf <- PRIME::configure_plc_python(python_path = python_path,
                                       log_target = log_target)

# ---- Logging and Execution ----
plc_log("\n\n\n 🚀 Running PRIMEloci: Prediction using PRIMEloci model", log_target)

prediction_cmd <- c(
  python_path, predict_script_path,
  "--script_dir", dirname(predict_script_path),
  "--profile_main_dir", profile_main_dir,
  "--combined_outdir", dirname(profile_main_dir),
  "--model_path", model_path,
  "--log_file", if (is.character(log_target)) log_target else "stdout",
  "--name_prefix", name_prefix
)

if (!is.null(num_cores)) {
  assert_that(is.numeric(num_cores), num_cores > 0)
  prediction_cmd <- c(prediction_cmd, "--num_core", as.character(num_cores))
}

plc_log(paste("🔧 Python command:",
              paste(shQuote(prediction_cmd), collapse = " ")),
        log_target, print_console = FALSE)

plc_log("🔹 Running Python prediction script...", log_target)
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
    plc_log(msg, log_target, level = "❌ ERROR")
    attr(msg, "status") <- 1
    msg
  }
)

if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
  stop("❌ Prediction script failed. Check the log for details.")
} else {
  plc_log("✅ DONE :: Prediction script executed successfully.", log_target)
}
