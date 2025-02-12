#!/bin/bash

# Stop the script if any command fails
set -e

# Default configuration file as an empty string (user must provide one)
CONFIG_FILE=""

# Default value for EXTENSION_DISTANCE
EXTENSION_DISTANCE=200

# Usage message
usage() {
    echo "Usage: $0 --config <config_file> [-1] [-2] [-3] [-4] [-5] [-6] [--all] [--pred] [--keeptmp]"
    echo "  --config   Specify the bash configuration file to use (required for all commands)"
    echo "  -1         Run step 1: _1_get_ctss_from_bw.r"
    echo "  -2         Run step 2: _2_get_tc_from_ctss.r"
    echo "  -3         Run step 3: _3_get_sld_window_from_tc.r"
    echo "  -4         Run step 4: _4_get_profile.r"
    echo "  -5         Run step 5: _5_predict_profile_probability.py"
    echo "  -6         Run step 6: _6_apply_post_processing_coreovlwith-d.r"
    echo "  --all      Run all steps"
    echo "  --pred     Run only steps 3-6 as ctss and tc rds objects are provided"
    echo "  --keeptmp  Keep temporary files created in steps 3-5"
    exit 1
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

# Initialize flags
all=false
pred=false
keeptmp=false
steps=()

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --config)
            shift
            CONFIG_FILE="$1"
            ;;
        -1) steps+=("1") ;;
        -2) steps+=("2") ;;
        -3) steps+=("3") ;;
        -4) steps+=("4") ;;
        -5) steps+=("5") ;;
        -6) steps+=("6") ;;
        --all) all=true ;;
        --pred) pred=true ;;
        --keeptmp) keeptmp=true ;;
        *) echo "Invalid option: $1" >&2; usage ;;
    esac
    shift
done

# Check if the config file is provided
if [[ -z "$CONFIG_FILE" ]]; then
    echo "Error: The --config option is required for all commands."
    exit 1
fi

# Check if the configuration file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found."
    exit 1
fi

# Load the specified configuration file
source "$CONFIG_FILE"

# Check if output directory is provided, if not, create one
if [[ ! -d $OUTPUT_DIR ]]; then
    echo "Creating output directory..."
    mkdir -p $OUTPUT_DIR
fi

# Function to create /PRIMEloci_tmp directory if needed
create_tmp_dir() {
    if [[ ! -d "$OUTPUT_DIR/PRIMEloci_tmp" ]]; then
        echo "Creating /PRIMEloci_tmp directory inside $OUTPUT_DIR..."
        mkdir -p "$OUTPUT_DIR/PRIMEloci_tmp"
    fi
    TMP_DIR="$OUTPUT_DIR/PRIMEloci_tmp"
}

# Function to combine BED files based on prefix before "_chr"
# Takes two arguments: input directory and output directory
combine_bed_files() {
  local input_dir="$1"
  local output_dir="$2"

  # Ensure the output directory exists
  mkdir -p "$output_dir"

  # Step 1: Find all unique prefixes before "_chr"
  prefixes=$(ls "$input_dir"/*.bed | sed -E 's/(.*)(_chr[^_]+).*/\1/' | sort -u)

  # Step 2: For each unique prefix, concatenate all matching files
  for prefix in $prefixes; do
    # Find all files that match the prefix (with any chr)
    matching_files=$(ls "$input_dir"/$(basename "$prefix")_chr*.bed)

    # Combine those files into one file in the output directory
    combined_file="$output_dir/$(basename "$prefix")_combined.bed"

    # Handle headers: Include header from the first file, skip for others
    head -n 1 $(echo $matching_files | cut -d ' ' -f1) > "$combined_file"  # Extract header from the first file
    for file in $matching_files; do
      tail -n +2 "$file" >> "$combined_file"  # Skip header (start from line 2)
    done

    echo "Combined files into $combined_file"
  done

}

# If --all is specified, add all steps 1-6
if $all; then
    steps=("1" "2" "3" "4" "5" "6")
fi

# If --pred is specified, add only steps 3-6
if $pred; then
    steps=("3" "4" "5" "6")
fi

# Run the specified steps
for step in "${steps[@]}"; do
    case $step in
        1)
            echo -e "\nRunning _1_get_ctss_from_bw.r"
            Rscript _1_get_ctss_from_bw.r \
                -i $CAGE_DIR \
                -m $DESIGN_MATRIX \
                -o $OUTPUT_DIR \
                -n $CTSS_RSE_NAME \
                -k
            ;;
        2)
            echo -e "\nRunning _2_get_tc_from_ctss.r"
            Rscript _2_get_tc_from_ctss.r \
                -i $OUTPUT_DIR/$CTSS_RSE_NAME \
                -o $OUTPUT_DIR \
                -n $TC_GRL_NAME \
                -e $EXTENSION_DISTANCE
            ;;
        3)
            # Create /PRIMEloci_tmp directory before running step 3
            create_tmp_dir
            echo -e "\nRunning _3_get_sld_window_from_tc.r"
            Rscript _3_get_sld_window_from_tc.r \
                -i $OUTPUT_DIR/$TC_GRL_NAME \
                -o $TMP_DIR \
                -n $SLD_TC_GRL_NAME \
                -s $SLD_WINDOW \
                -e $EXTENSION_DISTANCE
            ;;
        4)
            echo -e "\nRunning _4_get_tc_profile.r"
            Rscript _4_get_profile.r \
                -c $OUTPUT_DIR/$CTSS_RSE_NAME \
                -t $OUTPUT_DIR/PRIMEloci_tmp/$SLD_TC_GRL_NAME \
                -o $OUTPUT_DIR/PRIMEloci_tmp \
                -n $PROFILE_MAIN_DIR \
                -f $PROFILE_FILE_TYPE
            ;;
        5)
            echo -e "\nRunning _5_predict_profile_probability.py"
            python3 _5_predict_profile_probability.py \
                -w $SCRIPT_DIR \
                -m $MODEL_PATH \
                -p $OUTPUT_DIR/PRIMEloci_tmp/$PROFILE_MAIN_DIR \
                -n $PREFIX_OUT_NAME \
                -f $PROFILE_FILE_TYPE $CALIBRATION_FLAG
            combine_bed_files $OUTPUT_DIR/PRIMEloci_tmp/$PROFILE_MAIN_DIR/predictions $OUTPUT_DIR 
            ;;
        6)
            echo -e "\nRunning _6_apply_post_processing_coreovlwith-d.r"
            for FILE in $(find "$OUTPUT_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE ..."
                Rscript _6_apply_post_processing_coreovlwith-d.r \
                    -i "$FILE" \
                    -o $OUTPUT_DIR \
                    -t $THRESHOLD \
                    -d $SCORE_DIFF \
                    -m
                done
            ;;
    esac
done

# Function to clean up the temporary directory if needed
cleanup() {
    if ! $keeptmp && [[ -d "$TMP_DIR" ]]; then
        echo "Cleaning up temporary directory: $TMP_DIR"
        rm -rf "$TMP_DIR"
    fi
}

# Clean up the temporary directory only if running steps 3-5
if ($all || $pred) && ! $keeptmp; then
    echo "Cleaning up temporary files from steps 3-5"
    rm -rf $TMP_DIR
fi

# Trap to clean up if script is interrupted
# trap cleanup EXIT
