#!/bin/bash

# Stop the script if any command fails
set -e

# Default configuration file as an empty string (user must provide one)
CONFIG_FILE=""

# Usage message
usage() {
    echo "Usage: $0 --config <config_file> [-1] [-2] [-3] [-4] [-5] [-6] [--PRIMEloci] [--PRIMEloci_facet] [--keeptmp] [--num_cores <num_cores>]"
    echo "  --config            Specify the bash configuration file to use (required for all commands)"
    echo "  -0                  Run the check for input files for PRIMEloci_facet"
    echo "  -1                  Run step 1: _1_get_ctss_from_bw.r"
    echo "  -2                  Run step 2: _2_get_tc_from_ctss.r"
    echo "  -3                  Run step 3: _3_get_sld_window_from_tc.r"
    echo "  -4                  Run step 4: _4_get_profile.r"
    echo "  -5                  Run step 5: _5_predict_profile_probability.py"
    echo "  -6                  Run step 6: _6_apply_post_processing_coreovlwith-d.r"
    echo "  --PRIMEloci         Run all steps"
    echo "  --PRIMEloci_facet   Run only steps 0, 4, and 5 as ctss and tc rds objects are provided"
    echo "  --keeptmp           Keep temporary files created in steps 3-5, only if running PRIMEloci and PRIMEloci_facet"
    echo "  --num_cores         Number of cores to use for parallel processing"
    exit 1
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

# Initialize flags
PRIMEloci=false
PRIMEloci_facet=false
keeptmp=false
steps=()

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --config)
            shift
            CONFIG_FILE="$1"
            ;;
        -0) steps+=("0") ;;
        -1) steps+=("1") ;;
        -2) steps+=("2") ;;
        -3) steps+=("3") ;;
        -4) steps+=("4") ;;
        -5) steps+=("5") ;;
        -6) steps+=("6") ;;
        --PRIMEloci) PRIMEloci=true ;;
        --PRIMEloci_facet) PRIMEloci_facet=true ;;
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
#if [[ ! -d $OUTPUT_DIR ]]; then
#    echo "Creating output directory..."
#    mkdir -p $OUTPUT_DIR
#fi

# Function to create /PRIMEloci_tmp directory if needed
#create_tmp_dir() {
#    if [[ ! -d "$OUTPUT_DIR/PRIMEloci_tmp" ]]; then
#        echo "Creating /PRIMEloci_tmp directory inside $OUTPUT_DIR..."
#        mkdir -p "$OUTPUT_DIR/PRIMEloci_tmp"
#    fi
#    TMP_DIR="$OUTPUT_DIR/PRIMEloci_tmp"
#}
TMP_DIR="$OUTPUT_DIR/PRIMEloci_tmp"

ARGS=""
if [ -n "$NUM_CORES" ]; then
  ARGS="$ARGS -p $NUM_CORES"
fi

# Function to combine BED files based on prefix before "_chr"
# Takes two arguments: input directory and output directory
#combine_bed_files() {
#  local input_dir="$1"
#  local output_dir="$2"
#
#  # Ensure the output directory exists
#  mkdir -p "$output_dir"
#
#  # Step 1: Find all unique prefixes before "_chr"
#  prefixes=$(ls "$input_dir"/*.bed | sed -E 's/(.*)(_chr[^_]+).*/\1/' | sort -u)
#
#  # Step 2: For each unique prefix, concatenate all matching files
#  for prefix in $prefixes; do
#    # Find all files that match the prefix (with any chr)
#    matching_files=$(ls "$input_dir"/$(basename "$prefix")_chr*.bed)
#
#    # Combine those files into one file in the output directory
#    combined_file="$output_dir/$(basename "$prefix")_combined.bed"
#
#    # Handle headers: Include header from the first file, skip for others
#    head -n 1 $(echo $matching_files | cut -d ' ' -f1) > "$combined_file"  # Extract header from the first file
#    for file in $matching_files; do
#      tail -n +2 "$file" >> "$combined_file"  # Skip header (start from line 2)
#    done
#
#    echo "Combined files into $combined_file"
#  done
#
#}

# If --PRIMEloci is specified, add all steps 1-6
if $PRIMEloci; then
    steps=("1" "2" "3" "4" "5" "6")
fi

# If --PRIMEloci_facet is specified, add only steps 0, 3-5
if $PRIMEloci_facet; then
    steps=("0" "4" "5")
fi

# Run the specified steps
for step in "${steps[@]}"; do
    case $step in
        0)
            echo -e "\nRunning check for input files for PRIMEloci_facet"
            Rscript _0_validate_ctss_and_region.r \
                --ctss_rse $CTSS_RSE_RDS \
                --region $REGION_RDS \
                -o $OUTPUT_DIR \
            ;;
        1)
            echo -e "\nRunning _1_get_ctss_from_bw.r"
            Rscript _1_get_ctss_from_bw.r \
                -i $CAGE_BW_DIR \
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
                -n $TC_GRL_NAME
            ;;
        3)
            # Create /PRIMEloci_tmp directory before running step 3
            #create_tmp_dir
            echo -e "\nRunning _3_get_sld_window_from_tc.r"
            Rscript _3_get_sld_window_from_tc.r \
                -i $OUTPUT_DIR/$TC_GRL_NAME \
                -o $TMP_DIR \
                -n $SLD_TC_GRL_NAME \
                -s $SLD_WINDOW \
                $ARGS \
            ;;
        4)
            echo -e "\nRunning _4_get_tc_profile.r"
            if [ "$PRIMEloci_facet" = true ]; then
                # If running PRIMEloci_facet, use the provided ctss_rse and region
                CTSS_RSE="$CTSS_RSE_RDS"
                REGION="$REGION_RDS"
            elif [ "$PRIMEloci" = true ]; then
                # If running PRIMEloci, use the files from the previous steps
                CTSS_RSE="$OUTPUT_DIR/$CTSS_RSE_NAME"
                REGION="$TMP_DIR/$SLD_TC_GRL_NAME"
            else
                # -4 was specified directly; fallback to config values if available
                if [ -n "$CTSS_RSE_RDS" ] && [ -n "$REGION_RDS" ]; then
                    CTSS_RSE="$CTSS_RSE_RDS"
                    REGION="$REGION_RDS"
                else
                    # Fallback to previous step outputs
                    CTSS_RSE="$OUTPUT_DIR/$CTSS_RSE_NAME"
                    REGION="$TMP_DIR/$SLD_TC_GRL_NAME"
                fi
            fi
            Rscript _4_get_profile.r \
                --ctss_rse $CTSS_RSE \
                --region $REGION \
                --python_path $PYTHON_PATH \
                -o $OUTPUT_DIR \
                --profile_dir_name $PROFILE_MAIN_DIR \
                -f $PROFILE_FORMAT \
                $ARGS \
            ;;
        5)
            echo -e "\nRunning _5_predict_profile_probability.r"

            #if [ "$PRIMEloci" = true ]; then
            #    # If running PRIMEloci, use the files from the previous steps
            #    IN="$TMP_DIR"
            #else
            #    IN="$OUTPUT_DIR"
            #fi

            Rscript _5_predict_profile_probability.r \
                -i $TMP_DIR/$PROFILE_MAIN_DIR \
                --python_path $PYTHON_PATH \
                -m $MODEL_PATH \
                --name_prefix $PREFIX_OUT_NAME \
                $ARGS \
            ;;
        6)
            echo -e "\nRunning _6_apply_post_processing_coreovlwith-d.r"
            Rscript _6_apply_post_processing_coreovlwith-d.r \
                -i $TMP_DIR \
                --partial_name $PARTIAL_NAME \
                -o $OUTPUT_DIR \
                -t $THRESHOLD \
                -d $SCORE_DIFF \
                --core_width $CORE_WIDTH \
                $ARGS \
            ;;
    esac
done

# Function to clean up the temporary directory if needed
#cleanup() {
#    if ! $keeptmp && [[ -d "$TMP_DIR" ]]; then
#        echo "Cleaning up temporary directory: $TMP_DIR"
#        rm -rf "$TMP_DIR"
#    fi
#}

# Clean up the temporary directory only if running steps 3-5
if ($PRIMEloci || $PRIMEloci_facet) && ! $keeptmp; then
    echo "Cleaning up temporary directory: $TMP_DIR"
    rm -rf $TMP_DIR
fi
