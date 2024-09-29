#!/bin/bash

# Stop the script if any command fails
set -e

# Load the configuration
source ./bash_config_PRIMEloci.sh

# Usage message
usage() {
    echo "Usage: $0 [--tc] [--sld] [--all --tc] [--all --sld]"
    echo "  --tc      Run TC steps: 4, 5, 6 (input for 4 from step 2)"
    echo "  --sld     Run Sliding window steps: 3, 4, 5, 7 (input for 4 from step 3)"
    echo "  --all --tc   Run all steps: 1, 2, then TC steps (4, 5, 6)"
    echo "  --all --sld  Run all steps: 1, 2, then Sliding window steps (3, 4, 5, 7)"
    exit 1
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

run_all=false
use_tc=false
use_sld=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tc)
            use_tc=true
            ;;
        --sld)
            use_sld=true
            ;;
        --all)
            run_all=true
            ;;
        *)
            echo "Invalid option: $1" >&2
            usage
            ;;
    esac
    shift
done

# If --all is specified, add steps 1 and 2 first
if $run_all; then
    steps=("1" "2")
    if $use_tc; then
        steps+=("4" "5" "6")
    elif $use_sld; then
        steps+=("3" "4" "5" "7")
    else
        echo "Error: You must specify either --tc or --sld with --all."
        exit 1
    fi
# If --tc or --sld is specified without --all
elif $use_tc; then
    steps=("4" "5" "6")
elif $use_sld; then
    steps=("3" "4" "5" "7")
fi

# Run the specified steps
for step in "${steps[@]}"; do
    case $step in
        1)
            echo -e "\nRunning _get_ctss_from_bw.R"
            Rscript _get_ctss_from_bw.r -i $CAGE_DIR -m $DESIGN_MATRIX -o $OUTPUT_DIR -c $CTSS_RSE_NAME -k
            ;;
        2)
            echo -e "\nRunning _get_tc_grl.R"
            Rscript _get_tc_from_ctss.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -o $OUTPUT_DIR -t $TC_GRL_NAME -e $EXTENSION_DISTANCE
            ;;
        3)
            echo -e "\nRunning sliding window on tc_grl.R"
            Rscript sliding_window_script.R --input $OUTPUT_DIR/$TC_GRL_NAME --output $OUTPUT_DIR/tc_grl_sld.rds --slide_by 20 --expand_by 200
            ;;
        4)
            if $use_tc; then
                echo -e "\nRunning _get_tc_profiles.R (TC input)"
                Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -f $PROFILE_FILE_TYPE
            elif $use_sld; then
                echo -e "\nRunning _get_tc_profiles.R (SLD input)"
                Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/tc_grl_sld.rds -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -f $PROFILE_FILE_TYPE
            fi
            ;;
        5)
            echo -e "\nRunning _predict_profile_probabilities.py"
            python3 _predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_PATH -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -n $PREFIX_OUT_NAME -t $THRESHOLD -f $PROFILE_FILE_TYPE
            ;;
        6)
            echo -e "\nRunning _filter_bed_to_reduced_gr.R (Method 1)"
            for FILE in $(find "$PREDICTION_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE..."
                Rscript _filter_bed_to_reduced_gr.r -i "$FILE" -o $PREDICTION_DIR
            done
            ;;
        7)
            echo -e "\nRunning _filter_bed_to_reduced_gr.R (Method 2)"
            for FILE in $(find "$PREDICTION_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE with Method 2..."
                Rscript _filter_bed_to_reduced_gr_method2.r -i "$FILE" -o $PREDICTION_DIR
            done
            ;;
    esac
done
