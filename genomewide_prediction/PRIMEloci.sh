#!/bin/bash

# Stop the script if any command fails
set -e

# Load the configuration
source ./bash_config_PRIMEloci.sh

# Usage message
usage() {
    echo "Usage: $0 [-1] [-2] [-3] [-4] [-5] [-6] [-7] [--tc] [--sld] [--all --tc] [--all --sld]"
    echo "  -1       Run step 1: _get_ctss_from_bw.r"
    echo "  -2       Run step 2: _get_tc_from_ctss.r"
    echo "  -3       Run step 3: sliding window on _get_sld_windows.r"
    echo "  -4       Run step 4: _get_tc_profiles.r"
    echo "  -5       Run step 5: _predict_profile_probabilities.py"
    echo "  -6       Run step 6: _filter_max_nonoverlapping.r"
    echo "  -7       Run step 7: _filter_core_overlapping.r"
    echo "  --tc     Run TC steps: 4, 5, 6 (input for 4 from step 2)"
    echo "  --sld    Run Sliding window steps: 3, 4, 5, 7 (input for 4 from step 3)"
    echo "  --all --tc   Run all steps: 1, 2, then TC steps (4, 5, 6)"
    echo "  --all --sld  Run all steps: 1, 2, then Sliding window steps (3, 4, 5, 7)"
    exit 1
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

# Initialize flags
run_all=false
use_tc=false
use_sld=false
steps=()

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -1) steps+=("1") ;;
        -2) steps+=("2") ;;
        -3) steps+=("3") ;;
        -4) steps+=("4") ;;
        -5) steps+=("5") ;;
        -6) steps+=("6") ;;
        -7) steps+=("7") ;;
        --tc) use_tc=true ;;
        --sld) use_sld=true ;;
        --all) run_all=true ;;
        *) echo "Invalid option: $1" >&2; usage ;;
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
fi

# Run the specified steps
for step in "${steps[@]}"; do
    case $step in
        1)
            echo -e "\nRunning _get_ctss_from_bw.r"
            Rscript _get_ctss_from_bw.r -i $CAGE_DIR -m $DESIGN_MATRIX -o $OUTPUT_DIR -c $CTSS_RSE_NAME -k
            ;;
        2)
            echo -e "\nRunning _get_tc_from_ctss.r"
            Rscript _get_tc_from_ctss.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -o $OUTPUT_DIR -t $TC_GRL_NAME -e $EXTENSION_DISTANCE
            ;;
        3)
            echo -e "\nRunning sliding window on _get_sld_windows.r"
            Rscript _get_sld_windows.r -i $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR/$SLD_TC_GRL_NAME -s $SLD_WINDOW -e $EXTENSION_DISTANCE
            ;;
        4)
            if $use_tc; then
                echo -e "\nRunning _get_tc_profiles.r (TC input)"
                Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -f $PROFILE_FILE_TYPE
            elif $use_sld; then
                echo -e "\nRunning _get_tc_profiles.r (SLD input)"
                Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$SLD_TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -f $PROFILE_FILE_TYPE
            fi
            ;;
        5)
            echo -e "\nRunning _predict_profile_probabilities.py"
            python3 _predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_PATH -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -n $PREFIX_OUT_NAME -f $PROFILE_FILE_TYPE
            ;;
        6)
            echo -e "\nRunning _filter_max_nonoverlapping.r (Method 1)"
            for FILE in $(find "$PREDICTION_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE..."
                Rscript _filter_max_nonoverlapping.r -i "$FILE" -o $PREDICTION_DIR -t $THRESHOLD
            done
            ;;
        7)
            echo -e "\nRunning _filter_core_overlapping.r (Method 2)"
            for FILE in $(find "$PREDICTION_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE with Method 2..."
                Rscript _filter_core_overlapping.r -i "$FILE" -o $PREDICTION_DIR -t $THRESHOLD
            done
            ;;
    esac
done
