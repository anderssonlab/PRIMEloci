#!/bin/bash

# Stop the script if any command fails
set -e

# Load the configuration
source ./bash_config_PRIMEloci.sh

# Usage message
usage() {
    echo "Usage: $0 [-1] [-2] [-3] [-4] [-5] [-p] [--all]"
    echo "  -1        Run _get_ctss_from_bw.R"
    echo "  -2        Run _get_tc_grl.R"
    echo "  -3        Run _get_tc_profiles.R"
    echo "  -4        Run _predict_profile_probabilities.py"
    echo "  -5        Run _filter_bed_to_reduced_gr.R"
    echo "  -p        Run steps 3-5 together"
    echo "  --all     Run all steps (1-5)"
    exit 1
}
# To run all steps: ./script.sh --all
# To run steps 1, 2, and group 3-5 together: ./script.sh -1 -2 -p
# To run specific steps individually, e.g., 1, 3, and 4: ./script.sh -1 -3 -4

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

run_all=false
group_3_5=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -1)
            steps+=("1")
            ;;
        -2)
            steps+=("2")
            ;;
        -3)
            steps+=("3")
            ;;
        -4)
            steps+=("4")
            ;;
        -5)
            steps+=("5")
            ;;
        -p)
            group_3_5=true
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

if $run_all; then
    steps=("1" "2" "3" "4" "5")
elif $group_3_5; then
    steps+=("3" "4" "5")
fi

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
            echo -e "\nRunning _get_tc_profiles.R"
            Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -f $PROFILE_FILE_TYPE
            ;;
        4)
            echo -e "\nRunning _predict_profile_probabilities.py"
            python3 _predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_PATH -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -n $PREFIX_OUT_NAME -t $THRESHOLD -f $PROFILE_FILE_TYPE
            ;;
        5)
            echo -e "\nRunning _filter_bed_to_reduced_gr.R" 
            for FILE in $(find "$PREDICTION_DIR" -type f -name $PARTIAL_NAME); do
                echo "Processing $FILE..."
                Rscript _filter_bed_to_reduced_gr.r -i "$FILE" -o $PREDICTION_DIR
            done
            ;;
    esac
done