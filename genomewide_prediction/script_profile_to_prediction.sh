#!/bin/bash

# Set envirioment variables


### 4 ### predict_profile_probabilities.py

CELL_TYPE="HCT116"

MODEL_NAME="model_Meena_v3_k562_fullModel_C.sav"
MODEL_DIR="../data/resources/models"
THRESHOLD=0.2

SCRIPT_DIR="."
OUTPUT_DIR="../data/results"
PROFILE_MAIN_DIR=$CELL_TYPE"_tc_profiles"
PROFILE_SUB_DIR="tcs"


# Run scripts
python3 _predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_DIR/$MODEL_NAME -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -t $THRESHOLD
