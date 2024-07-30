#!/bin/bash

# Set envirioment variables


### 4 ### predict_profile_probabilities.py

CELL_TYPE="epc_poolSubK562"
MODEL_PREFIX="GM12878_wt10M"
MODEL_NAME="model_full_GM12878_wt10M.sav"
PROFILE_MAIN_DIR="epc_poolSub_profile"


MODEL_DIR="/home/zmk214/zmk214/model_PRIMEloci/2_train_model"
PREFIX_OUT_NAME=$CELL_TYPE"-on-"$MODEL_PREFIX"-model"

SCRIPT_DIR="."
OUTPUT_DIR="/home/zmk214/zmk214/model_PRIMEloci/5_eval_epcrispr_benchmark"


# Run scripts
python3 ../genomewide_prediction/_predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_DIR/$MODEL_NAME -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -n $PREFIX_OUT_NAME -r "pos" 
python3 ../genomewide_prediction/_predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_DIR/$MODEL_NAME -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -n $PREFIX_OUT_NAME -r "neg" 
