#!/bin/bash

# Set envirioment variables


### 1 ### get_ctss_from_bw.R

CAGE_DIR="../example/resources/cage_bw"
DESIGN_MATRIX="../example/resources/design_matrix_k562.tsv"

OUTPUT_DIR="../example/results"
CTSS_RSE_NAME="ctss_rse.rds"
# add -k for keeping only standard chromosomes


### 2 ### get_tc_grl.R

# OUTPUT_DIR and CTSS_RSE_NAME from above
TC_GRL_NAME="tc_grl.rds"
EXTENSION_DISTANCE=200


### 3 ### get_tc_profiles.R

# CTSS_RSE_NAME, TC_GRL_NAME, OUTPUT_DIR from above
PROFILE_MAIN_DIR="example_tc_profiles"
PROFILE_SUB_DIR="tcs"
# add -s if you want to save count profiles


### 4 ### predict_profile_probabilities.py

# OUTPUT_DIR, PROFILE_MAIN_DIR, and PROFILE_SUB_DIR from above
SCRIPT_DIR="."
MODEL_PATH="../example/resources/model_Meena_v3_gm12878_fullModel_C.sav"

THRESHOLD=0.2


# Run scripts
Rscript _get_ctss_from_bw.r -i $CAGE_DIR -m $DESIGN_MATRIX -o $OUTPUT_DIR -c $CTSS_RSE_NAME -k
Rscript _get_tc_from_ctss.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -o $OUTPUT_DIR -t $TC_GRL_NAME -e $EXTENSION_DISTANCE
Rscript _get_tc_profiles.r -c $OUTPUT_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR
python3 _predict_profile_probabilities.py -w $SCRIPT_DIR -m $MODEL_PATH -p $OUTPUT_DIR/$PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR -t $THRESHOLD
