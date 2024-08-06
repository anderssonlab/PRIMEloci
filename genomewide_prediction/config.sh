# config.sh

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
PREFIX_OUT_NAME="K562-on-GM12878-C-model"
THRESHOLD=0.2


### 5 ### predict_profile_probabilities.py

# OUTPUT_DIR, PROFILE_MAIN_DIR, PROFILE_SUB_DIR from above
PREDICTION_DIR=$OUTPUT_DIR/$PROFILE_MAIN_DIR/predictions/$PROFILE_SUB_DIR
PARTIAL_NAME="*pred_slt*.bed"