# config.sh

### 1 ### _1_get_ctss_from_bw.r
CAGE_DIR="../example/resources/cage_bw"
DESIGN_MATRIX="../example/resources/design_matrix_k562.tsv"
OUTPUT_DIR="../example/results"
CTSS_RSE_NAME="ctss_rse.rds"
# -k was set for keeping only standard chromosomes

### 2 ### _2_get_tc_from_ctss.r
# EXTENSION_DISTANCE=200 was fixed in the script
# OUTPUT_DIR and CTSS_RSE_NAME from above
TC_GRL_NAME="tc_grl.rds"

### 3 ### _3_get_sld_window_from_tc.r
# OUTPUT_DIR and TC_GRL_NAME from above
SLD_TC_GRL_NAME="tc_grl_sld.rds"
SLD_WINDOW=20

### 4 ### _4_get_profile.r
# CTSS_RSE_NAME, TC_GRL_NAME, OUTPUT_DIR from above
PROFILE_MAIN_DIR="example_tc_profiles"
PROFILE_FILE_TYPE="parquet"
# add -s if you want to save count profiles

### 5 ### _5_predict_profile_probability.py
# OUTPUT_DIR, PROFILE_MAIN_DIR, PROFILE_SUB_DIR, and PROFILE_FILE_TYPE from above
SCRIPT_DIR="."
MODEL_PATH="../model/PRIMEloci_GM12878_model_1.0.sav"
PREFIX_OUT_NAME="K562-on-PRIMEloci"

## 6 ### _6_apply_post_processing_coreovlwith-d.r
# OUTPUT_DIR from above
PARTIAL_NAME="*pred_all*_combined.bed"
THRESHOLD=0.75
SCORE_DIFF=0.10
