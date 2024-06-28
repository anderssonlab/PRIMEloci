#!/bin/bash

# Set envirioment variables


### 2 ### get_tc_grl.R

# OUTPUT_DIR and CTSS_RSE_NAME from above
CELL_TYPE="HCT116"

RESOURCES_DIR="../data/resources"
CTSS_RSE_NAME="poolSub."$CELL_TYPE".CTSSs.rmSingletons.rds"

OUTPUT_DIR="../data/results"
TC_GRL_NAME=$CELL_TYPE"_tc_grl.rds"
EXTENSION_DISTANCE=200


### 3 ### get_tc_profiles.R

# CTSS_RSE_NAME, TC_GRL_NAME, OUTPUT_DIR from above
PROFILE_MAIN_DIR=$CELL_TYPE"_tc_profiles"
PROFILE_SUB_DIR="tcs"
# add -s if you want to save count profiles


# Run scripts
Rscript _get_tc_from_ctss.r -c $RESOURCES_DIR/$CTSS_RSE_NAME -o $OUTPUT_DIR -t $TC_GRL_NAME -e $EXTENSION_DISTANCE
Rscript _get_tc_profiles.r -c $RESOURCES_DIR/$CTSS_RSE_NAME -t $OUTPUT_DIR/$TC_GRL_NAME -o $OUTPUT_DIR -n $PROFILE_MAIN_DIR -r $PROFILE_SUB_DIR
