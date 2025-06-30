#!/bin/bash

CELL_TYPE=GM12878_wt10M
INPUT_FILE_DHS=/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/E116-DNase.macs2.hg38.narrowPeak
INPUT_DIR_CAGE=/projects/ralab/data/CAGE/nucleiCAGE/bw_files/
INPUT_FILE_DESIGN=design_matrix_${CELL_TYPE}.tsv
INPUT_FILE_BLACKLIST=/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/hg38-blacklist.v2.bed

nohup Rscript 1_dhs.r -i ${INPUT_FILE_DHS} -b $INPUT_FILE_BLACKLIST -n $CELL_TYPE -d 200 > nohup_${CELL_TYPE}_1

nohup Rscript 2_cage_ocrlikeneg.r -i ${INPUT_DIR_CAGE} -m ${INPUT_FILE_DESIGN} -r 1_${CELL_TYPE}_ocr.RData -n $CELL_TYPE -d 200 >> nohup_${CELL_TYPE}_2

nohup Rscript 3_augmentation.r -r 2_${CELL_TYPE}_ocr_ocrlikeneg.RData -n $CELL_TYPE -d 200 >> nohup_${CELL_TYPE}_3

nohup Rscript 4_annotation.r -r 3_${CELL_TYPE}_augmentation.RData -c 2_${CELL_TYPE}_orig_ctss.RData -n $CELL_TYPE -u TRUE >> nohup_${CELL_TYPE}_4

nohup Rscript 5_validation_and_sltnegative.r -r 4_${CELL_TYPE}_augmentation_with_annotation.RData -c 2_${CELL_TYPE}_orig_ctss.RData -n $CELL_TYPE >> nohup_${CELL_TYPE}_5

nohup Rscript 6_create_profiles.r -r 5_${CELL_TYPE}_tr_vl_te.RData -c 2_${CELL_TYPE}_orig_ctss.RData -d 200 -O ${CELL_TYPE}_profiles > nohup_${CELL_TYPE}_6
