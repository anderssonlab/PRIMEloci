# Example run of genome-wide prediction:

module load gcc/11.2.0 R/4.2.2 R-packages/4.2.2
module load python/3.9.9 python-packages/3.9

Example model is model_Meena_v3_gm12878_fullModel_C.sav [trained on GM12878 whole cell CAGE due to the small size]


cd genomewide_prediction
./script_bw_to_prediction.sh
