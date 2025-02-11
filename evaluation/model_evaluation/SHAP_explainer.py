import shap
from lightgbm import LGBMClassifier
import os
import pickle
import joblib

# directories
#working_dir = "/home/zmk214/zmk214/PRIMEloci/evaluation/model_evaluation"
#model_dir = "/home/zmk214/zmk214/PRIMEloci_model/2_training/PRIMEloci_GM12878_model_1.0_trvl.sav"
working_dir = "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation"
model_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/model_PRIMEloci/PRIMEloci_GM12878_model_1.0.sav"

os.chdir(working_dir)

# Load the trained model
with open(model_dir, 'rb') as file:
    model = pickle.load(file)

explainer = shap.TreeExplainer(model)
joblib.dump(explainer, "shap_explainer_PRIMEloci_GM12878_model_1.0.pkl")
