import pyarrow.parquet as pq
from lightgbm import LGBMClassifier
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
from sklearn.metrics import roc_auc_score
from sklearn.utils import shuffle

#from python.extraction import extract_filenames
#from python.extraction import extract_ranges
#from python.extraction import extract_profiles
#from python.helpfn_for_prediction import set_index_if_exists
#from python.helpfn_for_prediction import adjust_genomic_positions_for_bed
#from python.helpfn_for_prediction import move_metadata_columns_to_ranges

import os
import pickle

# directories
working_dir = "/home/zmk214/zmk214/PRIMEloci/evaluation/model_evaluation"
model_dir = "/home/zmk214/zmk214/PRIMEloci_model/2_training/PRIMEloci_GM12878_model_1.0.sav"
#profile_dir = "/home/zmk214/zmk214/PRIMEloci_model/GM12878_wt10M_profiles_te/profiles_subtnorm"

os.chdir(working_dir)

# Load the trained model
with open(model_dir, 'rb') as file:
    model = pickle.load(file)

# Test data
#profile_file_ls = extract_filenames(profile_dir)
#
#pos_te_list = [element for element in profile_file_ls if "_pos_" in element]
#neg_te_list = [element for element in profile_file_ls if "_neg_" in element]
#
#pos_profiles = []
#
#for profile_filename in pos_te_list:
#    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
#    subtnorm_df = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=None)
#    subtnorm_df = set_index_if_exists(subtnorm_df)
#    subtnorm_np = extract_profiles(subtnorm_df)
#    pos_profiles.append(subtnorm_np)
#
#if pos_profiles:
#    pos_concatenated_profiles = np.concatenate(pos_profiles, axis=0)
#    pos_labels = [1] * len(pos_concatenated_profiles)
#    print("Final concatenated array shape:", pos_concatenated_profiles.shape)
#else:
#    print("No profiles were processed.")
#
#neg_profiles = []
#
#for profile_filename in neg_te_list:
#    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
#    subtnorm_df = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=None)
#    subtnorm_df = set_index_if_exists(subtnorm_df)
#    subtnorm_np = extract_profiles(subtnorm_df)
#    neg_profiles.append(subtnorm_np)
#
#if neg_profiles:
#    neg_concatenated_profiles = np.concatenate(neg_profiles, axis=0)
#    neg_labels = [0] * len(neg_concatenated_profiles)
#    print("Final concatenated array shape:", neg_concatenated_profiles.shape)
#else:
#    print("No profiles were processed.")
#
#all_concatenated_profiles = np.concatenate((pos_concatenated_profiles, neg_concatenated_profiles), axis=0)
#print("Final concatenated array shape:", all_concatenated_profiles.shape)
#all_labels = pos_labels + neg_labels
#
#test_X = all_concatenated_profiles
#test_y = np.array(all_labels)
#
## Shuffle the dataset
#shuffled_indices = np.random.permutation(len(test_y))
#test_X = test_X[shuffled_indices]
#test_y = test_y[shuffled_indices]

# Calculate or load SHAP explainer
# Uncomment below to load precomputed SHAP explainer
# with open("shap_explainer.pkl", "rb") as file:
#     explainer = pickle.load(file)
# else calculate and save
explainer = shap.TreeExplainer(model)
with open("shap_explainer_PRIMEloci_GM12878_model_1.0.pkl", "wb") as file:
    pickle.dump(explainer, file)

# Calculate SHAP values
#shap_values = explainer.shap_values(test_X)
#
## Save SHAP values for later analysis
#with open("shap_values_positive.pkl", "wb") as file:
#    pickle.dump(shap_values[1], file)
#
#with open("shap_values_negative.pkl", "wb") as file:
#    pickle.dump(shap_values[0], file)
#
## Save SHAP summary to CSV
#shap_summary_positive = pd.DataFrame(shap_values[1], columns=[f"Feature_{i}" for i in range(test_X.shape[1])])
#shap_summary_positive["Mean"] = shap_summary_positive.mean(axis=1)
#shap_summary_positive.to_csv("shap_summary_positive.csv", index=False)
#
#shap_summary_negative = pd.DataFrame(shap_values[0], columns=[f"Feature_{i}" for i in range(test_X.shape[1])])
#shap_summary_negative["Mean"] = shap_summary_negative.mean(axis=1)
#shap_summary_negative.to_csv("shap_summary_negative.csv", index=False)
#
## Plot SHAP summary for positive predictions
#plt.figure(figsize=(12, 8))
#shap.summary_plot(shap_values[1], test_X, plot_type="dot", color_bar_label="CAGE signal (processed)", show=False)
#plt.title("SHAP Summary Plot (Positive Predictions)", fontsize=16, color="brown")
#plt.xlabel("SHAP value (impact on model output)", fontsize=12)
#plt.savefig("shap_summary_positive.pdf")
#plt.show()
#
## Optionally, for negative predictions
#plt.figure(figsize=(12, 8))
#shap.summary_plot(shap_values[0], test_X, plot_type="dot", color_bar_label="CAGE signal (processed)", show=False)
#plt.title("SHAP Summary Plot (Negative Predictions)", fontsize=16, color="brown")
#plt.xlabel("SHAP value (impact on model output)", fontsize=12)
#plt.savefig("shap_summary_negative.pdf")
#plt.show()
