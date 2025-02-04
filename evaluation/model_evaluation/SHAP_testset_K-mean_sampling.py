#import pyarrow.parquet as pq
#from lightgbm import LGBMClassifier
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import shap
#from shap.explainers._tree import Tree
import os
import pickle

#from sklearn.metrics import roc_auc_score
#from sklearn.utils import shuffle

from python.extraction import extract_filenames
#from python.extraction import extract_ranges
from python.extraction import extract_profiles
from python.helpfn_for_prediction import set_index_if_exists
#from python.helpfn_for_prediction import adjust_genomic_positions_for_bed
#from python.helpfn_for_prediction import move_metadata_columns_to_ranges

from sklearn.cluster import KMeans
import joblib

# directories
working_dir = "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation"
profile_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/GM12878_wt10M_profiles_te/profiles_subtnorm"
#explainer_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/shap_explainer_PRIMEloci_GM12878_model_1.0_trvl.pkl"
output_sir = "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/results_SHAP"
explainer_dir = "shap_explainer_PRIMEloci_GM12878_model_1.0_trvl.pkl"

#os.chdir(working_dir)


# Step 1: Load SHAP Explainer from Pickle File
#with open(explainer_dir, "rb") as file:
#    explainer = pickle.load(file)
explainer = joblib.load(explainer_dir)
# 📌 Step 2: Verify the Explainer Object
print(type(explainer))  # Should print <class 'shap.explainers.tree.Tree'>

'''
total_clusters = 500


# Test data
profile_file_ls = extract_filenames(profile_dir)

pos_te_list = [element for element in profile_file_ls if "_pos_" in element]
neg_te_list = [element for element in profile_file_ls if "_neg_" in element]

pos_profiles = []

for profile_filename in pos_te_list:
    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
    subtnorm_df = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=None)
    subtnorm_df = set_index_if_exists(subtnorm_df)
    subtnorm_np = extract_profiles(subtnorm_df)
    pos_profiles.append(subtnorm_np)

if pos_profiles:
    pos_concatenated_profiles = np.concatenate(pos_profiles, axis=0)
    pos_labels = [1] * len(pos_concatenated_profiles)
    print("Final concatenated array shape:", pos_concatenated_profiles.shape)
else:
    print("No profiles were processed.")

neg_profiles = []

for profile_filename in neg_te_list:
    input_profiles_subtnorm_path = os.path.join(profile_dir, profile_filename)
    subtnorm_df = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=None)
    subtnorm_df = set_index_if_exists(subtnorm_df)
    subtnorm_np = extract_profiles(subtnorm_df)
    neg_profiles.append(subtnorm_np)

if neg_profiles:
    neg_concatenated_profiles = np.concatenate(neg_profiles, axis=0)
    neg_labels = [0] * len(neg_concatenated_profiles)
    print("Final concatenated array shape:", neg_concatenated_profiles.shape)
else:
    print("No profiles were processed.")

all_concatenated_profiles = np.concatenate((pos_concatenated_profiles, neg_concatenated_profiles), axis=0)
print("Final concatenated array shape:", all_concatenated_profiles.shape)
all_labels = pos_labels + neg_labels

test_X = all_concatenated_profiles
test_y = np.array(all_labels)

print(test_X[0])


# Step 3: Stratified K-Means Sampling on X_test

# Count instances per class in full test set
class_counts = np.bincount(test_y)  
total_samples = len(test_y)

# Calculate proportional clusters per class
class_proportions = class_counts / total_samples  
class_clusters = (class_proportions * total_clusters).astype(int)  # Assign clusters per class

# Initialize list to store sampled data
X_sampled_list, y_sampled_list = [], []

# Apply K-Means within each class
for class_label in np.unique(test_y):
    # Extract samples for the current class
    X_class = test_X[test_y == class_label]
    
    # Define number of clusters for this class
    n_clusters = max(1, class_clusters[class_label])  # Ensure at least 1 cluster

    # Apply K-Means to the subset of the class
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    kmeans.fit(X_class)

    # Collect the cluster centers as representative samples
    X_sampled_list.append(pd.DataFrame(kmeans.cluster_centers_, columns=test_X.columns))
    y_sampled_list.append([class_label] * n_clusters)  # Assign class labels

# Combine all stratified sampled data
X_sampled = pd.concat(X_sampled_list, ignore_index=True)
y_sampled = np.concatenate(y_sampled_list)


# 📌 Step 4: Compute SHAP Values for Sampled Data
shap_values = explainer(X_sampled)

# 📌 Step 5: Visualize SHAP Values with Beeswarm Plot
shap.plots.beeswarm(shap_values)

'''