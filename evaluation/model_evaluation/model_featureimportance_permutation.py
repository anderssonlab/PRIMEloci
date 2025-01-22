import pyarrow.parquet as pq
from lightgbm import LGBMClassifier
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
from sklearn.metrics import roc_auc_score
from sklearn.utils import shuffle

from python.extraction import extract_filenames
from python.extraction import extract_ranges
from python.extraction import extract_profiles
from python.helpfn_for_prediction import set_index_if_exists
from python.helpfn_for_prediction import adjust_genomic_positions_for_bed
from python.helpfn_for_prediction import move_metadata_columns_to_ranges

import os
import pickle

# directories
working_dir = "/home/zmk214/zmk214/PRIMEloci/evaluation/model_evaluation"
model_dir = "/home/zmk214/zmk214/PRIMEloci_model/2_training/PRIMEloci_GM12878_model_1.0_trvl.sav"
profile_dir = "/home/zmk214/zmk214/PRIMEloci_model/GM12878_wt10M_profiles_te/profiles_subtnorm"

os.chdir(working_dir)

# Load the trained model
with open(model_dir, 'rb') as file:
    model = pickle.load(file)

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

# Uncomment the following lines to test with a balanced subset of data
pos_indices_sample = np.random.choice(np.where(test_y == 1)[0], 100, replace=False)  # 100 positive samples
neg_indices_sample = np.random.choice(np.where(test_y == 0)[0], 100, replace=False)  # 100 negative samples
selected_indices = np.concatenate([pos_indices_sample, neg_indices_sample])
test_X_sample = test_X[selected_indices]  # Subset of features
test_y_sample = test_y[selected_indices]  # Subset of labels

# Shuffle the dataset
shuffled_indices = np.random.permutation(len(test_y_sample))
test_X = test_X_sample[shuffled_indices]
test_y = test_y_sample[shuffled_indices]

# Shuffle the dataset
#shuffled_indices = np.random.permutation(len(test_y))
#test_X = test_X[shuffled_indices]
#test_y = test_y[shuffled_indices]




# Function to calculate Permutation Feature Importance
def permutation_feature_importance(model, X, y, metric=roc_auc_score, n_repeats=5):
    baseline_score = metric(y, model.predict_proba(X)[:, 1])  # Base performance
    feature_names = [f"Feature_{i}" for i in range(X.shape[1])]  # Map to original indices
    X = pd.DataFrame(X, columns=feature_names)
    importance = pd.DataFrame(index=X.columns, columns=range(n_repeats))

    for col in X.columns:
        for i in range(n_repeats):
            X_permuted = X.copy()
            X_permuted[col] = shuffle(X[col])  # Shuffle feature values
            permuted_score = metric(y, model.predict_proba(X_permuted)[:, 1])
            importance.loc[col, i] = baseline_score - permuted_score  # Drop in performance

    return importance.mean(axis=1).sort_values(ascending=False), importance

# Corrected mapping function
def map_feature_indices(features):
    """Map feature indices from 0-400 to -200 to 200"""
    return [int(feature.replace("Feature_", "")) - 200 for feature in features]

# Split data by predictions and true labels
print("Split data by predictions and true labels")
pred_probs = model.predict_proba(test_X)[:, 1]
pos_indices_pred = pred_probs >= 0.5
neg_indices_pred = pred_probs < 0.5

# Permutation Feature Importance
print("Calculating Permutation Feature Importance for positive predictions...")
importance_pos_pred_mean, importance_pos_pred_repeats = permutation_feature_importance(
    model, test_X[pos_indices_pred], test_y[pos_indices_pred]
)
print("Calculating Permutation Feature Importance for negative predictions...")
importance_neg_pred_mean, importance_neg_pred_repeats = permutation_feature_importance(
    model, test_X[neg_indices_pred], test_y[neg_indices_pred]
)

# Save importance objects
with open("importance_pos_pred_mean.pkl", "wb") as file:
    pickle.dump(importance_pos_pred_mean, file)

with open("importance_neg_pred_mean.pkl", "wb") as file:
    pickle.dump(importance_neg_pred_mean, file)

with open("importance_pos_pred_repeats.pkl", "wb") as file:
    pickle.dump(importance_pos_pred_repeats, file)

with open("importance_neg_pred_repeats.pkl", "wb") as file:
    pickle.dump(importance_neg_pred_repeats, file)

# Uncomment below lines to load saved importance objects if needed
# with open("importance_pos_pred_mean.pkl", "rb") as file:
#     importance_pos_pred_mean = pickle.load(file)

# with open("importance_neg_pred_mean.pkl", "rb") as file:
#     importance_neg_pred_mean = pickle.load(file)

# with open("importance_pos_pred_repeats.pkl", "rb") as file:
#     importance_pos_pred_repeats = pickle.load(file)

# with open("importance_neg_pred_repeats.pkl", "rb") as file:
#     importance_neg_pred_repeats = pickle.load(file)

# Map feature names to -200 to 200
mapped_features = map_feature_indices(importance_pos_pred_mean.index)

# Update index to new mapped features
importance_pos_pred_mean.index = mapped_features
importance_neg_pred_mean.index = mapped_features

# Sort by feature index
importance_pos_pred_mean = importance_pos_pred_mean.sort_index()
importance_neg_pred_mean = importance_neg_pred_mean.sort_index()

# Save plot for positive predictions
plt.figure(figsize=(12, 6))
plt.bar(importance_pos_pred_mean.index, importance_pos_pred_mean.values, color='green', alpha=0.7)
plt.title('Permutation Feature Importance (Positive Predictions)')
plt.xlabel('Feature Index (-200 to 200)')
plt.ylabel('Mean Importance (Decrease in AUC)')
plt.grid()
plt.savefig("positive_permutation_importance.pdf")
plt.show()

# Save plot for negative predictions
plt.figure(figsize=(12, 6))
plt.bar(importance_neg_pred_mean.index, importance_neg_pred_mean.values, color='red', alpha=0.7)
plt.title('Permutation Feature Importance (Negative Predictions)')
plt.xlabel('Feature Index (-200 to 200)')
plt.ylabel('Mean Importance (Decrease in AUC)')
plt.grid()
plt.savefig("negative_permutation_importance.pdf")
plt.show()