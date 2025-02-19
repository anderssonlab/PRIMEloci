import numpy as np
import joblib
import shap
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list
import pickle
import pandas as pd
import seaborn as sns

# directories
#working_dir = "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation"
#profile_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_C_endres.npz"
#explainer_dir = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/shap_explainer_PRIMEloci_GM12878_model_1.0.pkl"
#output_dir = "/Users/natsudanav/Desktop/PRIMEloci/evaluation/model_evaluation/results_SHAP"

working_dir = "/home/zmk214/zmk214/PRIMEloci/evaluation/model_evaluation"
profile_dir = "/home/zmk214/zmk214/K562_C_endres.npz"
explainer_dir = "/home/zmk214/zmk214/shap_explainer_PRIMEloci_GM12878_model_1.0.pkl"

outname = "K562_C_endres"
os.chdir(working_dir)

# Step 1: Load SHAP Explainer
explainer = joblib.load(explainer_dir)
print(type(explainer))

# Step 2: Load profiles and extract the stored array
data = np.load(profile_dir)
profiles = data["test_X"]
print(f"Loaded filtered profiles: {profiles.shape}")

# Step 3: Compute SHAP values (if not already computed)
shap_values = explainer(profiles)
# Define the output file path
shap_pkl_path = "SHAP_value_"+outname+".pkl"

# Save SHAP values to Pickle
with open(shap_pkl_path, "wb") as file:
    pickle.dump(shap_values, file)
print(f"SHAP values saved to {shap_pkl_path}")

# # Load SHAP values from Pickle
# with open(shap_pkl_path, "rb") as file:
#     loaded_shap_values = pickle.load(file)
# 
# print("SHAP values loaded successfully!")

# Step 4: Visualize SHAP Values with Summary Plot

## Generate the summary plot
# plt.figure(figsize=(200, 20))
# shap.summary_plot(shap_values, profiles, show=False, max_display=401)
# plt.savefig("SHAP_beeswarm_"+outname+".pdf", format="pdf", dpi=300, bbox_inches="tight")
# 
# 
# ## individual bar plots
# # Define the output PDF file
# output_bar_path = "SHAP_individual_barplot_" + outname + ".pdf"
# num = 10
# # Open a multi-page PDF
# with PdfPages(output_bar_path) as pdf:
# 
#     shap.plots.bar(shap_values, max_display=401)
#     plt.gcf().set_size_inches(200, 20)
#     pdf.savefig(plt.gcf())
#     plt.close()
# 
#     for i in range(num):
#         shap.plots.bar(shap_values[i], max_display=401)
#         plt.gcf().set_size_inches(200, 20)
#         pdf.savefig(plt.gcf())
#         plt.close()
# 
# print(f"SHAP bar plots saved to {output_bar_path}")
# 
# # CLUSTERING HEATMAP
# shap_df = pd.DataFrame(shap_values.values, columns=shap_values.feature_names)
# # Perform hierarchical clustering on **samples (rows)** instead of features
# linkage_matrix = linkage(shap_df, method="ward")  # Cluster samples (rows)
# ordered_rows = leaves_list(linkage_matrix)  # Get ordered row indices
# # Reorder the DataFrame based on sample clustering
# shap_df_reordered = shap_df.iloc[ordered_rows, :]  # Reorder rows (samples)
# # Set up figure size
# plt.figure(figsize=(70, 50))
# sns.heatmap(shap_df_reordered.T, cmap="coolwarm", center=0, robust=True, linewidths=0.5)
# plt.xlabel("Clustered Samples", fontsize=14)
# plt.ylabel("Features", fontsize=14)
# plt.title("SHAP Heatmap (Samples Clustered, No Dendrogram)", fontsize=16)
# plt.savefig("SHAP_clustering_heatmap_"+outname+".pdf", format="pdf", dpi=300, bbox_inches="tight")


'''
# Stack SHAP values and assign cohort labels
combined_shap_values = shap.Explanation(
    values=np.vstack([shap_values_pos.values, shap_values_neg.values]),  # Combine SHAP values
    base_values=np.hstack([shap_values_pos.base_values, shap_values_neg.base_values]),  # Base values
    data=np.vstack([shap_values_pos.data, shap_values_neg.data]),  # Feature values
    feature_names=shap_values_pos.feature_names  # Keep feature names
)
# Create cohort labels for + and - groups
cohort_labels = ["Positive"] * len(shap_values_pos.values) + ["Negative"] * len(shap_values_neg.values)
plt.figure(figsize=(200, 10))  # Adjust width and height as needed
# Compute and plot mean absolute SHAP values per cohort
shap.plots.bar(combined_shap_values.cohorts(cohort_labels).abs.mean(0), max_display=401)
# Save the plot as high-resolution PDF
plt.savefig("shap_bar_cohorts.pdf", format="pdf", dpi=300, bbox_inches="tight")
'''
