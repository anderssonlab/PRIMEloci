import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


shap_pkl_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/SHAP_value_K562_C_endres.pkl"
outname = "pooled_K562_C_1"

# Load SHAP values from Pickle file
with open(shap_pkl_path, "rb") as file:
    shap_values = pickle.load(file)

# Check the type of the SHAP object
print(f" SHAP object type: {type(shap_values)}")

# Check available attributes in the SHAP object
print(f" SHAP object attributes: {dir(shap_values)}")

# Check the shape of SHAP values
print(f" SHAP values shape: {shap_values.values.shape}")

# Preview the first few SHAP values
print(" First few SHAP values:\n", shap_values.values[:5])

# Check the base values (expected model output)
print(" Base values (expected value per prediction):\n", shap_values.base_values[:5])

# Check feature names (if available)
if hasattr(shap_values, "feature_names") and shap_values.feature_names is not None:
    print(" Feature names:\n", shap_values.feature_names)
else:
    print("⚠️ No feature names found in SHAP object.")


# Convert SHAP values to a DataFrame
manual_colnames = range(-200, 201)
shap_df = pd.DataFrame(shap_values.values[:, :, 1], columns=manual_colnames)
#shap_df = pd.DataFrame(shap_values.values, columns=shap_values.feature_names)
print(shap_df.head())

# Box Plot with Rotated X-Labels
plt.figure(figsize=(24, 6))
sns.boxplot(data=shap_df, palette="coolwarm", linewidth=1, showfliers=False)

# Formatting
plt.xlabel("Feature Index (-200 to 200)", fontsize=12)
plt.ylabel("SHAP Value", fontsize=14)
plt.title("Box Plot of SHAP Values per Feature", fontsize=16)

# Rotate x-axis labels for better readability
plt.xticks(rotation=90, fontsize=3)

# Save the plot
plt.savefig(outname + "_shap_boxplot.pdf", format="pdf", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()


# Compute SHAP statistics per feature
shap_stats_df = pd.DataFrame({
    "Feature": manual_colnames,
    "Mean_Abs_SHAP_Value": shap_df.abs().mean(axis=0),  # Mean absolute SHAP values
    "Sum_SHAP_Value": shap_df.sum(axis=0),  # Sum of SHAP values
    "Mean_SHAP_Value": shap_df.mean(axis=0)  # Mean signed SHAP values
})
#shap_stats_df = pd.DataFrame({
#    "Feature": shap_values.feature_names,
#    "Mean_Abs_SHAP_Value": shap_df.abs().mean(axis=0),  # Mean absolute SHAP values
#    "Sum_SHAP_Value": shap_df.sum(axis=0),  # Sum of SHAP values
#    "Mean_SHAP_Value": shap_df.mean(axis=0)  # Mean signed SHAP values
#})

# Save SHAP statistics to CSV
shap_stats_csv_path = outname + "_shap_statistics_per_feature.csv"
shap_stats_df.to_csv(shap_stats_csv_path, index=False)

# Bar Plots for SHAP Statistics
# Mean Absolute SHAP Value
plt.figure(figsize=(12, 6))
plt.bar(shap_stats_df["Feature"], shap_stats_df["Mean_Abs_SHAP_Value"], color="purple")
plt.xlabel("Feature Position (-200 to 200)", fontsize=14)
plt.ylabel("Mean Absolute SHAP Value", fontsize=14)
plt.title("Mean Absolute SHAP Value per Feature", fontsize=16)
plt.xticks(rotation=90, fontsize=8)
plt.savefig(outname + "_mean_abs_shap_value_bar_plot.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.close()

# Sum of SHAP Values
plt.figure(figsize=(12, 6))
plt.bar(shap_stats_df["Feature"], shap_stats_df["Sum_SHAP_Value"], color="blue")
plt.xlabel("Feature Position (-200 to 200)", fontsize=14)
plt.ylabel("Sum of SHAP Values", fontsize=14)
plt.title("Sum of SHAP Values per Feature", fontsize=16)
plt.xticks(rotation=90, fontsize=8)
plt.savefig(outname + "_sum_shap_value_bar_plot.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.close()

# Mean SHAP Value
plt.figure(figsize=(12, 6))
plt.bar(shap_stats_df["Feature"], shap_stats_df["Mean_SHAP_Value"], color="green")
#plt.bar(range(len(shap_stats_df["Mean_SHAP_Value"])), shap_stats_df["Mean_SHAP_Value"], color="green")
plt.xlabel("Feature Position (-200 to 200)", fontsize=14)
plt.ylabel("Mean SHAP Value", fontsize=14)
plt.title("Mean SHAP Value per Feature", fontsize=16)
plt.xticks(rotation=90, fontsize=8)
plt.savefig(outname + "_mean_shap_value_bar_plot.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.close()

print("SHAP feature statistics saved to:", shap_stats_csv_path)
print("Violin plot saved to: shap_violin_plot.pdf")
print("Mean Absolute SHAP Value bar plot saved to: mean_abs_shap_value_bar_plot.pdf")
print("Sum of SHAP Values bar plot saved to: sum_shap_value_bar_plot.pdf")
print("Mean SHAP Value bar plot saved to: mean_shap_value_bar_plot.pdf")
