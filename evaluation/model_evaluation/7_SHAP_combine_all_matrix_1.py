import pandas as pd
import pickle
import pyreadr

shap_pkl_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/SHAP_value_K562_N_endres_v2.pkl"
shap_index_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/SHAP/K562_endres_v2/K562_N_endres_v2.index"
subtnorm_profile_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles_subtnorm/K562_N_combined.csv"
count_profile_path = "/Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles/K562_N_combined.csv"

outname = "K562_N_endres_v2"

manual_colnames = range(-200, 201)

# SHAP_values

with open(shap_pkl_path, "rb") as file:
    shap_values = pickle.load(file)
shap_df = pd.DataFrame(shap_values.values[:, :, 1], columns=manual_colnames)
print(shap_df.head())
print(shap_df.shape)

shap_index_df = pd.read_csv(shap_index_path, header=None, sep='\t')
print(shap_index_df.head())
print(shap_index_df.shape)

shap_df.index = shap_index_df[0].values
print(shap_df.head())

# profile subtnorm
subtnorm_df = pd.read_csv(subtnorm_profile_path, header=0, index_col="rownames")
subtnorm_df.index.name = None
subtnorm_df.columns = manual_colnames
print(subtnorm_df.head())
print(subtnorm_df.shape)

# profile count
count_df = pd.read_csv(count_profile_path, header=0, index_col="rownames")
count_df.index.name = None
print(count_df.head())
print(count_df.shape)

# Split the DataFrame into two based on column prefixes
plus_count_df = count_df.loc[:, count_df.columns.str.startswith("Plus")]
plus_count_df.columns = manual_colnames
minus_count_df = count_df.loc[:, count_df.columns.str.startswith("Minus")]
minus_count_df.columns = manual_colnames
print(plus_count_df.head())
print(plus_count_df.shape)
print(minus_count_df.head())
print(minus_count_df.shape)


# Check if all DataFrames have the same index, column names (order matters), and shape
dataframes = [shap_df, subtnorm_df, plus_count_df, minus_count_df]
index_check = all(dataframes[0].index.equals(df.index) for df in dataframes)
columns_check = all(dataframes[0].columns.equals(df.columns) for df in dataframes)
shape_check = all(dataframes[0].shape == df.shape for df in dataframes)

print(f"All DataFrames have the same index: {index_check}")
print(f"All DataFrames have the same column names: {columns_check}")
print(f"All DataFrames have the same shape: {shape_check}")


# Save each DataFrame as an RDS file
pyreadr.write_rds(outname + "_shap_value_df.rds", shap_df)
pyreadr.write_rds(outname + "_profile_subtnorm_df.rds", subtnorm_df)
pyreadr.write_rds(outname + "_profile_count_plus_df.rds", plus_count_df)
pyreadr.write_rds(outname + "_profile_count_minus_df.rds", minus_count_df)

index_df = pd.DataFrame(minus_count_df.index)
index_df.columns = ["rownames"]
index_df.index.name = None
print(index_df.head())
pyreadr.write_rds(outname + "_index.rds", index_df)


print("All DataFrames saved as .rds files.")
