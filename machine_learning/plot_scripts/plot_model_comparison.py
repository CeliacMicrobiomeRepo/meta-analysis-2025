"""

Automatically plot the results of the best models in one dataset.

- plots a box and whisker plot of AUCs across all folds/replicates of the best model hyperparameter for each model in one dataset
- writes a 3:4 png and an svg file

"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Inputs ------------------------------------------
# VAL_METHOD = "logo"
# DATASET_DIR_PATH = "./datasets/duodenal_noGFD_Genus_after_log10"
# VAL_METHOD = "kfold"
# DATASET_DIR_PATH = "./datasets/duodenal_noGFD_Genus_after_log10"
# VAL_METHOD = "logo"
# DATASET_DIR_PATH = "./datasets/stool_noGFD_Genus_after_log10"
# VAL_METHOD = "kfold"
# DATASET_DIR_PATH = "./datasets/stool_noGFD_Genus_after_log10"

# VAL_METHOD = "logo"
# DATASET_DIR_PATH = "./datasets/duodenal_noGFD_ASV_after_log10"
# VAL_METHOD = "kfold"
# DATASET_DIR_PATH = "./datasets/duodenal_noGFD_ASV_after_log10"
# VAL_METHOD = "logo"
# DATASET_DIR_PATH = "./datasets/stool_noGFD_ASV_after_log10"
VAL_METHOD = "kfold"
DATASET_DIR_PATH = "./datasets/stool_noGFD_ASV_after_log10"

RESULTS_SUBDIR_FILE_NAME = VAL_METHOD + "_results/best_models_results.tsv"
RESULTS_PATH = os.path.join(DATASET_DIR_PATH, RESULTS_SUBDIR_FILE_NAME)
DATASET_NAME = os.path.basename(DATASET_DIR_PATH)

# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "./plots/model_comparison"
TITLE = "Different Model Performances on " + VAL_METHOD.upper() + " Within The " + DATASET_NAME + " Dataset"
FILE_NAME = "model_comparison_" + DATASET_NAME + "_" + VAL_METHOD


# Plot ------------------------------------------
# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR_PATH):
    os.makedirs(OUTPUT_DIR_PATH)

# Read and process results
df = pd.read_csv(RESULTS_PATH, sep='\t')

# Dynamically get fold names based on validation method
if VAL_METHOD == "kfold":
    fold_names = [f'Fold_{i}' for i in range(1, 11)]
else:  # logo
    # Get all column names that aren't 'Replicate' or 'Model'
    fold_names = [col for col in df.columns if col not in ['Replicate', 'Model', 'Mean']]

# Filter out the 'Average' rows and get AUC values across folds
df_values = df[df['Replicate'] != 'Average'].melt(
    id_vars=['Replicate', 'Model'],
    value_vars=fold_names,
    var_name='Fold',
    value_name='AUC'
)

# Average the AUCs across the folds, so there is one value per replicate per model
df_values = df_values.groupby(['Replicate', 'Model'])['AUC'].mean().reset_index()

# Sort models alphabetically and capitalize model names
unique_models = [model.upper() for model in df_values['Model'].unique()]  # Convert to uppercase first
df_values['Model'] = pd.Categorical(
    df_values['Model'].str.upper(),  # Capitalize model names
    categories=sorted(unique_models)  # Use the pre-processed uppercase models
)

# Create the plot
plt.figure(figsize=(12, 9))
sns.set_style("whitegrid")
sns.boxplot(
    data=df_values,
    x='Model',
    y='AUC',
    width=0.7,
    showfliers=False
)

# Customize the plot
plt.xticks(rotation=45, ha='right')
plt.xlabel('')
plt.ylabel('AUC')
plt.title(TITLE)

# Set y-axis limits and ticks
plt.ylim(0.5, 1.0)  # Set y-axis limits
plt.yticks(np.arange(0.5, 1.05, 0.05))


# Export ------------------------------------------
# Save as PNG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.png'), dpi=300, bbox_inches='tight')
# Save as SVG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.svg'), format='svg', bbox_inches='tight')

# Show plot
plt.show()

plt.close()
