"""

Automatically plot the results the best model between some selected datasets.

- plots a box and whisker plot of (mean of all folds) AUCs for each replicate of the best model for each dataset
- writes a png and an svg file

"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



# Inputs ------------------------------------------
# Input datasets to plot
# DATASET_DIR_PATHS_AND_NAMES = [
#     ("./datasets/duodenal_noGFD_ASV_after_log10/", "noV3,ASV"),
#     ("./datasets/duodenal_noGFD_Genus_noV3_after_log10/", "noV3,Genus"),
#     ("./datasets/duodenal_noGFD_Genus_after_log10/", "wV3,Genus")
# ]
DATASET_DIR_PATHS_AND_NAMES = [
    ("./datasets/stool_wGFD_ASV_after_log10/", "wGFD,ASV"),
    ("./datasets/stool_wGFD_Genus_noV3_after_log10/", "wGFD-noV3,Genus"),
    ("./datasets/stool_wGFD_Genus_after_log10/", "wGFD,Genus"),
    ("./datasets/stool_noGFD_ASV_after_log10/", "noGFD,ASV"),
    ("./datasets/stool_noGFD_Genus_noV3_after_log10/", "noGFD-noV3,Genus"),
    ("./datasets/stool_noGFD_Genus_after_log10/", "noGFD,Genus")
]

# Which results to plot
# RESULTS_SUBDIR_FILE_NAME = "kfold_results/best_models_results.tsv"
RESULTS_SUBDIR_FILE_NAME = "logo_results/best_models_results.tsv"



# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "./plots/main_comparisons"

# Name of plot files to write
# FILE_NAME = "main_comparison_duodenal_kfold"
# FILE_NAME = "main_comparison_duodenal_logo"
# FILE_NAME = "main_comparison_stool_kfold"
FILE_NAME = "main_comparison_stool_logo"

# Title of plot
# TITLE = "Best Model Performance Across Different Duodenal Datasets (K-fold)"
# TITLE = "Best Model Performance Across Different Duodenal Datasets (logo)"
# TITLE = "Best Model Performance Across Different Stool Datasets (K-fold)"
TITLE = "Best Model Performance Across Different Stool Datasets (logo)"

# Create a color mapping for transformations
TRANSFORMATION_COLORS = {
    'ASV': '#1f77b4',  # blue
    'Genus': '#ff7f0e'  # orange
}



# Plot ------------------------------------------
# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR_PATH):
    os.makedirs(OUTPUT_DIR_PATH)

# Read and process results from each dataset
dataset_results = []
for path, name in DATASET_DIR_PATHS_AND_NAMES:
    results_path = os.path.join(path, RESULTS_SUBDIR_FILE_NAME)
    if os.path.exists(results_path):
        df = pd.read_csv(results_path, sep='\t')
        # Identify the best model
        # Get only the average rows and find the best model based on mean AUC
        best_model = df[df['Replicate'] == 'Average'].loc[df[df['Replicate'] == 'Average']['Mean'].idxmax(), 'Model']
        df = df[df['Model'] == best_model]  # Keep only rows for best model
        # Remove the 'Average' rows
        df = df[df['Replicate'] != 'Average']
        # Create a DataFrame with the AUC values
        df_values = pd.DataFrame({
            'AUC': df['Mean'].values,  # Use the Mean column for AUC values
            'Dataset': [name] * len(df)  # Repeat dataset name for each row
        })
        dataset_results.append(df_values)
        print(name, "Best Model:", best_model)

# Combine all results
all_results = pd.concat(dataset_results, ignore_index=True)

# Create the plot
plt.figure(figsize=(12, 9))
sns.set_style("whitegrid")

# Extract transformation type from dataset name (second part after comma)
all_results['Transformation'] = all_results['Dataset'].apply(lambda x: x.split(',')[1])

# Create boxplot with colors based on transformation
sns.boxplot(
    data=all_results,
    x='Dataset',
    y='AUC',
    hue='Transformation',
    width=0.7,
    showfliers=False,
    palette=TRANSFORMATION_COLORS,  # Simply pass the color dictionary
    legend=False
)

# Customize the plot
plt.xticks(rotation=45, ha='right')
plt.xlabel('')
plt.ylabel('AUC')
plt.title(TITLE)

# Set y-axis limits and ticks
plt.ylim(0.5, 1.0)  # Set y-axis limits
plt.yticks(np.arange(0.5, 1.05, 0.05))

# Adjust layout to prevent label cutoff
plt.tight_layout()


# Export ------------------------------------------
# Save as PNG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.png'), dpi=300, bbox_inches='tight')
# Save as SVG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.svg'), format='svg', bbox_inches='tight')

# Show plot
plt.show()

plt.close()
