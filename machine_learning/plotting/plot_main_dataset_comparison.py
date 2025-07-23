"""
Automatically plot the results the best model between some selected datasets.

- plots a box and whisker plot of mean AUCs for each replicate of the best model for each dataset
- writes a png and an svg file
"""

# Imports ------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns




# Evalulation method
EVALUATION_METHOD = "logo"  # logo or kfold

# Inputs ------------------------------------------

# Input datasets to plot
DATASET_DIR_PATHS_AND_NAMES = [
    ("~/Repos/meta-analysis/machine_learning/datasets/duodenum_active_log10_after/", "Duodenum Active"),
    ("~/Repos/meta-analysis/machine_learning/datasets/stool_active_log10_after/", "Stool Active"),
    ("~/Repos/meta-analysis/machine_learning/datasets/stool_prospective_log10_after/", "Stool Prospective"),
    ("~/Repos/meta-analysis/machine_learning/datasets/stool_treated_log10_after/", "Stool Treated"),
]

# Which results to plot
RESULTS_SUBDIR_FILE_NAME = EVALUATION_METHOD + "_results/best_models_results.tsv"



# Outputs ------------------------------------------
OUTPUT_DIR_PATH = os.path.expanduser("~/Repos/meta-analysis/machine_learning/results")

# Name of plot files to write
FILE_NAME = "main_comparison_" + EVALUATION_METHOD

# Title of plot
TITLE = "Best Model Performance Across Different Datasets (" + EVALUATION_METHOD + ")"



# Plot ------------------------------------------
# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR_PATH):
    os.makedirs(OUTPUT_DIR_PATH)

# Read and process results from each dataset
dataset_results = []
for path, name in DATASET_DIR_PATHS_AND_NAMES:
    results_path = os.path.join(os.path.expanduser(path), RESULTS_SUBDIR_FILE_NAME)
    print(results_path)
    print(os.path.exists(results_path))
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

# Create boxplot with colors based on transformation
sns.boxplot(
    data=all_results,
    x='Dataset',
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
