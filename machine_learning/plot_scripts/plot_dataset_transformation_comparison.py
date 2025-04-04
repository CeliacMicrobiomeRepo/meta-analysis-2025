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


FEATURE = "ASV"   # "ASV" or "Genus"


# Inputs ------------------------------------------
DATASET_DIR_PATHS_AND_NAMES = [
    ("./datasets/stool_wGFD_" + FEATURE + "_before_tss_zscore/", "Before,TSS,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_before_tss/", "Before,TSS"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_tss_zscore/", "After,TSS,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_tss/", "After,TSS"),
    ("./datasets/stool_wGFD_" + FEATURE + "_before_log10_zscore/", "Before,Log10,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_before_log10/", "Before,Log10"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_log10_zscore/", "After,Log10,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_log10/", "After,Log10"),
    ("./datasets/stool_wGFD_" + FEATURE + "_before_clr_zscore/", "Before,CLR,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_before_clr/", "Before,CLR"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_clr_zscore/", "After,CLR,Z"),
    ("./datasets/stool_wGFD_" + FEATURE + "_after_clr/", "After,CLR")
]
RESULTS_SUBDIR_FILE_NAME = "kfold_results/best_models_results.tsv"

# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "./plots/transformation_comparisons"

# Create a color mapping for transformations
TRANSFORMATION_COLORS = {
    'TSS': '#1f77b4',  # blue
    'Log10': '#2ca02c',  # green
    'CLR': '#ff7f0e'  # orange
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
        print(results_path)
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
plt.title('Best Model Performance Across Transformations for ' + FEATURE + ' Features')

# Set y-axis limits and ticks
plt.ylim(0.7, 0.95)  # Set y-axis limits
plt.yticks(np.arange(0.7, 1.0, 0.05))

# Adjust layout to prevent label cutoff
plt.tight_layout()


# Export ------------------------------------------
# Save as PNG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, 'transformations_' + FEATURE + '.png'), dpi=300, bbox_inches='tight')
# Save as SVG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, 'transformations_' + FEATURE + '.svg'), format='svg', bbox_inches='tight')

# Show plot 
plt.show()

plt.close()
