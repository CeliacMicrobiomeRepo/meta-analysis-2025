"""
Automatically plot the results of the best model between some selected datasets.

- plots a box and whisker plot of mean AUCs for each replicate of the best model for each dataset
- writes a png file
"""

# Imports ------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



# Inputs ------------------------------------------
# Analysis group name
ANALYSIS_GROUP_NAME = "stool_treated"
# Eval method
EVAL_METHOD = "logo"

# Dataset directory paths and names
DATASET_DIR_PATHS_AND_NAMES = [
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_after/", "After,CLR"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_before/", "Before,CLR"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_zscore_after/", "After,CLR,Z"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_zscore_before/", "Before,CLR,Z"),

    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_after/", "After,TSS"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_before/", "Before,TSS"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_zscore_after/", "After,TSS,Z"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_zscore_before/", "Before,TSS,Z"),

    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_after/", "After,Log10-sum"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_before/", "Before,Log10-sum"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_zscore_after/", "After,Log10-sum,Z"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_zscore_before/", "Before,Log10-sum,Z"),

    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_after/", "After,Log10"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_before/", "Before,Log10"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_zscore_after/", "After,Log10,Z"),
    ("/home/haig/Repos/meta-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_zscore_before/", "Before,Log10,Z"),
]
RESULTS_SUBDIR_FILE_NAME = EVAL_METHOD + "_results/best_models_results.tsv"
# Draws a red line at the highest median out of all datasets
TOP_MEDIAN_LINE_COLOR = True
# Maximum value for the y-axis
Y_AXIS_MAX = 1.0   # <- Make divisible by 0.05
# Minimum value for the y-axis
Y_AXIS_MIN = 0.5   # <- Make divisible by 0.05


# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "/home/haig/Repos/meta-analysis/machine_learning/results/transformation_comparison/"

# Create a color mapping for transformations
TRANSFORMATION_COLORS = {
    'TSS': '#1f77b4',  # blue
    'Log10': '#2ca02c',  # green
    'Log10-sum': '#d62728',  # red
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
    print(results_path)
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

if TOP_MEDIAN_LINE_COLOR:
    medians = all_results.groupby('Dataset')['AUC'].median()
    top_median = medians.max()
    plt.axhline(y=top_median, color='red', linestyle='--', linewidth=1)

# Customize the plot
plt.xticks(rotation=45, ha='right')
plt.xlabel('')
plt.ylabel('AUC')
plt.title('Best Model Performance Across Transformations')

# Set y-axis limits and ticks
plt.ylim(Y_AXIS_MIN, Y_AXIS_MAX)  # Set y-axis limits
plt.yticks(np.arange(Y_AXIS_MIN, Y_AXIS_MAX + 0.001, 0.05))

# Adjust layout to prevent label cutoff
plt.tight_layout()


# Export ------------------------------------------
# Save as PNG
plt.savefig(os.path.join(OUTPUT_DIR_PATH, EVAL_METHOD + '_' + ANALYSIS_GROUP_NAME + '.png'), dpi=300, bbox_inches='tight')

# Show plot 
plt.show()

plt.close()
