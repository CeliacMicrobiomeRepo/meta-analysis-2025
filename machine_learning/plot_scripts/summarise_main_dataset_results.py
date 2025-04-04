"""
Combine dataset results into a single TSV file.
The columns of the final tsv will be: 
Dataset, Body Site, Filter, Features, 
# Datasets, # Samples, # Positive Samples, # Negative Samples, 
LODO AUC, KFOLD AUC, Mean XSET AUC
"""

import os
import pandas as pd

# Inputs ------------------------------------------
# List of dataset directories with metadata
DATASET_DIR_PATHS_AND_NAMES = [
    ("./datasets/duodenal_noGFD_ASV_after_log10/", "Duodenal", "noGFD", "ASV"),
    ("./datasets/duodenal_noGFD_Genus_noV3_after_log10/", "Duodenal", "noGFD-noV3", "Genus"),
    ("./datasets/duodenal_noGFD_Genus_after_log10/", "Duodenal", "noGFD", "Genus"),
    ("./datasets/stool_wGFD_ASV_after_log10/", "Stool", "", "ASV"),
    ("./datasets/stool_wGFD_Genus_noV3_after_log10/", "Stool", "noV3", "Genus"),
    ("./datasets/stool_wGFD_Genus_after_log10/", "Stool", "", "Genus"),
    ("./datasets/stool_noGFD_ASV_after_log10/", "Stool", "noGFD", "ASV"),
    ("./datasets/stool_noGFD_Genus_noV3_after_log10/", "Stool", "noGFD-noV3", "Genus"),
    ("./datasets/stool_noGFD_Genus_after_log10/", "Stool", "noGFD", "Genus")
]

# Output ------------------------------------------
# Path to save the combined TSV file
OUTPUT_FILE_PATH = "./results_summary_main_datasets.tsv"

# Initialize list to collect all dataset information
combined_data = []

# Process each dataset directory
for dataset_dir, body_site, filter_, features in DATASET_DIR_PATHS_AND_NAMES:
    # Read summary.tsv
    summary_path = os.path.join(dataset_dir, "summary.tsv")
    summary_df = pd.read_csv(summary_path, sep='\t')
    
    # Extract required summary information
    dataset_name = summary_df.loc[0, 'Dataset Name']
    number_of_datasets = summary_df.loc[0, 'Number of Datasets Included']
    positive_samples = summary_df.loc[0, 'Positive Samples Included']
    negative_samples = summary_df.loc[0, 'Negative Samples Included']
    total_samples = summary_df.loc[0, 'Total Included Samples']
    
    # Read mean_xset_results.tsv and get maximum Mean AUC
    mean_xset_path = os.path.join(dataset_dir, "logo_results", "mean_xset_results.tsv")
    mean_xset_df = pd.read_csv(mean_xset_path, sep='\t')
    mean_xset_auc = mean_xset_df['Mean AUC'].max()
    
    # Read summary_best_model.tsv for LODO AUC
    lodo_summary_path = os.path.join(dataset_dir, "logo_results", "summary_best_model.tsv")
    lodo_df = pd.read_csv(lodo_summary_path, sep='\t', header=None)
    lodo_auc = float(lodo_df.iloc[0,0])
    
    # Read kfold_results/summary_best_model.tsv for KFOLD AUC
    kfold_summary_path = os.path.join(dataset_dir, "kfold_results", "summary_best_model.tsv")
    kfold_df = pd.read_csv(kfold_summary_path, sep='\t', header=None)
    kfold_auc = float(kfold_df.iloc[0,0])
    
    # Append the data to combined_data list
    combined_data.append({
        'Dataset': dataset_name,
        'Body Site': body_site,
        'Filter': filter_,
        'Features': features,
        '# Datasets': number_of_datasets,
        '# Samples': total_samples,
        '# Positive Samples': positive_samples,
        '# Negative Samples': negative_samples,
        'LODO AUC': lodo_auc,
        'KFOLD AUC': kfold_auc,
        'Mean XSET AUC': mean_xset_auc
    })

# Create a DataFrame from the combined data
combined_df = pd.DataFrame(combined_data)

# Reorder columns as specified
combined_df = combined_df[['Dataset', 'Body Site', 'Filter', 'Features',
                           '# Datasets', '# Samples', '# Positive Samples', '# Negative Samples',
                           'LODO AUC', 'KFOLD AUC', 'Mean XSET AUC']]

# Export the combined DataFrame to a TSV file
combined_df.to_csv(OUTPUT_FILE_PATH, sep='\t', index=False)