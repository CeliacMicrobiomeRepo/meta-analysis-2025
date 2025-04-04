"""
A script for construction of a dataset suitable for machine learning.

Uses the abundance_table.tsv and sample_data.tsv from the preprocessing of the truncated ASVs (output from the phylogseq_to_tsv.R script)

abundance_table.tsv:
 - Columns are ASVs (with headers the ASV sequences)
 - Rows are samples (with first column being the sample IDs)

sample_data.tsv
 - Columns are metadata (e.g. Dataset_ID, SRA_Run_ID, Group_Present_Study)
 - Rows are samples (with first column being the sample IDs)

Steps:
 - Filters samples according to options
 - Transforms abundances
 - Reduces metadata to essential components
 - Exports labels and features

"""


# SET UP ============================================

# Imports -------------------------------------------
import pandas as pd
import os
from composition_stats import ilr, clr, multiplicative_replacement
import numpy as np

# Main options ------------------------------------------
# Body site (duodenal or stool)
BODY_SITE = "duodenal"
# Include GFD samples?
WITH_GFD = False
# Transform with 'tss', 'log10', 'ilr', 'clr'
TRANFORMATION = 'log10'  
# Apply Z-score to transformed abundances (to each bacteria)
APPLY_ZSCORE = False             # False or True
# Filter before or after transformation?
WHEN_FILTERING = 'after'    # 'before' or 'after'

# Other options ------------------------------------------
# Pseudo count
PSEUDO_COUNT = 1e-6

# Input data ------------------------------------------
# TSV file containing metadata for each sample
SAMPLE_DATA_TSV_PATH = "./input_data/" + BODY_SITE + "_v4/unfiltered_sample_data.tsv"
# TSV file containing relative abundances of ASVs in samples
ASV_TABLE_TSV_PATH = "./input_data/" + BODY_SITE + "_v4/unfiltered_asv_table.tsv"    
# Column containing sample IDs
SAMPLE_ID_COLUMN = "Sample_ID"
# Column containing dataset IDs
DATASET_ID_COLUMN = "Dataset_ID"
# Columns to use as labels
LABEL_COLUMNS = ["Diagnosed_Celiac"]
# Extra features
INCLUDE_COLUMNS_AS_FEATURES = []


# Filters ------------------------------------------
# Exclude rows if
EXCLUDE_ROWS_WITH_VALUES = {
    SAMPLE_ID_COLUMN: ["-", "", "NA"],
    DATASET_ID_COLUMN: ["-", "", "NA"],
    "Gluten_Free_Diet": [] if WITH_GFD else [True]
}
# Only include rows if
ONLY_INCLUDE_ROWS_WITH_VALUES = {
    "Diagnosed_Celiac": [True, False],
    "Any_Significant_Factor": [False]
}
# Exclude specific datasets
EXCLUDE_DATASETS = [] # 16S_27_Federica 16S_80_Garcia
# Keep only specific columns
KEEP_ONLY_COLUMNS = [SAMPLE_ID_COLUMN, DATASET_ID_COLUMN] + LABEL_COLUMNS + INCLUDE_COLUMNS_AS_FEATURES
# Exclude samples
EXCLUDE_SAMPLES = ["SRR1107516", "ERR1551255", "ERR1551306", "SRR18231165", "SRR6885558"]
# Minimum number sample per class (prevent imbalance in classes)
MIN_SAMPLES_PER_CLASS = 5




# ASV filter ------------------------------
# Minimum average abundance across a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_AVERAGE_ABUNDANCE_IN_DATASET = 0.001
# Minimum proportion of all samples in a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_PREVALENCE_IN_SAMPLES_IN_DATASET = 0.1
# ...X proportion of the datasets:
IN_N_PROPORTION_DATASETS = 0.3




# Output options ------------------------------------------
# Output directory path
# OUTPUT_DIR_PATH = "./datasets/stool_wGFD_ASV_after_tss_zscore/"  
OUTPUT_DIR_PATH = f"./datasets/{BODY_SITE}_{'wGFD' if WITH_GFD else 'noGFD'}_ASV_{WHEN_FILTERING}_{TRANFORMATION}{'_zscore' if APPLY_ZSCORE else ''}/"
# Filtered samples metadata file path
FILTERED_SAMPLES_TSV_OUTPUT_PATH = os.path.join(OUTPUT_DIR_PATH, "sample_labels.tsv")
# ASV abundance matrix file path
ASV_ABUNDANCE_MATRIX_TSV_OUTPUT_PATH = os.path.join(OUTPUT_DIR_PATH, "sample_asv_abundances.tsv")
# About file location
ABOUT_PATH = os.path.join(OUTPUT_DIR_PATH, "about_dataset.txt")



def apply_total_sum_scaling(abundances_df):
    """
    Normalise/transform the abundances of an abundance matrix so that: 
        - The all abundances in a sample sum to 1.0

    - Called total sum scaling in this paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160169
    
    (rows are samples, columns are taxonomic units)
    """
    return abundances_df.div(abundances_df.sum(axis=1), axis=0)

def apply_zscore(abundances_df):
    """
    Normalise/transform the abundances of an abundance matrix so that: 
        - The features have a mean of zero and a unit variance across all samples

    - Subtracts the means of features, then divides by SDs of features
    
    (rows are samples, columns are taxonomic units)
    """
    # Calculate the mean and standard deviation for each column (feature)
    feature_means = abundances_df.mean(axis=0)
    feature_stds = abundances_df.std(axis=0, ddof=0)  # ddof=0 for population standard deviation

    # Avoid division by zero for features with zero standard deviation
    feature_stds = feature_stds.replace(0, 1e-6)

    # Normalize the data (z-score transformation)
    normalized_df = (abundances_df - feature_means) / feature_stds

    return normalized_df

def transform_log(abundances_df):
    """
    Apply log10 transformation to the abundances in an abundance matrix.
    1. Log10 transformation to all values (with a small constant to avoid log(0))
    2. Divide all values by the sum of the values in the sample

    Args:
        abundances_df (pd.DataFrame): A DataFrame where each row corresponds 
                                      to a sample and each column corresponds 
                                      to a taxonomic unit's abundance.

    Returns:
        pd.DataFrame: A DataFrame with the same shape as the input, where each 
                      value is log-transformed.
    """
    # Add a small constant to avoid log(0)
    log_transformed_df = abundances_df.apply(lambda col: np.log(col + 1e-6))

    # Scale by the sample total
    log_transformed_df = log_transformed_df.div(log_transformed_df.sum(axis=1), axis=0)
    
    return log_transformed_df

def transform_clr(abundances_df):
    """
    Apply the Centered Log-Ratio (CLR) transformation to the abundances in an abundance matrix.

    The transformation is applied to each row of the input DataFrame, which 
    represents the abundances of different taxonomic units in a sample.

    Prior to transformation, zero abundances are replaced with a small value using multiplicative_replacement.
    """

    # Replace zero abundances with a small value using multiplicative_replacement
    abundances_array = multiplicative_replacement(abundances_df.values, PSEUDO_COUNT)

    # Convert back to DataFrame
    abundances_df = pd.DataFrame(abundances_array, index=abundances_df.index, columns=abundances_df.columns)

    # Apply the clr transformation to each row
    transformed_df = abundances_df.apply(lambda row: pd.Series(clr(row.values), index=abundances_df.columns), axis=1)
    
    return transformed_df

def transform_ilr(abundances_df):
    """
    Apply the Isometric Log-Ratio (ILR) transformation to the abundances in an abundance matrix.

    The transformation is applied to each row of the input DataFrame, which 
    represents the abundances of different taxonomic units in a sample.

    Prior to transformation, zero abundances are replaced with a small value using multiplicative_replacement.
    """

    # Replace zero abundances with a small value using multiplicative_replacement
    abundances_array = multiplicative_replacement(abundances_df.values, PSEUDO_COUNT)

    # Apply the ILR transformation to each row and store the results
    transformed_array = np.apply_along_axis(ilr, axis=1, arr=abundances_array)

    # Create a DataFrame with one less column than the input
    transformed_df = pd.DataFrame(
        transformed_array,
        index=abundances_df.index,
        columns=[f"FEATURE_{i}" for i in range(transformed_array.shape[1])]
    )

    return transformed_df



# MAIN ================================================
# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR_PATH):
    os.makedirs(OUTPUT_DIR_PATH)



# Read data ------------------------------------
# This should be read counts

# Read samples metadata
samples_df = pd.read_csv(SAMPLE_DATA_TSV_PATH, sep="\t")
# Read asv table
abundance_table_df = pd.read_csv(ASV_TABLE_TSV_PATH, sep="\t")



# Filter samples ------------------------------------

print(f"\nInitial number of samples: {len(samples_df)}")

samples_to_remove = []

# Exclude rows with specific values
for column, excluded_values in EXCLUDE_ROWS_WITH_VALUES.items():
    excluded_samples = samples_df[samples_df[column].isin(excluded_values)][SAMPLE_ID_COLUMN].tolist()
    samples_to_remove.extend(excluded_samples)
    if excluded_samples:
        print(f"Excluding {len(excluded_samples)} samples due to excluded values in '{column}'.")

# Only include rows with specific values 
for column, required_values in ONLY_INCLUDE_ROWS_WITH_VALUES.items():
    excluded_samples = samples_df[~samples_df[column].isin(required_values)][SAMPLE_ID_COLUMN].tolist()
    samples_to_remove.extend(excluded_samples)
    if excluded_samples:
        print(f"Excluding {len(excluded_samples)} samples not matching required values in '{column}'.")

# Exclude specific datasets
if EXCLUDE_DATASETS:
    excluded_dataset_samples = samples_df[samples_df[DATASET_ID_COLUMN].isin(EXCLUDE_DATASETS)][SAMPLE_ID_COLUMN].tolist()
    samples_to_remove.extend(excluded_dataset_samples)
    print(f"Excluding {len(excluded_dataset_samples)} samples from excluded datasets.")

# Exclude specific samples
if EXCLUDE_SAMPLES:
    samples_to_remove.extend(EXCLUDE_SAMPLES)
    print(f"Excluding {len(EXCLUDE_SAMPLES)} specific samples (if present).")

# Remove duplicates in samples_to_remove
samples_to_remove = list(set(samples_to_remove))

# Keep only specific columns
samples_df = samples_df[KEEP_ONLY_COLUMNS]

# Remove samples in samples_to_remove from both dataframes
samples_df = samples_df[~samples_df[SAMPLE_ID_COLUMN].isin(samples_to_remove)].reset_index(drop=True)
abundance_table_df = abundance_table_df[~abundance_table_df[SAMPLE_ID_COLUMN].isin(samples_to_remove)].reset_index(drop=True)

# Exclude datasets with less than MIN_SAMPLES_PER_CLASS samples (!! ONLY DOES LABEL_COLUMNS[0])
class_counts_per_dataset = samples_df.groupby(DATASET_ID_COLUMN)[LABEL_COLUMNS[0]].value_counts().unstack(fill_value=0)
datasets_to_exclude = class_counts_per_dataset[
    (class_counts_per_dataset < MIN_SAMPLES_PER_CLASS).any(axis=1)
].index.tolist()

if datasets_to_exclude:
    # Get samples to remove
    samples_to_remove_min_class = samples_df[
        samples_df[DATASET_ID_COLUMN].isin(datasets_to_exclude)
    ][SAMPLE_ID_COLUMN].tolist()
    
    # Remove samples from both dataframes
    samples_df = samples_df[~samples_df[DATASET_ID_COLUMN].isin(datasets_to_exclude)].reset_index(drop=True)
    abundance_table_df = abundance_table_df[~abundance_table_df[SAMPLE_ID_COLUMN].isin(samples_to_remove_min_class)].reset_index(drop=True)
    
    print(f"Excluding {len(samples_to_remove_min_class)} samples from {len(datasets_to_exclude)} datasets with less than {MIN_SAMPLES_PER_CLASS} samples per class.")

print(f"Final number of samples: {len(samples_df)}\n")


asv_abundance_matrix = abundance_table_df.set_index(SAMPLE_ID_COLUMN)



# Determine which ASVs should be filtered -----------------------------------------
# This simply decides which ASVs to remove based on abundances and prevalences
# Actual removal is done before/after transformation depending on FILTER_BEFORE_TRANSFORMATION
tss_asv_abundance_matrix = apply_total_sum_scaling(asv_abundance_matrix.copy())

# Calculate the threshold number of datasets
num_datasets = len(samples_df[DATASET_ID_COLUMN].unique())
dataset_threshold = round(IN_N_PROPORTION_DATASETS * num_datasets)
print(f"Filtering {len(asv_abundance_matrix.columns)} ASVs that do not meet the abundance and prevalence conditions in at least {IN_N_PROPORTION_DATASETS * 100}% of datasets.")
print(f"This corresponds to a minimum of {dataset_threshold} datasets out of {num_datasets}.")

# Calculate the average abundance for each ASV in each dataset
dataset_avg_abundances = tss_asv_abundance_matrix.reset_index().merge(
    samples_df[[SAMPLE_ID_COLUMN, DATASET_ID_COLUMN]], 
    on=SAMPLE_ID_COLUMN
).drop(SAMPLE_ID_COLUMN, axis=1).groupby(DATASET_ID_COLUMN).mean()

# Get indices of ASVs to remove based on abundance
asvs_meeting_abundance_threshold = (dataset_avg_abundances >= MIN_AVERAGE_ABUNDANCE_IN_DATASET).sum() >= dataset_threshold
asvs_to_remove_abundance = dataset_avg_abundances.columns[~asvs_meeting_abundance_threshold].tolist()

# Print number of ASVs to remove from abundance filter
print(f"Number of ASVs to remove based on abundance threshold: {len(asvs_to_remove_abundance)}")

# Calculate the prevalence for each ASV in each dataset
dataset_prevalences = tss_asv_abundance_matrix.reset_index().merge(
    samples_df[[SAMPLE_ID_COLUMN, DATASET_ID_COLUMN]], 
    on=SAMPLE_ID_COLUMN
).drop(SAMPLE_ID_COLUMN, axis=1).groupby(DATASET_ID_COLUMN).apply(
    lambda x: (x > 0).sum() / len(x), 
    include_groups=False  # Add this parameter to address the warning
)

# Get indices of ASVs to remove based on prevalence
asvs_meeting_prevalence_threshold = (dataset_prevalences >= MIN_PREVALENCE_IN_SAMPLES_IN_DATASET).sum() >= dataset_threshold
asvs_to_remove_prevalence = dataset_prevalences.columns[~asvs_meeting_prevalence_threshold].tolist()

# Print number of ASVs to remove from prevalence filter
print(f"Number of ASVs to remove based on prevalence threshold: {len(asvs_to_remove_prevalence)}")

# Print the number of ASVs to be removed from both (excluding duplicates)
all_asvs_to_remove = list(set(asvs_to_remove_abundance + asvs_to_remove_prevalence))
if 'Dataset_ID' in all_asvs_to_remove:
    all_asvs_to_remove.remove('Dataset_ID')
print(f"Number of ASVs to remove based on abundance OR prevalence threshold: {len(all_asvs_to_remove)}")
print(f"...Leaving {len(asv_abundance_matrix.columns) - len(all_asvs_to_remove)} ASVs in the dataset.")




# Filter ASVs before transformation -----------------------------------------
# This removes ASVs, the data is kept as untransformed read counts

if WHEN_FILTERING == 'before': 

    # Remove ASVs
    asv_abundance_matrix = asv_abundance_matrix.drop(columns=all_asvs_to_remove)


# Normalisation and transformation ------------------------------------------
# These are applied to raw read counts

# Transformation
if TRANFORMATION == 'tss':
    asv_abundance_matrix = apply_total_sum_scaling(asv_abundance_matrix)

elif TRANFORMATION == 'log10':
    asv_abundance_matrix = transform_log(asv_abundance_matrix)

elif TRANFORMATION == 'ilr':
    asv_abundance_matrix = transform_ilr(asv_abundance_matrix)

elif TRANFORMATION == 'clr':
    asv_abundance_matrix = transform_clr(asv_abundance_matrix)


# Z-score bacteria
if APPLY_ZSCORE:
    asv_abundance_matrix = apply_zscore(asv_abundance_matrix)



# Filter ASVs after transformation -----------------------------------------
# This removes ASVs, the data is kept as it is

if WHEN_FILTERING == 'after': 

    # Remove ASVs
    asv_abundance_matrix = asv_abundance_matrix.drop(columns=all_asvs_to_remove)




# Add extra features to the abundance matrix dataframe ----------------------------------
if INCLUDE_COLUMNS_AS_FEATURES:
    # Encode each feature column as an integer starting at 0
    for feature in INCLUDE_COLUMNS_AS_FEATURES:
        samples_df[feature] = samples_df[feature].astype('category').cat.codes

    # Merge the encoded features into the taxonomic abundance matrix
    features_df = samples_df.set_index(SAMPLE_ID_COLUMN)[INCLUDE_COLUMNS_AS_FEATURES]
    asv_abundance_matrix = features_df.join(asv_abundance_matrix)
    print(f"Added extra features ({INCLUDE_COLUMNS_AS_FEATURES}) to the taxonomic abundance matrix")




# Summary of labels ------------------------------------------

# Print the number of samples with each label in each dataset
print("\nNumber of samples with each label in each dataset:")
for dataset_id, group in samples_df.groupby(DATASET_ID_COLUMN):
    print(dataset_id + ":")
    label_counts = group[LABEL_COLUMNS[0]].value_counts()
    for label, count in label_counts.items():
        print(f"  {label}: {count}")

# Print the number of samples with each label overall
print("\nTotal:")
overall_label_counts = samples_df[LABEL_COLUMNS[0]].value_counts()
for label, count in overall_label_counts.items():
    print(f"  {label}: {count}")

print()


# Write dataset to file ------------------------------------------
# Features
asv_abundance_matrix.to_csv(ASV_ABUNDANCE_MATRIX_TSV_OUTPUT_PATH, sep="\t")
# Labels
samples_df.to_csv(FILTERED_SAMPLES_TSV_OUTPUT_PATH, sep="\t", index=False)

# Write about.txt
with open(ABOUT_PATH, "w") as about_file:
    about_file.write("Dataset Construction Summary\n")
    about_file.write("Input Files:\n")
    about_file.write(f"  - Sample Data: {SAMPLE_DATA_TSV_PATH}\n")
    about_file.write(f"  - ASV Table: {ASV_TABLE_TSV_PATH}\n\n")
    about_file.write("Output Files:\n")
    about_file.write(f"  - Filtered Samples Metadata: {FILTERED_SAMPLES_TSV_OUTPUT_PATH}\n")
    about_file.write(f"  - ASV Abundance Matrix: {ASV_ABUNDANCE_MATRIX_TSV_OUTPUT_PATH}\n\n")
    about_file.write("Filters Applied:\n")
    about_file.write("  - Exclude rows if:\n")
    for column, values in EXCLUDE_ROWS_WITH_VALUES.items():
        about_file.write(f"    - {column}: {values}\n")
    about_file.write("  - Only include rows if:\n")
    for column, values in ONLY_INCLUDE_ROWS_WITH_VALUES.items():
        about_file.write(f"    - {column}: {values}\n")
    about_file.write(f"  - Exclude datasets: {EXCLUDE_DATASETS}\n")
    about_file.write(f"  - Exclude samples: {EXCLUDE_SAMPLES}\n\n")
    about_file.write("Transformation Options:\n")
    about_file.write(f"  - Pseudo count: {PSEUDO_COUNT}\n")
    about_file.write(f"  - Transformation: {TRANFORMATION}\n\n")
    about_file.write("Additional Options:\n")
    about_file.write(f"  - Additional Features: {INCLUDE_COLUMNS_AS_FEATURES}\n\n")
    about_file.write("Summary:\n")
    about_file.write(f"  - Initial number of samples: {len(samples_df) + len(samples_to_remove)}\n")
    about_file.write(f"  - Final number of samples: {len(samples_df)}\n")
    about_file.write(f"  - Samples removed: {len(samples_to_remove)}\n")
    about_file.write("Sample Label Details:\n")
    about_file.write("Overall label distribution:\n")
    for label, count in overall_label_counts.items():
        about_file.write(f"  - {label}: {count}\n")
    about_file.write("\nLabel distribution by dataset:\n")
    for dataset_id, group in samples_df.groupby(DATASET_ID_COLUMN):
        about_file.write(f"  Dataset {dataset_id}:\n")
        label_counts = group[LABEL_COLUMNS[0]].value_counts()
        for label, count in label_counts.items():
            about_file.write(f"    - {label}: {count}\n")
# Summary TSV
summary_file_path = os.path.join(OUTPUT_DIR_PATH, "summary.tsv")
summary_data = {
    "Dataset Name": [os.path.basename(OUTPUT_DIR_PATH.rstrip('/'))],
    "Number of Datasets Included": [samples_df[DATASET_ID_COLUMN].nunique()],
    "Positive Samples Included": [overall_label_counts.get(True, 0)],
    "Negative Samples Included": [overall_label_counts.get(False, 0)],
    "Total Included Samples": [len(samples_df)],
    "Total Excluded Samples": [len(samples_to_remove)]
}
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file_path, sep="\t", index=False)



print("Done!")
