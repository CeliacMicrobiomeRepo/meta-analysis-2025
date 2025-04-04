"""
A script for construction of a dataset suitable for machine learning.

Uses assigned taxonomies from full length ASVs (output from a dada2_pipeline_X.R script)

It assumes a specific directory structure and file format.
 - Takes a file of samples and filters it according to the options below
 - Reads all the taxonomic data for the samples
 - Constructs a taxonomic abundance matrix
 - Filters the taxonomic abundance matrix according to the options below
 - Writes the filtered taxonomic abundance matrix and filtered samples metadata to file
"""

# SET UP ============================================

# Imports -------------------------------------------
import pandas as pd
import os
from composition_stats import ilr, clr, multiplicative_replacement
import numpy as np
from datetime import datetime



# Main options ------------------------------------------
# Body site (duodenal or stool)
BODY_SITE = "stool"
# Include GFD samples?
WITH_GFD = True
# Transform with 'tss', 'log10', 'ilr', 'clr'
TRANFORMATION = 'log10'  
# Apply Z-score to transformed abundances (to each bacteria)
APPLY_ZSCORE = False             # False or True
# Filter before or after transformation?
WHEN_FILTERING = 'after'    # 'before' or 'after'
# Rank at which to give labels
LABEL_RANK = 5      # 0 = Kingdom, 1 = Phylum, 2 = Class, 3 = Order, 4 = Family, 5 = Genus, 6 = Species
# Invalid taxonomic labels
INVALID_TAX_LABELS = [None, "NA", ""]
# Invalid strings for taxonomic labels to contain
INVALID_TAX_LABELS_TO_CONTAIN = ["unknown"]
# Exclude V3 amplicons
NO_V3 = True

# Other options ------------------------------------------
# Pseudo count
PSEUDO_COUNT = 1e-6

# Input data ------------------------------------------
# TSV file containing all samples with whatever metadata
SAMPLES_TSV_PATH = "./input_data/all_samples.tsv"
# Column containing sample IDs
SAMPLE_ID_COLUMN = "Sample_ID"
# Column containing dataset IDs
DATASET_ID_COLUMN = "Dataset_ID"
# Columns to use as labels
LABEL_COLUMNS = ["Diagnosed_Celiac"]
# Extra features
INCLUDE_COLUMNS_AS_FEATURES = []
# Data directory path
DATA_DIR_PATH = "D:/microbiome_sequencing_datasets/celiac_16s_datasets/"

# Filters ------------------------------------------
# Exclude rows if
EXCLUDE_ROWS_WITH_VALUES = {
    SAMPLE_ID_COLUMN: ["-", "", "NA"],
    DATASET_ID_COLUMN: ["-", "", "NA"],
    "Diagnosed_Celiac": ["-", "", "NA"],
    "Gluten_Free_Diet": ["-", "", "NA"],
    "Amplicon_Region": ["V3"] if NO_V3 else [],
}
# Only include rows if
ONLY_INCLUDE_ROWS_WITH_VALUES = {
    "Sample_Site": [BODY_SITE], 
    "Sequencing_Type": ["16S"],
    "Any_Significant_Factor": [False],
    "Gluten_Free_Diet": [False, True] if WITH_GFD else [False]
}
# Exclude specific datasets
EXCLUDE_DATASETS = []
# Keep only specific columns
KEEP_ONLY_COLUMNS = [SAMPLE_ID_COLUMN, DATASET_ID_COLUMN] + LABEL_COLUMNS + INCLUDE_COLUMNS_AS_FEATURES
# Exclude samples
EXCLUDE_SAMPLES = ["SRR1107516", "ERR1551255", "ERR1551306", "SRR18231165", "SRR6885558"]
# Minimum number sample per class (prevent imbalance in classes)
MIN_SAMPLES_PER_CLASS = 5

# Taxa filter ------------------------------
# Minimum average abundance across a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_AVERAGE_ABUNDANCE_IN_DATASET = 0.001
# Minimum proportion of all samples in a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_PREVALENCE_IN_SAMPLES_IN_DATASET = 0.1
# ...X proportion of the datasets:
IN_N_PROPORTION_DATASETS = 0.3



# Output options ------------------------------------------
# Output directory path
# OUTPUT_DIR_PATH = "./datasets/stool_wGFD_Genus_after_tss_zscore/"  
OUTPUT_DIR_PATH = f"./datasets/{BODY_SITE}_{'wGFD' if WITH_GFD else 'noGFD'}_Genus{'_noV3' if NO_V3 else ''}_{WHEN_FILTERING}_{TRANFORMATION}{'_zscore' if APPLY_ZSCORE else ''}/"
# Filtered samples metadata file path
FILTERED_SAMPLES_TSV_PATH = os.path.join(OUTPUT_DIR_PATH, "sample_labels.tsv")
# Taxonomic abundance matrix file path
TAXONOMIC_ABUNDANCE_MATRIX_TSV_PATH = os.path.join(OUTPUT_DIR_PATH, "sample_taxa_abundances.tsv")
# About file location
ABOUT_PATH = os.path.join(OUTPUT_DIR_PATH, "about_dataset.txt")



# FUNCTIONS ===========================================

def do_samples_exist(samples, data_dir_path, transposed=False):
    """
    Check that every sample exists in the asv_abundances.tsv file in data_dir_path.
    If transposed=True, assume that the samples are in the column headers and the asvs are in the rows.
    """
    # Get asv abundances file path
    asv_abundances_file_path = os.path.join(data_dir_path, "asv_abundances.tsv")
    # Read asv abundances
    asv_abundances_df = pd.read_csv(asv_abundances_file_path, sep="\t", index_col=0)
    # Get samples in asv abundances file
    if transposed:
        found_samples = asv_abundances_df.columns.tolist()
    else:
        found_samples = asv_abundances_df.index.tolist()
    # Check that every sample exists in the asv abundances file
    all_exist = all([sample in found_samples for sample in samples])
    return all_exist

def read_all_taxa_as_dictionary(dataset_dirs):
    """
    Read taxonomic abundances for all samples and combine into a single dictionary.
        Assumes there is a taxonomy.tsv file in each dataset directory.
        A snippet of the taxonomy.tsv files:
                Kingdom	Phylum	Class	Order	Family	Genus	Species
            CCTACGGGTGGC	Bacteria	Firmicutes	Clostridia	Lachnospirales	Lachnospiraceae	Blautia	NA
            CCTACGGGTGGCCTAG	Bacteria	Firmicutes	Clostridia	Lachnospirales	Lachnospiraceae	Agathobacter	NA
            CCTACGGGTGGGCCG	Bacteria	Firmicutes	Clostridia	Oscillospirales	Ruminococcaceae	Subdoligranulum	NA

    Args:
        dataset_dirs (list): List of all dataset directories

    Returns:
        dict: A dictionary mapping ASVs to their taxonomic labels.
            ASVs are represented as tuples like (dataset_dir, asv_sequence)
            Taxonomic labels are represented as tuples with labels in this order:
                ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    """
    taxa_dict = {}
    for dataset_dir in dataset_dirs:
        # Construct the path to the taxonomy.tsv file
        taxonomy_file_path = os.path.join(DATA_DIR_PATH, dataset_dir, "taxonomy.tsv")
        
        # Read the taxonomy.tsv file
        if os.path.exists(taxonomy_file_path):
            taxonomy_df = pd.read_csv(taxonomy_file_path, sep="\t", index_col=0)
            
            # Iterate over each ASV in the taxonomy file
            for asv_sequence, row in taxonomy_df.iterrows():
                # Create a tuple for the taxonomic labels
                taxonomic_labels = tuple(row)

                # Replace nan with None
                taxonomic_labels = tuple(label if pd.notna(label) else None for label in taxonomic_labels)
                
                # Ensure a length of 7 by padding with Nones
                taxonomic_labels = taxonomic_labels + (None,) * (7 - len(taxonomic_labels))

                # Map the ASV to its taxonomic labels in the dictionary
                taxa_dict[(dataset_dir, asv_sequence)] = taxonomic_labels
        else:
            print(f"Warning: {taxonomy_file_path} does not exist.")
    
    return taxa_dict

def is_tax_label_valid(tax_label):
    """
    Check if a taxonomic label is valid.
    """
    # Check that the taxonomic label is not invalid
    if tax_label in INVALID_TAX_LABELS:
        return False
    # Check that the taxonomic label does not contain invalid strings
    if any(invalid_string in tax_label.lower() for invalid_string in INVALID_TAX_LABELS_TO_CONTAIN):
        return False
    # Valid!
    return True

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

def construct_taxonomic_abundance_matrix(all_samples, all_taxa_dictionary, rank):
    """
    Construct a taxonomic abundance matrix.

    Args:
        all_samples (list): A list of tuples, where each tuple contains a 
                            dataset directory and a sample ID.
        all_taxa_dictionary (dict): A dictionary mapping ASVs to their 
                                    taxonomic labels. ASVs are represented 
                                    as tuples like (dataset_dir, asv_sequence), 
                                    and taxonomic labels are tuples with 
                                    labels in this order: 
                                    ("Kingdom", "Phylum", "Class", "Order", 
                                    "Family", "Genus", "Species").

    Returns:
        pd.DataFrame: A DataFrame representing the taxonomic abundance matrix, 
                      with samples as rows and taxonomic units as columns.
                      - The column headers are the tuple format of taxonomic labels.
                      - The index is the sample IDs (strings)
    """

    # Reduce the taxonomic labels to the desired rank (may include Nones !!!)
    all_taxa_dictionary = {asv: taxonomic_labels[rank] for asv, taxonomic_labels in all_taxa_dictionary.items()}
    
    # Initialize an empty DataFrame with samples as the index
    taxonomic_abundance_matrix = pd.DataFrame(index=[sample_id for _, sample_id in all_samples], dtype='float64')

    # Set the columns of the empty matrix to be the unique taxonomic labels from the dictionary
    unique_taxonomic_labels = list(set([taxonomic_label for taxonomic_label in all_taxa_dictionary.values()]))
    taxonomic_abundance_matrix = taxonomic_abundance_matrix.reindex(columns=unique_taxonomic_labels, fill_value=0.0)

    # For each sample (dataset by dataset), get the ASVs and add their abundances to the correct place in the matrix
    dataset_dir_of_loaded_abundances_df = None
    for i, (dataset_dir, sample_id) in enumerate(sorted(all_samples)):
        # Print progress every 100 samples
        if i % 50 == 0:
            print(f"Processed {i}/{len(all_samples)} samples")

        # Load taxonomic abundances for this sample
        if dataset_dir_of_loaded_abundances_df != dataset_dir:
            # Next dataset
            dataset_dir_of_loaded_abundances_df = dataset_dir
            # The columns are the ASV sequences and the rows are the sample IDs (the first column is the sample ID)
            abundances_df = pd.read_csv(os.path.join(DATA_DIR_PATH, dataset_dir, "asv_abundances.tsv"), sep="\t", index_col=0)
            # Replace the ASV sequences with the taxonomic labels
            abundances_df.columns = abundances_df.columns.map(lambda asv_seq: all_taxa_dictionary[(dataset_dir, asv_seq)])

        # For every ASV in this sample, add its abundance to the matrix
        # (group by the taxonomy and sum the abundances to handle duplicate columns)
        summed_abundances = abundances_df.loc[sample_id].groupby(level=0).sum()
        
        # Update the taxonomic abundance matrix with the summed abundances
        taxonomic_abundance_matrix.loc[sample_id, summed_abundances.index] += summed_abundances.values

    # Sort columns by total in descending order
    taxonomic_abundance_matrix = taxonomic_abundance_matrix.reindex(columns=taxonomic_abundance_matrix.sum(axis=0).sort_values(ascending=False).index)

    # Remove columns with 0 total abundance
    taxonomic_abundance_matrix = taxonomic_abundance_matrix.loc[:, (taxonomic_abundance_matrix.sum(axis=0) != 0)]

    return taxonomic_abundance_matrix


# MAIN ================================================
# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR_PATH):
    os.makedirs(OUTPUT_DIR_PATH)



# Read and filter samples ------------------------------------
# Read samples metadata
try:
    samples_df = pd.read_csv(SAMPLES_TSV_PATH, sep="\t", encoding='utf-8')
except UnicodeDecodeError:
    # Try different encodings if UTF-8 fails
    samples_df = pd.read_csv(SAMPLES_TSV_PATH, sep="\t", encoding='latin1')


# Filter samples ------------------------------------
# Filter samples
num_inital_samples = len(samples_df)
print(f"\nInitial number of samples: {num_inital_samples}")

# Exclude rows with specific values
for col, values in EXCLUDE_ROWS_WITH_VALUES.items():
    samples_df = samples_df[~samples_df[col].isin(values)]
    print(f"After excluding {col} with values {values}: {len(samples_df)} samples")

# Only include rows with specific values 
for col, values in ONLY_INCLUDE_ROWS_WITH_VALUES.items():
    samples_df = samples_df[samples_df[col].isin(values)]
    print(f"After including only {col} with values {values}: {len(samples_df)} samples")

# Exclude specific datasets
samples_df = samples_df[~samples_df[DATASET_ID_COLUMN].isin(EXCLUDE_DATASETS)]
print(f"After excluding datasets {EXCLUDE_DATASETS}: {len(samples_df)} samples")

# Exclude specific samples
samples_df = samples_df[~samples_df[SAMPLE_ID_COLUMN].isin(EXCLUDE_SAMPLES)]
print(f"After excluding samples {EXCLUDE_SAMPLES}: {len(samples_df)} samples")

# Exclude datasets with less than MIN_SAMPLES_PER_CLASS samples (!! ONLY DOES LABEL_COLUMNS[0])
class_counts_per_dataset = samples_df.groupby(DATASET_ID_COLUMN)[LABEL_COLUMNS[0]].value_counts().unstack(fill_value=0)
datasets_to_exclude = class_counts_per_dataset[
    (class_counts_per_dataset < MIN_SAMPLES_PER_CLASS).any(axis=1)
].index.tolist()
samples_df = samples_df[~samples_df[DATASET_ID_COLUMN].isin(datasets_to_exclude)]
print(f"After excluding datasets with less than {MIN_SAMPLES_PER_CLASS} samples: {len(samples_df)} samples")

num_final_samples = len(samples_df)
print(f"Final number of samples: {num_final_samples}\n")

# Keep only specific columns
samples_df = samples_df[KEEP_ONLY_COLUMNS]


# Checks ------------------------------------
# Check that SAMPLE_ID_COLUMN is unique for every row
if not samples_df[SAMPLE_ID_COLUMN].is_unique:
    raise ValueError(f"{SAMPLE_ID_COLUMN} is not unique for every row")
# Check that every unique DATASET_ID_COLUMN is a real directory in DATA_DIR_PATH
if not samples_df[DATASET_ID_COLUMN].isin(os.listdir(DATA_DIR_PATH)).all():
    raise ValueError(f"Every unique {DATASET_ID_COLUMN} must be a real directory in {DATA_DIR_PATH}")
# Check that every row has a value in LABEL_COLUMNS
for col in LABEL_COLUMNS:
    if samples_df[col].isnull().all().any():
        raise ValueError(f"Every row must have a value in {col}")
# Check that every row has a value in INCLUDE_COLUMNS_AS_FEATURES
if INCLUDE_COLUMNS_AS_FEATURES:
    for col in INCLUDE_COLUMNS_AS_FEATURES:
        if samples_df[col].isnull().all().any():
            raise ValueError(f"Every row must have a value in {col}")
# Check that every sample exists in DATA_DIR_PATH
for dataset_id in samples_df[DATASET_ID_COLUMN].unique():
    # Get samples for this dataset
    this_dataset_df = samples_df[samples_df[DATASET_ID_COLUMN] == dataset_id]
    sample_ids_list = this_dataset_df[SAMPLE_ID_COLUMN].tolist()
    # Check that every sample exists in DATA_DIR_PATH
    if not do_samples_exist(sample_ids_list, DATA_DIR_PATH + dataset_id):
        raise ValueError(f"Every sample must exist in {DATA_DIR_PATH}")
print("All checks passed.\n")


# Read all taxa ------------------------------------------
# Get all dataset directories
all_dataset_dirs = [os.path.join(DATA_DIR_PATH, dataset_id) for dataset_id in samples_df[DATASET_ID_COLUMN].unique()]
# Read all taxa as a dictionary (maps ASVs to their taxonomic labels - like: (dataset_dir, asv_sequence) -> ('Bacteria', 'Firmicutes', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 'Blautia', None))
all_taxa_dictionary = read_all_taxa_as_dictionary(all_dataset_dirs)


# Construct taxanomic abundance matrix ------------------------------------------
# Get all samples as tuples like (dataset_dir, sample_id)
all_samples = [(DATA_DIR_PATH + dataset_id, sample_id) for dataset_id in samples_df[DATASET_ID_COLUMN].unique() for sample_id in samples_df[samples_df[DATASET_ID_COLUMN] == dataset_id][SAMPLE_ID_COLUMN].unique()]
# Construct taxanomic abundance matrix (column headers are the taxonomic labels at rank LABEL_RANK)
taxonomic_abundance_matrix = construct_taxonomic_abundance_matrix(all_samples, all_taxa_dictionary, LABEL_RANK)

# Remove columns that are invalid according to is_tax_label_valid
taxonomic_abundance_matrix = taxonomic_abundance_matrix.loc[:, taxonomic_abundance_matrix.columns.map(is_tax_label_valid)]



# Determine which taxa should be filtered -----------------------------------------
# This simply decides which taxa to remove based on abundances and prevalences
# Actual removal is done before/after transformation depending on WHEN_FILTERING
print(f"\nInitial number of taxa: {taxonomic_abundance_matrix.shape[1]}")

# Create TSS version for filtering calculations
tss_taxonomic_abundance_matrix = apply_total_sum_scaling(taxonomic_abundance_matrix.copy())

# Make the index into a column
tss_taxonomic_abundance_matrix = tss_taxonomic_abundance_matrix.reset_index().rename(columns={'index': SAMPLE_ID_COLUMN})


# Calculate the threshold number of datasets
num_datasets = len(samples_df[DATASET_ID_COLUMN].unique())
dataset_threshold = round(IN_N_PROPORTION_DATASETS * num_datasets)
print(f"\nFiltering taxa that do not meet the abundance and prevalence conditions in at least {IN_N_PROPORTION_DATASETS * 100}% of datasets.")
print(f"This corresponds to a minimum of {dataset_threshold} datasets out of {num_datasets}.")

# Calculate the average abundance for each taxon in each dataset
dataset_avg_abundances = tss_taxonomic_abundance_matrix.reset_index().merge(
    samples_df[[SAMPLE_ID_COLUMN, DATASET_ID_COLUMN]], 
    on=SAMPLE_ID_COLUMN
).drop(SAMPLE_ID_COLUMN, axis=1).groupby(DATASET_ID_COLUMN).mean()

# Get taxa to remove based on abundance threshold
taxa_meeting_abundance_threshold = (dataset_avg_abundances >= MIN_AVERAGE_ABUNDANCE_IN_DATASET).sum() >= dataset_threshold
taxa_to_remove_abundance = dataset_avg_abundances.columns[~taxa_meeting_abundance_threshold].tolist()

# Print number of taxa to remove from abundance filter
print(f"Number of taxa to remove based on abundance threshold: {len(taxa_to_remove_abundance)}")

# Calculate the prevalence for each taxon in each dataset
dataset_prevalences = tss_taxonomic_abundance_matrix.reset_index().merge(
    samples_df[[SAMPLE_ID_COLUMN, DATASET_ID_COLUMN]], 
    on=SAMPLE_ID_COLUMN
).drop(SAMPLE_ID_COLUMN, axis=1).groupby(DATASET_ID_COLUMN).apply(
    lambda x: (x > 0).sum() / len(x),
    include_groups=False
)

# Get taxa to remove based on prevalence threshold
taxa_meeting_prevalence_threshold = (dataset_prevalences >= MIN_PREVALENCE_IN_SAMPLES_IN_DATASET).sum() >= dataset_threshold
taxa_to_remove_prevalence = dataset_prevalences.columns[~taxa_meeting_prevalence_threshold].tolist()

# Print number of taxa to remove from prevalence filter
print(f"Number of taxa to remove based on prevalence threshold: {len(taxa_to_remove_prevalence)}")

# Get all taxa to remove (excluding duplicates)
all_taxa_to_remove = list(set(taxa_to_remove_abundance + taxa_to_remove_prevalence))
if 'Dataset_ID' in all_taxa_to_remove:
    all_taxa_to_remove.remove('Dataset_ID')
print(f"Number of taxa to remove based on abundance OR prevalence threshold: {len(all_taxa_to_remove)}")
print(f"...Leaving {taxonomic_abundance_matrix.shape[1] - len(all_taxa_to_remove)} taxa in the dataset.\n")



# Filter taxa before transformation -----------------------------------------
# This removes taxa, the data is kept as untransformed read counts

if WHEN_FILTERING == 'before': 

    # Remove taxa
    taxonomic_abundance_matrix = taxonomic_abundance_matrix.drop(columns=all_taxa_to_remove)



# Normalisation and transformation ------------------------------------------
# These are applied to raw read counts

# Transformation
if TRANFORMATION == 'tss':
    taxonomic_abundance_matrix = apply_total_sum_scaling(taxonomic_abundance_matrix)

elif TRANFORMATION == 'log10':
    taxonomic_abundance_matrix = transform_log(taxonomic_abundance_matrix)

elif TRANFORMATION == 'ilr':
    taxonomic_abundance_matrix = transform_ilr(taxonomic_abundance_matrix)

elif TRANFORMATION == 'clr':
    taxonomic_abundance_matrix = transform_clr(taxonomic_abundance_matrix)


# Z-score bacteria
if APPLY_ZSCORE:
    taxonomic_abundance_matrix = apply_zscore(taxonomic_abundance_matrix)





# Filter taxa after transformation -----------------------------------------
# This removes taxa, the data is kept as it is

if WHEN_FILTERING == 'after': 

    # Remove taxa
    taxonomic_abundance_matrix = taxonomic_abundance_matrix.drop(columns=all_taxa_to_remove)






# Add extra features to the abundance matrix dataframe ----------------------------------
if INCLUDE_COLUMNS_AS_FEATURES:
    # Encode each feature column as an integer starting at 0
    for feature in INCLUDE_COLUMNS_AS_FEATURES:
        samples_df[feature] = samples_df[feature].astype('category').cat.codes

    # Merge the encoded features into the taxonomic abundance matrix
    features_df = samples_df.set_index(SAMPLE_ID_COLUMN)[INCLUDE_COLUMNS_AS_FEATURES]
    taxonomic_abundance_matrix = features_df.join(taxonomic_abundance_matrix)
    print(f"Added extra features ({INCLUDE_COLUMNS_AS_FEATURES}) to the taxonomic abundance matrix")




# Summary of labels ------------------------------------------

# Print the number of samples with each label overall
print()
print("\nNumber of samples with each label overall:")
overall_label_counts = samples_df[LABEL_COLUMNS[0]].value_counts()
for label, count in overall_label_counts.items():
    print(f"  {label}: {count}")

# Print the number of samples with each label in each dataset
print("\nNumber of samples with each label in each dataset:")
for dataset_id, group in samples_df.groupby(DATASET_ID_COLUMN):
    print(f"Dataset {dataset_id}:")
    label_counts = group[LABEL_COLUMNS[0]].value_counts()
    for label, count in label_counts.items():
        print(f"  {label}: {count}")

print()

# Write dataset to file ------------------------------------------
# Features
taxonomic_abundance_matrix.to_csv(TAXONOMIC_ABUNDANCE_MATRIX_TSV_PATH, sep="\t")
# Labels
samples_df.to_csv(FILTERED_SAMPLES_TSV_PATH, sep="\t", index=False)

print("Done!")

# Write about.txt
with open(ABOUT_PATH, "w") as about_file:
    about_file.write("Dataset Construction Summary\n")
    about_file.write("Input Files:\n")
    
    about_file.write("Output Files:\n")
    
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
    for label, count in overall_label_counts.items():
        about_file.write(f"  - {label}: {count}\n")
    about_file.write("\nLabel distribution by dataset:\n")
    for dataset_id, group in samples_df.groupby(DATASET_ID_COLUMN):
        about_file.write(f"  Dataset {dataset_id}:\n")
        label_counts = group[LABEL_COLUMNS[0]].value_counts()
        for label, count in label_counts.items():
            about_file.write(f"    - {label}: {count}\n")

    # Filters:
    about_file.write(f"\nFiltering taxa by average abundance:\n")
    about_file.write(f"  - Minimum average abundance: {MIN_AVERAGE_ABUNDANCE_IN_DATASET * 100}%\n")
    about_file.write(f"  - Required in at least {IN_N_PROPORTION_DATASETS * 100}% of datasets ({dataset_threshold} datasets)\n")
    about_file.write(f"\nFiltering taxa by prevalence:\n")
    about_file.write(f"  - Minimum prevalence: {MIN_PREVALENCE_IN_SAMPLES_IN_DATASET * 100}%\n")
    about_file.write(f"  - Required in at least {IN_N_PROPORTION_DATASETS * 100}% of datasets ({dataset_threshold} datasets)\n")


    # Add timestamp
    about_file.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    # Add constants for clarity
    about_file.write("Constants Used:\n")
    about_file.write(f"  - Pseudo count: {PSEUDO_COUNT}\n")
    about_file.write(f"  - Minimum average abundance across all samples for a taxonomic unit to be included: {MIN_AVERAGE_ABUNDANCE_IN_DATASET}\n")
    about_file.write(f"  - Minimum proportion of all samples for a taxonomic unit to be included: {MIN_PREVALENCE_IN_SAMPLES_IN_DATASET}\n")
    about_file.write(f"  - Transformation Applied: {TRANFORMATION}\n\n")
    about_file.write("Filters Applied:\n")
    about_file.write("  - Exclude rows with specific values:\n")
    for column, values in EXCLUDE_ROWS_WITH_VALUES.items():
        about_file.write(f"    - {column}: {values}\n")
    about_file.write("  - Include rows only with specific values:\n")
    for column, values in ONLY_INCLUDE_ROWS_WITH_VALUES.items():
        about_file.write(f"    - {column}: {values}\n")
    about_file.write(f"  - Exclude datasets: {EXCLUDE_DATASETS}\n")
    about_file.write(f"  - Exclude samples: {EXCLUDE_SAMPLES}\n\n")

    # Transformation-specific details
    if TRANFORMATION in ['clr', 'ilr', 'clr-z', 'ilr-z']:
        about_file.write("Transformation Details:\n")
        about_file.write(f"  - Applied centered log-ratio (CLR): {'clr' in TRANFORMATION}\n")
        about_file.write(f"  - Applied isometric log-ratio (ILR): {'ilr' in TRANFORMATION}\n")
        about_file.write(f"  - Applied Z-score normalization: {'-z' in TRANFORMATION}\n\n")

    # Label-specific information
    about_file.write("Sample Label Summary:\n")
    overall_label_counts = samples_df[LABEL_COLUMNS[0]].value_counts()
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
    "Total Excluded Samples": [num_inital_samples - num_final_samples]
}
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file_path, sep="\t", index=False)

print(f"Summary file written to: {summary_file_path}")




print("Done!")
