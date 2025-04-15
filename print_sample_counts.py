"""
You can run this script to get a look at the samples used in analysis.

For analysis, we focused on only:
 - 16S sequencing data
 - Stool or Duodenal samples
 - Samples from Non-prospective Studies

For specific analyses we did further filtering not seen here.
"""

# Imports
import os
import urllib.request
import pandas as pd

# Download the data file if it does not exist
DATA_URL = 'https://raw.githubusercontent.com/CeliacMicrobiomeRepo/celiac-repository/main/all_samples.tsv'
DATA_FILE = 'all_samples.tsv'

if not os.path.exists(DATA_FILE):
    print(f"Downloading {DATA_FILE} from {DATA_URL}...")
    urllib.request.urlretrieve(DATA_URL, DATA_FILE)
    print(f"Downloaded {DATA_FILE}.")
else:
    print(f"{DATA_FILE} already exists. Using local copy.")

# ======================================================================================
# STOOL SAMPLES
# ======================================================================================
print("\n------------------------------------------------------------------------------ \nSTOOL SAMPLES \n------------------------------------------------------------------------------ \n")
# Load data
df = pd.read_csv('all_samples.tsv', sep='\t', dtype=str, keep_default_na=False, na_values=[])
print("Original:", len(df))

# Apply filters
f1 = df[df['Sequencing_Type'] != 'SG']
print("After Sequencing_Type != 'SG':", len(f1))
f2 = f1[f1['Group_Prospective_Study'] == 'NA']
print("After Group_Prospective_Study == 'NA':", len(f2))
f3 = f2[f2['Sample_Site'].isin(['stool'])]
print("After Sample_Site in ['stool']:", len(f3))

# Print 
print("\n----------\nTotal samples:", len(f3))
print('\nCounts by Group:')
print(f3['Group'].value_counts(dropna=False))
print('\nCounts by Dataset_ID:')
print(f3['Dataset_ID'].value_counts(dropna=False))
print("\n----------\nTotal Datasets:", len(f3['Dataset_ID'].unique()))
print("\nCounts by Dataset_ID and Diagnosed_Celiac:")
stool_ct = pd.crosstab(f3['Dataset_ID'], f3['Diagnosed_Celiac'])
stool_ct.loc['Total'] = stool_ct.sum()
print(stool_ct)
# ======================================================================================




# ======================================================================================
# DUODENAL SAMPLES
# ======================================================================================
print("\n\n\n------------------------------------------------------------------------------ \nDUODENAL SAMPLES \n------------------------------------------------------------------------------ \n")
# Load data
df = pd.read_csv('all_samples.tsv', sep='\t', dtype=str, keep_default_na=False, na_values=[])
print("Original:", len(df))

# Apply filters
f1 = df[df['Sequencing_Type'] != 'SG']
print("After Sequencing_Type != 'SG':", len(f1))
f2 = f1[f1['Group_Prospective_Study'] == 'NA']
print("After Group_Prospective_Study == 'NA':", len(f2))
f3 = f2[f2['Sample_Site'].isin(['duodenal'])]
print("After Sample_Site in ['duodenal']:", len(f3))

# Print 
print("\n----------\nTotal samples:", len(f3))
print('\nCounts by Group:')
print(f3['Group'].value_counts(dropna=False))
print('\nCounts by Dataset_ID:')
print(f3['Dataset_ID'].value_counts(dropna=False))
print("\n----------\nTotal Datasets:", len(f3['Dataset_ID'].unique()))
print("\nCounts by Dataset_ID and Diagnosed_Celiac:")
duod_ct = pd.crosstab(f3['Dataset_ID'], f3['Diagnosed_Celiac'])
duod_ct.loc['Total'] = duod_ct.sum()
print(duod_ct)
# ======================================================================================





# ======================================================================================
# STOOL AND DUODENAL SAMPLES
# ======================================================================================
print("\n\n\n------------------------------------------------------------------------------ \nSTOOL AND DUODENAL SAMPLES \n------------------------------------------------------------------------------ \n")
# Load data
df = pd.read_csv('all_samples.tsv', sep='\t', dtype=str, keep_default_na=False, na_values=[])
print("Original:", len(df))

# Apply filters
f1 = df[df['Sequencing_Type'] != 'SG']
print("After Sequencing_Type != 'SG':", len(f1))
f2 = f1[f1['Group_Prospective_Study'] == 'NA']
print("After Group_Prospective_Study == 'NA':", len(f2))
f3 = f2[f2['Sample_Site'].isin(['duodenal', 'stool'])]
print("After Sample_Site in ['duodenal', 'stool']:", len(f3))

# Print 
print("\n----------\nTotal samples:", len(f3))
print('\nCounts by Group:')
print(f3['Group'].value_counts(dropna=False))
print('\nCounts by Dataset_ID:')
print(f3['Dataset_ID'].value_counts(dropna=False))
print("\n----------\nTotal Datasets:", len(f3['Dataset_ID'].unique()))
print("\nCounts by Dataset_ID and Diagnosed_Celiac:")
combo_ct = pd.crosstab(f3['Dataset_ID'], f3['Diagnosed_Celiac'])
combo_ct.loc['Total'] = combo_ct.sum()
print(combo_ct)
# ======================================================================================

