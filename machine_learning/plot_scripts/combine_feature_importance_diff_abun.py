"""
This script combines the feature importance and differential abundance files
- creates a single file 
- creates plots
"""

import os
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats


MODELS = ['rf', 'svm', 'lr', 'xgboost', 'mlp']
PLOT_MODEL = 'rf'           # rf or xgboost
BODY_SITE = "stool"     # stool or duodenal

OUTPUT_DIR = "./feature_importance_diff_abun/"
OUTPUT_COMBINED_FILE = OUTPUT_DIR + BODY_SITE + "_combined_feature_importance.tsv"

# The diff abun files
DIFF_ABUN_ANCOMBC2_FILE = "./feature_importance_diff_abun/daa/" + BODY_SITE + "/ancombc2_results_all.csv"
DIFF_ABUN_DESEQ2_FILE = "./feature_importance_diff_abun/daa/" + BODY_SITE + "/deseq2_results_all.csv"
DIFF_ABUN_LEFSE_FILE = "./feature_importance_diff_abun/daa/" + BODY_SITE + "/lefse_results_all.csv"
# Significant ASVs (These the the actual ones we are using for DAA results)
SIGNIFICANT_ASVS_FILE = "./feature_importance_diff_abun/daa/" + BODY_SITE + "/significant_asvs.tsv"

# The feature importance file always has columns:       ASV	mlp	xgboost	svm	lr	rf	Taxonomy
FEATURE_IMPORTANCE_FILE = "./datasets/" + BODY_SITE + "_noGFD_ASV_after_log10/kfold_results/feature_importance_best_models.tsv"


PADJ_THRESHOLD = None


# Read feature importance file
feature_importance_df = pd.read_csv(FEATURE_IMPORTANCE_FILE, sep='\t')

# Read diff abun files
diff_abun_ancombc2_df = pd.read_csv(DIFF_ABUN_ANCOMBC2_FILE)
diff_abun_deseq2_df = pd.read_csv(DIFF_ABUN_DESEQ2_FILE)
diff_abun_lefse_df = pd.read_csv(DIFF_ABUN_LEFSE_FILE)


# How many ASVs are in each of the 4 files?
feature_importance_asvs = set(feature_importance_df['ASV'])
diff_abun_ancombc2_asvs = set(diff_abun_ancombc2_df['ASV'])
diff_abun_deseq2_asvs = set(diff_abun_deseq2_df['ASV'])
diff_abun_lefse_asvs = set(diff_abun_lefse_df['ASV'])
print(f"Number of ASVs in feature importance file: {len(feature_importance_asvs)}")
print(f"Number of ASVs in ancombc2 file: {len(diff_abun_ancombc2_asvs)}")
print(f"Number of ASVs in deseq2 file: {len(diff_abun_deseq2_asvs)}")
print(f"Number of ASVs in lefse file: {len(diff_abun_lefse_asvs)}")


# Which ASVs are in each of the diff abun files but not in the feature importance file?
diff_abun_ancombc2_asvs_not_in_feature_importance = diff_abun_ancombc2_asvs - feature_importance_asvs
diff_abun_deseq2_asvs_not_in_feature_importance = diff_abun_deseq2_asvs - feature_importance_asvs
diff_abun_lefse_asvs_not_in_feature_importance = diff_abun_lefse_asvs - feature_importance_asvs
print(f"Number of ASVs in ancombc2 file but not in feature importance file: {len(diff_abun_ancombc2_asvs_not_in_feature_importance)}")
print(f"Number of ASVs in deseq2 file but not in feature importance file: {len(diff_abun_deseq2_asvs_not_in_feature_importance)}")
print(f"Number of ASVs in lefse file but not in feature importance file: {len(diff_abun_lefse_asvs_not_in_feature_importance)}")

# Remove rows in the diff abun df where the ASV is not in the feature_importance_asvs
diff_abun_ancombc2_df = diff_abun_ancombc2_df[diff_abun_ancombc2_df['ASV'].isin(feature_importance_asvs)]
diff_abun_deseq2_df = diff_abun_deseq2_df[diff_abun_deseq2_df['ASV'].isin(feature_importance_asvs)]
diff_abun_lefse_df = diff_abun_lefse_df[diff_abun_lefse_df['ASV'].isin(feature_importance_asvs)]
print(f"Removed {len(diff_abun_ancombc2_asvs_not_in_feature_importance)} ASVs not in feature importance file")

# Which ASVs are in the feature importance file, but not in any of the diff abun files?
feature_importance_asvs_not_in_diff_abun = feature_importance_asvs - diff_abun_ancombc2_asvs - diff_abun_deseq2_asvs - diff_abun_lefse_asvs
print(f"Number of ASVs in feature importance file but not in any diff abun file: {len(feature_importance_asvs_not_in_diff_abun)}")

# Remove rows in the feature importance df where the ASV is not in any of the diff abun files
feature_importance_df = feature_importance_df[feature_importance_df['ASV'].isin(diff_abun_ancombc2_asvs | diff_abun_deseq2_asvs | diff_abun_lefse_asvs)]
print(f"Removed {len(feature_importance_asvs_not_in_diff_abun)} ASVs not in diff abun files")

# How many ASVs are in each of the 4 files? and how many overall?
feature_importance_asvs = set(feature_importance_df['ASV'])
diff_abun_ancombc2_asvs = set(diff_abun_ancombc2_df['ASV'])
diff_abun_deseq2_asvs = set(diff_abun_deseq2_df['ASV'])
diff_abun_lefse_asvs = set(diff_abun_lefse_df['ASV'])
all_asvs = feature_importance_asvs | diff_abun_ancombc2_asvs | diff_abun_deseq2_asvs | diff_abun_lefse_asvs
print()
print(f"Number of ASVs in feature importance file (after filtering ASVs): {len(feature_importance_asvs)}")
print(f"Number of ASVs in ancombc2 file (after filtering ASVs): {len(diff_abun_ancombc2_asvs)}")
print(f"Number of ASVs in deseq2 file (after filtering ASVs): {len(diff_abun_deseq2_asvs)}")
print(f"Number of ASVs in lefse file (after filtering ASVs): {len(diff_abun_lefse_asvs)}")
print(f"Number of ASVs in all files (after filtering ASVs): {len(all_asvs)}")

# Get ASVs missing from the the diff abun files
missing_asvs_ancombc2 = all_asvs - diff_abun_ancombc2_asvs
missing_asvs_deseq2 = all_asvs - diff_abun_deseq2_asvs
missing_asvs_lefse = all_asvs - diff_abun_lefse_asvs

print(len(missing_asvs_ancombc2))
print(len(missing_asvs_deseq2))
print(len(missing_asvs_lefse))

# Filter the ancombc2 df to only include the columns "ASV" "Log2FC" "padj" "pval" 
# + rename "Log2FC" to "ANCOMBC2_Log2FC" and "padj" to "ANCOMBC2_padj" and "pval" to "ANCOMBC2_pval"
diff_abun_ancombc2_df = diff_abun_ancombc2_df[['ASV', 'Log2FC', 'padj', 'pval']]
diff_abun_ancombc2_df.rename(columns={'Log2FC': 'ANCOMBC2_Log2FC', 'padj': 'ANCOMBC2_padj', 'pval': 'ANCOMBC2_pval'}, inplace=True)
# Filter the deseq2 df to only include the columns "ASV" "log2FoldChange" "padj" "pvalue" 
# + rename "log2FoldChange" to "DESeq2_Log2FC" and "padj" to "DESeq2_padj" and "pvalue" to "DESeq2_pval"
diff_abun_deseq2_df = diff_abun_deseq2_df[['ASV', 'log2FoldChange', 'padj', 'pvalue']]
diff_abun_deseq2_df.rename(columns={'log2FoldChange': 'DESeq2_Log2FC', 'padj': 'DESeq2_padj', 'pvalue': 'DESeq2_pval'}, inplace=True)
# Filter the lefse df to only include the columns "ASV" "scores" 
# + rename "scores" to "LEfSe_LDA_score"
diff_abun_lefse_df = diff_abun_lefse_df[['ASV', 'scores']]
diff_abun_lefse_df.rename(columns={'scores': 'LEfSe_LDA_score'}, inplace=True)

# Add the missing ASVs as rows to the diff abun dfs (with all values 0)
if len(missing_asvs_ancombc2) > 0:
    cols = ['ASV', 'ANCOMBC2_Log2FC', 'ANCOMBC2_padj', 'ANCOMBC2_pval']
    data = [[asv, 99999, 99999, 99999] for asv in missing_asvs_ancombc2]
    diff_abun_ancombc2_df = pd.concat([diff_abun_ancombc2_df, pd.DataFrame(columns=cols, data=data)])
if len(missing_asvs_deseq2) > 0:
    cols = ['ASV', 'DESeq2_Log2FC', 'DESeq2_padj', 'DESeq2_pval']
    data = [[asv, 99999, 99999, 99999] for asv in missing_asvs_deseq2]
    diff_abun_deseq2_df = pd.concat([diff_abun_deseq2_df, pd.DataFrame(columns=cols, data=data)])
if len(missing_asvs_lefse) > 0:
    cols = ['ASV', 'LEfSe_LDA_score']
    data = [[asv, 99999] for asv in missing_asvs_lefse]
    diff_abun_lefse_df = pd.concat([diff_abun_lefse_df, pd.DataFrame(columns=cols, data=data)])
# Fill empty values with 0
diff_abun_ancombc2_df = diff_abun_ancombc2_df.fillna(0)
diff_abun_deseq2_df = diff_abun_deseq2_df.fillna(0)
diff_abun_lefse_df = diff_abun_lefse_df.fillna(0)

# Combine the 4 dfs on the ASV column (but still keep rows where the ASV is not in all dfs)
combined_df = pd.merge(feature_importance_df, diff_abun_ancombc2_df, on='ASV', how='outer')
combined_df = pd.merge(combined_df, diff_abun_deseq2_df, on='ASV', how='outer')
combined_df = pd.merge(combined_df, diff_abun_lefse_df, on='ASV', how='outer')

# Read significant ASVs file
significant_asvs_df = pd.read_csv(SIGNIFICANT_ASVS_FILE, sep='\t')
# Make a dictionary which maps the ASV to the "Num Sig" column (possible values for num sig are 1, 2 or 3)
significant_asvs_dict = dict(zip(significant_asvs_df['ASV'], significant_asvs_df['Num Sig']))

# Write the combined df to a file
combined_df.to_csv(OUTPUT_COMBINED_FILE, sep='\t', index=False)
print(f"Combined file saved to {OUTPUT_COMBINED_FILE}")

# Get absolute values of the columns 'ANCOMBC2_Log2FC', 'DESeq2_Log2FC' and 'LEfSe_LDA_score'
combined_df['ANCOMBC2_Log2FC'] = combined_df['ANCOMBC2_Log2FC'].abs()
combined_df['DESeq2_Log2FC'] = combined_df['DESeq2_Log2FC'].abs()
combined_df['LEfSe_LDA_score'] = combined_df['LEfSe_LDA_score'].abs()

# Make a scatter plot of the feature importance vs the diff abun Log2FC and LDAs
diff_abun_columns = [('ANCOMBC2_Log2FC', 'ANCOMBC2_padj'), ('DESeq2_Log2FC', 'DESeq2_padj'), ('LEfSe_LDA_score', None)]
for diff_abun_column, colour_column in diff_abun_columns:
    plot_combined_df = combined_df.copy()
    # Do not show ASVs that have values of 99999
    plot_combined_df = plot_combined_df[plot_combined_df[diff_abun_column] != 99999]

    # If padj in colour_column, do not show ASVs that have padj > PADJ_THRESHOLD
    if colour_column is not None and PADJ_THRESHOLD is not None and 'padj' in colour_column:
        plot_combined_df = plot_combined_df[plot_combined_df[colour_column] <= PADJ_THRESHOLD]

    if colour_column is None:
        plt.scatter(plot_combined_df[diff_abun_column], plot_combined_df[PLOT_MODEL])
    else:
        # Use colour
        scatter = plt.scatter(plot_combined_df[diff_abun_column], plot_combined_df[PLOT_MODEL], c=plot_combined_df[colour_column])
        # Add a legend for the colour
        plt.colorbar(scatter, label=colour_column)
    
    # Run a linear regression on the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(plot_combined_df[diff_abun_column], plot_combined_df[PLOT_MODEL])
    print(f"Slope: {slope}, Intercept: {intercept}, R-value: {r_value}, P-value: {p_value}, Std Error: {std_err}")

    # Add the regression line to the plot
    plt.plot(plot_combined_df[diff_abun_column], slope * plot_combined_df[diff_abun_column] + intercept, color='red')
    # Add the equation to the plot
    plt.text(0.05, 0.95, f'y = {slope:.2f}x + {intercept:.2f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Add the R-value to the plot
    plt.text(0.05, 0.90, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')


    plt.xlabel(diff_abun_column)
    plt.ylabel(f'Feature Importance ({PLOT_MODEL.upper()})')
    plt.ylim(0, max(plot_combined_df[PLOT_MODEL]) * 1.1)
    plt.show()




# Make a scatter plot of the feature importance vs the diff abun padj values
diff_abun_columns = [('ANCOMBC2_padj', None), ('DESeq2_padj', None)]
for x_column, colour_column in diff_abun_columns:
    plot_combined_df = combined_df.copy()
    # Do not show ASVs that have values of 99999
    plot_combined_df = plot_combined_df[plot_combined_df[x_column] != 99999]

    # If padj in x_column, do not show ASVs that have padj > PADJ_THRESHOLD
    if PADJ_THRESHOLD is not None and 'padj' in x_column:
        plot_combined_df = plot_combined_df[plot_combined_df[x_column] <= PADJ_THRESHOLD]

    if colour_column is None:
        plt.scatter(plot_combined_df[x_column], plot_combined_df[PLOT_MODEL])
    else:
        # Use colour
        scatter = plt.scatter(plot_combined_df[x_column], plot_combined_df[PLOT_MODEL], c=plot_combined_df[colour_column])
        # Add a legend for the colour
        plt.colorbar(scatter, label=colour_column)
    
    # Run a linear regression on the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(plot_combined_df[x_column], plot_combined_df[PLOT_MODEL])

    # Add the regression line to the plot
    plt.plot(plot_combined_df[x_column], slope * plot_combined_df[x_column] + intercept, color='red')
    # Add the equation to the plot
    plt.text(0.05, 0.95, f'y = {slope:.2f}x + {intercept:.2f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    # Add the R-value to the plot
    plt.text(0.05, 0.90, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')


    plt.xlabel(x_column)
    plt.ylabel(f'Feature Importance ({PLOT_MODEL.upper()})')
    plt.ylim(0, max(plot_combined_df[PLOT_MODEL]) * 1.1)
    plt.show()




# Make a boxplot of the feature importance of "All ASVs", "Sig ASVs (>0)", "Sig ASVs (>1)", "Sig ASVs (>2)"
# Get the ASVs that have a num sig of 1 or more, 2 or more, 3 or more
significant_asvs_1_or_more = [asv for asv, num_sig in significant_asvs_dict.items() if num_sig >= 1]
significant_asvs_2_or_more = [asv for asv, num_sig in significant_asvs_dict.items() if num_sig >= 2]
significant_asvs_3_or_more = [asv for asv, num_sig in significant_asvs_dict.items() if num_sig >= 3]
print(len(significant_asvs_1_or_more))
print(len(significant_asvs_2_or_more))
print(len(significant_asvs_3_or_more))

# Subset feature importance values for each group
all_asvs_importance = combined_df[PLOT_MODEL]
sig_asvs_1_or_more_importance = combined_df[combined_df['ASV'].isin(significant_asvs_1_or_more)][PLOT_MODEL]
sig_asvs_2_or_more_importance = combined_df[combined_df['ASV'].isin(significant_asvs_2_or_more)][PLOT_MODEL]
sig_asvs_3_or_more_importance = combined_df[combined_df['ASV'].isin(significant_asvs_3_or_more)][PLOT_MODEL]

# Make labels for each group
labels = [f'All ASVs (n={len(all_asvs_importance)})', f'Sig ASVs (1+) (n={len(sig_asvs_1_or_more_importance)})', f'Sig ASVs (2+) (n={len(sig_asvs_2_or_more_importance)})', f'Sig ASVs (3+) (n={len(sig_asvs_3_or_more_importance)})']

# Make title
title = f'Feature Importance Boxplot ({PLOT_MODEL.upper()})'

# Make plot
plt.boxplot(
    [all_asvs_importance, sig_asvs_1_or_more_importance, sig_asvs_2_or_more_importance, sig_asvs_3_or_more_importance],
    labels=labels
)
plt.title(title)
plt.ylabel('Feature Importance')
plt.show()


# Make a scatter plot of the feature importance vs the diff abun Log2FC and LDAs (SIGNIFICANT ASVS COLOURED)
diff_abun_columns = [('ANCOMBC2_Log2FC', None), ('DESeq2_Log2FC', None), ('LEfSe_LDA_score', None)]
for diff_abun_column, colour_column in diff_abun_columns:
    plot_combined_df = combined_df.copy()
    # Do not show ASVs that have values of 99999
    plot_combined_df = plot_combined_df[plot_combined_df[diff_abun_column] != 99999]

    # If padj in colour_column, do not show ASVs that have padj > PADJ_THRESHOLD
    if colour_column is not None and PADJ_THRESHOLD is not None and 'padj' in colour_column:
        plot_combined_df = plot_combined_df[plot_combined_df[colour_column] <= PADJ_THRESHOLD]

    # Add a colour column based on the ASV's presence in sig_asvs_1_or_more_importance, sig_asvs_2_or_more_importance, sig_asvs_3_or_more_importance
    plot_combined_df['colour'] = ['grey'] * len(plot_combined_df['ASV'])
    # Set the colour to yellow where the ASV is in significant_asvs_1_or_more
    plot_combined_df.loc[plot_combined_df['ASV'].isin(significant_asvs_1_or_more), 'colour'] = 'yellow'
    # Set the colour to orange where the ASV is in significant_asvs_2_or_more
    plot_combined_df.loc[plot_combined_df['ASV'].isin(significant_asvs_2_or_more), 'colour'] = 'orange'
    # Set the colour to red where the ASV is in significant_asvs_3_or_more
    plot_combined_df.loc[plot_combined_df['ASV'].isin(significant_asvs_3_or_more), 'colour'] = 'red'
    
    # Plot the data with colours based on if ASV
    plt.scatter(plot_combined_df[diff_abun_column], plot_combined_df[PLOT_MODEL], c=plot_combined_df['colour'])

    # Select colours for labels (for legend)
    labels = [f'All ASVs (n={len(all_asvs_importance)})', f'Sig ASVs (1+) (n={len(sig_asvs_1_or_more_importance)})', f'Sig ASVs (2+) (n={len(sig_asvs_2_or_more_importance)})', f'Sig ASVs (3+) (n={len(sig_asvs_3_or_more_importance)})']
    colours = ['grey', 'yellow', 'orange', 'red']

    # Add a legend for the colours of labels
    plt.legend(labels, colours)

    plt.xlabel(diff_abun_column)
    plt.ylabel(f'Feature Importance ({PLOT_MODEL.upper()})')
    plt.ylim(0, max(plot_combined_df[PLOT_MODEL]) * 1.1)
    plt.show()
