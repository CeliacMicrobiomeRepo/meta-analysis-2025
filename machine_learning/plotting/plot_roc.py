"""
Plot the ROC curve for a given model and dataset.

- reads the ROC curve data from a TSV file
- plots a ROC curve
- writes a png 
"""

# Imports ------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os

# Options ------------------------------------------
ANALYSIS_GROUP = "duodenum_active"
EVALUATION_METHOD = "kfold"
BEST_MODEL = True
MODEL_NAME = "rf" # <- used only if BEST_MODEL is False

if BEST_MODEL:
    # Read summary_best_model.tsv to get the model name
    summary_best_model_path = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/summary_best_model.tsv'
    summary_best_model_df = pd.read_csv(summary_best_model_path, sep='\t')
    # Get 2nd column name
    MODEL_NAME = summary_best_model_df.columns[1].lower()

TSV_FILE = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/best_models/roc_' + MODEL_NAME + '_data.tsv'
TITLE = ANALYSIS_GROUP + ' ' + MODEL_NAME.upper() + ' ' + EVALUATION_METHOD.upper() + ' ROC Curve'
OUT_FILE = TITLE.lower().replace(' ', '_').replace('-', '_') + '.png'

OUTPUT_DIR = '/home/haig/Repos/meta-analysis/machine_learning/results/roc_plots/'

# Main ------------------------------------------
# Read the data
data = pd.read_csv(TSV_FILE, sep='	')
fpr = data['fpr']
tpr = data['tpr']

# Calculate AUC
roc_auc = auc(fpr, tpr)

# Create plot
plt.figure(figsize=(6, 6))

# Plot the ROC curve
plt.plot(fpr, tpr, color='darkorange', lw=2)

# Plot the diagonal line for reference
plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')

# Set labels and title
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(TITLE)

# Set axis limits
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])

# Make it square
plt.gca().set_aspect('equal', adjustable='box')

# Display AUC value bottom right
plt.text(0.95, 0.05, f'AUC = {roc_auc:.3f}',
         verticalalignment='bottom', horizontalalignment='right',
         transform=plt.gca().transAxes,
         fontsize=12)

# Save the figure
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = os.path.join(OUTPUT_DIR, OUT_FILE)
plt.savefig(output_path, bbox_inches='tight')

print(f"Plot saved to {output_path}")
