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
TSV_FILE = 'machine_learning/datasets/stool_treated_log10_after/kfold_results/best_models/roc_mlp_data.tsv'
TITLE = 'Stool Treated MLP K-fold ROC Curve'
OUT_FILE = TITLE.lower().replace(' ', '_').replace('-', '_') + '.png'

OUTPUT_DIR = '~/Repos/meta-analysis/machine_learning/results/roc_curves/'

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
