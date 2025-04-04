# This script takes the results of ASV truncation across multiple datasets
# It combines them into one phyloseq object and filters according to presence and prevalence

# Outputs:
# ps0 - Phyloseq object with no filtering of ASVs and untransformed or normalised read counts
# ps1 - Phyloseq object after filtering of ASVs, but before transformation or normalisation of read counts
# ps2 - Phyloseq object after filtering of ASVs and transformation/normalisation of read counts



# Set up ----------------------------------------

# Import
library(dada2)
library(phyloseq)
library(ggplot2)


BASE_DIR <- "C:/Users/haig/UC Enterprise Dropbox/Haig Bishop/celiac_microbiome_data/16S_pipeline/"

# Paths to databases
TAXONOMY_TRAIN_SET <- "D:/16S_databases/silva_nr99_v138.1_train_set.fa.gz"
SPECIES_ASSIGNMENT_SET <- "D:/16S_databases/silva_species_assignment_v138.1.fa.gz"

# Path to samples.tsv file (rows are samples and columns are metadata)
ALL_SAMPLES_TSV <- file.path(BASE_DIR, "all_samples.tsv")

# Directory in which all datasets are contained
DATASETS_DIR <- "D:/microbiome_sequencing_datasets/celiac_16s_datasets"
# Directories which are in the DATASETS_DIR
DATASET_NAMES <- c(
  "16S_136_Nobel",           # Stool
  "16S_20_Rawson",           # Stool
  "16S_60_Shi",              # Stool
  "16S_96_Quagliariello",    # Stool
  "16S_27_Fornasaro",        # Stool
  "16S_49_Turjeman",         # Stool
  "16S_102_Bodkhe",          # Stool and Duodenum
  "16S_80_Garcia",           # Duodenum (ignore stool samples from this)
  "16S_119_Salamon"          # Duodenum
)
# Path to sequence table within in dataset directory ^
TRUNC_SUBDIR_NAME <- "core_v4_truncation/seqs_table.rds"

# Output paths
PS0_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps0.rds")
PS1_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps1.rds")
PS2_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps2.rds")
PS0_taxa_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps0_taxa.csv")
PS1_taxa_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps1_taxa.csv")
PS2_taxa_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps2_taxa.csv")
# Make the directory if it doesnt exist (for PS0_output and PS0_taxa_output)
dir.create(dirname(PS0_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(PS0_taxa_output), recursive = TRUE, showWarnings = FALSE)
# How must of the samples are made up by these core ASVs?
TOTAL_ABUNDANCE_DATA_PATH <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_core_asv_sample_abundance.csv")
TOTAL_ABUNDANCE_PLOT_PATH <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_core_asv_abundance_boxplot.png")
# TSV output paths
PS2_otu_tsv_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps2_asv_table.tsv")
PS2_sample_tsv_output <- file.path(DATASETS_DIR, "core_v4_truncation", "combined_phyloseq_objects", "stool_ps2_sample_data.tsv")


# Sample filtering options --------
# Column containing sample IDs
SAMPLE_ID_COLUMN = "Sample_ID"
# Column containing dataset IDs
DATASET_ID_COLUMN = "Dataset_ID"
# Columns to use as labels
LABEL_COLUMN = "Diagnosed_Celiac"
# Exclude rows if
EXCLUDE_ROWS_WITH_VALUES = list(
  "Dataset_ID" = c("-", "", NA),
  "Gluten_Free_Diet" = c()
)
# Only include rows if
ONLY_INCLUDE_ROWS_WITH_VALUES = list(
  "Diagnosed_Celiac" = c(TRUE, FALSE),
  "Any_Significant_Factor" = c(FALSE),
  "Sample_Site" = c("stool")
)
# Exclude specific datasets
EXCLUDE_DATASETS = c("16S_80_Garcia")
# Keep only specific columns
KEEP_ONLY_COLUMNS = c(DATASET_ID_COLUMN, LABEL_COLUMN, "Sample_Site", "Any_Significant_Factor", "Gluten_Free_Diet")
# Exclude samples
EXCLUDE_SAMPLES <- c("SRR1107516", "ERR1551255", "ERR1551306", "SRR18231165", "SRR6885558")


# ASV filtering options --------------
# Minimum average abundance across a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_AVERAGE_ABUNDANCE_IN_DATASET <- 0.001
# Minimum proportion of all samples in a dataset (in X% of datasets) for a taxonomic unit to be included
MIN_PREVALENCE_IN_SAMPLES_IN_DATASET <- 0.1
# ...X proportion of the datasets:
IN_N_PROPORTION_DATASETS <- 0.3

# A final filter to Samples ------------------
# Must have at least this proportion of their abundance made of the core ASVs
MIN_TOTAL_ABUNDANCE_CORE_ASVS <- 0.25


# Load sequence table files ----------------------------------------

# Get a list of all the sequence table files (for each dataset)
seqtab_files <- file.path(DATASETS_DIR, DATASET_NAMES, TRUNC_SUBDIR_NAME)

# Check the list of files to ensure correctness
print(seqtab_files)
file.exists(seqtab_files)

# Load each sequence table into a list
seqtab_list <- lapply(seqtab_files, readRDS)
seqtab_list <- lapply(seqtab_list, as.matrix)



# Load metadata -----------------------------------

# Get the metadata for all samples (rows are samples and columns are metadata)
all_samples_table <- read.delim(ALL_SAMPLES_TSV, header = TRUE, row.names= 1)



# Filter samples ----------------------------------------
# Make a deep copy
filtered_sample_table <- all_samples_table[, ]

# Print number of samples to start
cat("Number of samples initially:", nrow(filtered_sample_table), "\n")

for (colname in names(EXCLUDE_ROWS_WITH_VALUES)) {
  vals <- EXCLUDE_ROWS_WITH_VALUES[[colname]]
  if (length(vals) > 0) {
    filtered_sample_table <- filtered_sample_table[!(filtered_sample_table[[colname]] %in% vals), ]
  }
}

for (colname in names(ONLY_INCLUDE_ROWS_WITH_VALUES)) {
  vals <- ONLY_INCLUDE_ROWS_WITH_VALUES[[colname]]
  filtered_sample_table <- filtered_sample_table[filtered_sample_table[[colname]] %in% vals, ]
}

# Exclude specific datasets
if (length(EXCLUDE_DATASETS) > 0) {
  filtered_sample_table <- filtered_sample_table[!(filtered_sample_table[[DATASET_ID_COLUMN]] %in% EXCLUDE_DATASETS), ]
}

# Exclude specific samples
filtered_sample_table <- filtered_sample_table[!(rownames(filtered_sample_table) %in% EXCLUDE_SAMPLES), ]

# Keep only specific columns
filtered_sample_table <- filtered_sample_table[, KEEP_ONLY_COLUMNS, drop=FALSE]

# Keep only samples found in seq_table
seq_table <- Reduce(mergeSequenceTables, seqtab_list)
filtered_sample_table <- filtered_sample_table[rownames(filtered_sample_table) %in% rownames(seq_table), ]

# Print number of samples to end
cat("Number of samples after filtering:", nrow(filtered_sample_table), "\n")

# Print the number of samples with each label in each dataset
samples_overview <- table(
  filtered_sample_table[[DATASET_ID_COLUMN]],
  filtered_sample_table[[LABEL_COLUMN]]
)
print(samples_overview)



# Combine ----------------------------------------

# Merge all the loaded sequence tables
seq_table <- Reduce(mergeSequenceTables, seqtab_list)

# Assign taxonomy
taxa <- assignTaxonomy(seq_table, TAXONOMY_TRAIN_SET, multithread=TRUE)
taxa <- addSpecies(taxa, SPECIES_ASSIGNMENT_SET, verbose=TRUE)
# Convert taxonomy to a format compatible with phyloseq
tax_table_obj <- tax_table(taxa)

# Construct the overall phyloseq object
ps_unfiltered <- phyloseq(otu_table(seq_table, taxa_are_rows=FALSE), 
                          tax_table_obj, sample_data(filtered_sample_table))

# Check phyloseq object before phylum filtering
ps_unfiltered

# Filter unwanted taxa (Cyanobacteria, Chloroflexi, uncharacterized/NA)
ps0 <- subset_taxa(ps_unfiltered, !is.na(Phylum) & !Phylum %in% c("Cyanobacteria", "Chloroflexi"))

# Check phyloseq object after phylum filtering
ps0



# Determine which ASVs should be filtered out ----------------------------

sample_data_df <- as.data.frame(sample_data(ps0))
datasets <- unique(sample_data_df[[DATASET_ID_COLUMN]])
tss_counts <- sweep(otu_table(ps0), 1, sample_sums(ps0), "/")

# Abundance filtering
asv_means_by_dataset <- sapply(datasets, function(ds) {
  ds_samples <- rownames(sample_data_df[sample_data_df[[DATASET_ID_COLUMN]] == ds, , drop=FALSE])
  colMeans(tss_counts[ds_samples, , drop=FALSE])
})
abundance_threshold_datasets <- round(IN_N_PROPORTION_DATASETS * length(datasets))
asv_abundance_count <- rowSums(asv_means_by_dataset >= MIN_AVERAGE_ABUNDANCE_IN_DATASET)
asv_to_keep_abundance <- asv_abundance_count >= abundance_threshold_datasets

cat("ASVs failing abundance filter:", sum(!asv_to_keep_abundance), "\n")

# Prevalence filtering
asv_prevalence_by_dataset <- sapply(datasets, function(ds) {
  ds_samples <- rownames(sample_data_df[sample_data_df[[DATASET_ID_COLUMN]] == ds, , drop=FALSE])
  colMeans((otu_table(ps0)[ds_samples, , drop=FALSE]) > 0)
})
prevalence_threshold_datasets <- round(IN_N_PROPORTION_DATASETS * length(datasets))
asv_prevalence_count <- rowSums(asv_prevalence_by_dataset >= MIN_PREVALENCE_IN_SAMPLES_IN_DATASET)
asv_to_keep_prevalence <- asv_prevalence_count >= prevalence_threshold_datasets

cat("ASVs failing prevalence filter:", sum(!asv_to_keep_prevalence), "\n")

asv_to_remove_abundance <- which(!asv_to_keep_abundance)
asv_to_remove_prevalence <- which(!asv_to_keep_prevalence)
cat("ASVs removed total (unique):", length(unique(c(asv_to_remove_abundance, asv_to_remove_prevalence))), "\n")

asv_to_keep <- asv_to_keep_abundance & asv_to_keep_prevalence

cat("ASVs after both filters:", sum(asv_to_keep), "\n")




# Present Total Percentage Abundance of Core ASVs -----------------------------------------

# Calculate mean percentage abundance of core ASVs for each dataset
core_asv_abundance <- sapply(datasets, function(ds) {
    ds_samples <- rownames(sample_data_df[sample_data_df[[DATASET_ID_COLUMN]] == ds, , drop=FALSE])
    mean(rowSums(tss_counts[ds_samples, asv_to_keep, drop=FALSE])) * 100
})
core_asv_df <- data.frame("Total Percentage Abundance of Core ASVs" = core_asv_abundance, row.names = datasets)
print(core_asv_df)

# Calculate percentage abundance of core ASVs for each sample
core_asv_sample_abundance <- data.frame(
    Dataset = rep(datasets, sapply(datasets, function(ds) {
        sum(sample_data_df[[DATASET_ID_COLUMN]] == ds)
    })),
    Sample = rownames(sample_data_df),
    Core_ASV_Abundance = unlist(lapply(datasets, function(ds) {
        ds_samples <- rownames(sample_data_df[sample_data_df[[DATASET_ID_COLUMN]] == ds, , drop=FALSE])
        rowSums(tss_counts[ds_samples, asv_to_keep, drop=FALSE]) * 100
    }))
)
raw_data_df <- core_asv_sample_abundance
write.csv(raw_data_df, TOTAL_ABUNDANCE_DATA_PATH, row.names = FALSE)


# Create a box-and-whisker plot
boxplot <- ggplot(raw_data_df, aes(x = Dataset, y = Core_ASV_Abundance)) +
    geom_boxplot() +
    labs(
        title = "Distribution of Core ASV Abundance Across Datasets",
        x = "Dataset",
        y = "Percentage Abundance of Core ASVs"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
    filename = TOTAL_ABUNDANCE_PLOT_PATH,
    plot = boxplot,
    width = 10,
    height = 6
)
print(boxplot)



# Filter samples based on total abundance of core ASVs -------------------------

# Identify samples to keep
samples_to_keep <- core_asv_sample_abundance$Core_ASV_Abundance / 100 >= MIN_TOTAL_ABUNDANCE_CORE_ASVS

# Filter the phyloseq object to include only samples meeting the criteria
ps0 <- prune_samples(
    rownames(sample_data(ps0)) %in% core_asv_sample_abundance$Sample[samples_to_keep],
    ps0
)

# Print updated sample count
cat("Number of samples after filtering by core ASV abundance:", nsamples(ps0), "\n")






# Filter ASVs -----------------------------------------

# Remove ASVs determined by the filtering section
ps1 <- prune_taxa(asv_to_keep, ps0)

# Check phyloseq object
ps1



# Transform abundances -----------------------------------

# Apply total sum scaling
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))

# Check phyloseq object
ps2



# Overview resulting data --------------------------------

# Print the number of samples with each label in each dataset
cat("Samples per label per dataset:\n")
print(table(sample_data(ps2)[[DATASET_ID_COLUMN]], sample_data(ps2)[[LABEL_COLUMN]]))

# Print the number of samples with each label overall
cat("Overall label counts:\n")
print(table(sample_data(ps2)[[LABEL_COLUMN]]))



# Export --------------------------------------------


# ps0 - Phyloseq object with no filtering of ASVs and untransformed or normalised read counts 
# ps1 - Phyloseq object after filtering of ASVs, but before transformation or normalisation of read counts
# ps2 - Phyloseq object after filtering of ASVs and transformation/normalisation of read counts

# Save phyloseq objects
saveRDS(ps0, PS0_output)
saveRDS(ps1, PS1_output)
saveRDS(ps2, PS2_output)


# Save taxonomy
write.csv(as.data.frame(tax_table(ps0)), PS0_taxa_output)
write.csv(as.data.frame(tax_table(ps1)), PS1_taxa_output)
write.csv(as.data.frame(tax_table(ps2)), PS2_taxa_output)

# Save OTU table and sample data from ps2 as TSV files
ps2_otu <- as.data.frame(as.matrix(otu_table(ps2)))
ps2_sample <- as.data.frame(as.matrix(sample_data(ps2)))

ps2_otu <- cbind(Sample.ID = rownames(ps2_otu), ps2_otu)
ps2_sample <- cbind(Sample.ID = rownames(ps2_sample), ps2_sample)

write.table(ps2_otu, PS2_otu_tsv_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(ps2_sample, PS2_sample_tsv_output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

