
# This script takes the results of the ASV truncation (from truncate_asvs.R) and combines selected datasets into single 
# phyloseq objects and filters ASVs according to presence and prevalence, ready for analysis.

# Outputs:
# ps0 -    Phyloseq object with no filtering of ASVs and untransformed or normalised read counts
# ps1 -    Phyloseq object after filtering of ASVs, but before transformation or normalisation of read counts
# ps2 -    Phyloseq object after filtering of ASVs and transformation/normalisation of read counts
# ps0_taxa.csv                    -    Taxonomy tables for the ps0 phyloseq object
# ps1_taxa.csv                    -    Taxonomy tables for the ps1 phyloseq object
# ps2_taxa.csv                    -    Taxonomy tables for the ps2 phyloseq object
# core_asv_sample_abundance.csv   -    Percentage abundance of core ASVs for each sample
# core_asv_abundance_boxplot.png  -    Boxplot showing distribution of core ASV abundance across datasets
# ps2_asv_table.tsv               -    ASV table from ps2
# ps2_sample_data.tsv             -    Sample data from ps2
# filtered_samples.tsv            -    Filtered sample table (samples used in the analysis)


# Requirements:
#  - R 4.5.1
#  - R packages: dada2, phyloseq, ggplot2, Biostrings
#  - SILVA databases (DADA2): silva_nr99_v138.1_train_set.fa.gz, silva_species_assignment_v138.1.fa.gz



# SET UP ====================================
# Load necessary packages
library(dada2)
library(phyloseq)
library(ggplot2)
library(Biostrings)



# Path to samples.tsv and low_read_samples.tsv files (from Celiac Microbiome Repository https://github.com/CeliacMicrobiomeRepo/celiac-repository/)
ALL_SAMPLES_TSV <- "/home/haig/Repos/celiac-repository/all_samples.tsv"
LOW_READ_SAMPLES_TSV <- "/home/haig/Repos/celiac-repository/low_read_samples.tsv"

# OPTIONS ===================================
# Directories of datasets to include
# Assumed to each include the following files:
#   - seqs.fna: FASTA file with all ASV sequences (output from DADA2)
#   - asv_abundances_transposed.tsv: TSV file with ASV abundances
DATASET_DIRS <- c(
    # V4 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_27_Fornasaro",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_49_Turjeman",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_102_Bodkhe",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_50_Bibbo",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_80_Garcia",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_68_Girdhar",

    # V3-V4 datasets ------------
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_136_Nobel",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_60_Shi",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_96_Quagliariello",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_227_Oscarsson",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_99_Lionetti",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_5_Senicar",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_27_Federica",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_20_Rawson",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_24_Giacomin",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_62_Tian",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_1211_Milletich",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_119_Salamon"

    # V4-V6 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_56_Laffaldano",

    # non-V4 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_20_Caminero",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_179_Verdu",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_25_Francavilla",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_39_Olivares",
    )
# Each dataset directory should contain the following 2 files:
# e.g.
#   v4_truncation_stool_prospective/*
#   v4_truncation_stool_active/*
#   v4_truncation_stool_treated/*
#   v4_truncation_duodenal_active/*
SEQS_FNA_FILE_PATH <- "v4_truncation_stool_active/seqs.fna"   # <--- [!!!] Change per analysis
ASV_ABUNDANCES_FILE_PATH <- "v4_truncation_stool_active/asv_abundances_transposed.tsv"   # <--- [!!!] Change per analysis


# Paths to databases
TAXONOMY_TRAIN_SET <- "/mnt/secondary/16S_databases/silva_nr99_v138.1_train_set.fa.gz"
SPECIES_ASSIGNMENT_SET <- "/mnt/secondary/16S_databases/silva_species_assignment_v138.1.fa.gz"


# Output paths
OUT_DIR <- "/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects"
ANALYSIS_DIR_NAME <- "stool_active_phyloseq_objects"   # <--- [!!!] Change per analysis
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
#
PS0_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps0.rds")
PS1_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps1.rds")
PS2_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps2.rds")
PS0_taxa_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps0_taxa.csv")
PS1_taxa_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps1_taxa.csv")
PS2_taxa_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps2_taxa.csv")
# Make the directory if it doesnt exist (for PS0_output and PS0_taxa_output)
dir.create(dirname(PS0_output), recursive = TRUE)
dir.create(dirname(PS0_taxa_output), recursive = TRUE)
# How much of the samples are made up by these core ASVs?
TOTAL_ABUNDANCE_DATA_PATH <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "core_asv_sample_abundance.csv")
TOTAL_ABUNDANCE_PLOT_PATH <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "core_asv_abundance_boxplot.png")
# TSV output paths
PS2_otu_tsv_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps2_asv_table.tsv")
PS2_sample_tsv_output <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "ps2_sample_data.tsv")
# Filtered sample table
FILTERED_SAMPLES_TSV_OUTPUT <- file.path(OUT_DIR, ANALYSIS_DIR_NAME, "filtered_samples.tsv")

# Record all console output to a log file
sink(file.path(OUT_DIR, ANALYSIS_DIR_NAME, "combination_console_output.log"), split = TRUE)


# Sample filtering options --------
# Column containing dataset IDs
DATASET_ID_COLUMN = "Dataset_ID"
# Columns to use as labels
LABEL_COLUMN = "Diagnosed_Celiac"   # <--- [!!!] Change per analysis ("Diagnosed_Celiac" or "Will_Develop_Celiac")
# Exclude data
EXCLUDE_ROWS_WITH_VALUES = list(
  # Exclude incomplete data
  "Dataset_ID" = c("-", "", NA),
  # Exclude shotgun datasets
  "Sequencing_Type" = c("shotgun", "", NA),
  # Exclude any data with significant factors
  "Any_Significant_Factor" = c("-", "", NA, TRUE)
)
# Include data. Each element in this list is a logical expression string that will be used to filter the samples.
# Expressions are combined with AND.
FILTERING_INCLUSION_RULES = list(               # <--- [!!!] Change per analysis

  # Stool Prospecitve ---
  # "Will_Develop_Celiac %in% c(TRUE, FALSE)",
  # "Sample_Site == 'stool'"


  # Stool Active ---
  "(Diagnosed_Celiac == TRUE & Gluten_Free_Diet == FALSE) | (Diagnosed_Celiac == FALSE & Gluten_Free_Diet == FALSE)",
  "Sample_Site == 'stool'"


  # Stool Treated ---
  # "(Diagnosed_Celiac == TRUE & Gluten_Free_Diet == TRUE) | (Diagnosed_Celiac == FALSE & Gluten_Free_Diet == FALSE)",
  # "Sample_Site == 'stool'"


  # Duodenum Active ---
  # "(Diagnosed_Celiac == TRUE & Gluten_Free_Diet == FALSE) | (Diagnosed_Celiac == FALSE & Gluten_Free_Diet == FALSE)",
  # "Sample_Site == 'duodenal'"

)
# Downsample Milletich dataset to 26 control samples
DOWNSAMPLE_MILLETICH <- FALSE   # <--- [!!!] Change per analysis (use for prospective analysis)
# Exclude non-illumina datasets
EXCLUDE_NON_ILLUMINA_DATA <- TRUE
# Exclude specific datasets
EXCLUDE_DATASETS = c()
# Keep only specific columns
KEEP_ONLY_COLUMNS = unique(c(DATASET_ID_COLUMN, LABEL_COLUMN, "Sample_Site", "Any_Significant_Factor", "Gluten_Free_Diet", "Seq_Tech", "Will_Develop_Celiac", "Diagnosed_Celiac", "Sex"))
# Exclude samples in LOW_READ_SAMPLES_TSV
EXCLUDE_SAMPLES <- unique(read.delim(LOW_READ_SAMPLES_TSV, header = TRUE)[, 2])


# ASV filtering options --------------
# Each ASV must occur in at least this proportion of samples at a non-zero abundance across all samples
MIN_PREVALENCE_OVER_ALL_SAMPLES <- 0.1

# A final filter to Samples ------------------
# Must have at least this proportion of their abundance made of the core ASVs
MIN_TOTAL_ABUNDANCE_CORE_ASVS <- 0.01 # (this should never happen anyway)


# Load sequence table files ----------------------------------------

# Get paths for fasta and abundance files
seqs_fna_files <- file.path(DATASET_DIRS, SEQS_FNA_FILE_PATH)
asv_abundances_files <- file.path(DATASET_DIRS, ASV_ABUNDANCES_FILE_PATH)

# Check the list of files to ensure correctness
cat("\n--- Checking Input Files ---\n")
cat("Checking for seqs.fna files:\n")
print(data.frame(path=seqs_fna_files, exists=file.exists(seqs_fna_files)))
cat("\nChecking for asv_abundances_transposed.tsv files:\n")
print(data.frame(path=asv_abundances_files, exists=file.exists(asv_abundances_files)))
cat("\n")


# Load each pair of files into a sequence table and add to list
seqtab_list <- lapply(seq_along(DATASET_DIRS), function(i) {
    dn <- DATASET_DIRS[i]
    seqs_fna_path <- seqs_fna_files[i]
    abund_path <- asv_abundances_files[i]

    # Check if files exist
    if (!file.exists(seqs_fna_path) || !file.exists(abund_path)) {
        warning(paste("Skipping dataset", dn, "due to missing files."))
        return(NULL)
    }

    # Read sequences and abundance table
    seqs <- readDNAStringSet(seqs_fna_path)
    # Using check.names=FALSE to prevent R from converting special characters in sample names
    abund_table <- read.delim(abund_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

    # Transpose abundance table to have samples as rows
    abund_table_t <- t(abund_table)

    # Create a map from ASV ID to sequence
    seq_map <- as.character(seqs)
    names(seq_map) <- names(seqs)

    # Match ASV IDs in abundance table with sequences
    asv_ids <- colnames(abund_table_t)

    # Ensure all ASV IDs from abundance table are in the FASTA file
    if (!all(asv_ids %in% names(seq_map))) {
        missing_asvs <- asv_ids[!(asv_ids %in% names(seq_map))]
        stop(paste("For dataset", dn, "not all ASV IDs from abundance table found in FASTA file. Missing ASVs:", paste(missing_asvs, collapse=", ")))
    }

    new_colnames <- seq_map[asv_ids]

    if (any(is.na(new_colnames)) || any(new_colnames == "")) {
        stop(paste("Some ASV IDs from abundance table could not be mapped to sequences for dataset", dn))
    }

    colnames(abund_table_t) <- unname(new_colnames)

    return(as.matrix(abund_table_t))
})

# Remove NULL elements from list (from skipped datasets)
seqtab_list <- seqtab_list[!sapply(seqtab_list, is.null)]

if (length(seqtab_list) == 0) {
    stop("No sequence tables could be loaded. Please check file paths and contents.")
}





# Load metadata -----------------------------------

# Get the metadata for all samples (rows are samples and columns are metadata)
all_samples_table <- read.delim(ALL_SAMPLES_TSV, header = TRUE, row.names= 2)

# Convert columns with only "True"/"False" strings to logical type
for (colname in names(all_samples_table)) {
  # Check if the column is of character type
  if (is.character(all_samples_table[[colname]])) {
    unique_vals <- unique(na.omit(all_samples_table[[colname]]))
    # Check if all unique non-NA values are some form of "true" or "false"
    if (length(unique_vals) > 0 && all(tolower(unique_vals) %in% c("true", "false"))) {
      cat("Converting column '", colname, "' to logical.\n", sep="")
      all_samples_table[[colname]] <- as.logical(all_samples_table[[colname]])
    }
  }
}
cat("\n")

# Filter samples ----------------------------------------
# Make a deep copy
filtered_sample_table <- all_samples_table[, ]

# Print number of samples to start
cat("\n--- Sample Filtering ---\n")
cat("Initial number of samples:", nrow(filtered_sample_table), "\n")

for (colname in names(EXCLUDE_ROWS_WITH_VALUES)) {
  vals <- EXCLUDE_ROWS_WITH_VALUES[[colname]]
  if (length(vals) > 0) {
    initial_count <- nrow(filtered_sample_table)
    filtered_sample_table <- filtered_sample_table[!(filtered_sample_table[[colname]] %in% vals), ]
    cat("Samples after excluding rows where '", colname, "' is one of [", paste(vals, collapse=", "), "]: ", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n", sep="")
  }
}

cat("\n")

for (rule in FILTERING_INCLUSION_RULES) {
  initial_count <- nrow(filtered_sample_table)
  # The rule is a string expression. We evaluate it in the context of the filtered_sample_table
  keep_indices <- which(with(filtered_sample_table, eval(parse(text = rule))))
  filtered_sample_table <- filtered_sample_table[keep_indices, ]
  cat("Samples after applying filter rule '", rule, "': ", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n", sep="")
}

cat("\n")

# Exclude specific datasets
if (length(EXCLUDE_DATASETS) > 0) {
  initial_count <- nrow(filtered_sample_table)
  filtered_sample_table <- filtered_sample_table[!(filtered_sample_table[[DATASET_ID_COLUMN]] %in% EXCLUDE_DATASETS), ]
  cat("Samples after excluding datasets [", paste(EXCLUDE_DATASETS, collapse=", "), "]: ", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n\n", sep="")
}

# Exclude specific samples
initial_count <- nrow(filtered_sample_table)
filtered_sample_table <- filtered_sample_table[!(rownames(filtered_sample_table) %in% EXCLUDE_SAMPLES), ]
cat("Samples after excluding specific samples listed in LOW_READ_SAMPLES_TSV:", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n\n")

# Keep only specific columns
initial_count <- nrow(filtered_sample_table)
filtered_sample_table <- filtered_sample_table[, KEEP_ONLY_COLUMNS, drop=FALSE]
cat("Applied column filter to keep only:", paste(KEEP_ONLY_COLUMNS, collapse=", "), "\n")
cat("Number of samples remains:", nrow(filtered_sample_table), "\n\n")

# Keep only samples found in seq_table (often removed non-V4 datasets)
seq_table <- Reduce(mergeSequenceTables, seqtab_list)
initial_count <- nrow(filtered_sample_table)
filtered_sample_table <- filtered_sample_table[rownames(filtered_sample_table) %in% rownames(seq_table), ]
cat("Samples after keeping only those present in the sequence tables:", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n\n")

# Optionally downsample (ONLY FOR Milletich dataset)
if (DOWNSAMPLE_MILLETICH) {
  set.seed(42)
  N_SAMPLES_TO_KEEP <- 26
  GROUP_DATASET_ID <- "16S_1211_Milletich"
  GROUP_LABEL <- FALSE
  if (N_SAMPLES_TO_KEEP > 0) {
    # Indices of samples in the group to downsample
    group_indices <- which(
      filtered_sample_table[[DATASET_ID_COLUMN]] == GROUP_DATASET_ID &
      filtered_sample_table[[LABEL_COLUMN]] == GROUP_LABEL
    )
    if (length(group_indices) > N_SAMPLES_TO_KEEP) {
      # Randomly select indices to keep from the group
      indices_to_keep_from_group <- sample(group_indices, N_SAMPLES_TO_KEEP)
      # Indices of all other samples
      other_indices <- setdiff(seq_len(nrow(filtered_sample_table)), group_indices)
      # Combine indices and subset the table, preserving original order
      final_indices <- sort(c(other_indices, indices_to_keep_from_group))
      initial_count <- nrow(filtered_sample_table)
      filtered_sample_table <- filtered_sample_table[final_indices, ]
      cat("Downsampling group '", GROUP_DATASET_ID, "' where '", LABEL_COLUMN, "' is '", GROUP_LABEL, "' to ", N_SAMPLES_TO_KEEP, " samples.\n", sep="")
      cat("Number of samples after downsampling:", nrow(filtered_sample_table), " (removed ", initial_count - nrow(filtered_sample_table), ")\n\n")
    }
  }
}

# Optionally remove non-illumina datasets
if (EXCLUDE_NON_ILLUMINA_DATA) {
  initial_count <- nrow(filtered_sample_table)
  illumina_mask <- grepl("Illumina", filtered_sample_table$Seq_Tech, ignore.case = TRUE)
  filtered_sample_table <- filtered_sample_table[illumina_mask, ]
  removed_count <- initial_count - nrow(filtered_sample_table)
  cat("After removing non-Illumina datasets:", nrow(filtered_sample_table), " (removed ", removed_count, ")\n\n", sep="")
}


# Remove datasets with an imbalance of case-control samples
dataset_is_balanced <- function(n_treated_celiac, n_active_celiac, n_healthy) {
  return(n_healthy >= 7 & (n_treated_celiac >= 7 | n_active_celiac >= 7))
}

initial_count <- nrow(filtered_sample_table)
dataset_ids <- unique(filtered_sample_table[[DATASET_ID_COLUMN]])
balanced_dataset_ids <- c()
excluded_datasets <- c()

for (dataset_id in dataset_ids) {
    dataset_df <- filtered_sample_table[filtered_sample_table[[DATASET_ID_COLUMN]] == dataset_id, ]
    
    # Identify celiac samples (diagnosed or prospective)
    celiac_mask <- (dataset_df$Will_Develop_Celiac %in% TRUE) | (dataset_df$Diagnosed_Celiac %in% TRUE)
    n_celiac <- sum(celiac_mask, na.rm = TRUE)

    # Split celiac samples into treated (on a gluten-free diet) and active
    if (n_celiac > 0) {
      celiac_df <- dataset_df[which(celiac_mask), ]
      n_treated_celiac <- sum(celiac_df$Gluten_Free_Diet %in% TRUE, na.rm = TRUE)
    } else {
      n_treated_celiac <- 0
    }
    n_active_celiac <- n_celiac - n_treated_celiac

    # Healthy controls are the remaining samples
    n_healthy <- nrow(dataset_df) - n_celiac

    if (dataset_is_balanced(n_treated_celiac, n_active_celiac, n_healthy)) {
        balanced_dataset_ids <- c(balanced_dataset_ids, dataset_id)
    } else {
        cat("Dataset '", dataset_id, "' is imbalanced (", n_treated_celiac, " treated celiacs, ", n_active_celiac, " active celiacs, ", n_healthy, " healthy controls)\n", sep="")
        excluded_datasets <- c(excluded_datasets, dataset_id)
    }
}

if (length(excluded_datasets) > 0) {
  cat("Excluding imbalanced datasets: [", paste(excluded_datasets, collapse=", "), "]\n", sep="")
}

filtered_sample_table <- filtered_sample_table[filtered_sample_table[[DATASET_ID_COLUMN]] %in% balanced_dataset_ids, ]
removed_count <- initial_count - nrow(filtered_sample_table)
cat("After removing datasets with imbalanced celiac vs healthy controls:", nrow(filtered_sample_table), " (removed ", removed_count, " samples from ", length(excluded_datasets), " datasets)\n\n", sep="")


# Print number of samples after downsampling
cat("--- Final Sample Table ---", "\n")
cat("Total samples remaining after all filtering:", nrow(filtered_sample_table), "\n\n")

# Print the number of samples with each label in each dataset
cat("Sample counts per dataset and label:\n")
samples_overview <- table(
  filtered_sample_table[[DATASET_ID_COLUMN]],
  filtered_sample_table[[LABEL_COLUMN]]
)
print(samples_overview)

# Print overall label counts
cat("\nOverall label counts:\n")
print(table(filtered_sample_table[[LABEL_COLUMN]]))
cat("\n")

# Save the filtered sample table
write.table(filtered_sample_table, FILTERED_SAMPLES_TSV_OUTPUT, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)






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

# An ASV must occur in at least 10% of samples (rounded down) at a non-zero abundance.
otu <- as.matrix(otu_table(ps0))
min_samples <- floor(MIN_PREVALENCE_OVER_ALL_SAMPLES * nsamples(ps0))
asv_to_keep <- colSums(otu > 0) >= min_samples

cat("\n--- ASV Filtering ---\n")
cat("Total number of samples:", nsamples(ps0), "\n")
cat("Minimum prevalence for ASVs:", MIN_PREVALENCE_OVER_ALL_SAMPLES, "\n")
cat("Minimum number of samples an ASV must be present in:", min_samples, "\n")
cat("Total ASVs before filtering:", ntaxa(ps0), "\n")
cat("ASVs to keep based on prevalence:", sum(asv_to_keep), "\n")
cat("ASVs to remove:", ntaxa(ps0) - sum(asv_to_keep), "\n\n")




# Present Total Percentage Abundance of Core ASVs -----------------------------------------

# Create a data.frame from sample_data
sample_data_df <- as(sample_data(ps0), "data.frame")

# Apply total sum scaling to the phyloseq object
ps0_tss <- transform_sample_counts(ps0, function(x) x / sum(x))

# Get the TSS-transformed counts
tss_counts <- as(otu_table(ps0_tss), "matrix")

# Get the dataset IDs
datasets <- unique(sample_data_df[[DATASET_ID_COLUMN]])

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
ggsave(filename = TOTAL_ABUNDANCE_PLOT_PATH, plot = boxplot, width = 10, height = 6, bg = "white")
print(boxplot)





# Check samples have a bare minimum % total abundance of core ASVs -------------------------

# Throw an error if any samples have less than MIN_TOTAL_ABUNDANCE_CORE_ASVS of their abundance made up by the core ASVs
if (any(core_asv_sample_abundance$Core_ASV_Abundance / 100 < MIN_TOTAL_ABUNDANCE_CORE_ASVS)) {
    stop("Some samples have less than MIN_TOTAL_ABUNDANCE_CORE_ASVS of their abundance made up by the core ASVs")
}





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

# Stop writing to console 
sink()
