
# Script for combination of 16S datasets according to their ASV alignments

# Steps:
# 1. Writes a single FASTA file with all ASV sequences combined from all datasets
#   - ASV IDs are updated by adding their dataset ID (e.g. "ASV_0028" -> "16S_24_Giacomin;ASV_0028")
# 2. Run mothur on the FASTA file to align to the reference alignment
#   - An existing alignment of full length 16S genes is used...
#   - Like silva.nr_v138_1.align (from: https://mothur.org/wiki/silva_reference_files/)
# 3. Truncates the alignment so that the overlap is the best
#   - Ensuring that all ASVs overlap with no overhangs the best you can
# 4. For each dataset, replaces the ASVs with truncated versions and export to a subdirectory
#   - ASVs which failed to align well are discarded
#   - Truncated ASVs which are identical are aggregated, summing abundances
#   - ASVs IDs are changed and are still specific to each dataset




# SET UP ====================================
# Load necessary packages
library(dada2)
library(phyloseq)
library(DECIPHER)
library(Biostrings)
library(ggplot2)

# Paths to databases
TAXONOMY_TRAIN_SET <- "D:/16S_databases/silva_nr99_v138.1_train_set.fa.gz"
SPECIES_ASSIGNMENT_SET <- "D:/16S_databases/silva_species_assignment_v138.1.fa.gz"


# OPTIONS ===================================
# Directories of datasets to include
DATASET_DIRS <- c("D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_80_Garcia",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_119_Salamon",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_136_Nobel",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_20_Rawson",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_60_Shi",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_96_Quagliariello",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_27_Fornasaro",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_49_Turjeman",
                  "D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_102_Bodkhe"
                  )

# Intermediate file for mothur alignment
ALL_ASVS_FASTA_PATH <- "all_core_v4_asvs.fasta"

# Truncated mothur alignment
ALL_ASVS_TRUNC_FASTA_PATH <- "all_core_v4_asvs_trunc.fasta"

# Subdirectory to output combination into for each dataset
OUTPUT_SUBDIR_NAME <- "core_v4_truncation"

# Option to skip taxonomic and phyloseq stuff (slow)
SAVE_PHYLOSEQ_AND_TAXONOMY <- TRUE

# About the output data
ABOUT_TXT_STR <- "This truncation is of all 'core V4 datasets' which are all datasets which cover V4, are duodenal/stool, aren't prospective, and are not imbalanced."

# Set the working directory (for mothur)
setwd("C:/Users/Haig/Documents/tent")





# Combine all ASVs into one FASTA file ===========================

# Open connection to FASTA file in "append" mode
con <- file(ALL_ASVS_FASTA_PATH, open = "a")

# For all datasets
for (dataset_dir in DATASET_DIRS) {
    # Read all ASVs from the dataset
    seqs_path <- file.path(dataset_dir, "seqs.fna")
    seqs <- readLines(seqs_path)
    
    # Extract the directory name
    dataset_name <- basename(dataset_dir)

    # Give ASVs new IDs (descriptions)
    # (e.g. "ASV_0028" -> "16S_24_Giacomin;ASV_0028")
    seqs <- gsub("^>(ASV_\\d+)", paste0(">", dataset_name, ";\\1"), seqs)
    
    # Write to the output FASTA file
    writeLines(seqs, con)
}

# Close the connection
close(con)






# Run mothur to align ASVs ==========================

# Run mothur (in terminal)
#   mothur

# Run alignment
#   e.g. 
#   align.seqs(fasta=all_core_v4_asvs.fasta, reference=D:/16S_databases/silva.nr_v138_1/silva.nr_v138_1.align)
#   which outputs -> all_core_v4_asvs.align
alignment_path <- sub("\\.fasta$", ".align", ALL_ASVS_FASTA_PATH)

# Filter empty columns in alignment (optional - if not enough RAM)
#   e.g. 
#   filter.seqs(fasta=all_core_v4_asvs.align)
#   which outputs -> all_core_v4_asvs.filter.fasta
alignment_path <- sub("\\.fasta$", ".filter.fasta", ALL_ASVS_FASTA_PATH)

# Read output alignment file
alignment <- readDNAStringSet(alignment_path, format="fasta")
alignment_names <- names(alignment)

# Convert to a matrix
seqs <- as.character(alignment)
seqs_split <- strsplit(seqs, '')
alignment_mat <- do.call(rbind, seqs_split)
rownames(alignment_mat) <- alignment_names

# Remove objects from memory
rm(alignment)
rm(seqs_split)
gc()

# Print before
cat("Dimensions of whole alignment matrix:\n")
print(dim(alignment_mat))





# Decide truncation positions ============================

# For every sequence, get starting and ending positions relative to the alignment
start_positions <- apply(alignment_mat, 1, function(seq) min(which(!seq %in% c(".", "-"))))
end_positions <- apply(alignment_mat, 1, function(seq) max(which(!seq %in% c(".", "-"))))

# Print min, max, and median for starting and ending positions
print("Start Positions:")
print(summary(start_positions))
print("End Positions:")
print(summary(end_positions))

# Print dataframe of frequencies of starting positions
start_freq <- as.data.frame(table(start_positions))
colnames(start_freq) <- c("Position", "Frequency")
print("Frequencies of starting positions:")
print(start_freq)

# Print dataframe of frequencies of ending positions
end_freq <- as.data.frame(table(end_positions))
colnames(end_freq) <- c("Position", "Frequency")
print("Frequencies of ending positions:")
print(end_freq)

# Define truncation positions 
trunc_start <- 14944
trunc_end <- 22587






# Truncate ASVs ============================
cat("Truncating sequences from position", trunc_start, "to", trunc_end, "\n")

# Truncate the whole alignment from trunc_start to trunc_end
alignment_mat_truncated <- alignment_mat[, trunc_start:trunc_end]

# Reconstruct sequences
truncated_sequences <- apply(alignment_mat_truncated, 1, paste0, collapse='')
names(truncated_sequences) <- rownames(alignment_mat_truncated)

# Remove gaps from truncated sequences for mapping
truncated_sequences_no_gaps <- gsub("[.-]", "", truncated_sequences)

# Write to file (ALL_ASVS_TRUNC_FASTA_PATH)
trunc_seqs_no_gaps_dna <- DNAStringSet(truncated_sequences_no_gaps)
names(trunc_seqs_no_gaps_dna) <- alignment_names
writeXStringSet(trunc_seqs_no_gaps_dna, filepath=ALL_ASVS_TRUNC_FASTA_PATH, format="fasta")

# Create mapping of all original ASV sequences to truncated sequences
sequence_mapping <- truncated_sequences_no_gaps
names(sequence_mapping) <- names(truncated_sequences)  # names are 'dataset_name;ASV_ID'

# Set poorly aligned sequences in the mapping to NA
# (sequences which do not cover this truncation range - i.e. they start after trunc_start or end before trunc_end)
sequence_mapping <- as.list(sequence_mapping)
sequence_mapping[start_positions > trunc_start | end_positions < trunc_end] <- NA

# Print min, max and median lengths of successfully truncated sequences
print("Truncated sequences:")
sequence_lengths <- sapply(sequence_mapping, function(x) if(!is.null(x)) nchar(x) else NA)
sequence_lengths <- sequence_lengths[!is.na(sequence_lengths)]
cat("Min length:", min(sequence_lengths), "\n")
cat("Median length:", median(sequence_lengths), "\n")
cat("Max length:", max(sequence_lengths), "\n")
cat("Std length:", sd(sequence_lengths), "\n")



# Export truncated ASVs ========================

# Initialize a list to store truncated lengths for each dataset
truncated_lengths_list <- list()

# For all datasets
for (dataset_dir in DATASET_DIRS) {
    
    # Extract the directory name (last part of the path)
    dataset_name <- basename(dataset_dir)
    cat("\n\n---\nStarting dataset:", dataset_name, "\n")
    
    # Load the ASV abundance table (matrix where column headers are original sequences)
    dataset_sequence_table <- readRDS(file.path(dataset_dir, "seqs_table.rds"))
    
    # Load the sequences and get mapping between ASV IDs and sequences
    seqs_fna_path <- file.path(dataset_dir, "seqs.fna")
    asv_seqs <- readDNAStringSet(seqs_fna_path, format="fasta")
    
    # Create mapping from ASV IDs to sequences
    asv_seq_strings <- as.character(asv_seqs)
    names(asv_seq_strings) <- names(asv_seqs)  # ASV IDs
    
    # Create full IDs used in the alignment
    full_ids <- paste0(dataset_name, ";", names(asv_seqs))
    
    # Get truncated sequences corresponding to each ASV
    truncated_seqs_for_dataset <- sequence_mapping[full_ids]
    
    # Calculate and print the percentage of ASVs with NULL headers (poorly aligned)
    null_percentage <- sum(is.na(truncated_seqs_for_dataset)) / length(truncated_seqs_for_dataset) * 100
    cat("Percentage of ASVs with NULL truncated sequences due to poor alignment:", round(null_percentage, 4), "%  (", sum(is.na(truncated_seqs_for_dataset)), "/", length(truncated_seqs_for_dataset), ")\n")
    
    # If no sequences worked
    if (sum(is.na(truncated_seqs_for_dataset)) == length(truncated_seqs_for_dataset)) {
        cat("WARNING: 0 sequences were able to be truncated for dataset ", dataset_name, "\n")
        cat("Skipping processing ", dataset_name, "\n")
        next
    }
    
    # Create mapping from original sequences to truncated sequences
    original_to_truncated <- truncated_seqs_for_dataset
    names(original_to_truncated) <- asv_seq_strings  # original sequences
    
    # Update column headers with truncated sequences using the sequence mapping
    # Replace original column names with truncated sequences if they exist in the mapping
    colnames(dataset_sequence_table) <- sapply(colnames(dataset_sequence_table), function(seq) {
        if (!is.null(original_to_truncated[[seq]])) {
            # cat("Original:", seq, "-> Truncated:", original_to_truncated[[seq]], "\n")
            original_to_truncated[[seq]]
        } else {
            #cat("Original:", seq, "-> Truncated: NA (no mapping available)\n")
            NA  # Mark as NA if there's no truncated sequence available
        }
    })
    # Drop columns with NA values (where truncated sequences were not available)
    dataset_sequence_table <- dataset_sequence_table[, !is.na(colnames(dataset_sequence_table))]
    
    # Aggregate by duplicate column names, summing abundances of identical truncated sequences
    truncated_aggregated_sequence_table <- t(rowsum(t(dataset_sequence_table), group=colnames(dataset_sequence_table)))
    
    # Print before duplicate aggregation
    cat("Number of ASVs before aggregation of duplicates after truncation:", ncol(dataset_sequence_table), "\n")
    cat("Number of ASVs after aggregation of duplicates after truncation:", ncol(truncated_aggregated_sequence_table), "\n")
    
    # Get truncated sequences for length calculation
    truncated_sequences_dataset <- colnames(truncated_aggregated_sequence_table)
    truncated_lengths <- nchar(truncated_sequences_dataset)
    
    # Store the truncated lengths in the list
    truncated_lengths_list[[dataset_name]] <- truncated_lengths
    
    # Print to user
    min_length <- min(truncated_lengths)
    median_length <- median(truncated_lengths)
    max_length <- max(truncated_lengths)
    cat("\nMinimum length of the truncated ASVs:", min_length, "\n")
    cat("Median length of the truncated ASVs:", median_length, "\n")
    cat("Maximum length of the truncated ASVs:", max_length, "\n\n")
    
    if (SAVE_PHYLOSEQ_AND_TAXONOMY) {
        # Assign taxonomy and species
        cat("Assigning taxonomy to truncated sequences...\n")
        truncated_sequence_strings <- colnames(truncated_aggregated_sequence_table)
        truncated_sequence_set <- DNAStringSet(truncated_sequence_strings)
        taxonomy <- assignTaxonomy(truncated_sequence_strings, TAXONOMY_TRAIN_SET, multithread=TRUE)
        taxonomy <- addSpecies(taxonomy, SPECIES_ASSIGNMENT_SET)
        
        # Sample names
        sample_names <- rownames(truncated_aggregated_sequence_table)
        
        # Create a metadata data frame (with no metadata)
        sample_metadata <- data.frame(sample_names = sample_names)
        rownames(sample_metadata) <- sample_metadata$sample_names
        sample_metadata <- sample_data(sample_metadata)
        
        # Create phyloseq object
        phyloseq_obj <- phyloseq(
            otu_table(truncated_aggregated_sequence_table, taxa_are_rows=FALSE),
            sample_metadata,
            tax_table(taxonomy)
        )
        
        # Add DNA sequences to phyloseq object
        dna_sequences <- Biostrings::DNAStringSet(taxa_names(phyloseq_obj))
        names(dna_sequences) <- taxa_names(phyloseq_obj)
        phyloseq_obj <- merge_phyloseq(phyloseq_obj, dna_sequences)
    }
    
    
    # Write outputs ---------------------------
    
    cat("Writing outputs for this dataset\n")
    output_dir <- file.path(dataset_dir, OUTPUT_SUBDIR_NAME)
    if (!dir.exists(output_dir)) dir.create(output_dir)
    
    # Export ASV table as a tsv file (asv_abundances.tsv into dataset_dir)
    asv_abundances_path <- file.path(output_dir, "asv_abundances.tsv")
    asv_abundances <- as.data.frame(otu_table(phyloseq_obj))
    write.table(asv_abundances, file=asv_abundances_path, sep="\t", quote=FALSE, col.names=NA)
    
    # Export ASV sequences as a fasta file (seqs.fna into dataset_dir)
    asv_seqs_path <- file.path(output_dir, "seqs.fna")
    asv_seqs <- colnames(truncated_aggregated_sequence_table)
    names(asv_seqs) <- paste0(">ASV_", seq_along(asv_seqs))
    asv_seqs_fasta <- paste(names(asv_seqs), asv_seqs, sep="\n")
    write(asv_seqs_fasta, file=asv_seqs_path)
    
    if (SAVE_PHYLOSEQ_AND_TAXONOMY) {
        # Export taxonomy as a tsv file (taxonomy.tsv into dataset_dir)
        taxonomy_path <- file.path(output_dir, "taxonomy.tsv")
        taxonomy_df <- as.data.frame(tax_table(phyloseq_obj))
        write.table(taxonomy_df, file=taxonomy_path, sep="\t", quote=FALSE, col.names=NA)
        
        # Export the phyoseq object as a RDS file (phyloseq_obj.rds into dataset_dir)
        phyloseq_obj_path <- file.path(output_dir, "phyloseq_obj.rds")
        saveRDS(phyloseq_obj, file = phyloseq_obj_path)
    }
    
    # Export truncated_aggregated_sequence_table as a RDS file
    truncated_aggregated_sequence_table_path <- file.path(output_dir, "seqs_table.rds")
    saveRDS(truncated_aggregated_sequence_table, file = truncated_aggregated_sequence_table_path)
    
    # Transpose ASV abundances table and rename sequences with ASV IDs and "ASV_ID" in the first corner cell
    asv_abundances_transposed <- t(asv_abundances)
    asv_ids <- paste0("ASV_", seq_len(nrow(asv_abundances_transposed)))
    rownames(asv_abundances_transposed) <- asv_ids
    asv_abundances_transposed <- as.data.frame(asv_abundances_transposed)
    asv_abundances_transposed <- cbind(ASV_ID = rownames(asv_abundances_transposed), asv_abundances_transposed)
    rownames(asv_abundances_transposed) <- NULL
    asv_abundances_transposed_path <- file.path(output_dir, "asv_abundances_transposed.tsv")
    write.table(asv_abundances_transposed, file=asv_abundances_transposed_path, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Write about text
    about_txt_path <- file.path(output_dir, "about.txt")
    about_txt_content <- paste0(
        ABOUT_TXT_STR, "\n\n",
        "Date of Combination: ", Sys.Date(), "\n",
        "Dataset: ", dataset_name, "\n",
        "Taxonomy Training Set Path: ", TAXONOMY_TRAIN_SET, "\n",
        "Species Assignment Set Path: ", SPECIES_ASSIGNMENT_SET, "\n",
        "Sequences were truncated based on positions from ", trunc_start, " to ", trunc_end, " in the alignment."
    )
    write(about_txt_content, file = about_txt_path)
    
    
    cat(paste("Dataset done:", dataset_dir, "\n"))
}
cat("Done!\n")


# Convert the list to a data frame for ggplot
truncated_lengths_df <- do.call(rbind, lapply(names(truncated_lengths_list), function(name) {
    data.frame(Dataset = name, Truncated_Length = truncated_lengths_list[[name]])
}))


# Box-and-whisker plot for all datasets' truncated lengths
boxplot <- ggplot(truncated_lengths_df, aes(x = Dataset, y = Truncated_Length)) +
    geom_boxplot() +
    labs(title = "Truncated ASV Lengths per Dataset", x = "Dataset", y = "Truncated Length") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(boxplot)


# Prepare a data frame with normalized densities for each dataset
density_data <- do.call(rbind, lapply(names(truncated_lengths_list), function(dataset_name) {
    lengths <- truncated_lengths_list[[dataset_name]]
    density_obj <- density(lengths)  # Calculate density
    
    # Normalize the density by dividing by the maximum value
    data.frame(
        Truncated_Length = density_obj$x,
        Density = density_obj$y / max(density_obj$y),  # Normalize to max 1
        Dataset = dataset_name
    )
}))

# Frequencies of lengths between datasets (density plot)
density_plot <- ggplot(density_data, aes(x = Truncated_Length, y = Density, color = Dataset)) +
    geom_line() +
    labs(title = "Normalized Density of Truncated ASV Lengths Across Datasets",
         x = "Truncated Length",
         y = "Normalized Density") +
    theme_minimal() +
    theme(legend.position = "right")
print(density_plot)


# Make a directory in the working directory called OUTPUT_SUBDIR_NAME and put in the plots and also write the dataframes end_freq and start_freq
# Define the output directory name
output_main_dir <- file.path(getwd(), OUTPUT_SUBDIR_NAME)
if (!dir.exists(output_main_dir)) dir.create(output_main_dir)

# Save start_freq and end_freq data frames as CSV files
write.csv(start_freq, file = file.path(output_main_dir, "start_freq.csv"), row.names = FALSE)
write.csv(end_freq, file = file.path(output_main_dir, "end_freq.csv"), row.names = FALSE)

# Save plots to the output directory
# Save the boxplot as a PNG
boxplot_path <- file.path(output_main_dir, "truncated_asv_lengths_boxplot.png")
ggsave(boxplot_path, plot = boxplot, width = 8, height = 6, dpi = 300)

# Save the density plot as a PNG
density_plot_path <- file.path(output_main_dir, "truncated_asv_lengths_density_plot.png")
ggsave(density_plot_path, plot = density_plot, width = 8, height = 6, dpi = 300)

cat("Data frames and plots have been saved to the output directory:", output_main_dir, "\n")


