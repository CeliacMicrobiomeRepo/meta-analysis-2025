# Pipeline for processing gut metagenomic raw reads into abundance tables + more
# This code is largely taken from  https://benjjneb.github.io/dada2/tutorial.html

# This script is used for paired end Illumina sequencing


# SET UP =============================
# Load necessary packages
library(ggplot2)
library(dada2)
library(Biostrings)
library(phyloseq)

# Set seed for reproducibility
set.seed(24)

# Path to dataset
DATASET_DIR <- "D:/metagenome_datasets/unprocessed_datasets_final/16S_119_Salamon"
# Path to raw data
RAW_DATA_DIR <- "D:/metagenome_datasets/unprocessed_datasets_final/16S_119_Salamon/fastqs"
list.files(RAW_DATA_DIR)
length(list.files(RAW_DATA_DIR))



# If a plots directory doesn't exist, create it
dir_path <- file.path(DATASET_DIR, "plots")
if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
}

# Define constants...
# Foward and reverse fastq file suffixes
FASTQ_SUFFIX_1 <- "_L001_R1.gz"
FASTQ_SUFFIX_2 <- "_L001_R2.gz"
# File suffixes for filtered read files
FILTERED_SUFFIX_F <- "_F_filtered.fastq.gz"
FILTERED_SUFFIX_R <- "_R_filtered.fastq.gz"
# Paths to databases
TAXONOMY_TRAIN_SET <- "D:/16S_databases/silva_nr99_v138.1_train_set.fa.gz"
SPECIES_ASSIGNMENT_SET <- "D:/16S_databases/silva_species_assignment_v138.1.fa.gz"


# Get file & sample names =============================
# Split forward and reverse fastq files + sort 
forward_reads <- sort(list.files(RAW_DATA_DIR, pattern=FASTQ_SUFFIX_1, full.names = TRUE))
reverse_reads <- sort(list.files(RAW_DATA_DIR, pattern=FASTQ_SUFFIX_2, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq OR SAMPLENAME.fastq
sample_names <- sapply(strsplit(gsub("\\.fastq\\.gz$", "", basename(forward_reads)), "_"), `[`, 1)

# Check if they look good
sample_names
print(length(sample_names))
print(length(forward_reads))



# Plot quality =================================
# Plot quality profiles of the forward and reverse reads
fw_quality_plots <- plotQualityProfile(forward_reads[1:2])
rv_quality_plots <- plotQualityProfile(reverse_reads[1:2])

# Show the plots
print(fw_quality_plots)
print(rv_quality_plots)



# Save the plots to DATASET_DIR
fw_plot_filename <- file.path(DATASET_DIR, "plots/forward_reads_quality_profile.png")
rv_plot_filename <- file.path(DATASET_DIR, "plots/reverse_reads_quality_profile.png")
ggsave(fw_plot_filename, plot = fw_quality_plots, width = 10, height = 8)
ggsave(rv_plot_filename, plot = rv_quality_plots, width = 10, height = 8)




# Filtering and trimming reads =============================
# Define file paths for filtered reads
filtered_forward_reads <- file.path(RAW_DATA_DIR, "filtered", paste0(sample_names, FILTERED_SUFFIX_F))
filtered_reverse_reads <- file.path(RAW_DATA_DIR, "filtered", paste0(sample_names, FILTERED_SUFFIX_R))
names(filtered_forward_reads) <- sample_names
names(filtered_reverse_reads) <- sample_names

# Filter and trim the reads
TRUNC_LENGTHS <- c(295, 210) # Optional and differs between datasets (!!!)
filter_summary <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads,
                                truncLen=TRUNC_LENGTHS, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                                compress=TRUE, multithread=FALSE, verbose=TRUE)

# Show the summary of the top few reads
head(filter_summary)



# Error Rates ========================
# Learn error rates
error_rates_forward <- learnErrors(filtered_forward_reads, multithread=TRUE)
error_rates_reverse <- learnErrors(filtered_reverse_reads, multithread=TRUE)

# Plot error rates
fw_error_plots <- plotErrors(error_rates_forward, nominalQ=TRUE)
rv_error_plots <- plotErrors(error_rates_reverse, nominalQ=TRUE)

# Show the plots
print(fw_error_plots)
print(rv_error_plots)

# Save the plots to DATASET_DIR
fw_err_plot_filename <- file.path(DATASET_DIR, "plots/forward_reads_error_rates.png")
rv_err_plot_filename <- file.path(DATASET_DIR, "plots/reverse_reads_error_rates.png")
ggsave(fw_err_plot_filename, plot = fw_error_plots, width = 10, height = 8)
ggsave(rv_err_plot_filename, plot = rv_error_plots, width = 10, height = 8)



# Sample Inference ========================
# Use DADA! (dereplicate sequences and infer sample composition)
dada_forward <- dada(filtered_forward_reads, err=error_rates_forward, multithread=TRUE, verbose=TRUE)
dada_reverse <- dada(filtered_reverse_reads, err=error_rates_reverse, multithread=TRUE, verbose=TRUE)
# Inspect the returned dada-class objects
dada_forward[[1]]
dada_reverse[[1]]



# Merge paired reads ===========================
# Merge paired-end reads
merged_pairs <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse, filtered_reverse_reads, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(merged_pairs[[1]])



# ASV Table construction =======================
# Construct sequence table
sequence_table <- makeSequenceTable(merged_pairs)
dim(sequence_table)

# Inspect distribution of sequence lengths
sequence_lengths <- table(nchar(getSequences(sequence_table)))
sequence_lengths

# Threshold with max and min sequence lengths
# sequence_table <- sequence_table[,nchar(colnames(sequence_table)) %in% 438:467] # Optional and differs between datasets (!!!)



# Chimeras ===============================================
# Remove chimeric sequences
sequence_table_no_chimera <- removeBimeraDenovo(sequence_table, method="consensus", multithread=TRUE, verbose=TRUE)
dim(sequence_table_no_chimera)
sum(sequence_table_no_chimera) / sum(sequence_table)



# Track Reads ==============================================
# Track reads through the pipeline
get_unique_count <- function(x) sum(getUniques(x))
tracking_reads <- cbind(filter_summary, sapply(dada_forward, get_unique_count), sapply(dada_reverse, get_unique_count), 
                        sapply(merged_pairs, get_unique_count), rowSums(sequence_table_no_chimera))
colnames(tracking_reads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(tracking_reads) <- sample_names
# Show tracking of the top few reads
head(tracking_reads)
# Show tracking of all reads
tracking_reads





# Make ASVs Reverse Complement (OPTIONAL) ==============================================
# (used for 16S_20_Caminero and 16S_179_Verdu)
colnames(sequence_table_no_chimera) <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(colnames(sequence_table_no_chimera))))




# Assign Taxonomy ==============================================
# Assign taxonomy and species
taxonomy <- assignTaxonomy(sequence_table_no_chimera, TAXONOMY_TRAIN_SET, multithread=TRUE)
taxonomy <- addSpecies(taxonomy, SPECIES_ASSIGNMENT_SET)

# Print top few taxonomies
taxonomy_display <- taxonomy
rownames(taxonomy_display) <- NULL
head(taxonomy_display)



# Phyloseq ==================================
# Create a metadata data frame (with no metadata)
sample_metadata <- data.frame(sample_names = sample_names)
# Ensure the row names are set to the sample names
rownames(sample_metadata) <- sample_metadata$sample_names
sample_metadata <- sample_data(sample_metadata)

# Create phyloseq object
phyloseq_obj <- phyloseq(otu_table(sequence_table_no_chimera, taxa_are_rows=FALSE), 
                         sample_metadata, tax_table(taxonomy))

# Add DNA sequences to phyloseq object
dna_sequences <- Biostrings::DNAStringSet(taxa_names(phyloseq_obj))
names(dna_sequences) <- taxa_names(phyloseq_obj)
phyloseq_obj <- merge_phyloseq(phyloseq_obj, dna_sequences)

# Display phyloseq object
phyloseq_obj



# Export ==========================================
# Export ASV table as a tsv file (asv_abundances.tsv into DATASET_DIR)
asv_abundances_path <- file.path(DATASET_DIR, "asv_abundances.tsv")
asv_abundances <- as.data.frame(otu_table(phyloseq_obj))
write.table(asv_abundances, file=asv_abundances_path, sep="\t", quote=FALSE, col.names=NA)

# Export ASV sequences as a fasta file (seqs.fna into DATASET_DIR)
asv_seqs_path <- file.path(DATASET_DIR, "seqs.fna")
asv_seqs <- getSequences(sequence_table_no_chimera)
names(asv_seqs) <- paste0(">ASV_", seq_along(asv_seqs))
asv_seqs_fasta <- paste(names(asv_seqs), asv_seqs, sep="\n")
write(asv_seqs_fasta, file=asv_seqs_path)

# Export taxonomy as a tsv file (taxonomy.tsv into DATASET_DIR)
taxonomy_path <- file.path(DATASET_DIR, "taxonomy.tsv")
taxonomy_df <- as.data.frame(tax_table(phyloseq_obj))
write.table(taxonomy_df, file=taxonomy_path, sep="\t", quote=FALSE, col.names=NA)

# Export tracking reads as a tsv file (tracking_reads.tsv into DATASET_DIR)
tracking_reads_path <- file.path(DATASET_DIR, "tracking_reads.tsv")
write.table(tracking_reads, file=tracking_reads_path, sep="\t", quote=FALSE, col.names=NA)

# Export sequence lengths (before any thresholding) as a tsv file (sequence_lengths.tsv into DATASET_DIR)
sequence_lengths_path <- file.path(DATASET_DIR, "sequence_lengths.tsv")
write.table(sequence_lengths, file=sequence_lengths_path, sep="\t", quote=FALSE, col.names=NA)

# Export the phyoseq object as a RDS file (phyloseq_obj.rds into DATASET_DIR)
phyloseq_obj_path <- file.path(DATASET_DIR, "phyloseq_obj.rds")
saveRDS(phyloseq_obj, file = phyloseq_obj_path)
# To read this use: phyloseq_obj <- readRDS(file.path(DATASET_DIR, "phyloseq_obj.rds"))

# Export sequence_table_no_chimera as a RDS file
sequence_table_no_chimera_path <- file.path(DATASET_DIR, "seqs_table.rds")
saveRDS(sequence_table_no_chimera, file = sequence_table_no_chimera_path)
# To read this use: sequence_table_no_chimera <- readRDS(file.path(DATASET_DIR, "sequence_table_no_chimera.rds"))

# Export filter summary as a tsv file
filter_summary_path <- file.path(DATASET_DIR, "filter_summary.tsv")
write.table(filter_summary, file=filter_summary_path, sep="\t", quote=FALSE, col.names=NA)

# Save merged pairs as RDS file
merged_pairs_path <- file.path(DATASET_DIR, "merged_pairs.rds")
saveRDS(merged_pairs, file = merged_pairs_path)
# To read this use: merged_pairs <- readRDS(file.path(DATASET_DIR, "merged_pairs.rds"))

# Save dada objects for forward and reverse reads as RDS files
dada_forward_path <- file.path(DATASET_DIR, "dada_forward.rds")
saveRDS(dada_forward, file = dada_forward_path)
# To read this use: dada_forward <- readRDS(file.path(DATASET_DIR, "dada_forward.rds"))

dada_reverse_path <- file.path(DATASET_DIR, "dada_reverse.rds")
saveRDS(dada_reverse, file = dada_reverse_path)
# To read this use: dada_reverse <- readRDS(file.path(DATASET_DIR, "dada_reverse.rds"))

# Transpose ASV abundances table and rename sequences with ASV IDs and "ASV_ID" in the first corner cell
asv_abundances_transposed <- t(asv_abundances)
asv_ids <- paste0("ASV_", seq_len(nrow(asv_abundances_transposed)))
rownames(asv_abundances_transposed) <- asv_ids
asv_abundances_transposed <- as.data.frame(asv_abundances_transposed)
asv_abundances_transposed <- cbind(ASV_ID = rownames(asv_abundances_transposed), asv_abundances_transposed)
rownames(asv_abundances_transposed) <- NULL
asv_abundances_transposed_path <- file.path(DATASET_DIR, "asv_abundances_transposed.tsv")
write.table(asv_abundances_transposed, file=asv_abundances_transposed_path, sep="\t", quote=FALSE, row.names=FALSE)

