# Differential Abundance Code Using ANCOM-BC2

# Takes a phyloseq object with non-transformed counts and filters ASVs to remove zero variance ASVs
# This phyloseq object likely originated from the combine_trunc_asvs_to_phyloseq.R script

# Requirements:
#  - R 4.4.0
#  - R packages: phyloseq, ANCOMBC, dplyr, tidyr, metafor, tibble, readr, stringr, purrr



# ── 0) Set Up ────────────────────────────

# Imports
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(metafor)
library(tibble)
library(readr)
library(stringr)
library(purrr)

# Primary biological factor
#   ("Diagnosed_Celiac" or "Will_Develop_Celiac")
GROUPING_VAR    <- "Diagnosed_Celiac"

# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "~/Repos/meta-analysis/analysis/differential_abundance_results/stool_active"
dir.create(out_dir, recursive = TRUE)



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1 object should be used (filtered ASVs, but non-transformed counts)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps <- readRDS("~/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps1.rds")

# 1. Coerce your sample_data to a data.frame
sd_df <- as(sample_data(ps), "data.frame")

# 2. Find all Dataset_IDs that have at least one sample
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))

# 3. Split into a named list, pruning both samples and zero‐count taxa
ps_list <- lapply(dataset_ids, function(d) {
  # get the sample names for this dataset
  samps <- rownames(sd_df)[sd_df$Dataset_ID == d]
  # prune to only those samples
  ps_sub <- prune_samples(samps, ps)                        
  # drop any taxa that now have zero counts
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids

# 4. Quick check
str(ps_list)

# Note some ASVs will have zero variance and will need to be removed before ANCOM
bad_asvs <- c("")

# Iterate over each element of ps_list, removing the unwanted ASVs
ps_list_no_bad <- lapply(
    ps_list,
    function(ps) {
        keep <- !taxa_names(ps) %in% bad_asvs   # logical vector: TRUE = keep this ASV
        prune_taxa(keep, ps)
    })
ps_list <- ps_list_no_bad




# ── 2) Clean up each phyloseq: drop zero‐var ASVs ────────────────────────────
ps_list_clean <- lapply(ps_list, function(ps) {
  filter_taxa(ps, function(x) var(x) > 0, prune = TRUE)
})



# ── 3) Functions ────────────────────────────
# helper utilities 
get_slope_cols <- function(df, prefix) {
    cols <- grep(paste0("^", prefix), names(df), value = TRUE)
    cols[!grepl("Intercept", cols)]
}

get_pq_col <- function(df, type = c("p", "q")) {
    type <- match.arg(type)
    # Match both "p_Diagnosed_CeliacTRUE" and "p_val_Diagnosed_CeliacTRUE" (same for q)
    pat <- paste0("^", type, "(_val)?_.*", GROUPING_VAR)
    grep(pat, names(df), value = TRUE)[1]
}



# ── 4) Run ANCOM‑BC² per dataset ───────────────────────────────────────────
ancom_res_list <- list()
per_dataset_sig <- list()

for (d in names(ps_list_clean)) {
    res <- ancombc2(
        data         = ps_list_clean[[d]],
        fix_formula  = GROUPING_VAR,
        group        = GROUPING_VAR,
        p_adj_method = "fdr",
        struc_zero   = TRUE,
        pseudo_sens  = TRUE,
        prv_cut      = 0.10,
        lib_cut      = 1000,
        alpha        = 0.05,
        n_cl         = 4,
        verbose      = TRUE
    )$res
    
    # add significance column (p < .05 & q < .05)
    p_col <- get_pq_col(res, "p")
    q_col <- get_pq_col(res, "q")
    if (!is.na(p_col) && !is.na(q_col)) {
        res <- res %>%
            mutate(sig_pq = !is.na(.data[[p_col]]) & !is.na(.data[[q_col]]) &
                       .data[[p_col]] < 0.05 & .data[[q_col]] < 0.05)
    } else {
        res$sig_pq <- FALSE
    }
    
    write.csv(res, file = file.path(out_dir, sprintf("ancombc2_%s.csv", d)), row.names = FALSE)
    
    ancom_res_list[[d]]  <- res
    per_dataset_sig[[d]] <- tibble(taxon = res$taxon, sig = res$sig_pq)
}



# ── 5) Wide beta / se table for meta‑analysis ──────────────────────────────
meta_df <- bind_rows(lapply(names(ancom_res_list), function(d) {
    df <- ancom_res_list[[d]]
    lfc <- get_slope_cols(df, paste0("lfc_", GROUPING_VAR))
    se  <- get_slope_cols(df, paste0("se_", GROUPING_VAR))
    if (length(lfc) != 1 || length(se) != 1) return(NULL)
    df %>%
        select(taxon, all_of(lfc), all_of(se)) %>%
        rename_with(~"beta", all_of(lfc)) %>%
        rename_with(~"se",   all_of(se)) %>%
        mutate(dataset = d)
})) %>%
    pivot_wider(id_cols = taxon, names_from = dataset, values_from = c(beta, se), names_sep = "__")

write.csv(meta_df, file = file.path(out_dir, "meta_input_beta_se.csv"), row.names = FALSE)



# ── 6) Random‑effects REML meta‑analysis ───────────────────────────────────
beta_cols <- grep("^beta__", names(meta_df), value = TRUE)
se_cols   <- sub("^beta__", "se__", beta_cols)

pooled_df <- bind_rows(lapply(seq_len(nrow(meta_df)), function(i) {
    row  <- meta_df[i, ]
    yi   <- as.numeric(row[beta_cols])
    sei  <- as.numeric(row[se_cols])
    keep <- is.finite(yi) & is.finite(sei) & (sei > 0)
    if (sum(keep) < 2) return(NULL)
    m <- rma(yi = yi[keep], sei = sei[keep], method = "REML")
    data.frame(taxon = row$taxon,
               beta_pooled = m$b,
               se_pooled   = m$se,
               ci_lb       = m$ci.lb,
               ci_ub       = m$ci.ub,
               pval        = m$pval,
               p_adj       = p.adjust(m$pval, method = "BH"),
               tau2        = m$tau2,
               Q           = m$QE,
               df_Q        = m$k - 1,
               Q_p         = m$QEp,
               I2          = m$I2)
})) %>%
    mutate(sig_pq_meta = pval < 0.05 & p_adj < 0.05)

write.csv(pooled_df, file = file.path(out_dir, "pooled_asv_results.csv"), row.names = FALSE)



# ── 7) Combined (all‑samples) ANCOM‑BC² run ────────────────────────────────
comb_res <- ancombc2(
    data         = ps,
    fix_formula  = GROUPING_VAR,
    group        = GROUPING_VAR,
    p_adj_method = "fdr",
    struc_zero   = TRUE,
    pseudo_sens  = TRUE,
    prv_cut      = 0.10,
    lib_cut      = 1000,
    alpha        = 0.05,
    n_cl         = 4,
    verbose      = TRUE
)$res

p_comb <- get_pq_col(comb_res, "p")
q_comb <- get_pq_col(comb_res, "q")
if (!is.na(p_comb) && !is.na(q_comb)) {
    comb_res <- comb_res %>%
        mutate(sig_pq_comb = !is.na(.data[[p_comb]]) & !is.na(.data[[q_comb]]) &
                   .data[[p_comb]] < 0.05 & .data[[q_comb]] < 0.05)
} else {
    comb_res$sig_pq_comb <- FALSE
}

write.csv(comb_res, file = file.path(out_dir, "ancombc2_combined.csv"), row.names = FALSE)



# ── 8) Heterogeneity table ─────────────────────────────────────────────────
het_df <- pooled_df %>% select(taxon, tau2, Q, df_Q, Q_p, I2)
write.csv(het_df, file = file.path(out_dir, "heterogeneity_per_asv.csv"), row.names = FALSE)

message("Per‑dataset, pooled, and combined results written to: ", out_dir)



# ── 9) Confidence Scoring ─────────────────────────────────────────────────
# (the +5 / +2.5 / +1.25 scheme)


# ── 9.1)  meta-analysis (pooled REML)  
pooled_meta <- read_csv(
  file.path(out_dir, "pooled_asv_results.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    taxon,
    beta_pooled,
    sig_pq_meta = sig_pq_meta,
    sign_meta   = sign(beta_pooled),
    score_pool  = ifelse(sig_pq_meta, 5, 0)
  )

# ── 9.2)  combined run (all samples together)  
comb <- read_csv(
  file.path(out_dir, "ancombc2_combined.csv"),
  show_col_types = FALSE
)

lfc_comb <- grep(paste0("^lfc_.*", GROUPING_VAR), names(comb), value = TRUE)[1]

comb_sig <- comb %>%
  transmute(
    taxon,
    sig_pq_comb = sig_pq_comb,
    sign_comb   = sign(.data[[lfc_comb]])
  )

# ── 9.3)  per-dataset significance & direction  
dataset_files <- list.files(
  out_dir,
  pattern = "^ancombc2_.*\\.csv$",
  full.names = TRUE
) |>
  # drop the combined file
  keep(~ !str_detect(.x, "combined"))

per_dataset <- lapply(dataset_files, function(f) {
  dname <- str_remove(basename(f), "^ancombc2_|\\.csv$")
  df    <- read_csv(f, show_col_types = FALSE)

  lfc_col <- grep(paste0("^lfc_.*", GROUPING_VAR), names(df), value = TRUE)[1]

  df %>%
    transmute(
      taxon,
      !!paste0("sig_", dname)  := sig_pq,
      !!paste0("dir_", dname)  := sign(.data[[lfc_col]])
    )
})

per_ds_wide <- purrr::reduce(
  per_dataset,                 # <- the list of data frames
  dplyr::full_join,            # <- the reducer function
  by = "taxon"
) %>%
  mutate(
    across(starts_with("sig_"), ~replace_na(., FALSE)),
    across(starts_with("dir_"), ~replace_na(.,  0))
  )


# ── 9.4)  calculate cross-dataset score (needs *matching sign*)  
sig_cols <- grep("^sig_", names(per_ds_wide), value = TRUE)
dir_cols <- grep("^dir_", names(per_ds_wide), value = TRUE)

indiv_score_df <- per_ds_wide %>%
  left_join(select(pooled_meta, taxon, sign_meta), by = "taxon") %>%
  rowwise() %>%
  mutate(
    # datasets where ASV was tested (dir ≠ 0 or sig flag present)
    n_tested = sum(abs(c_across(all_of(dir_cols))) > 0),
    # datasets where it's significant *and* direction matches pooled sign
    n_sig_ok = sum(c_across(all_of(sig_cols)) &
                   (c_across(all_of(dir_cols)) == sign_meta)),
    score_indiv = case_when(
      n_tested > 0 & n_sig_ok == n_tested              ~ 2.5,   # all match
      n_tested > 0 & n_sig_ok >= ceiling(n_tested / 2) ~ 1.25,  # ≥ half match
      TRUE                                             ~ 0
    )
  ) %>%
  ungroup() %>%
  select(taxon, score_indiv)

# ── 9.5)  merge & final score 
confidence_df <- pooled_meta %>%
  full_join(comb_sig,        by = "taxon") %>%
  full_join(indiv_score_df,  by = "taxon") %>%
  mutate(
    # +2.5 only if combined run is significant *and* sign matches meta
    score_comb = ifelse(sig_pq_comb & sign_comb == sign_meta, 2.5, 0),
    # replace any missing score with 0
    across(starts_with("score_"), ~replace_na(., 0)),
    conf_score = score_pool + score_indiv + score_comb
  ) %>%
  arrange(desc(conf_score), taxon)

write_csv(confidence_df,
          file.path(out_dir, "asv_confidence_scores.csv"))

message("✓  asv_confidence_scores.csv written to ", out_dir)



# ── 10) Attach taxonomy from phyloseq object (optional) ──────────────────────


tax_df <- tax_table(ps) |>
  as.data.frame() |>
  rownames_to_column("taxon")          # keeps ASV/OTU ID in a column

confidence_tax <- confidence_df |>
  left_join(tax_df, by = "taxon")      # keeps all score columns + taxonomy

write_csv(confidence_tax,
          file.path(out_dir, "asv_confidence_scores_w_tax.csv"))

message("✓  taxonomy-annotated confidence file written to ",
        file.path(out_dir, "asv_confidence_scores_w_tax.csv"))



# ── 11) Creating a csv with all the key results ──────────────────────────────

# 11.1) master confidence table (scores already computed earlier)
conf_df <- read_csv(file.path(out_dir,
                              "asv_confidence_scores_w_tax.csv"),
                    show_col_types = FALSE)

# 11.2) pooled CI only (conf_df already has beta_pooled)  ---------------------
pool_df <- read_csv(file.path(out_dir, "pooled_asv_results.csv"),
                    show_col_types = FALSE) %>%
  select(taxon, ci_lb, ci_ub)      # <-- beta_pooled left out


# 11.3) mean TSS relative abundance from combined phyloseq
ra_vec <- taxa_sums(ps) /
          sum(taxa_sums(ps))
ra_tbl <- tibble(taxon = names(ra_vec),
                 mean_RA = ra_vec)

rank_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")

final_df <- conf_df %>%
  left_join(pool_df, by = "taxon") %>%   # adds ci_lb / ci_ub, no clash
  left_join(ra_tbl,  by = "taxon") %>%
  filter(conf_score >= 5) %>%
  mutate(
    source  = str_trim(paste0(
                if_else(score_pool  > 0, "Meta; ",        ""),
                if_else(score_comb  > 0, "Combined; ",    ""),
                if_else(score_indiv > 0, "Per-dataset; ", "")
              )) |> stringr::str_remove(";\\s*$"),
    Display = coalesce(!!!syms(rank_order), taxon)
  ) %>%
  select(taxon, Display, all_of(rank_order),
         conf_score, source,
         beta_pooled, ci_lb, ci_ub,
         mean_RA,
         score_pool, score_comb, score_indiv)

write_csv(final_df, file.path(out_dir, "asv_high_confidence.csv"))

