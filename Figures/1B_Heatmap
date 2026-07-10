############################################
# 1. Load libraries
############################################

library(DESeq2)
library(dplyr)
library(openxlsx)
library(rtracklayer)
library(tibble)
library(stringr)

############################################
# 2. Read featureCounts matrices
############################################

se <- read.delim("../Processing/single_end/counts_unique/gene_counts.txt", comment.char = "#")
pe <- read.delim("../Processing/paired_end/counts_unique/gene_counts.txt", comment.char = "#")

se_counts <- se[, c("Geneid", colnames(se)[7:ncol(se)])]
pe_counts <- pe[, c("Geneid", colnames(pe)[7:ncol(pe)])]

combined <- full_join(se_counts, pe_counts, by = "Geneid")
combined[is.na(combined)] <- 0

rownames(combined) <- combined$Geneid
combined$Geneid <- NULL

counts_matrix <- as.matrix(combined)
mode(counts_matrix) <- "numeric"

############################################
# 3. Load TPM (SE + PE)
############################################

tpm_se <- read.delim("../Processing/single_end/gene_vs_sample_matrix_TPM.tsv")
tpm_pe <- read.delim("../Processing/paired_end/gene_vs_sample_matrix_TPM.tsv")

tpm_combined <- full_join(tpm_se, tpm_pe, by = "Gene")
tpm_combined[is.na(tpm_combined)] <- 0

rownames(tpm_combined) <- tpm_combined$Gene
tpm_combined$Gene <- NULL

tpm_matrix <- as.matrix(tpm_combined)
mode(tpm_matrix) <- "numeric"

############################################
# 5. Clean gene IDs
############################################

clean_ids <- function(x) {
  x <- sub("\\..*", "", x)
  x <- sub("::.*", "", x)
  return(x)
}

rownames(counts_matrix) <- clean_ids(rownames(counts_matrix)) 
rownames(tpm_matrix)   <- clean_ids(rownames(tpm_matrix))

############################################
# 6. Clean sample names
############################################

colnames(counts_matrix) <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts_matrix))

tpm_matrix   <- tpm_matrix[, colnames(counts_matrix)]

############################################
# 7. DO NOT ALIGN GENES (process separately)
############################################

# Keep matrices independent:
# counts_matrix → DESeq2
# tpm_matrix → TPM analysis
# count_salmon → CPM

############################################
# 9. TE / gene split (TPM + CPM separately)
############################################

# TPM
gene_tpm <- tpm_matrix[!grepl("^TE_", rownames(tpm_matrix)), ]
te_tpm   <- tpm_matrix[grepl("^TE_", rownames(tpm_matrix)), ]

te_family_tpm <- rowsum(te_tpm, rownames(te_tpm))

final_tpm <- rbind(gene_tpm, te_family_tpm)

 


############################################
# 10. Load metadata
############################################

meta1 <- read.csv("PRJNA606047_SraRunTable.csv")
meta2 <- read.csv("PRJNA282396_SraRunTable.csv")
meta3 <- read.csv("PRJNA773654_SraRunTable.csv")

meta1_small <- meta1 %>% select(Run, BioProject, cell_line, treatment)

meta2_small <- meta2 %>%
  dplyr::select(Run, BioProject, cell_line, source_name) %>%
  dplyr::rename(treatment = source_name)

meta3_small <- meta3 %>% select(Run, BioProject, cell_line, treatment)

all_meta <- bind_rows(meta1_small, meta2_small, meta3_small)

############################################
# 11. Clean metadata
############################################

meta_clean <- all_meta %>%
  mutate(
    Condition = case_when(
      str_detect(treatment, regex("hypoxia|oxygen", TRUE)) ~ "hypoxia",
      str_detect(treatment, regex("normoxia|control", TRUE)) ~ "normoxia",
      TRUE ~ NA_character_
    ),
    cell_line = case_when(
      cell_line == "A549" & BioProject == "PRJNA606047" ~ "A549_1",
      cell_line == "A549" & BioProject == "PRJNA773654" ~ "A549_2",
      str_detect(treatment, "HIF1Anull") ~ "HIF1AKO",
      TRUE ~ cell_line
    )
  ) %>%
  filter(!is.na(Condition))

############################################
# 12. Remove unwanted samples
############################################

remove_runs <- c("SRR16531801", "SRR16531802", "SRR16531803")

meta_clean <- meta_clean %>%
  filter(!Run %in% remove_runs)

counts_matrix <- counts_matrix[, !(colnames(counts_matrix) %in% remove_runs)]
tpm_matrix    <- tpm_matrix[, !(colnames(tpm_matrix) %in% remove_runs)]


############################################
# 13. Align metadata
############################################

valid_runs <- intersect(colnames(counts_matrix), meta_clean$Run)

counts_matrix <- counts_matrix[, valid_runs]
tpm_matrix    <- tpm_matrix[, valid_runs]

meta_aligned <- meta_clean %>%
  dplyr::filter(Run %in% valid_runs) %>%
  dplyr::slice(match(valid_runs, Run))

############################################
# 14. Gene annotation
############################################

gtf <- import("genes_plus_TE.sorted_hg38.gtf.gz")

gene_map <- gtf %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::mutate(gene_id = sub("\\..*", "", gene_id)) %>%   
  dplyr::distinct() %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(gene_name = dplyr::first(gene_name))

############################################
# 15. Excel output (DE + TPM + CPM)
############################################

wb <- createWorkbook()

cell_lines <- unique(meta_aligned$cell_line)
cell_lines <- setdiff(unique(meta_aligned$cell_line), "HCT116")
for (cell in cell_lines) {
  set.seed(0)
  meta_sub <- meta_aligned %>% filter(cell_line == cell)
  counts_sub <- counts_matrix[, meta_sub$Run]
  tpm_sub    <- final_tpm[, meta_sub$Run]
  
  meta_sub <- meta_sub[match(colnames(counts_sub), meta_sub$Run), ]
  
  if (length(unique(meta_sub$Condition)) < 2) next
  
  keep_genes <- rownames(counts_sub)[rowSums(counts_sub) >= 10]
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  tpm_sub    <- tpm_sub[rownames(tpm_sub) %in% keep_genes, , drop = FALSE]
  
  ##########################################
  # DESeq2
  ##########################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = ~ Condition
  )
  
  dds$Condition <- relevel(dds$Condition, ref = "normoxia")
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  res_df <- res_df %>%
    left_join(gene_map, by = "gene_id") %>%
    filter(!is.na(padj)) %>%
    arrange(padj)
  
  ##########################################
  # Write DE
  ##########################################
  
  addWorksheet(wb, paste0(cell, "_DE"))
  writeData(wb, paste0(cell, "_DE"), res_df)
  
  ##########################################
  # Write TPM
  ##########################################
  
  tpm_out <- tpm_sub %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_map, by = "gene_id")
  
  addWorksheet(wb, paste0(cell, "_TPM"))
  writeData(wb, paste0(cell, "_TPM"), tpm_out)
  
}

############################################
# 16. Save
############################################

saveWorkbook(wb, "DESeq2_TPM_results.xlsx", overwrite = TRUE)






############################################
# 1. Load libraries
############################################

library(openxlsx)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

############################################
# 2. Read Excel file
############################################

file <- "DESeq2_TPM_results.xlsx"
sheets <- getSheetNames(file)

de_sheets  <- sheets[grepl("_DE$", sheets)]
tpm_sheets <- sheets[grepl("_TPM$", sheets)]
de_sheets  <- de_sheets[!grepl("HIF1AKO", de_sheets)]
tpm_sheets <- tpm_sheets[!grepl("HIF1AKO", tpm_sheets)]

############################################
# 3. Collect significant features (union)
#    Criteria: padj < 0.01 and |log2FC| > 1
############################################

sig_list <- list()

for (s in de_sheets) {
  df <- read.xlsx(file, sheet = s)
  
  sig_features <- df %>%
    filter(!is.na(padj),
           !is.na(log2FoldChange),
           padj < 0.01,
           abs(log2FoldChange) > 1) %>%
    pull(gene_id)
  
  sig_list[[s]] <- sig_features
}

all_sig_features <- unique(unlist(sig_list))

############################################
# 4. Process each TPM sheet separately
#    Scale within each cell line
############################################

scaled_list <- list()
sample_info_list <- list()

for (s in tpm_sheets) {
  
  df <- read.xlsx(file, sheet = s)
  
  df_mat <- df %>%
    select(-gene_name) %>%
    column_to_rownames("gene_id")
  
  df_mat <- as.matrix(df_mat)
  mode(df_mat) <- "numeric"
  
  # keep only significant features present in this sheet
  common_features <- intersect(all_sig_features, rownames(df_mat))
  df_mat <- df_mat[common_features, , drop = FALSE]
  
  # infer cell line from sheet name
  cell <- sub("_TPM$", "", s)
  
  # get condition from metadata
  cond <- meta_aligned$Condition[
    match(colnames(df_mat), meta_aligned$Run)
  ]
  
  # log2(TPM + 1)
  log_tpm <- log2(df_mat + 1)
  
  # row-center within this cell line only
  row_means <- rowMeans(log_tpm, na.rm = TRUE)
  scaled_mat <- sweep(log_tpm, 1, row_means, "-")
  
  scaled_list[[cell]] <- scaled_mat
  
  info <- data.frame(
    Sample = colnames(scaled_mat),
    CellLine = cell,
    Condition = cond,
    stringsAsFactors = FALSE
  )
  rownames(info) <- info$Sample
  sample_info_list[[cell]] <- info
}

############################################
# 5. Merge scaled matrices across cell lines
############################################

# Collect all features
all_features <- unique(unlist(lapply(scaled_list, rownames)))

# Expand matrices to same feature space
scaled_list_full <- lapply(scaled_list, function(mat) {
  missing <- setdiff(all_features, rownames(mat))
  
  if (length(missing) > 0) {
    add_mat <- matrix(NA_real_, nrow = length(missing), ncol = ncol(mat))
    rownames(add_mat) <- missing
    colnames(add_mat) <- colnames(mat)
    mat <- rbind(mat, add_mat)
  }
  
  mat[all_features, , drop = FALSE]
})

############################################
# Merge all matrices
############################################

scaled_matrix <- do.call(cbind, scaled_list_full)

############################################
# Align annotation
############################################

annotation_df <- bind_rows(sample_info_list)

# ensure correct order matching matrix columns
annotation_df <- annotation_df[match(colnames(scaled_matrix), annotation_df$Sample), ]

############################################
# Rename samples → CellLine_condition.rep
############################################

annotation_df <- annotation_df %>%
  group_by(CellLine, Condition) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

annotation_df$new_label <- paste0(
  annotation_df$CellLine, "_",
  annotation_df$Condition, ".",
  annotation_df$rep_id
)

# apply new names
colnames(scaled_matrix) <- annotation_df$new_label

############################################
# REMOVE specific samples
############################################

remove_samples <- c("A549_1_hypoxia.1", "A549_1_normoxia.1")

keep_cols <- !colnames(scaled_matrix) %in% remove_samples

scaled_matrix <- scaled_matrix[, keep_cols, drop = FALSE]
annotation_df <- annotation_df[keep_cols, , drop = FALSE]

############################################
# Recreate condition vector (IMPORTANT)
############################################

cond <- annotation_df$Condition
names(cond) <- colnames(scaled_matrix)

############################################
# 6. Optional variance filtering
############################################

row_var <- apply(scaled_matrix, 1, var, na.rm = TRUE)
scaled_matrix <- scaled_matrix[row_var > 0.1, , drop = FALSE]

############################################
# Consistency filtering
############################################

threshold <- 11

keep_rows <- apply(scaled_matrix, 1, function(x) {
  
  hyp <- x[cond == "hypoxia"]
  nor <- x[cond == "normoxia"]
  
  hyp_up   <- sum(hyp > 0, na.rm = TRUE)
  hyp_down <- sum(hyp < 0, na.rm = TRUE)
  
  nor_up   <- sum(nor > 0, na.rm = TRUE)
  nor_down <- sum(nor < 0, na.rm = TRUE)
  
  (hyp_up >= threshold) |
    (hyp_down >= threshold) |
    (nor_up >= threshold) |
    (nor_down >= threshold)
})

scaled_matrix <- scaled_matrix[keep_rows, , drop = FALSE]

############################################
# Remove NA / Inf rows BEFORE clustering
############################################

good_rows <- apply(scaled_matrix, 1, function(x) all(is.finite(x)))
scaled_matrix <- scaled_matrix[good_rows, , drop = FALSE]

############################################
# 7. Unbiased clustering
############################################

library(dendextend)

row_hclust <- hclust(dist(scaled_matrix))
row_dend <- as.dendrogram(row_hclust)
row_dend <- rev(row_dend)





col_hclust <- hclust(dist(t(scaled_matrix)))
col_dend <- as.dendrogram(col_hclust)
desired_order <- colnames(scaled_matrix)[
  order(annotation_df$Condition)  # normoxia first, hypoxia second
]
library(dendextend)
col_dend <- rotate(col_dend, rev(desired_order))

############################################
# Optional: 4 row clusters
############################################

row_groups <- cutree(row_hclust, k = 2)

row_order <- order.dendrogram(as.dendrogram(row_hclust))
ordered_groups <- row_groups[row_order]
group_levels <- unique(ordered_groups)
group_map <- setNames(seq_along(group_levels), group_levels)

row_split <- factor(group_map[as.character(row_groups)], levels = 1:2)

############################################
# 8. Annotations
############################################

############################################
# 8. Annotations
############################################

cond_cols <- c(
  "hypoxia" = "orangered",
  "normoxia" = "turquoise3"
)

cell_cols <- c(
  "A549_1" = "plum2",
  "A549_2" = "plum2",
  "HeLa"   = "grey70",
  "H460"   = "cadetblue3",
  "RKO"    = "sandybrown"
)

ha_cond <- HeatmapAnnotation(
  Condition = annotation_df$Condition,
  col = list(Condition = cond_cols),
  show_annotation_name = TRUE,
  annotation_height = unit(3, "mm")
)

top_ha <- ha_cond

ha_violin <- HeatmapAnnotation(
  violin = anno_density(
    scaled_matrix,
    which = "column",
    type = "violin",
    gp = gpar(fill = NA, col = "black"),
    axis_param = list(side = "left")
  ),
  annotation_name_side = "right",
  annotation_name_rot = 0,
  annotation_height = unit(3.8, "cm")
)

top_ha <- c(ha_violin, ha_top)
ha_top -> top_ha

############################################
# 9. Color function
############################################

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("gold", "white", "blueviolet")
)

# Create gene_id → gene_name mapping
gene_map_vec <- setNames(gene_map$gene_name, gene_map$gene_id)

# Replace rownames ONLY for visualization
gene_labels <- gene_map_vec[rownames(scaled_matrix)]

# fallback to gene_id if missing
gene_labels[is.na(gene_labels)] <- rownames(scaled_matrix)[is.na(gene_labels)]

# make unique (avoid duplicate gene names crashing heatmap)
gene_labels <- make.unique(gene_labels)

rownames(scaled_matrix) <- gene_labels





############################################
# 10. Heatmap
############################################

set.seed(6)
ht <- Heatmap(
  scaled_matrix,
  name = "scaled exp",
  col = col_fun,
  na_col = "grey95",
  use_raster = FALSE,
  top_annotation = top_ha,
  
  #cluster_rows = row_dend,
  cluster_columns = col_dend,
  
  show_row_names = FALSE,
  column_names_rot = 90,
  
  km = 4,   # ✅ FIX
  row_title = NULL,
  
  column_title = NULL,
  
  heatmap_legend_param = list(
    title = "scaled exp",
    at = c(-2, -1, 0, 1, 2)
  )
)


############################################
# Draw
############################################

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)




library(ComplexHeatmap)
library(grid)

############################################
# 1. K-means clustering
############################################

set.seed(6)
km <- kmeans(scaled_matrix, centers = 4)

row_split <- factor(km$cluster, levels = c(2,3,1,4))

############################################
# 2. Heatmap (NO rect_gp here)
############################################
 
ht <- Heatmap(
  scaled_matrix,
  name = "scaled exp",
  col = col_fun,
  na_col = "grey95",
  use_raster = FALSE,
  top_annotation = top_ha,
  
  cluster_columns = col_dend,
  cluster_rows = TRUE,
  
  show_row_names = FALSE,
  column_names_rot = 90,
  
  row_split = row_split,
  row_title = NULL,
  column_title = NULL,
  
  heatmap_legend_param = list(
    title = "scaled exp",
    at = c(-2, -1, 0, 1, 2)
  )
)

############################################
# 3. Draw heatmap
############################################

ht_drawn <- draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

############################################
# 4. Add ONLY vertical lines
############################################

for (s in seq_len(length(unique(row_split)))) {
  
  decorate_heatmap_body("scaled exp", {
    
    n <- ncol(scaled_matrix)
    
    for (i in seq_len(n - 1)) {
      
      x <- unit(i / n, "npc")
      
      grid.lines(
        x = unit.c(x, x),
        y = unit.c(unit(0, "npc"), unit(1, "npc")),
        gp = gpar(col = "black", lwd = 0.5)
      )
    }
    
  }, slice = s)
}

library(grid)

png("heatmap_with_vertical_lines.png", width = 2800, height = 3600, res = 300)

# 1. draw heatmap
ht_drawn <- draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# 2. add vertical lines (ALL slices)
n_slices <- length(levels(row_split))

for (s in seq_len(n_slices)) {
  
  decorate_heatmap_body("scaled exp", {
    
    n <- ncol(scaled_matrix)
    
    for (i in seq_len(n - 1)) {
      
      x <- unit(i / n, "npc")
      
      grid.lines(
        x = unit.c(x, x),
        y = unit.c(unit(0, "npc"), unit(1, "npc")),
        gp = gpar(col = "black", lwd = 0.5)
      )
    }
    
  }, slice = s)
}

dev.off()
