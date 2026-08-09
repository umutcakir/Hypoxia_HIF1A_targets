set.seed(0)
library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)

## ---- Load Seurat object ----
sc <- readRDS("92_Schwarz_HS_Seurat_with_TEs.rds")

## ---- Remove multiplets and undetermined (chained correctly) ----
sc_singlets <- subset(sc, subset = Sample_Tag != "Multiplet")
sc_singlets <- subset(sc_singlets, subset = Sample_Tag != "Undetermined")

## ---- Read sample annotation from Excel ----
raw <- read_excel("Sample annotation.xlsx", sheet = "Sheet1", col_names = FALSE)
colnames(raw) <- c("Sample_Tag", "Field", "Value")

sample_annotation <- raw %>%
  fill(Sample_Tag, .direction = "down") %>%
  mutate(Sample_Tag = as.integer(Sample_Tag)) %>%
  pivot_wider(names_from = Field, values_from = Value)


print(sample_annotation)

## ---- Assign Donor and Treatment to Seurat metadata ----
tag_num <- as.integer(gsub("SampleTag(\\d+)_hs", "\\1", as.character(sc_singlets$Sample_Tag)))

sc_singlets$Donor     <- sample_annotation$Donor[match(tag_num, sample_annotation$Sample_Tag)]
sc_singlets$Treatment <- sample_annotation$Treatment[match(tag_num, sample_annotation$Sample_Tag)]

## ---- Sanity checks ----
table(sc_singlets$Sample_Tag, sc_singlets$Donor)
table(sc_singlets$Sample_Tag, sc_singlets$Treatment)

## ---- Plots ----
DimPlot(sc_singlets, group.by = "Cell_Type_Experimental")
DimPlot(sc_singlets, group.by = "Donor")
DimPlot(sc_singlets, group.by = "Treatment")


## ===============================
## Seurat + Harmony pipeline
## BD Rhapsody – singlets only
## ===============================

library(Seurat)
library(harmony)
library(ggplot2)

## -------------------------------
## 1. Set identities (optional)
## -------------------------------
Idents(sc_singlets) <- sc_singlets$Cell_Type_Experimental
sc_singlets[["percent.mt"]] <- PercentageFeatureSet(sc_singlets, pattern = "^MT-")

VlnPlot(sc_singlets, features = "percent.mt", group.by = "Sample_Tag", pt.size = 0)

cat("Cells before filtering:", ncol(sc_singlets), "\n")
# Apply filter (adjust threshold based on what the distribution shows)
sc_singlets <- subset(sc_singlets, subset = percent.mt < 25)
cat("Cells after filtering:", ncol(sc_singlets), "\n")



## -------------------------------
## 2. RNA preprocessing
## -------------------------------
sc_singlets <- NormalizeData(
  sc_singlets,
  normalization.method = "LogNormalize"
)

sc_singlets <- FindVariableFeatures(
  sc_singlets,
  selection.method = "vst"
)

sc_singlets <- ScaleData(
  sc_singlets,
  features = VariableFeatures(sc_singlets)
)

## -------------------------------
## 3. PCA
## -------------------------------
sc_singlets <- RunPCA(
  sc_singlets,
  features = VariableFeatures(sc_singlets),
  npcs = 50
)

ElbowPlot(sc_singlets, ndims = 50)

## -------------------------------
## 4. Harmony integration
##    (batch = Sample_Tag)
## -------------------------------
sc_singlets <- RunHarmony(
  object = sc_singlets,
  group.by.vars = "Sample_Tag",
  reduction.use = "pca",
  dims.use = 1:15,
  assay.use = "RNA"
)

## -------------------------------
## 5. UMAP on Harmony
## -------------------------------
sc_singlets <- RunUMAP(
  sc_singlets,
  reduction = "harmony",
  dims = 1:15
)

## -------------------------------
## 6. Neighbors & clustering
## -------------------------------
sc_singlets <- FindNeighbors(
  sc_singlets,
  reduction = "harmony",
  dims = 1:15
)

sc_singlets <- FindClusters(
  sc_singlets,
  resolution = 0.6
)

## -------------------------------
## 7. Visualization
## -------------------------------
DimPlot(
  sc_singlets,
  reduction = "umap",
  group.by = "Sample_Tag"
)

DimPlot(
  sc_singlets,
  reduction = "umap",
  group.by = "Cell_Type_Experimental",
  label = TRUE
)

DimPlot(
  sc_singlets,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)



saveRDS(sc_singlets, file = "sc_singlets.RDS")



set.seed(0)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---- Output directory ----
outdir <- "cluster_id_figures"
dir.create(outdir, showWarnings = FALSE)

## =========================================================
## 1. ADT normalization (CLR) + average expression per cluster
## =========================================================
DefaultAssay(sc_singlets) <- "ADT"
sc_singlets <- NormalizeData(
  sc_singlets,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2
)

avg_adt <- AverageExpression(sc_singlets, assay = "ADT", group.by = "seurat_clusters")$ADT
cat("\n==== Average ADT expression per cluster ====\n")
print(round(avg_adt, 2))
write.csv(round(avg_adt, 2), file.path(outdir, "avg_adt_per_cluster.csv"))

## =========================================================
## 2. RNA marker panel
## =========================================================
DefaultAssay(sc_singlets) <- "RNA"

rna_markers <- c(
  "CD8A", "CD8B",
  "TRDC", "TRGC1", "TRGC2",
  "FOXP3", "IL2RA",
  "ITGAX", "CLEC9A", "FCER1A", "CD1C",
  "LYZ", "FCGR3A",
  "MKI67", "TOP2A", "PCNA"
)

p_dotplot <- DotPlot(sc_singlets, features = rna_markers, group.by = "seurat_clusters") +
  RotatedAxis()
ggsave(file.path(outdir, "rna_marker_dotplot.png"), p_dotplot, width = 14, height = 8, dpi = 300)

## =========================================================
## 3. QC per cluster
## =========================================================
p_qc <- VlnPlot(
  sc_singlets,
  features = c("nCount_RNA", "nFeature_RNA", "nCount_ADT"),
  group.by = "seurat_clusters",
  pt.size = 0
)
ggsave(file.path(outdir, "qc_violin_per_cluster.png"), p_qc, width = 16, height = 6, dpi = 300)

## =========================================================
## 4. ADT FeaturePlots
## =========================================================
adt_markers <- rownames(sc_singlets[["ADT"]])
DefaultAssay(sc_singlets) <- "ADT"
p_adt_feat <- FeaturePlot(
  sc_singlets,
  features = adt_markers,
  reduction = "umap",
  min.cutoff = "q05", max.cutoff = "q95",
  ncol = 4
)
ggsave(file.path(outdir, "adt_featureplots.png"), p_adt_feat, width = 16, height = 8, dpi = 300)
DefaultAssay(sc_singlets) <- "RNA"

## =========================================================
## 5. Cluster UMAP + sizes
## =========================================================
p_clusters <- DimPlot(sc_singlets, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "umap_seurat_clusters.png"), p_clusters, width = 8, height = 7, dpi = 300)

cat("\n==== Cluster sizes ====\n")
print(table(sc_singlets$seurat_clusters))
write.csv(as.data.frame(table(sc_singlets$seurat_clusters)),
          file.path(outdir, "cluster_sizes.csv"), row.names = FALSE)

cat("\nAll figures and tables saved to:", normalizePath(outdir), "\n")


sc_singlets_clean <- subset(sc_singlets, subset = seurat_clusters != "5")

# Re-cluster using the SAME existing harmony embedding (already computed, not recalculated)
sc_singlets_clean <- FindNeighbors(sc_singlets_clean, reduction = "harmony", dims = 1:15)
sc_singlets_clean <- FindClusters(sc_singlets_clean, resolution = 0.6)

# UMAP coordinates are untouched — same embedding as before, just new cluster labels/colors
DimPlot(sc_singlets_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

table(sc_singlets_clean$seurat_clusters)

# Optional: keep the old UMAP layout under a different name before overwriting
sc_singlets_clean[["umap_old"]] <- sc_singlets_clean[["umap"]]

ElbowPlot(sc_singlets_clean, ndims = 50)
set.seed(0)
sc_singlets_clean <- RunUMAP(sc_singlets_clean, reduction = "harmony", dims = 1:15)
DimPlot(sc_singlets_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# You can always go back to comparing both:
DimPlot(sc_singlets_clean, reduction = "umap_old", group.by = "seurat_clusters", label = TRUE)

DimPlot(sc_singlets_clean, reduction = "umap", group.by = "Cell_Type_Experimental", label = TRUE)

DimPlot(sc_singlets_clean, reduction = "umap", group.by = "Sample_Tag", label = TRUE)



set.seed(0)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

outdir <- "cluster_id_figures_v2"
dir.create(outdir, showWarnings = FALSE)

## =========================================================
## 0. QC violin (sanity check, mt% and complexity already filtered)
## =========================================================
p_qc <- VlnPlot(sc_singlets_clean, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                group.by = "seurat_clusters", pt.size = 0)
ggsave(file.path(outdir, "qc_violin_per_cluster.png"), p_qc, width = 16, height = 6, dpi = 300)

## =========================================================
## 1. ADT normalization (CLR) + average expression per cluster
## =========================================================
DefaultAssay(sc_singlets_clean) <- "ADT"
sc_singlets_clean <- NormalizeData(sc_singlets_clean, assay = "ADT", normalization.method = "CLR", margin = 2)

avg_adt <- AverageExpression(sc_singlets_clean, assay = "ADT", group.by = "seurat_clusters")$ADT
write.csv(round(avg_adt, 2), file.path(outdir, "avg_adt_per_cluster.csv"))

## =========================================================
## 2. RNA marker panel
## =========================================================
DefaultAssay(sc_singlets_clean) <- "RNA"
rna_markers <- c(
  "CD8A", "CD8B", "TRDC", "TRGC1", "TRGC2",
  "FOXP3", "IL2RA",
  "ITGAX", "CLEC9A", "FCER1A", "CD1C",
  "LYZ", "FCGR3A",
  "MKI67", "TOP2A", "PCNA"
)
p_dotplot <- DotPlot(sc_singlets_clean, features = rna_markers, group.by = "seurat_clusters") + RotatedAxis()
ggsave(file.path(outdir, "rna_marker_dotplot.png"), p_dotplot, width = 14, height = 8, dpi = 300)

## =========================================================
## 3. ADT FeaturePlots
## =========================================================
adt_markers <- rownames(sc_singlets_clean[["ADT"]])
DefaultAssay(sc_singlets_clean) <- "ADT"
p_adt_feat <- FeaturePlot(sc_singlets_clean, features = adt_markers, reduction = "umap",
                          min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
ggsave(file.path(outdir, "adt_featureplots.png"), p_adt_feat, width = 16, height = 8, dpi = 300)
DefaultAssay(sc_singlets_clean) <- "RNA"

## =========================================================
## 4. Cluster UMAP + sizes
## =========================================================
p_clusters <- DimPlot(sc_singlets_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
ggsave(file.path(outdir, "umap_seurat_clusters.png"), p_clusters, width = 8, height = 7, dpi = 300)

write.csv(as.data.frame(table(sc_singlets_clean$seurat_clusters)),
          file.path(outdir, "cluster_sizes.csv"), row.names = FALSE)

cat("\nAll figures and tables saved to:", normalizePath(outdir), "\n")


cluster_annotation <- c(
  "0"  = "T_CD4",
  "1"  = "T_CD4",
  "2"  = "NK_CD8_gd_mixed",
  "3"  = "Monocyte_classical",
  "4"  = "T_CD4",
  "5"  = "T_CD4",
  "6"  = "B",
  "7"  = "Monocyte_classical",
  "8"  = "T_CD8_gd",
  "9"  = "NK_CD16bright",
  "10" = "T_CD4",
  "11" = "Dendritic",
  "12" = "Proliferating_lymphocyte",
  "13" = "Monocyte_nonclassical",
  "14" = "Rare_activated",
  "15" = "Monocyte_intermediate"
)

sc_singlets_clean$Cell_Type_Final <- unname(cluster_annotation[as.character(sc_singlets_clean$seurat_clusters)])

DimPlot(sc_singlets_clean, reduction = "umap", group.by = "Cell_Type_Final", label = TRUE, repel = TRUE) + NoLegend()
table(sc_singlets_clean$Cell_Type_Final)


cluster_annotation <- c(
  "0"  = "T_CD4_naive",
  "1"  = "T_CD4_naive",
  "2"  = "NK_mature_cytotoxic",
  "3"  = "Monocyte_classical",
  "4"  = "T_CD4_memory",
  "5"  = "T_ISG",
  "6"  = "B_naive",
  "7"  = "Monocyte_classical",
  "8"  = "MAIT",
  "9"  = "Neutrophil",
  "10" = "T_ISG",
  "11" = "Basophil",
  "12" = "Proliferating",
  "13" = "Neutrophil",
  "14" = "pDC",
  "15" = "Monocyte_classical"
)

sc_singlets_clean$Cell_Type_Final <- unname(cluster_annotation[as.character(sc_singlets_clean$seurat_clusters)])
DimPlot(sc_singlets_clean, reduction = "umap", group.by = "Cell_Type_Final", label = TRUE, repel = TRUE) 
table(sc_singlets_clean$Cell_Type_Final)


cluster_annotation_display_v2 <- c(
  "0"  = "CD4+ T cells",
  "1"  = "CD4+ T cells",
  "2"  = "NK cells",
  "3"  = "Monocytes",
  "4"  = "CD4+ T cells",
  "5"  = "CD4+ T cells",
  "6"  = "B cells",
  "7"  = "Monocytes",
  "8"  = "MAIT cells",
  "9"  = "Neutrophils",
  "10" = "CD4+ T cells",
  "11" = "Basophils",
  "12" = "Proliferating cells",
  "13" = "Neutrophils",
  "14" = "pDC",
  "15" = "Monocytes"
)

sc_singlets_clean$Cell_Type_Display <- unname(cluster_annotation_display_v2[as.character(sc_singlets_clean$seurat_clusters)])

celltype_colors_v2 <- c(
  "CD4+ T cells"        = "#E69F00",
  "NK cells"            = "#009E73",
  "MAIT cells"          = "#F0E442",
  "Monocytes"           = "#0072B2",
  "B cells"             = "#D55E00",
  "Neutrophils"         = "#CC79A7",
  "Basophils"           = "#999999",
  "pDC"                 = "#882255",
  "Proliferating cells" = "#8B4513"
)

p <- DimPlot(
  sc_singlets_clean,
  reduction = "umap",
  group.by = "Cell_Type_Display",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5,
  pt.size = 0.3
) +
  scale_color_manual(values = celltype_colors_v2) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  coord_fixed()

print(p)
ggsave("umap_celltype_final.png", p, width = 9, height = 7, dpi = 600)
saveRDS(sc_singlets_clean, file = "sc_singlets_clean_final_named.RDS")



sc_singlets_clean <- readRDS("sc_singlets_clean_final_named.RDS")


 


set.seed(0)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

outdir <- "DE_results_bimod_001pct"
dir.create(outdir, showWarnings = FALSE)

DefaultAssay(sc_singlets_clean) <- "RNA"

contrasts_to_run <- list(
  c("HOX", "NOX"),
  c("ROX", "NOX"),
  c("IFN", "NOX"),
  c("ROX", "HOX")
)

cell_types <- unique(sc_singlets_clean$Cell_Type_Display)
de_results <- list()

## minimum cells per group required to attempt a comparison
min_cells <- 10

"""
#pDC:
== NK cells == cell counts per treatment:
  
  HOX  IFN  NOX  ROX 
1266 1096 1256 1711 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03m 13s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03m 55s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 59s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03m 30s

== Monocytes == cell counts per treatment:
  
  HOX  IFN  NOX  ROX 
1277 1537 1468 1795 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=09m 39s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=11m 17s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10m 47s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10m 02s

== CD4+ T cells == cell counts per treatment:
  
  HOX  IFN  NOX  ROX 
4298 4336 4587 6229 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10m 28s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=12m 39s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10m 45s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13m 19s

== B cells == cell counts per treatment:
  
  HOX IFN NOX ROX 
282 335 325 425 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 24s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 42s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 27s
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 33s

== Neutrophils == cell counts per treatment:
  
  HOX IFN NOX ROX 
191 114  82 224 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=30s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=23s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=34s  

== Basophils == cell counts per treatment:
  
  HOX IFN NOX ROX 
40  39  48  62 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=16s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=20s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=15s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=18s  

== MAIT cells == cell counts per treatment:
  
  HOX IFN NOX ROX 
265 210 282 343 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=47s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=52s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=46s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=47s  

== pDC == cell counts per treatment:
  
  HOX IFN NOX ROX 
9  25  10  13 
Skipping pDC HOX vs NOX - insufficient cells (need >= 10 per group)
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=11s  
Skipping pDC ROX vs HOX - insufficient cells (need >= 10 per group)

== Proliferating cells == cell counts per treatment:
  
  HOX IFN NOX ROX 
18  28  24  14 
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=25s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=34s  
|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=22s  
"""
set.seed(0)
for (ct in cell_types) {
  
  sub_obj <- subset(sc_singlets_clean, subset = Cell_Type_Display == ct)
  Idents(sub_obj) <- sub_obj$Treatment
  
  group_sizes <- table(sub_obj$Treatment)
  cat("\n==", ct, "== cell counts per treatment:\n")
  print(group_sizes)
  
  for (contrast in contrasts_to_run) {
    trt1 <- contrast[1]; trt2 <- contrast[2]
    
    if (!(trt1 %in% names(group_sizes)) || !(trt2 %in% names(group_sizes))) next
    if (group_sizes[trt1] < min_cells || group_sizes[trt2] < min_cells) {
      cat("Skipping", ct, trt1, "vs", trt2, "- insufficient cells (need >=", min_cells, "per group)\n")
      next
    }
    
    markers <- tryCatch({
      FindMarkers(
        sub_obj,
        ident.1 = trt1,
        ident.2 = trt2,
        test.use = "bimod",
        logfc.threshold = 0,
        min.pct = 0.01
      )
    }, error = function(e) {
      cat("FindMarkers failed for", ct, trt1, "vs", trt2, ":", conditionMessage(e), "\n")
      NULL
    })
    if (is.null(markers) || nrow(markers) == 0) next
    
    markers_df <- markers %>%
      mutate(
        gene = rownames(markers),
        feature_type = ifelse(grepl("^TE-", gene), "TE", "gene"),
        cell_type = ct,
        comparison = paste0(trt1, "_vs_", trt2)
      ) %>%
      relocate(gene, feature_type, cell_type, comparison)
    
    de_results[[paste(ct, trt1, trt2, sep = "_")]] <- markers_df
    
    safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
    write.csv(markers_df,
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_all.csv")),
              row.names = FALSE)
    write.csv(markers_df %>% filter(feature_type == "gene"),
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_genes.csv")),
              row.names = FALSE)
    write.csv(markers_df %>% filter(feature_type == "TE"),
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_TEs.csv")),
              row.names = FALSE)
    
    ## volcano plot (bimod output uses avg_log2FC and p_val_adj)
    p_volcano <- markers_df %>%
      filter(!is.na(p_val_adj)) %>%
      ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = feature_type)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c(gene = "steelblue", TE = "firebrick")) +
      labs(title = paste0(ct, ": ", trt1, " vs ", trt2),
           x = "avg log2FC", y = "-log10(p_val_adj)") +
      theme_classic(base_size = 12)
    
    ggsave(file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_volcano.png")),
           p_volcano, width = 7, height = 6, dpi = 300)
  }
}

## =========================================================
## Combine + summary
## =========================================================
de_results_all <- bind_rows(de_results)
write.csv(de_results_all, file.path(outdir, "all_celltypes_all_comparisons_bimod.csv"), row.names = FALSE)

summary_table <- de_results_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0) %>%
  count(cell_type, comparison, feature_type) %>%
  pivot_wider(names_from = feature_type, values_from = n, values_fill = 0)

write.csv(summary_table, file.path(outdir, "significant_hits_summary_bimod.csv"), row.names = FALSE)
print(summary_table)


de_results_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0) -> de_results_significant

write.csv(de_results_significant, file.path(outdir, "de_results_significant.csv"), row.names = FALSE)


cat("\nAll bimod DE results, summary table, and volcano plots saved to:", normalizePath(outdir), "\n")







set.seed(0)

library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

outdir <- "DE_results_pseudobulk"
dir.create(outdir, showWarnings = FALSE)

DefaultAssay(sc_singlets_clean) <- "RNA"

## =========================================================
## 1. Build a clean grouping ID (avoids string-parsing issues
##    with spaces/"+" in cell type names like "CD4+ T cells")
## =========================================================
sc_singlets_clean$pseudobulk_id <- paste0(
  "S", as.integer(factor(paste(
    sc_singlets_clean$Cell_Type_Display,
    sc_singlets_clean$Donor,
    sc_singlets_clean$Treatment,
    sep = "___"
  )))
)

lookup <- sc_singlets_clean@meta.data %>%
  distinct(pseudobulk_id, Cell_Type_Display, Donor, Treatment)

cell_counts <- sc_singlets_clean@meta.data %>%
  dplyr::count(pseudobulk_id, name = "n_cells")
lookup <- left_join(lookup, cell_counts, by = "pseudobulk_id")

min_cells <- 10
lookup_valid <- lookup %>% filter(n_cells >= min_cells)

cat("Pseudobulk groups before filtering:", nrow(lookup), "\n")
cat("Pseudobulk groups after filtering (>=", min_cells, "cells):", nrow(lookup_valid), "\n")
write.csv(lookup, file.path(outdir, "pseudobulk_group_summary.csv"), row.names = FALSE)

## =========================================================
## 2. Aggregate counts
## =========================================================
pseudobulk <- AggregateExpression(
  sc_singlets_clean,
  assays = "RNA",
  group.by = "pseudobulk_id",
  return.seurat = FALSE
)$RNA

colnames(pseudobulk) <- gsub("^.*\\.(S[0-9]+)$", "\\1", colnames(pseudobulk))
pseudobulk <- pseudobulk[, colnames(pseudobulk) %in% lookup_valid$pseudobulk_id, drop = FALSE]

## =========================================================
## 3. Loop over cell types: run DESeq2, extract contrasts
## =========================================================
cell_types <- unique(lookup_valid$Cell_Type_Display)
de_results <- list()

contrasts_to_run <- list(
  c("HOX", "NOX"),
  c("ROX", "NOX"),
  c("IFN", "NOX"),
  c("ROX", "HOX")
)

for (ct in cell_types) {
  
  meta_sub <- lookup_valid %>% filter(Cell_Type_Display == ct)
  if (nrow(meta_sub) < 4) {
    cat("Skipping", ct, "- too few pseudobulk samples overall\n")
    next
  }
  
  cols_this_ct <- meta_sub$pseudobulk_id
  counts_sub <- pseudobulk[, cols_this_ct, drop = FALSE]
  
  coldata <- meta_sub
  rownames(coldata) <- coldata$pseudobulk_id
  coldata <- coldata[colnames(counts_sub), ]
  coldata$Treatment <- factor(coldata$Treatment, levels = c("NOX", "HOX", "ROX", "IFN"))
  coldata$Donor <- factor(coldata$Donor)
  
  if (any(table(coldata$Treatment) < 2)) {
    cat("Skipping", ct, "- insufficient donor replication per treatment\n")
    next
  }
  
  dds <- tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = round(counts_sub), colData = coldata, design = ~ Donor + Treatment)
    DESeq(dds)
  }, error = function(e) {
    cat("DESeq2 failed for", ct, ":", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(dds)) next
  
  for (contrast in contrasts_to_run) {
    trt1 <- contrast[1]; trt2 <- contrast[2]
    if (!(trt1 %in% coldata$Treatment) || !(trt2 %in% coldata$Treatment)) next
    
    res <- results(dds, contrast = c("Treatment", trt1, trt2))
    res_df <- as.data.frame(res) %>%
      mutate(
        gene = rownames(res),
        feature_type = ifelse(grepl("^TE-", gene), "TE", "gene"),
        cell_type = ct,
        comparison = paste0(trt1, "_vs_", trt2)
      ) %>%
      relocate(gene, feature_type, cell_type, comparison)
    
    de_results[[paste(ct, trt1, trt2, sep = "_")]] <- res_df
    
    safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
    write.csv(res_df,
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_all.csv")),
              row.names = FALSE)
    write.csv(res_df %>% filter(feature_type == "gene"),
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_genes.csv")),
              row.names = FALSE)
    write.csv(res_df %>% filter(feature_type == "TE"),
              file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_TEs.csv")),
              row.names = FALSE)
    
    p_volcano <- res_df %>%
      filter(!is.na(padj)) %>%
      ggplot(aes(x = log2FoldChange, y = -log10(padj), color = feature_type)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c(gene = "steelblue", TE = "firebrick")) +
      labs(title = paste0(ct, ": ", trt1, " vs ", trt2),
           x = "log2 Fold Change", y = "-log10(padj)") +
      theme_classic(base_size = 12)
    
    ggsave(file.path(outdir, paste0(safe_ct, "_", trt1, "_vs_", trt2, "_volcano.png")),
           p_volcano, width = 7, height = 6, dpi = 300)
  }
}

## =========================================================
## 4. Combine everything + summary table
## =========================================================
de_results_all <- bind_rows(de_results)
write.csv(de_results_all, file.path(outdir, "all_celltypes_all_comparisons_pseudobulk.csv"), row.names = FALSE)

summary_table <- de_results_all %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0) %>%
  dplyr::count(cell_type, comparison, feature_type) %>%
  pivot_wider(names_from = feature_type, values_from = n, values_fill = 0)

write.csv(summary_table, file.path(outdir, "significant_hits_summary_pseudobulk.csv"), row.names = FALSE)
print(summary_table)

cat("\nAll pseudobulk DE results, summary table, and volcano plots saved to:", normalizePath(outdir), "\n")








 




 


set.seed(0)
library(Seurat)
library(dplyr)

DefaultAssay(sc_singlets_clean) <- "RNA"

outdir <- "DE_HOX_vs_NOX_by_celltype"
dir.create(outdir, showWarnings = FALSE)

cell_types <- unique(sc_singlets_clean$Cell_Type_Display)

de_results <- list()

for (ct in cell_types) {
  
  cat("\n==== Processing:", ct, "====\n")
  
  # Subset to this cell type AND only HOX/NOX cells
  subset_obj <- subset(
    sc_singlets_clean,
    subset = Cell_Type_Display == ct & Treatment %in% c("HOX", "NOX")
  )
  
  # Check there are enough cells in each group before running
  treat_counts <- table(subset_obj$Treatment)
  cat("Cell counts:\n")
  print(treat_counts)
  
  if (length(treat_counts) < 2 || any(treat_counts < 10)) {
    cat("Skipping", ct, "- insufficient cells in one or both groups\n")
    next
  }
  
  Idents(subset_obj) <- "Treatment"
  
  markers <- tryCatch({
    FindMarkers(
      subset_obj,
      ident.1 = "HOX",
      ident.2 = "NOX",
      test.use = "bimod",
      min.pct = 0.1,
      logfc.threshold = 0.25
    )
  }, error = function(e) {
    cat("Error for", ct, ":", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(markers)) next
  
  markers$gene <- rownames(markers)
  markers <- markers[order(markers$p_val_adj), ]
  
  de_results[[ct]] <- markers
  
  # Save per-cell-type CSV
  safe_name <- gsub("[^A-Za-z0-9]", "_", ct)
  write.csv(markers, file.path(outdir, paste0("DE_", safe_name, "_HOXvsNOX.csv")), row.names = FALSE)
}

cat("\nAll results saved to:", normalizePath(outdir), "\n")

# Quick summary: number of significant DE genes per cell type
summary_table <- sapply(de_results, function(x) sum(x$p_val_adj < 0.05))
print(summary_table)
write.csv(data.frame(CellType = names(summary_table), n_sig_genes = summary_table),
          file.path(outdir, "summary_n_sig_genes.csv"), row.names = FALSE)


 






set.seed(0)
library(Seurat)
library(dplyr)

DefaultAssay(sc_singlets_clean) <- "RNA"

outdir <- "cluster_markers_final"
dir.create(outdir, showWarnings = FALSE)

## Run FindAllMarkers: one-vs-rest for every cluster
all_markers <- FindAllMarkers(
  sc_singlets_clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

# Save the full table
write.csv(all_markers, file.path(outdir, "all_cluster_markers.csv"), row.names = FALSE)

# Save top 15 markers per cluster (ranked by avg_log2FC) for quick scanning
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15)
write.csv(top_markers, file.path(outdir, "top15_markers_per_cluster.csv"), row.names = FALSE)

print(top_markers, n = 300)
cat("\nSaved to:", normalizePath(outdir), "\n")






## 1. What genes actually define cluster 5 vs everything else?
cluster5_markers <- FindMarkers(
  sc_singlets,
  ident.1 = "5",
  group.by = "seurat_clusters",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)
head(cluster5_markers[order(-cluster5_markers$avg_log2FC), ], 20)

## 2. Check T-lineage RNA markers (in case ADT is just internalized)
FeaturePlot(sc_singlets, features = c("CD3D", "CD3E", "CD2", "IL7R", "CCR7", "SELL"), reduction = "umap")

## 3. Check for non-immune contaminants (RBC / platelet)
FeaturePlot(sc_singlets, features = c("HBB", "HBA1", "HBA2", "PPBP", "PF4"), reduction = "umap")

## 4. Mitochondrial content (true QC check for cell death/degradation)
sc_singlets[["percent.mt"]] <- PercentageFeatureSet(sc_singlets, pattern = "^MT-")
VlnPlot(sc_singlets, features = "percent.mt", group.by = "seurat_clusters")

## 5. Ribosomal gene content (another QC/stress signal)
sc_singlets[["percent.ribo"]] <- PercentageFeatureSet(sc_singlets, pattern = "^RPS|^RPL")
VlnPlot(sc_singlets, features = "percent.ribo", group.by = "seurat_clusters")







## ---- Install/load DoubletFinder ----
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}
library(DoubletFinder)
library(Seurat)
library(dplyr)

## ---- 1. pK identification (no ground-truth doublets, so use paramSweep) ----
# Use RNA assay, PCA reduction already computed (harmony-corrected PCs if you used those for clustering)
sweep_res <- paramSweep(sc_singlets, PCs = 1:15, sct = FALSE)  # adjust PCs to match what you used for clustering
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
bcmvn <- find.pK(sweep_stats)

# pick the pK with the max BCmetric (mean-variance normalized bimodality coefficient)
pK_optimal <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  pull(pK) %>%
  as.numeric() %>%
  as.character() %>%
  as.numeric()

print(pK_optimal)

## ---- 2. Estimate expected doublet rate ----
# Homotypic doublet proportion based on your cluster sizes
annotations <- sc_singlets$seurat_clusters
homotypic_prop <- modelHomotypic(annotations)

# Expected doublet rate: BD Rhapsody/10x-style estimate ~ 0.8% per 1000 cells loaded is a common rule of thumb,
# but adjust this based on your actual loading concentration if known.
n_cells <- ncol(sc_singlets)
doublet_rate <- 0.05   # fixed ~5% as a more realistic starting estimate; adjust if you have the actual loading-based rate
nExp_poi <- round(doublet_rate * n_cells)
nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))  # ADJUST based on your actual cell loading target

nExp_poi <- round(doublet_rate * n_cells)
nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

## ---- 3. Run DoubletFinder ----
sc_singlets <- doubletFinder(
  sc_singlets,
  PCs = 1:15,
  pN = 0.25,
  pK = pK_optimal,
  nExp = nExp_poi_adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

# the new metadata column name will look like "DF.classifications_0.25_<pK>_<nExp>"
doublet_col <- grep("DF.classifications", colnames(sc_singlets@meta.data), value = TRUE)
print(doublet_col)

## ---- 4. Check doublet calls against cluster 12 ----
table(sc_singlets$seurat_clusters, sc_singlets@meta.data[[doublet_col]])

DimPlot(sc_singlets, reduction = "umap", group.by = doublet_col)
 










## -------------------------------
## 8. Set sample identity for DE
## -------------------------------
Idents(sc_singlets) <- "Sample_Name"












## ============================================
## Automatic immune cell annotation with SingleR
## ============================================

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

## --------------------------------------------
## 1. Convert Seurat -> SingleCellExperiment
## --------------------------------------------
sce <- as.SingleCellExperiment(sc_singlets)

ref_hpca <- HumanPrimaryCellAtlasData()

pred_hpca <- SingleR(
  test = sce,
  ref = ref_hpca,
  labels = ref_hpca$label.main
)

sc_singlets$HPCA_label <- pred_hpca$labels
sc_singlets$HPCA_delta <- pred_hpca$delta.next

## --------------------------------------------
## 5. Visualization
## --------------------------------------------
DimPlot(
  sc_singlets,
  reduction = "umap",
  group.by = "HPCA_label",
  label = TRUE,
  repel = TRUE
)

## --------------------------------------------
## 6. Confidence check
## --------------------------------------------
FeaturePlot(
  sc_singlets,
  features = "SingleR_score"
)

## Optional: flag low-confidence cells
sc_singlets$SingleR_confidence <- ifelse(
  sc_singlets$SingleR_score < 0.3,
  "Low",
  "High"
)

table(sc_singlets$SingleR_confidence)

## --------------------------------------------
## 7. Marker validation (recommended)
## --------------------------------------------
FeaturePlot(
  sc_singlets,
  features = c(
    "CD3D", "CD3E",     # T cells
    "CD4", "CD8A",     # T subsets
    "MS4A1",           # B cells
    "NKG7", "GNLY",    # NK cells
    "LYZ", "S100A8", "S100A9"  # Monocytes
  ),
  ncol = 3
)

## --------------------------------------------
## 8. (Optional) Set identity to SingleR labels
## --------------------------------------------
Idents(sc_singlets) <- "SingleR_label"


## ============================================
## Cluster-level broad immune annotation (FIXED)
## ============================================

## 1. Convert clusters to character ONCE
clusters_char <- as.character(sc_singlets$seurat_clusters)

## 2. Define cluster → broad cell type mapping
cluster_to_celltype <- c(
  "0"  = "T cells",
  "1"  = "Monocytes",
  "2"  = "NK cells",
  "3"  = "T cells",
  "4"  = "T cells",
  "5"  = "T cells",
  "6"  = "T cells",
  "7"  = "B cells",
  "8"  = "Monocytes",
  "9"  = "T cells",
  "10" = "Non-classical monocytes",
  "11" = "NK cells",
  "12" = "Monocytes",
  "13" = "Monocytes",
  "14" = "B cells",
  "15" = "B cells",
  "16" = "Non-classical monocytes",
  "17" = "T cells",
  "18" = "Monocytes",
  "19" = "B cells",
  "20" = "B cells"
)

## 3. Sanity check: did we miss any clusters?
stopifnot(
  all(unique(clusters_char) %in% names(cluster_to_celltype))
)

## 4. Assign broad labels (SAFE indexing)
sc_singlets$Cluster_Broad <- unname(
  cluster_to_celltype[clusters_char]
)

## 5. Set identities
Idents(sc_singlets) <- "Cluster_Broad"

## 6. Check result
table(sc_singlets$Cluster_Broad)

## 7. Plot
DimPlot(
  sc_singlets,
  reduction = "umap",
  group.by = "Cluster_Broad",
  label = TRUE,
  repel = TRUE
)

FeaturePlot(sc_singlets, features = rownames(sc_singlets), max.cutoff = "q95")

FeaturePlot(
  sc_singlets,
  features = rownames(sc_singlets),
  max.cutoff = "q95",
  ncol = 2
)

saveRDS(sc_singlets, file = "sc_singlets_named.RDS")
 
