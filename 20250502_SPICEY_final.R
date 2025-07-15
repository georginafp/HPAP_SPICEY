## .......................................................................... ##
# Compute Shannon entropy of RETSI/GETSI values across cell types, per region
# this tells you how uniformly a region is active across cell types.
# For example, if RETSI is high in one cell type and near-zero in others,
# entropy is low (cell-type-specific).
# If RETSI is evenly spread across all cell types, entropy is high (ubiquitous).
## .......................................................................... ##

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tidyr)
library(regioneR)
library(purrr)
library(Signac)
library(purrr)
library(ggside)
library(viridis)
library(tibble)
library(ggsignif)
library(RColorBrewer)
library(ggplotify)
library(pheatmap)
library(patchwork)
library(scales)
library(data.table)
library(here)


# Function to extract top cell type per region for a given measure
get_top_celltype_prop <- function(df, measure) {
  df %>%
    data.frame() %>%
    filter(annotation == "Distal",
          (RETSI_entropy < 1 & GETSI_entropy < 1)) %>%
    group_by(region) %>%
    slice_max(order_by = .data[[measure]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    count(cell_type) %>%
    mutate(
      spicey_measure = measure,
      prop = n / sum(n)
    )
}

get_regions <- function(path, from_rds = FALSE) {
  gwas <- if (from_rds) {
    readRDS(path)
  } else {
    data.table::fread(path) %>%
      select("CHR_ID", "CHR_POS") %>%
      mutate(CHR_ID = paste0("chr", CHR_ID),
             CHR_POS = as.numeric(CHR_POS)) %>%
      makeGRangesFromDataFrame(seqnames.field = "CHR_ID",
                               start.field = "CHR_POS",
                               end.field = "CHR_POS",
                               keep.extra.columns = FALSE,
                               na.rm = TRUE)
  }

  subsetByOverlaps(retsi_gr_final, gwas) %>%
    data.frame() %>%
    mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
    pull(region) %>%
    unique()
}


make_zscore_heatmap <- function(gr, region_df, score_var) {
  gr %>%
    data.frame() %>%
    mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
    left_join(region_df, by = "region") %>%
    filter(
      !is.na(DISEASE),
      RETSI_entropy < 1,
      GETSI_entropy < 1,
      !is.na(RETSI),
      !is.na(GETSI)
    ) %>%
    group_by(DISEASE, cell_type) %>%
    summarise(score = mean(.data[[score_var]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = cell_type, values_from = score) %>%
    column_to_rownames("DISEASE") %>%
    as.matrix()
}



# Plotting function
plot_zscore_heatmap <- function(mat_wide, value_name, title, filename) {
  mat_df <- mat_wide %>%
    as.data.frame() %>%
    rownames_to_column("DISEASE") %>%
    pivot_longer(-DISEASE, names_to = "cell_type", values_to = "value") %>%
    mutate(
      DISEASE = factor(DISEASE, levels = disease_order),
      cell_type = factor(cell_type, levels = cell_type_order),
      z = scale(value)[, 1]
    )

  z_min <- min(mat_df$z, na.rm = TRUE)
  z_max <- max(mat_df$z, na.rm = TRUE)

  ggplot(mat_df, aes(x = cell_type, y = DISEASE, fill = z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(11, "RdBu")),
      values = scales::rescale(c(z_min, 0, z_max)),
      name = paste("Z-scored\n", value_name)
    ) +
    labs(x = "Cell type", y = "Disease") +
    theme_gray() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 35, hjust = 1),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12),
      plot.title = element_text(size = 16)
    ) +
    ggtitle(title) -> p

  ggsave(here("test", filename), plot = p, width = 6, height = 5, units = "in", dpi = 300)
}



# CTRL-HPAP ---------------------------------------------------------------------
re <- readRDS("~/scratch/projects/spicey/test/data/CTRL_HPAP.rds")
Idents(re) <- re$RNA_scSorter.Pred_Type.predicted.id

cell_colors <- c("Alpha" = "#008E80",
                 "Ductal" = "#5A7A98",
                 "Beta" = "#D37972",
                 "Delta" = "#fcbbb6",
                 "Acinar" = "#B4609B",
                 "Gamma" = "#C6E6FF",
                 "Activated-stellate" = "#fbb679",
                 "Quiescent-stellate" = "#ffdd6c",
                 "Unknown" = "#c7c7c7",
                 "Immune" = "#7E9F66",
                 "Epsilon" = "#A3D2CA",
                 "Endothelial" = "#87abbf")


DimPlot(re, reduction="umap.har", group.by="RNA_scSorter.Pred_Type.predicted.id",
                label=TRUE,
                raster=FALSE, shuffle = TRUE, pt.size=.01) +
  labs(x="UMAP 1", y="UMAP 2") +
  scale_color_manual(values = cell_colors) +
  theme_gray() +
  ggtitle("") +
  theme(legend.position = "none",
        ggside.panel.scale.x = 0.3,
        axis.title = element_text(size = 14),  # Axis titles (x and y)
        axis.text = element_text(size = 12),  # Axis text
        strip.text = element_text(size = 12),  # Facet labels
        plot.title = element_text(size = 16)) +
  NoLegend()

ggsave(paste0(here::here("test"), "/HPAP_CTRL_UMAP.png"),
       plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)



# ATAC -------------------------------------------------------------------------
merged_df <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA/ATAC_CTRL_SEURAT_WL.rds")
merged_df <- unlist(GRangesList(merged_df))

retsi_gr <- merged_df %>%
  data.frame() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand)) %>%
  regioneR::toGRanges()

anno <- ChIPseeker::annotatePeak(retsi_gr,
                                 TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                 verbose=FALSE)

retsi_gr$distanceToTSS <- data.frame(anno)$distanceToTSS
retsi_gr$annotation <- "Distal"
retsi_gr$annotation[abs(retsi_gr$distanceToTSS)<=2e3] <- "Promoter"

## annotate peaks
genes <- readRDS("/homes/users/gfuentes/scratch/projects/coculture/data/_targets/objects/GENES")$gr
retsi_gr <- retsi_gr %>%
  plyranges::mutate(nearestGeneSymbol = genes$external_gene_name[nearest(retsi_gr, promoters(genes, 1, 0))])
names(retsi_gr) <- NULL


## RETSI computation ----
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = 2^(avg_FC),
         p_val_adj = p.adjust(p_val, method = "fdr")  # Adjust p-value using FDR
  ) %>%            # remove negative sign
  group_by(region) %>%
  mutate(
    N = n(),
    max_log2FC = max(avg_FC, na.rm = TRUE),
    max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)  # Avoid zero division
  ) %>%
  ungroup() %>%
  mutate(p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
         weight = -log10(p_val_adj)) %>%
  group_by(cell_type, region) %>%
  mutate(
    norm_log2FC = (avg_FC / max_log2FC),
    RETSI = (norm_log2FC * weight) / (N - 1)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-10         # To avoid log(0)
near_zero_thresh <- 0.05 # Values below this are treated as near-zero

retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean RETSI values
    RETSI_clean = ifelse(is.na(RETSI) | RETSI < near_zero_thresh, epsilon, RETSI),

    # Check if all RETSI values are effectively 0
    all_near_zero = all(RETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all values are near-zero, use uniform; otherwise, normalize RETSI_clean
    RETSI_adj = ifelse(all_near_zero, 1 / n(), RETSI_clean),
    RETSI_sum = sum(RETSI_adj, na.rm = TRUE),
    RETSI_prob = RETSI_adj / RETSI_sum,

    # Final renormalization in case of numeric drift
    RETSI_prob = RETSI_prob / sum(RETSI_prob, na.rm = TRUE),

    # Entropy computation
    entropy_component = -RETSI_prob * log2(RETSI_prob),
    RETSI_prob_check = sum(RETSI_prob, na.rm = TRUE)  # Sanity check
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types),
    RETSI_prob_check = first(RETSI_prob_check)
  ) %>%
  ungroup()


retsi_atac <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "region") %>%
  mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()
saveRDS(retsi_atac, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_RETSI.rds")




# RNA --------------------------------------------------------------------------
gr_list <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA/RNA_CTRL_SEURAT_WL.rds")

gr_df_list <- map(gr_list, ~ as.data.frame(.x))

source("/homes/users/gfuentes/scratch/projects/spicey/R/utils.R")
gr_list <- lapply(names(gr_df_list), function(name) {

  df <- gr_df_list[[name]]
  df$symbol <- rownames(df)
  df$ensembl_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys=df$symbol,
                                         column="ENSEMBL",
                                         keytype="SYMBOL",
                                         multiVals="first")

  df_annot <- dplyr::left_join(df,
                               biomart_genes()$df %>% dplyr::select(chromosome_name,
                                                                    start_position,,
                                                                    end_position,
                                                                    strand,
                                                                    ensembl_gene_id,
                                                                    gene_biotype),
                               by=c(ensembl_id="ensembl_gene_id"))

  if (nrow(df_annot) == 0) {
    message("Missing entry annotation for ", name)
    return(NULL)
  }

  df_annot %>%
    dplyr::filter(!is.na(chromosome_name),
                  !is.na(start_position),
                  !is.na(end_position),
                  gene_biotype == "protein_coding") %>%
    dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
    dplyr::select(chromosome_name,
                  start_position,
                  end_position,
                  strand,
                  everything()) %>%
    dplyr::mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
    makeGRangesFromDataFrame(seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field=c("start", "start_position"),
                             end.field=c("end", "stop", "end_position"),
                             strand.field="strand",
                             keep.extra.columns = T) %>%
    sort()
})


# Rename list with cell types
names(gr_list) <- names(gr_df_list)
gr_list <- Filter(Negate(is.null), gr_list) # Remove NULLs


# Bind together in one GRanges
retsi_gr <- imap(gr_list, ~ {
  mcols(.x)$cell_type <- .y  # .x = GRanges, .y = name
  .x
}) %>%
  purrr::reduce(c) #214896


retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = 2^(avg_FC),
         # avg_FC = ifelse(avg_FC < 0, 0, avg_FC)  # Clip negative values
         p_val_adj = p.adjust(p_val, method = "fdr")  # Adjust p-value using FDR
  ) %>%            # remove negative sign
  group_by(symbol) %>%
  mutate(
    N = n(),
    max_log2FC = max(avg_FC, na.rm = TRUE),
    max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)  # Avoid zero division
  ) %>%
  ungroup() %>%
  mutate(p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
         weight = -log10(p_val_adj)) %>%
  group_by(cell_type, symbol) %>%
  mutate(
    norm_log2FC = (avg_FC / max_log2FC),
    GETSI = (norm_log2FC * weight) / (N - 1)
    # GETSI = (norm_log2FC * ifelse(p_val_adj <= 0.05, 1, 0)) / (N - 1)  # Set weight to 1 or 0 based on p_val_adj
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-10       # To avoid log(0)
near_zero_thresh <- 0.05 # Values below this are treated as near-zero

retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean GETSI values
    GETSI_clean = ifelse(is.na(GETSI) | GETSI < near_zero_thresh, epsilon, GETSI),

    # Check if all GETSI values are effectively 0
    all_near_zero = all(GETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all values are near-zero, use uniform; otherwise, normalize GETSI_clean
    GETSI_adj = ifelse(all_near_zero, 1 / n(), GETSI_clean),
    GETSI_sum = sum(GETSI_adj, na.rm = TRUE),
    GETSI_prob = GETSI_adj / GETSI_sum,

    # Final renormalization in case of numeric drift
    GETSI_prob = GETSI_prob / sum(GETSI_prob, na.rm = TRUE),

    # Entropy computation
    entropy_component = -GETSI_prob * log2(GETSI_prob),
    GETSI_prob_check = sum(GETSI_prob, na.rm = TRUE)  # Sanity check
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types),
    GETSI_prob_check = first(GETSI_prob_check)
  ) %>%
  ungroup()

retsi_rna <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "region") %>%
  mutate(GETSI = log1p(GETSI)) %>%
  regioneR::toGRanges()
saveRDS(retsi_rna, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_GETSI.rds")



# SPICEY bind ------------------------------------------------------------------
retsi_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_GETSI.rds")
retsi_rna <- retsi_rna %>%
  as.data.frame() %>%
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
  dplyr::rename(GETSI_entropy = norm_entropy) %>%
  distinct()

retsi_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_RETSI.rds")
retsi_atac <- retsi_atac %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, annotation, nearestGeneSymbol, cell_type, RETSI, norm_entropy) %>%
  dplyr::rename(symbol = nearestGeneSymbol,
                RETSI_entropy = norm_entropy)


retsi_gr_final <- retsi_atac %>%
  left_join(retsi_rna %>%
              dplyr::select(symbol, GETSI, cell_type, GETSI_entropy),  # Select relevant columns
            by = c("symbol", "cell_type")) %>%
  regioneR::toGRanges()




## markers .....................................................................
ubiquitous <- c("VAPA", "RPLP0", "HPRT1", "TBP",
                "MRPL19", "ARF1","RAB7A") # Ubiquitous

# ubiquitous <- data.table::fread("/homes/users/gfuentes/shared_data/refs/USTFs_Kribelbauer_2024.csv")
# ubiquitous <- ubiquitous$UST
ct_markers <- read.delim("/homes/users/gfuentes/shared_data/refs/hi_populations_markers.txt",
                         stringsAsFactors=F)
ct_markers <- ct_markers[ct_markers$CellType %in%
                           unique(retsi_gr_final$cell_type),] # cell type spz

# ct_markers <- data.table::fread("/homes/users/gfuentes/shared_data/refs/Islet_important_Motifs.csv")


# Entropy at control sites ----
retsi_labels <- bind_rows(
  # Tissue-specific data
  retsi_gr_final %>%
    data.frame() %>%
    dplyr::filter(genes_coacc %in% ct_markers$Gene) %>%
    mutate(type = "Tissue-specific"),

  # Ubiquitous data
  retsi_gr_final %>%
    data.frame() %>%
    dplyr::filter(genes_coacc %in% ubiquitous) %>%
    mutate(type = "Ubiquitous")) %>%
  mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific")),
         RETSI = as.numeric(RETSI),
         GETSI_coacc = as.numeric(GETSI_coacc)) %>%
  filter(
    !is.na(RETSI),
    !is.na(GETSI_coacc),
    !is.na(type)) %>%
  {rownames(.) <- NULL; .}



retsi_labels %>%
  filter(!is.na(GETSI_coacc),
         !is.na(RETSI)) %>%
  pivot_longer(cols = c(RETSI_entropy, GETSI_entropy),
               names_to = "spicey_measure",
               values_to = "value") %>%
  mutate(spicey_measure = gsub("_entropy", "", spicey_measure)) %>%
  ggplot(aes(x = type,
             y = value,
             fill = type)) +
  geom_violin() +
  geom_boxplot(alpha = 0.4) +  # Adjust outlier behavior if needed
  theme_gray() +
  labs(x = "", y = "Entropy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = c("Ubiquitous" = "#E69F00",
                               # "Others" = "#009E73",
                               "Tissue-specific" = "#009E73")) +
  # ylim(c(1,1.2))+
  geom_signif(comparisons = list(
    c("Tissue-specific", "Ubiquitous")),
    y_position = c(1, 1.1),  # Adjust y position for annotations
    map_signif_level = TRUE,  # Display significance level
    textsize = 2.3,
    vjust = 0.01,               # Moves text slightly above the line
    tip_length = 0.02         # Negative value to flip arrows **upwards**
  ) +
  facet_grid(~spicey_measure) +
  labs(title = "")

ggsave(paste0(here::here("HPAP_CTRL/figs"), "/SPICEY_ENTROPY_CT_MARKERS.png"),
       plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)




# SPICEY at control sites ----
# spicey_at_control <- bind_rows(
#   # Tissue-specific data
#   retsi_gr_final %>%
#     data.frame() %>%
#     inner_join(ct_markers %>% dplyr::select(-c("Source")),
#                by = c("symbol" = "Gene", "cell_type" = "CellType")) %>%
#     distinct(cell_type, symbol, .keep_all = TRUE) %>%
#     mutate(type = "Tissue-specific"),
#
#   # Ubiquitous data
#   retsi_gr_final %>%
#     data.frame() %>%
#     dplyr::filter(symbol %in% ubiquitous) %>%
#     mutate(type = "Ubiquitous")) %>%
#   mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific")),
#          RETSI = as.numeric(RETSI)) %>%
#   filter(!is.na(RETSI),
#          !is.na(GETSI),
#          !is.na(type)) %>%
#   {rownames(.) <- NULL; .}



spicey_at_control <- bind_rows(
  # Tissue-specific data
  retsi_gr_final %>%
    data.frame() %>%
    dplyr::filter(symbol %in% ct_markers$Gene) %>%
    mutate(type = "Tissue-specific"),

  # Ubiquitous data
  retsi_gr_final %>%
    data.frame() %>%
    dplyr::filter(symbol %in% ubiquitous) %>%
    mutate(type = "Ubiquitous")) %>%
  mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific")),
         RETSI = as.numeric(RETSI)) %>%
  filter(!is.na(RETSI),
         !is.na(GETSI),
         !is.na(type)) %>%
  {rownames(.) <- NULL; .}




## Promoters ----
spicey_at_control %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  filter(
    (type == "Ubiquitous" & (RETSI_entropy > 0.5 | GETSI_entropy > 0.5)) |
    (type == "Tissue-specific" & (RETSI_entropy < 0.25 | GETSI_entropy < 0.25))
  ) %>%
  filter(annotation == "Promoter") %>%
  group_by(region) %>%
  # filter(any(RETSI != 0, na.rm = TRUE)) %>%
  ungroup() %>%
  # mutate(RETSI = log1p(RETSI)) %>%
  pivot_longer(cols = c(RETSI, GETSI),
               names_to = "spicey_measure",
               values_to = "value") %>%
  ggplot(aes(x = value, color = type, fill = type)) +
  geom_density(alpha = .6) +
  ggside::geom_xsideboxplot(aes(y = type, fill = type), orientation = "y") +
  ggside::scale_xsidey_discrete() +
  scale_y_continuous(name = "Density") +
  scale_color_manual(values=c("grey40", "#c57c2a")) + # around
  scale_fill_manual(values=c("#bcbc77", "#e6ba88")) + # fig
  labs(x = "Spicey measure") +
  ggtitle("Promoter regions") +
  theme(legend.position = "none",
        ggside.panel.scale.x = 0.3,
        axis.title = element_text(size = 14),  # Axis titles (x and y)
        axis.text = element_text(size = 12),  # Axis text
        strip.text = element_text(size = 12),  # Facet labels
        plot.title = element_text(size = 16)) +
  facet_grid(~ spicey_measure)



ggsave(paste0(here::here("test"), "/SPICEY_PROM_CT_MARKERS.png"),
       plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)



## Enhancers ----
enh_datasets <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/test/data/TEST_ENHANCERS/islet_enhancers_datasets_hg38.rds")
enh_ds <- enh_datasets$hg38_cluster_enhancerClusters
# enh_ds <- subsetByOverlaps(enh_datasets$hg38_cluster_enhancerClusters, enh_datasets$hg38_cluster_stretchEnhancers)
spicey_at_control_enh <- subsetByOverlaps(retsi_gr_final, enh_ds)
spicey_at_control_enh <- spicey_at_control_enh %>%
  data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  regioneR::toGRanges()

# How many enhancers does each cell type dominate?
# Which cell types most frequently have the highest RETSI in tissue-specific enhancers?
cell_order <- c("Unknown",
                "Gamma",
                "Delta",
                "Epsilon",
                "Immune",
                "Quiescent-stellate",
                "Endothelial",
                "Activated-stellate",
                "Alpha",
                "Beta",
                "Ductal",
                "Acinar")

# Now using the `cell_order` and `cell_colors` in the ggplot code
spicey_at_control_enh_top <- bind_rows(
  get_top_celltype_prop(spicey_at_control_enh, "RETSI"),
  get_top_celltype_prop(spicey_at_control_enh, "GETSI")
) %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_order),
    spicey_measure = factor(spicey_measure, levels = c("RETSI", "GETSI"))
  )

spicey_at_control_enh_top %>%
  ggplot(aes(x = spicey_measure, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  # geom_label(
  #   aes(label = comma(n, accuracy = 1)),
  #   position = position_fill(vjust = 0.5),
  #   size = 3
  # ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "% enhancer SPICEY-dominant") +
  scale_fill_manual(values = cell_colors) +
  xlab("Spicey measure") +
  theme_gray(base_size = 14) +
  theme(ggside.panel.scale.x = 0.3,
        axis.title = element_text(size = 14),  # Axis titles (x and y)
        axis.text = element_text(size = 12),  # Axis text
        strip.text = element_text(size = 12),  # Facet labels
        plot.title = element_text(size = 16)) +
  guides(fill = guide_legend(title = "Cell type")) +
  ggtitle("Enhancer regions")

ggsave(paste0(here::here("test"), "/SPICEY_CLUST_ENH.png"),
       plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# Genetics ----
disease_paths <- list(
  T1D = list(path = "/homes/users/gfuentes/shared_data/projects-data/t1d-snps/_targets/objects/GWAS", from_rds = TRUE),
  T2D = list(path = "/homes/users/gfuentes/shared_data/refs/T2D_GWAS_ASSOCIATIONS.tsv"),
  # COPD = list(path = "/homes/users/gfuentes/shared_data/refs/COPD_ASSOCIATIONS.tsv"),
  LUNG_CANCER = list(path = "/homes/users/gfuentes/shared_data/refs/LUNG_CANCER_ASSOCIATIONS.tsv"),
  SCHIZO = list(path = "/homes/users/gfuentes/shared_data/refs/SCHIZOFRENIA_ASSOCIATIONS.tsv"),
  RHEUMATOID = list(path = "/homes/users/gfuentes/shared_data/refs/RHEUMATOID_ASSOCIATIONS.tsv"),
  PSORIASIS = list(path = "/homes/users/gfuentes/shared_data/refs/PSORIASIS_ASSOCIATIONS.tsv")
  # GLAUCOMA = list(path = "/homes/users/gfuentes/shared_data/refs/GLAUCOMA_ASSOCIATIONS.tsv")
)

region_lists <- lapply(disease_paths, function(x) {
  get_regions(x$path, from_rds = x$from_rds %||% FALSE)
})

region_df <- stack(region_lists) %>%
  rename(region = values, DISEASE = ind)


# Compute z-score for spicey
mat_wide_retsi <- make_zscore_heatmap(retsi_gr_final, region_df, "RETSI")
mat_wide_getsi <- make_zscore_heatmap(retsi_gr_final, region_df, "GETSI")

# Use RETSI clustering order for both
disease_order <- rownames(mat_wide_retsi)[hclust(dist(mat_wide_retsi))$order]
cell_type_order <- colnames(mat_wide_retsi)[hclust(dist(t(mat_wide_retsi)))$order]

# plot heatmaps
plot_zscore_heatmap(mat_wide_retsi, "RETSI", "RETSI enrichment of GWAS-SNPs", "RETSI_GWAS.png")
plot_zscore_heatmap(mat_wide_getsi, "GETSI", "GETSI enrichment of GWAS-SNPs", "GETSI_GWAS.png")
