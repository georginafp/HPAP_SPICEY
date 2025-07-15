## OTHER SPICEY PLOTS Ç# # UMAP: T1D and RETSI ----
# HPAP_T1D <- readRDS("/homes/users/gfuentes/scratch/projects/RETSI_OLD/T1D_HPAP.rds")
# Idents(HPAP_T1D) <- HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id
# DefaultAssay(HPAP_T1D) <- "ATAC"
#
#
# GWAS <- readRDS("/homes/users/gfuentes/shared_data/projects-data/t1d-snps/_targets/objects/GWAS")
# T1D_spicey <- subsetByOverlaps(retsi_final_all, GWAS) %>%
#   as.data.frame() %>%
#   mutate(region = paste0(seqnames, "-", start, "-", end)) %>%
#   regioneR::toGRanges()
#
#
# regions_to_remove <- T1D_spicey %>%
#   data.frame() %>%
#   group_by(region) %>%
#   summarise(all_zero_or_negative = all(RETSI < 0, na.rm = TRUE)) %>%
#   filter(all_zero_or_negative == TRUE) %>%
#   pull(region)  # Get the regions that need to be removed
# T1D_spicey_FILT_RETSI <- T1D_spicey[!T1D_spicey$region %in% regions_to_remove]
#
#
# retsi_mean_values <- T1D_spicey_FILT_RETSI %>%
#   data.frame() %>%
#   group_by(cell_type) %>%
#   summarise(mean_RETSI = mean(RETSI, na.rm = TRUE))
# retsi_mapping <- setNames(retsi_mean_values$mean_RETSI, retsi_mean_values$cell_type)
#
# missing_cell_types <- setdiff(unique(HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id), names(retsi_mapping))
# HPAP_T1D@meta.data$RETSI <- retsi_mapping[HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id]
# HPAP_T1D$RETSI <- as.numeric(HPAP_T1D$RETSI)
#
#
# # DimPlot(object = HPAP_T1D, label = TRUE, reduction = "umap.har",
# #         group.by = "RNA_scSorter.Pred_Type.predicted.id",
# #         label.size = 5, label.color = "grey10",
# #         cols = c(   "Alpha" = "#008E80",
# #                     "Ductal" = "#5A7A98",
# #                     "Beta" = "#D37972",
# #                     "Delta" = "#fcbbb6",
# #                     "Acinar" = "#B4609B",
# #                     "Gamma" = "#C6E6FF",
# #                     "Activated-stellate" = "#fbb679",
# #                     "Quiescent-stellate" = "#ffdd6c",
# #                     "Unknown" = "darkred",   # Light gray
# #                     "Immune" = "#6A0DAD",    # Deep purple
# #                     "Endothelial" = "#3D8A8D", # Teal
# #                     "Epsilon" = "red"    # Bright yellow-green
# #                     )) + NoLegend() +
# #   ggtitle("")
#
# # magma_colors <- viridis(100, option = "magma")
# # magma_colors_trimmed <- magma_colors[21:100]
#
# # magma_colors_trimmed <- brewer.pal(n = 10, name = 'RdBu')
# magma_colors_trimmed <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7",
#                           "#4393C3","#2166AC","#053061")
# retsi_umap_plot <- FeaturePlot(HPAP_T1D,
#                                features = "RETSI",
#                                reduction = "umap.har",
#                                pt.size = 0.5,
#                                label.color = "grey10") +
#   scale_color_gradientn(colors = magma_colors_trimmed, oob = scales::squish,
#                         limits=c(0,0.4)) +
#   ggtitle("RETSI at T1D-risk SNPs") +
#   theme_gray() +
#   theme(
#     ggside.panel.scale.x = 0.3,
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 12),
#     plot.title = element_text(size = 16)
#   ) +
#   labs(
#     x = "UMAP 1",  # Set x-axis label
#     y = "UMAP 2",  # Set y-axis label
#     color = "RETSI"  # Set legend title
#   )
#
# # Extract the UMAP coordinates and cell types for annotation
# umap_df <- as.data.frame(Embeddings(HPAP_T1D, reduction = "umap.har"))
# umap_df$cell_type <- HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id
# umap_df$cell_id <- rownames(umap_df)
#
# # Calculate the centroid for each cell type
# centroids <- umap_df %>%
#   group_by(cell_type) %>%
#   summarise(UMAP_1 = mean(umaphar_1), UMAP_2 = mean(umaphar_2)) %>%
#   left_join(retsi_mean_values, by = "cell_type")  # Add the mean RETSI values
#
# # Plot the UMAP with added cell type labels and mean RETSI values
# T1D_RETSI_UMAP <- retsi_umap_plot +
#   geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2,
#                                   label = cell_type),
#             size = 5, color = "grey10") +
#   ggtitle("RETSI at T1D-risk loci") +
#   theme_gray() +
#   theme(
#     ggside.panel.scale.x = 0.3,
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 12),
#     plot.title = element_text(size = 16)
#   ) +
#   labs(
#     x = "UMAP 1",  # Set x-axis label
#     y = "UMAP 2",  # Set y-axis label
#     color = "RETSI"  # Set legend title
#   )
#
# ggsave(paste0(here::here("test"), "/UMAP_RETSI_T1D.png"),
#        plot = T1D_RETSI_UMAP, width = 6, height = 5, units = "in", dpi = 300)
#
#
#
# # UMAP: T1D and GETSI ----------------------------------------------------------
# HPAP_T1D <- readRDS("/homes/users/gfuentes/scratch/projects/RETSI_OLD/T1D_HPAP.rds")
# Idents(HPAP_T1D) <- HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id
# DefaultAssay(HPAP_T1D) <- "RNA"
#
#
# GWAS <- readRDS("/homes/users/gfuentes/shared_data/projects-data/t1d-snps/_targets/objects/GWAS")
# T1D_spicey <- subsetByOverlaps(retsi_final_all, GWAS) %>%
#   as.data.frame() %>%
#   mutate(region = paste0(seqnames, "-", start, "-", end)) %>%
#   regioneR::toGRanges()
#
# regions_to_remove <- T1D_spicey %>%
#   data.frame() %>%
#   group_by(region) %>%
#   summarise(all_zero_or_negative = all(GETSI < 0, na.rm = TRUE)) %>%
#   filter(all_zero_or_negative == TRUE) %>%
#   pull(region)  # Get the regions that need to be removed
#
# T1D_spicey_FILT_GETSI <- T1D_spicey[!T1D_spicey$region %in% regions_to_remove]
#
# getsi_mean_values <- T1D_spicey_FILT_GETSI %>%
#   data.frame() %>%
#   group_by(cell_type) %>%
#   summarise(mean_GETSI = mean(GETSI, na.rm = TRUE))
# getsi_mapping <- setNames(getsi_mean_values$mean_GETSI, getsi_mean_values$cell_type)
#
# missing_cell_types <- setdiff(unique(HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id), names(getsi_mapping))
# HPAP_T1D@meta.data$GETSI <- getsi_mapping[HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id]
# HPAP_T1D$GETSI <- as.numeric(HPAP_T1D$GETSI)
#
#
# # DimPlot(object = HPAP_T1D, label = TRUE, reduction = "umap.har",
# #         group.by = "RNA_scSorter.Pred_Type.predicted.id",
# #         label.size = 5, label.color = "grey10",
# #         cols = c(   "Alpha" = "#008E80",
# #                     "Ductal" = "#5A7A98",
# #                     "Beta" = "#D37972",
# #                     "Delta" = "#fcbbb6",
# #                     "Acinar" = "#B4609B",
# #                     "Gamma" = "#C6E6FF",
# #                     "Activated-stellate" = "#fbb679",
# #                     "Quiescent-stellate" = "#ffdd6c",
# #                     "Unknown" = "darkred",   # Light gray
# #                     "Immune" = "#6A0DAD",    # Deep purple
# #                     "Endothelial" = "#3D8A8D", # Teal
# #                     "Epsilon" = "red"    # Bright yellow-green
# #                     )) + NoLegend() +
# #   ggtitle("")
#
#
# magma_colors_trimmed <- brewer.pal(n = 10, name = 'RdBu')
# # magma_colors_trimmed <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7",
# #                           "#4393C3","#2166AC","#053061")
#
# retsi_umap_plot <- FeaturePlot(HPAP_T1D,
#                                features = "GETSI",
#                                reduction = "umap.har",
#                                pt.size = 0.5,
#                                label.color = "grey10") +
#   scale_color_gradientn(colors = magma_colors_trimmed, oob = scales::squish,
#                         limits=c(0,0.4)) +
#   ggtitle("GETSI at T1D-risk SNPs") +
#   theme_gray() +
#   theme(
#     ggside.panel.scale.x = 0.3,
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 12),
#     plot.title = element_text(size = 16)
#   ) +
#   labs(
#     x = "UMAP 1",  # Set x-axis label
#     y = "UMAP 2",  # Set y-axis label
#     color = "GETSI"  # Set legend title
#   )
#
# # Extract the UMAP coordinates and cell types for annotation
# umap_df <- as.data.frame(Embeddings(HPAP_T1D, reduction = "umap.har"))
# umap_df$cell_type <- HPAP_T1D$RNA_scSorter.Pred_Type.predicted.id
# umap_df$cell_id <- rownames(umap_df)
#
# # Calculate the centroid for each cell type
# centroids <- umap_df %>%
#   group_by(cell_type) %>%
#   summarise(UMAP_1 = mean(umaphar_1), UMAP_2 = mean(umaphar_2)) %>%
#   left_join(getsi_mean_values, by = "cell_type")  # Add the mean GETSI values
#
# # Plot the UMAP with added cell type labels and mean GETSI values
# T1D_GETSI_UMAP <- retsi_umap_plot +
#   geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2,
#                                   label = cell_type),
#             size = 5, color = "grey10") +
#   ggtitle("GETSI at T1D-risk loci") +
#   theme_gray() +
#   theme(
#     ggside.panel.scale.x = 0.3,
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 12),
#     plot.title = element_text(size = 16)
#   ) +
#   labs(
#     x = "UMAP 1",  # Set x-axis label
#     y = "UMAP 2",  # Set y-axis label
#     color = "GETSI"  # Set legend title
#   )
#
# ggsave(paste0(here::here("test"), "/UMAP_GETSI_T1D.png"),
#        plot = T1D_GETSI_UMAP, width = 6, height = 5, units = "in", dpi = 300)
#

# Other plots ------------------------------------------------------------------
T1D_spicey_long <- T1D_spicey %>%
  data.frame() %>%
  gather(key = "Spicey", value = "value", RETSI, GETSI)  # Convert the data into long format

T1D_spicey_long %>%
  ggplot(aes(x = factor(cell_type), y = value, fill = cell_type)) +
  geom_violin(alpha = 0.7, width = 1.2) +  # Violin plot
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.2) +  # Boxplot without outliers
  geom_jitter(aes(color = cell_type), alpha = 0.6, size = 1, width = 0.15) +  # Jittered points
  scale_fill_brewer(palette = "Set3") +  # Color palette for fill
  scale_color_brewer(palette = "Set3") +  # Color palette for jittered points
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(
    x = "Cell Type",
    y = "Value",
    title = "Distribution of RETSI and GETSI Values by Cell Type"
  ) +
  theme(legend.position = "none") +  # Remove legend for color
  facet_grid(. ~ Spicey)  # Facet by Spicey (RETSI or GETSI)



library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(viridis)
library(RColorBrewer)


RETSI_mat <- T1D_spicey %>%
  as.data.frame() %>%
  select(region, cell_type, RETSI) %>%
  pivot_wider(names_from = cell_type, values_from = RETSI) %>%
  column_to_rownames("region") %>%
  as.matrix()

RETSI_mat[is.na(RETSI_mat)] <- 0
RETSI_mat <- RETSI_mat[matrixStats::rowSds(RETSI_mat) > 0, ]

pheatmap(RETSI_mat,
         color = magma_colors_trimmed,
         # color = colorRampPalette(rev(brewer.pal(n = 4, name = "RdYlBu")))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         # scale = "row",  # Normalize across regions
         show_rownames = FALSE,
         main = "RETSI at T1D-risk loci")



## CANDIDATES ----
regions <- c(
  "chr1:172746045-172746940",
  "chr1:212721251-212722145",
  "chr5:132463905-132465398",
  "chr22:30195897-30196863"
)

# Parsear a GRanges
gr <- do.call(c, lapply(regions, function(r) {
  parts <- strcapture("^(chr[^:]+):(\\d+)-(\\d+)$", r, data.frame(chr=character(), start=integer(), end=integer()))
  GRanges(seqnames = parts$chr, ranges = IRanges(start = parts$start, end = parts$end))
}))


candidate_t1d <- subsetByOverlaps(retsi_final_all, gr)
candidate_t1d <- candidate_t1d %>%
  data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  regioneR::toGRanges()


candidate_t1d %>%
  data.frame() %>%
  ggplot(aes(x = GETSI, y = cell_type)) +
  geom_point(size = 2, color = "#1f78b4") +
  facet_wrap(~ region, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(
    x = "RETSI",
    y = "Cell Type"
  ) +
  theme_gray() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



retsi_heatmap <- candidate_t1d %>%
  data.frame() %>%
  mutate(region = factor(region),
         cell_type = factor(cell_type)) %>%
  ggplot(aes(x = region, y = cell_type, fill = RETSI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "RETSI") +
  theme_minimal(base_size = 12) +
  labs(title = "RETSI",
       x = "Region",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())


getsi_heatmap <- candidate_t1d %>%
  data.frame() %>%
  mutate(region = factor(region),
         cell_type = factor(cell_type)) %>%
  ggplot(aes(x = region, y = cell_type, fill = GETSI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "GETSI") +
  theme_minimal(base_size = 12) +
  labs(title = "GETSI",
       x = "Region",
       y = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())


ggpubr::ggarrange(retsi_heatmap, getsi_heatmap)
