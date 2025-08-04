# .......................#
## QC HPAP-CTRL samples ##
# .......................#

library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggh4x)
library(ggsci)

ctrls <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/CTRL_HPAP.rds")

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


# Explore samples to keep ----
props <- ctrls@meta.data %>%
  data.frame() %>%
  count(orig.ident, cell_type = .data[["RNA_scSorter.Pred_Type.predicted.id"]]) %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

wide_props <- props %>%
  dplyr::select(orig.ident, cell_type, prop) %>%
  pivot_wider(names_from = cell_type, values_from = prop, values_fill = 0)

prop_matrix <- as.data.frame(wide_props[,-1])
rownames(prop_matrix) <- wide_props$orig.ident
dist_mat <- dist(prop_matrix)
hc <- hclust(dist_mat)

ordered_samples <- hc$order
wide_props$orig.ident <- factor(wide_props$orig.ident, levels = wide_props$orig.ident[ordered_samples])

plot_data <- wide_props %>%
  pivot_longer(cols = -orig.ident, names_to = "cell_type", values_to = "prop") %>%
  left_join(ctrls@meta.data %>%
              dplyr::select(orig.ident, RNA_scSorter.Pred_Type.predicted.id),
            by = "orig.ident")

total_cells <- ctrls@meta.data %>%
  count(orig.ident) %>%
  rename(total_cells = n)

plot_data <- plot_data %>%
  left_join(total_cells, by = "orig.ident")


# Plot 1 (bottom plot - proportions)
plot1 <- ggplot(plot_data, aes(x = orig.ident, fill = RNA_scSorter.Pred_Type.predicted.id)) +
  geom_bar(position = "fill", color = "black") +
  scale_x_dendrogram(hclust = hc, expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cell_colors) +  # ✅ Correct usage for fill colors
  labs(x = "HPAP sample", y = "% cells", fill = "Cell type") +
  theme_gray(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.length = unit(8, "pt"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.margin = margin(-4, 0, 0, 0)
  )

# Plot 2 (top plot - total cells)
plot2 <- ggplot(total_cells, aes(x = orig.ident, y = total_cells)) +
  geom_bar(stat = "identity", fill = "#99807f", color = "#755e5c", width = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)))+
  labs(x = NULL, y = "Total cells") +
  theme_gray(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(3, 0, 1, 0)  # ⬅️ NO margin
  )

# Combine with cowplot and remove vertical space
cowplot::plot_grid(
  plot2, plot1,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(0.8, 3),  # Adjust as needed
  greedy = TRUE
)


# Save the final plot
ggsave(paste0(here::here("HPAP_CTRL/figs"), "/HPAP_celltypes_clustering.png"),
       plot = last_plot(), width = 7, height = 6, units = "in", dpi = 300)


# Subset samples ----
samples_to_keep <- c("HPAP-053", "HPAP-075", "HPAP-077", "HPAP-049", "HPAP-063")
hpap_subset <- subset(ctrls, subset = orig.ident %in% samples_to_keep)
saveRDS(hpap_subset, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/SUBSET_CTRL_HPAP_2.rds")



# UMAP plot
hpap_subset$cell_type <- hpap_subset$integrated_clusters_names
hpap_subset$cell_type <- gsub("_A", "", hpap_subset$cell_type)
hpap_subset$cell_type <- gsub("_B", "", hpap_subset$cell_type)
hpap_subset$cell_type <- gsub("Act_stellate", "Stellate", hpap_subset$cell_type)
hpap_subset$cell_type <- gsub("Quies_stellate", "Stellate", hpap_subset$cell_type)

celltypes_to_keep <- names(which(table(hpap_subset$cell_type) > 100))
Idents(hpap_subset) <- hpap_subset$cell_type
hpap_subset <- subset(hpap_subset, idents = celltypes_to_keep)


cell_colors <- c("Alpha" = "#008E80",
                 "Ductal" = "#5A7A98",
                 "Beta" = "#D37972",
                 "Delta" = "#fcbbb6",
                 "Acinar" = "#B4609B",
                 "Gamma" = "#C6E6FF",
                 "Stellate" = "#fbb679",
                 "Unknown" = "#c7c7c7",
                 "Immune" = "#7E9F66",
                 "Epsilon" = "#A3D2CA",
                 "Endothelial" = "#87abbf")

DimPlot(hpap_subset, reduction = "umap.har", group.by = "cell_type",
        label = TRUE, pt.size = 0.5) +
  labs(title = "") +
  theme_gray() +
  NoLegend() +
  scale_color_manual(values = cell_colors)

ggsave(paste0(here::here("HPAP_CTRL/figs"), "/HPAP_UMAP_subset.png"),
       plot = last_plot(), width = 7, height = 6, units = "in", dpi = 300)


## all ---
ctrls$cell_type <- ctrls$integrated_clusters_names
ctrls$cell_type <- gsub("_A", "", ctrls$cell_type)
ctrls$cell_type <- gsub("_B", "", ctrls$cell_type)
ctrls$cell_type <- gsub("Act_stellate", "Stellate", ctrls$cell_type)
ctrls$cell_type <- gsub("Quies_stellate", "Stellate", ctrls$cell_type)

celltypes_to_keep <- names(which(table(ctrls$cell_type) > 100))
Idents(ctrls) <- ctrls$cell_type
ctrls <- subset(ctrls, idents = celltypes_to_keep)


cell_colors <- c("Alpha" = "#008E80",
                 "Ductal" = "#5A7A98",
                 "Beta" = "#D37972",
                 "Delta" = "#fcbbb6",
                 "Acinar" = "#B4609B",
                 "Gamma" = "#C6E6FF",
                 "Stellate" = "#fbb679",
                 "Unknown" = "#c7c7c7",
                 "Immune" = "#7E9F66",
                 "Epsilon" = "#A3D2CA",
                 "Endothelial" = "#87abbf")

DimPlot(ctrls, reduction = "umap.har", group.by = "cell_type",
        label = TRUE, pt.size = 0.5) +
  labs(title = "") +
  theme_gray() +
  NoLegend() +
  scale_color_manual(values = cell_colors)

ggsave(paste0(here::here("HPAP_CTRL/figs"), "/HPAP_UMAP.png"),
       plot = last_plot(), width = 7, height = 6, units = "in", dpi = 300)


(DimPlot(ctrls, reduction = "umap.har", group.by = "cell_type",
        label = TRUE, pt.size = 0.5) +
  labs(title = "") +
  theme_gray() +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_color_manual(values = cell_colors)) +
  ggtitle("HPAP samples") +
(DimPlot(hpap_subset, reduction = "umap.har", group.by = "cell_type",
          label = TRUE, pt.size = 0.5) +
  labs(title = "") +
  theme_gray() +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_color_manual(values = cell_colors) +
  ggtitle("Filtered HPAP samples"))

ggsave(paste0(here::here("HPAP_CTRL/figs"), "/HPAP_UMAP_FINAL.png"),
       plot = last_plot(), width = 8, height = 5, units = "in", dpi = 300)
