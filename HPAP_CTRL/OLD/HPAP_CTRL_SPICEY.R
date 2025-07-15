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


# ATAC -------------------------------------------------------------------------
retsi_gr <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/ATAC_retsi_gr.rds")

## RETSI computation ----
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = ifelse(avg_FC < 0, 0, avg_FC)  # Clip negative values to 0
         # avg_FC = (avg_FC)^2
  ) %>%            # Remove positive values
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
    # RETSI = log1p(RETSI)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-6  # To avoid log(0)
near_zero_thresh <- 1e-2  # values below this are treated as 0

retsi_entropy <- retsi_gr %>%
  data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean RETSI values: Replace NA or very small values with epsilon
    RETSI_clean = ifelse(is.na(RETSI) | RETSI < near_zero_thresh, epsilon, RETSI),

    # Check if all RETSI values are effectively 0
    all_near_zero = all(RETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all RETSI values are near-zero, assign uniform probability (1/n_cell_types),
    # otherwise, use the cleaned RETSI values
    RETSI_adj = ifelse(all_near_zero, 1 / n(), RETSI_clean),

    # Normalize the adjusted RETSI values
    RETSI_sum = sum(RETSI_adj, na.rm = TRUE),
    RETSI_prob = RETSI_adj / RETSI_sum,

    # Calculate the entropy component for each cell type
    entropy_component = -RETSI_prob * log2(RETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types)
  ) %>%
  ungroup()


retsi_gr_final <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = c("region")) %>%
  mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()

retsi_gr_final <- saveRDS("/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_RETSI.rds")

# RNA -------------------------------------------------------------------------
retsi_gr <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/RNA_retsi_gr.rds")


## GETSI computation ----
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = ifelse(avg_FC < 0, 0, avg_FC)  # Clip negative values to 0

         # avg_FC = (avg_FC)^2  # only positive values
  ) %>%
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
    GETSI = (norm_log2FC * weight) / (N - 1)
    # GETSI = log1p(GETSI)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-6  # To avoid log(0)
near_zero_thresh <- 1e-2  # values below this are treated as 0

retsi_entropy <- retsi_gr %>%
  data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean GETSI values: Replace NA or very small values with epsilon
    GETSI_clean = ifelse(is.na(GETSI) | GETSI < near_zero_thresh, epsilon, GETSI),

    # Check if all GETSI values are effectively 0
    all_near_zero = all(GETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all GETSI values are near-zero, assign uniform probability (1/n_cell_types),
    # otherwise, use the cleaned GETSI values
    GETSI_adj = ifelse(all_near_zero, 1 / n(), GETSI_clean),

    # Normalize the adjusted GETSI values
    GETSI_sum = sum(GETSI_adj, na.rm = TRUE),
    GETSI_prob = GETSI_adj / GETSI_sum,

    # Calculate the entropy component for each cell type
    entropy_component = -GETSI_prob * log2(GETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types)
  ) %>%
  ungroup()



retsi_gr_final <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = c("region")) %>%
  mutate(GETSI = log1p(GETSI)) %>%
  regioneR::toGRanges()

retsi_gr_final <- saveRDS("/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_GETSI.rds")



# SPICEY bind ------------------------------------------------------------------
retsi_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_GETSI.rds")
retsi_rna <- retsi_rna %>%
  as.data.frame() %>%
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
  dplyr::rename(GETSI_entropy = norm_entropy) %>%
  distinct()

retsi_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_RETSI.rds")
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






## test

spicey_at_control %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  # filter(
  #   (type == "Ubiquitous" & (RETSI_entropy > 0.25 & GETSI_entropy > 0.5)) |
  #     (type == "Tissue-specific" & (RETSI_entropy < 0.7 & GETSI_entropy < 0.7))
  # ) %>%
  filter(annotation == "Promoter") %>%
  group_by(region) %>%
  # filter(any(RETSI != 0, na.rm = TRUE)) %>%
  ungroup() %>%
  # mutate(RETSI = log1p(RETSI)) %>%
  pivot_longer(cols = c(RETSI),
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
