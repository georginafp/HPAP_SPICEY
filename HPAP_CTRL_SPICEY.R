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
         avg_FC = ifelse(avg_FC < 0, 0, avg_FC),  # Clip negative values to 0
         # avg_FC = (avg_FC)^2
         p_val_adj = p.adjust(p_val, method="BH")
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
    RETSI = (norm_log2FC * weight) / (N - 1),
    RETSI = log1p(RETSI)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-6  # To avoid log(0)
near_zero_thresh <- 1e-4  # Treat only truly tiny values as near-zero

retsi_entropy <- retsi_gr %>%
  data.frame() %>%
  # filter(!is.na(RETSI)) %>%
  group_by(region) %>%
  mutate(
    # Replace NA or extremely small RETSI with epsilon
    RETSI_clean = ifelse(is.na(RETSI) | RETSI < near_zero_thresh, epsilon, RETSI),

    # Sum over all RETSI_clean for the region
    RETSI_sum = sum(RETSI_clean, na.rm = TRUE),

    # Define probabilities for entropy (if sum is 0, use uniform)
    RETSI_prob = ifelse(RETSI_sum > 0, RETSI_clean / RETSI_sum, 1 / n()),

    # Entropy component
    entropy_component = -RETSI_prob * log2(RETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types)
  ) %>%
  ungroup()

retsi_atac <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = c("region")) %>%
  # mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()
# saveRDS(retsi_gr_final, "/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_RETSI.rds")



ubiquitous <- c("VAPA", "RPLP0", "HPRT1", "TBP",
                "MRPL19", "ARF1","RAB7A", "GAPDH") # Ubiquitous


ct_markers <- read.delim("/homes/users/gfuentes/shared_data/refs/hi_populations_markers.txt",
                         stringsAsFactors=F)
ct_markers <- ct_markers[ct_markers$CellType %in%
                           unique(retsi_atac$cell_type),] # cell type spz


# Entropy at control sites ----
retsi_labels <- bind_rows(
  # Tissue-specific data
  retsi_atac %>%
    data.frame() %>%
    dplyr::filter(nearestGeneSymbol %in% ct_markers$Gene) %>%
    mutate(type = "Tissue-specific"),

  # Ubiquitous data
  retsi_atac %>%
    data.frame() %>%
    dplyr::filter(nearestGeneSymbol %in% ubiquitous) %>%
    mutate(type = "Ubiquitous")) %>%
  mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific")),
         RETSI = as.numeric(RETSI)
         # GETSI = as.numeric(GETSI)
         ) %>%
  filter(
    !is.na(RETSI),
    # !is.na(GETSI),
    !is.na(type)) %>%
  {rownames(.) <- NULL; .}



retsi_labels %>%
  filter(
    # !is.na(GETSI),
    !is.na(RETSI)) %>%
  pivot_longer(cols = c(norm_entropy),
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





# RNA -------------------------------------------------------------------------
retsi_gr <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/RNA_retsi_gr.rds")

## GETSI computation ----
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = ifelse(avg_FC < 0, 0, avg_FC),  # Clip negative values to 0
         # avg_FC = (avg_FC)^2
         p_val_adj = p.adjust(p_val, method="BH")
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
    GETSI = (norm_log2FC * weight) / (N - 1),
    GETSI = log1p(GETSI)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


## Entropy computation ----
epsilon <- 1e-6  # To avoid log(0)
near_zero_thresh <- 1e-4  # Treat only truly tiny values as near-zero

retsi_entropy <- retsi_gr %>%
  data.frame() %>%
  # filter(!is.na(GETSI)) %>%
  group_by(region) %>%
  mutate(
    # Replace NA or extremely small GETSI with epsilon
    GETSI_clean = ifelse(is.na(GETSI) | GETSI < near_zero_thresh, epsilon, GETSI),

    # Sum over all GETSI_clean for the region
    GETSI_sum = sum(GETSI_clean, na.rm = TRUE),

    # Define probabilities for entropy (if sum is 0, use uniform)
    GETSI_prob = ifelse(GETSI_sum > 0, GETSI_clean / GETSI_sum, 1 / n()),

    # Entropy component
    entropy_component = -GETSI_prob * log2(GETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types)
  ) %>%
  ungroup()

retsi_rna <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = c("region")) %>%
  # mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()
# saveRDS(retsi_gr_final, "/homes/users/gfuentes/scratch/projects/spicey/CTRL_HPAP_RETSI.rds")



ubiquitous <- c("VAPA", "RPLP0", "HPRT1", "TBP",
                "MRPL19", "ARF1","RAB7A", "GAPDH") # Ubiquitous


ct_markers <- read.delim("/homes/users/gfuentes/shared_data/refs/hi_populations_markers.txt",
                         stringsAsFactors=F)
ct_markers <- ct_markers[ct_markers$CellType %in%
                           unique(retsi_rna$cell_type),] # cell type spz


# Entropy at control sites ----
retsi_labels <- bind_rows(
  # Tissue-specific data
  retsi_rna %>%
    data.frame() %>%
    dplyr::filter(symbol %in% ct_markers$Gene) %>%
    mutate(type = "Tissue-specific"),

  # Ubiquitous data
  retsi_rna %>%
    data.frame() %>%
    dplyr::filter(symbol %in% ubiquitous) %>%
    mutate(type = "Ubiquitous")) %>%
  mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific")),
         GETSI = as.numeric(GETSI)
         # GETSI = as.numeric(GETSI)
  ) %>%
  filter(
    !is.na(GETSI),
    # !is.na(GETSI),
    !is.na(type)) %>%
  {rownames(.) <- NULL; .}



retsi_labels %>%
  filter(
    # !is.na(GETSI),
    !is.na(GETSI)) %>%
  pivot_longer(cols = c(norm_entropy),
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


