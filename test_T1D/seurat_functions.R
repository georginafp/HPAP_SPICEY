library(Seurat)
library(dplyr)
library(tidyr)
library(regioneR)
library(purrr)
library(ggplot2)
library(Signac)
library(future)
library(DESeq2)

# Load dataset
re <- readRDS("/homes/users/gfuentes/scratch/projects/RETSI_OLD/T1D_HPAP.rds")
Idents(re) <- re$RNA_scSorter.Pred_Type.predicted.id
DefaultAssay(re) <- "ATAC"

# Set parallelization once at script startup
plan("multicore", workers = 4)
findmarkers <- function(re, celltype) {

  # Check if the celltype exists in the subset
  if (!(celltype %in% as.character(Idents(re)))) {
    message(paste("Skipping:", celltype, "because it is not present in the subset"))
    return(NULL)
  }

  # Skip if fewer than 3 cells
  if (sum(as.character(Idents(re)) == celltype) < 3) {
    message(paste("Skipping:", celltype, "due to insufficient cells"))
    return(NULL)
  }

  t <- FindMarkers(
    re,
    ident.1 = celltype,
    logfc.threshold = 0,
    min.pct = 0
  )

  t <- t %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "region") %>%
    separate(region, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
    mutate(
      cell_type = celltype
    ) %>%
    dplyr::select(chr, start, end, cell_type, p_val_adj, p_val, avg_log2FC) %>%
    regioneR::toGRanges()

  return(t)
}

# Define cell types
cell_types <- as.character(unique(Idents(re)))  # Convert Idents to character

valid_cell_types <- re$RNA_scSorter.Pred_Type.predicted.id %>%
  as.character() %>%
  table() %>%
  .[. >= 3] %>%
  names()


# Apply function only to valid cell types
gr_list <- expand.grid(celltype = valid_cell_types) %>%
  pmap(~ findmarkers(re, .x)) %>%
  compact()  # Remove NULL entries
names(gr_list) <- valid_cell_types

saveRDS(gr_list, "/homes/users/gfuentes/scratch/projects/spicey/test/data/RE_T1D_SEURAT_WILCOX.rds")






genes <- readRDS("/homes/users/gfuentes/scratch/projects/coculture/data/_targets/objects/GENES")$gr

retsi_gr <- lapply(gr_list, function(x) {
  anno <- ChIPseeker::annotatePeak(x,
                                   TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   verbose=FALSE)

  x$distanceToTSS <- data.frame(anno)$distanceToTSS
  x$annotation <- "Distal"
  x$annotation[abs(x$distanceToTSS)<=2e3] <- "Promoter"

  ## annotate peaks
  x <- x %>%
    plyranges::mutate(nearestGeneSymbol = genes$external_gene_name[nearest(x, promoters(genes, 1, 0))])
  return(x)
})

retsi_gr <- imap(retsi_gr, ~ {
  mcols(.x)$cell_type <- .y  # .x = GRanges, .y = name
  .x
}) %>%
  purrr::reduce(c) #4197843


# saveRDS(retsi_gr, "/homes/users/gfuentes/scratch/projects/spicey/test/data/RE_T1D_SEURAT_WILCOX.rds")


retsi_gr <- retsi_gr %>%
  data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_log2FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_log2FC = (avg_log2FC)^2) %>%      # avoid negative values
  group_by(region) %>%
  mutate(
    N = n(),
    max_log2FC = max(avg_log2FC, na.rm = TRUE),
    max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)  # Avoid zero division
  ) %>%
  ungroup() %>%
  mutate(p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
         weight = -log10(p_val_adj)) %>%
  group_by(cell_type, region) %>%
  mutate(
    norm_log2FC = (avg_log2FC / max_log2FC),
    RETSI = (norm_log2FC * weight) / (N - 1)
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()



# Scale RETSI ----
retsi_gr <- retsi_gr %>%
  mutate(RETSI = log1p(RETSI))    # rescale data


genes_int <- c("VAPA", "RPLP0", "HPRT1", "TBP",
               "MRPL19", "ARF1","RAB7A", "GAPDH")

retsi_long <- retsi_gr  %>%
  as.data.frame() %>%
  dplyr::filter(annotation == "Promoter") %>%
  distinct()

retsi_labels <- retsi_long %>%
  dplyr::filter(nearestGeneSymbol %in% genes_int) %>%
  group_by(nearestGeneSymbol, cell_type) %>%
  summarise(RETSI = mean(RETSI, na.rm = TRUE)) %>%
  ungroup()  # ungroup the data after summarizing


hkg <- retsi_long %>%
  data.frame() %>%
  dplyr::filter(!is.na(RETSI)) %>%
  ggplot(aes(factor(cell_type), RETSI)) +
  geom_boxplot(notch=TRUE) +
  geom_label(data = retsi_labels,
             aes(x = factor(cell_type),
                 y = RETSI,
                 label = nearestGeneSymbol,
                 color = nearestGeneSymbol),
             size = 3,
             label.padding = unit(0.1, "lines"),
             position = position_jitter(width = .3)) +
  scale_color_discrete(guide="none") +
  labs(x="cell type", y="RETSI")


markers <- read.delim("/homes/users/gfuentes/shared_data/refs/hi_populations_markers.txt",
                      stringsAsFactors=F)
markers <- markers[markers$CellType %in% unique(retsi_long$cell_type),]

retsi_labels <- retsi_long %>%
  inner_join(markers, by = c("nearestGeneSymbol" = "Gene","cell_type" = "CellType")) %>%
  distinct(cell_type, nearestGeneSymbol, .keep_all = TRUE) %>%
  group_by(nearestGeneSymbol, cell_type) %>%
  summarise(RETSI = mean(RETSI, na.rm = TRUE)) %>%  # Summarize only RETSI
  ungroup()  # ungroup the data after summarizing


celltype_spz <- retsi_long %>%
  data.frame() %>%
  dplyr::filter(!is.na(RETSI)) %>%
  ggplot(aes(factor(cell_type), RETSI)) +
  geom_boxplot(notch = TRUE) +
  geom_label(data = retsi_labels,
             aes(x = factor(cell_type),
                 y = RETSI,
                 label = nearestGeneSymbol,
                 color = nearestGeneSymbol), size = 3,
             label.padding = unit(0.1, "lines"),
             position = position_jitter(width = .3)) +
  scale_color_discrete(guide = "none") +
  labs(x = "Cell Type", y = "RETSI")

# plot
hkg + celltype_spz



saveRDS(retsi_gr, "/homes/users/gfuentes/scratch/projects/spicey/test/data/RE_T1D_SEURAT_RETSI.rds")
