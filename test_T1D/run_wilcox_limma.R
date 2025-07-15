.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

library(Seurat)
library(dplyr)
library(tidyr)
library(regioneR)
library(purrr)
library(furrr)
library(ggplot2)
library(Signac)
library(future)
library(DESeq2)

## Set parallel plan
plan("multicore", workers = 4)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
options(future.rng.onMisuse = "ignore")  # To avoid warning

## Load Seurat object
re <- readRDS("/homes/users/gfuentes/scratch/projects/RETSI_OLD/T1D_HPAP.rds")
Idents(re) <- re$RNA_scSorter.Pred_Type.predicted.id

## -------------------------------------------------------------------------
## ATAC analysis
## -------------------------------------------------------------------------
DefaultAssay(re) <- "ATAC"

findmarkers_atac <- function(celltype) {
  # Filter for current cell type
  subset_cells <- WhichCells(re, idents = celltype)

  if (length(subset_cells) < 3) {
    message(paste("Skipping:", celltype, "due to insufficient cells"))
    return(NULL)
  }

  # Try/catch in case FindMarkers fails
  tryCatch({
    t <- FindMarkers(
      re,
      ident.1 = celltype,
      test.use = "wilcox_limma",
      logfc.threshold = 0,
      min.pct = 0
    )

    t <- t %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "region") %>%
      separate(region, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
      mutate(cell_type = celltype) %>%
      dplyr::select(chr, start, end, cell_type, p_val_adj, p_val, avg_log2FC) %>%
      regioneR::toGRanges()

    return(t)
  }, error = function(e) {
    message(paste("Error in celltype:", celltype, "->", e$message))
    return(NULL)
  })
}

# Valid cell types (with ≥ 3 cells)
valid_cell_types <- re$RNA_scSorter.Pred_Type.predicted.id %>%
  as.character() %>%
  table() %>%
  .[. >= 3] %>%
  names()

# Run ATAC markers
gr_list_atac <- future_map(valid_cell_types, findmarkers_atac, .progress = TRUE)
gr_list_atac <- compact(gr_list_atac)  # Remove NULLs
names(gr_list_atac) <- valid_cell_types

saveRDS(gr_list_atac, "/homes/users/gfuentes/scratch/projects/spicey/test/data/RE_T1D_SEURAT_WILCOX_LIMMA.rds")


## -------------------------------------------------------------------------
## RNA analysis
## -------------------------------------------------------------------------
DefaultAssay(re) <- "RNA"

findmarkers_rna <- function(celltype) {
  subset_cells <- WhichCells(re, idents = celltype)

  if (length(subset_cells) < 3) {
    message(paste("Skipping:", celltype, "due to insufficient cells"))
    return(NULL)
  }

  tryCatch({
    t <- FindMarkers(
      re,
      ident.1 = celltype,
      test.use = "wilcox_limma",
      logfc.threshold = 0,
      min.pct = 0
    )

    t <- t %>%
      as.data.frame() %>%
      mutate(cell_type = celltype)

    return(t)
  }, error = function(e) {
    message(paste("Error in RNA celltype:", celltype, "->", e$message))
    return(NULL)
  })
}

# Run RNA markers
gr_list_rna <- future_map(valid_cell_types, findmarkers_rna, .progress = TRUE)
gr_list_rna <- compact(gr_list_rna)
names(gr_list_rna) <- valid_cell_types

saveRDS(gr_list_rna, "/homes/users/gfuentes/scratch/projects/spicey/test/data/GEX_T1D_SEURAT_WILCOX_LIMMA.rds")
