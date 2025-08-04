## Differential analysis ----


### ATAC ----
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

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
re <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/SUBSET_CTRL_HPAP_2.rds")
re$cell_type <- re$integrated_clusters_names
re$cell_type <- gsub("_A", "", re$cell_type)
re$cell_type <- gsub("_B", "", re$cell_type)
re$cell_type <- gsub("Act_stellate", "Stellate", re$cell_type)
re$cell_type <- gsub("Quies_stellate", "Stellate", re$cell_type)

celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
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
    min.pct = 0.01
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

valid_cell_types <- re$cell_type %>%
  as.character() %>%
  table() %>%
  .[. >= 3] %>%
  names()


# Apply function only to valid cell types
gr_list <- expand.grid(celltype = valid_cell_types) %>%
  pmap(~ findmarkers(re, .x)) %>%
  compact()  # Remove NULL entries

saveRDS(gr_list, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_ATAC_HPAP-CTRL.rds")



### RNA ----
# Load dataset
re <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/SUBSET_CTRL_HPAP_2.rds")
re$cell_type <- re$integrated_clusters_names
re$cell_type <- gsub("_A", "", re$cell_type)
re$cell_type <- gsub("_B", "", re$cell_type)
re$cell_type <- gsub("Act_stellate", "Stellate", re$cell_type)
re$cell_type <- gsub("Quies_stellate", "Stellate", re$cell_type)

celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
DefaultAssay(re) <- "RNA"

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
    min.pct = 0.01
  )

  t <- t %>%
    as.data.frame() %>%
    mutate(
      cell_type = celltype
    )
  return(t)
}

# Define cell types
cell_types <- as.character(unique(Idents(re)))  # Convert Idents to character
valid_cell_types <- re$cell_type %>%
  as.character() %>%
  table() %>%
  .[. >= 3] %>%
  names()

# Apply function only to valid cell types
gr_list <- expand.grid(celltype = valid_cell_types) %>%
  pmap(~ findmarkers(re, .x)) %>%
  compact()  # Remove NULL entries
names(gr_list) <- valid_cell_types
saveRDS(gr_list, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_RNA_HPAP-CTRL.rds")
