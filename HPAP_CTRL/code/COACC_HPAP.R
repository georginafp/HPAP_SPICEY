# Load libraries
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

library(Seurat)
library(Signac)
library(SeuratWrappers)
library(cicero)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(tibble)
library(regioneR)
library(stringr)

# Load Seurat object
so <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/test_T1D/data/T1D_HPAP.rds")
DefaultAssay(so) <- "ATAC"

cds <- SeuratWrappers::as.cell_data_set(x = so)
cds.cicero <- make_cicero_cds(cds,
                              reduced_coordinates = reducedDims(cds)$UMAP.HAR)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[paste0("chr", 1:22)]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(cds.cicero, genomic_coords = genome.df, sample_num = 100)

saveRDS(conns, "/homes/users/gfuentes/scratch/projects/spicey_old/test_T1D/data/T1D_LINKS.rds")
