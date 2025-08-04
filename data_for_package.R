# params
library(usethis)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(GenomicRanges)


# functions ----
sample_regions <- function(data_list, max_rows = 500) {
  lapply(data_list, function(obj) {
    n <- if (inherits(obj, "GRanges")) length(obj) else nrow(obj)
    if (n > max_rows) {
      indices <- sample(n, max_rows)
      if (inherits(obj, "GRanges")) {
        obj <- obj[indices]
      } else {
        obj <- obj[indices, , drop = FALSE]
      }
    }
    return(obj)
  })
}


get_gsheet <- function(url, sheet) {
  googlesheets4::gs4_deauth()
  googlesheets4::gs4_auth()
  gsheet <- googlesheets4::read_sheet(url, sheet)
  return(gsheet)}



clean_cicero_links <- function(df, columns = c("Peak1", "Peak2", "coaccess")) {
  df <- as.data.frame(df)
  cols_present <- intersect(columns, colnames(df))
  df <- df[, cols_present, drop = FALSE]

  df[] <- lapply(df, function(col) {
    if (is.factor(col)) {
      col <- droplevels(col)
      col <- as.character(col)
    } else if (inherits(col, "character")) {
      col <- as.character(col)
    }
    return(col)
  })

  return(df)
}



# data ----
data(atac)
data(rna)
cicero_links <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/CTRL_LINKS_2.rds")
retsi <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/HPAP_CTRL_SPICEY_COACC.rds")
url <- "https://docs.google.com/spreadsheets/d/1bdGW69IXu4Iy9Uinnfjd42xDtR6dkCOmX0B-T0qVxrI/edit?gid=0#gid=0"
ct_markers <- get_gsheet(url=url, sheet="LP")
ct_markers <- ct_markers %>%
  data.frame() %>%
  dplyr::filter(CellType %in% unique(SPICEY_ANNOTATED_COACC$cell_type))



# 1. Filter cell types ----
keep_cell_types <- c("Alpha", "Beta", "Acinar", "Delta", "Ductal")

atac_filtered <- lapply(atac, function(gr) {
  if (!"cell_type" %in% colnames(mcols(gr))) return(NULL)
  gr_filtered <- gr[mcols(gr)$cell_type %in% keep_cell_types]
  if (length(gr_filtered) > 0) return(gr_filtered) else return(NULL)
})
atac_filtered <- Filter(Negate(is.null), atac_filtered)
atac_filtered <- lapply(atac_filtered, function(x) as.data.frame(x))

rna_filtered <- lapply(rna, function(df) {
  if (!"cell_type" %in% colnames(df)) return(NULL)
  df_filtered <- df[df$cell_type %in% keep_cell_types, ]
  if (nrow(df_filtered) > 0) return(df_filtered) else return(NULL)
})
rna_filtered <- Filter(Negate(is.null), rna_filtered)



# 2. Filter features ----
## 2.1. Peaks
important_regions <- retsi %>%
  filter(gene_coacc %in% ct_markers$Gene,
         !is.na(region)) %>%
  pull(region) %>%
  unique()

random_regions <- sample_regions(atac_filtered, max_rows = 500) %>%
  bind_rows() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  filter(!is.na(region)) %>%
  pull(region) %>%
  unique()


atac_filtered <- lapply(atac_filtered, function(x) {
  x <- x %>%
    data.frame() %>%
    mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
    filter(region %in% important_regions |
           region %in% random_regions) %>%
    dplyr::select(-c(region, width, strand)) %>%
    regioneR::toGRanges()
  return(x)
})


## 2.2. Genes
important_genes <- ct_markers$Gene

random_genes <- sample_regions(rna_filtered, max_rows = 500) %>%
  bind_rows() %>%
  filter(!is.na(gene_id)) %>%
  pull(gene_id) %>%
  unique()

rna_filtered <- lapply(rna_filtered, function(x) {
  x <- x %>%
    data.frame() %>%
    filter(gene_id %in% important_genes |
           gene_id %in% random_genes)
  return(x)
})



# 3. Filter links present in atac ----
regions2keep <- unlist(GRangesList(atac_filtered)) %>%
  reduce() %>%
  data.frame() %>%
  mutate(region = paste0(seqnames, "-", start, "-", end)) %>%
  filter(!is.na(region)) %>%
  pull(region) %>%
  unique()


cicero_links_filtered2 <- cicero_links %>%
  data.frame() %>%
  filter(Peak1 %in% regions2keep |
         Peak2 %in% regions2keep,
         coaccess > 0.6) %>%
  mutate(Peak1 = as.factor(Peak1),
         Peak2 = as.factor(Peak2)) %>%
  dplyr::select(c(Peak1, Peak2, coaccess)) %>%
  as.data.frame()


cicero_clean <- clean_cicero_links(cicero_links_filtered2)



# Save them
rm(atac)
rm(rna)
rm(cicero_links)

atac <- atac_filtered
rna <- rna_filtered
cicero_links <- cicero_clean

use_data(atac, overwrite = TRUE, compress = "xz")
use_data(rna, overwrite = TRUE, compress = "xz")
use_data(cicero_links, overwrite = TRUE, compress = "xz")
