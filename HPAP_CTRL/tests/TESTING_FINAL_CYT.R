library(cicero)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(ggside)
library(ggridges)
library(ggsignif)
library(ggplotify)
library(here)
library(InteractionSet)
library(patchwork)
library(pheatmap)
library(purrr)
library(regioneR)
library(RColorBrewer)
library(scales)
library(Seurat)
library(Signac)
library(stringr)
library(tibble)
library(tidyr)
library(tools)
library(tidyverse)
library(viridis)
library(S4Vectors)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("R/retsi.R")
source("R/getsi.R")
source("R/utils.R")


atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/CYT/data/CYT_FINAL_ATAC.rds")
rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/CYT/data/CYT_FINAL_RNA.rds")

retsi <- spicey_retsi(atac)
getsi <- spicey_getsi(rna)

# COACCESSIBILITY ----
# LINKS

##########
## NOTE ##
##########

# only D-D or P-D links!
# We removed P-P links!

links <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/CYT/data/CYT_LINKS.rds") |> dplyr::filter(coaccess > 0.3)

ccan <- generate_ccans(links,
                       coaccess_cutoff_override=0.25,
                       tolerance_digits = 2)
links <- links |>
  left_join(ccan |> dplyr::rename(CCAN1="CCAN"), by=c(Peak1="Peak")) |>
  left_join(ccan |> dplyr::rename(CCAN2="CCAN"), by=c(Peak2="Peak"))


annotate_unique_peaks <- function(peak_vec) {
  df_peaks <- tibble(peak = unique(peak_vec)) %>%
    tidyr::separate(peak,
                    into = c("seqnames", "start", "end"),
                    sep = "-",
                    convert = TRUE,
                    remove = FALSE)  %>%
    data.frame()

  gr <- regioneR::toGRanges(df_peaks %>% dplyr::select(seqnames, start, end))
  anno <- ChIPseeker::annotatePeak(gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, verbose = FALSE) %>%
    data.frame()

  df <- tibble(
    peak = df_peaks$peak,
    distanceToTSS = anno$distanceToTSS,
    annotation = ifelse(abs(anno$distanceToTSS) <= 2000, "Promoter", "Distal"))

  return(df)
}

peak1_anno <- annotate_unique_peaks(links$Peak1)
peak2_anno <- annotate_unique_peaks(links$Peak2)


# Join annotations back to original links data
links_annotated <- links %>%
  left_join(peak1_anno, by = c("Peak1" = "peak")) %>%
  rename(distanceToTSS1 = distanceToTSS, annotation1 = annotation) %>%
  left_join(peak2_anno, by = c("Peak2" = "peak")) %>%
  rename(distanceToTSS2 = distanceToTSS, annotation2 = annotation)

links_annotated <- links_annotated %>%
  filter(!(annotation1 == "Promoter" & annotation2 == "Promoter"),
         !(annotation1 == "Distal" & annotation2 == "Distal"))  ## 1160704


links_hpap <- GInteractions(
  StringToGRanges(links_annotated$Peak1),
  StringToGRanges(links_annotated$Peak2),
  coacc=links_annotated$coaccess,
  CCAN1=links_annotated$CCAN1,
  CCAN2=links_annotated$CCAN2
)


# Annotate with coaccessibility
get_ccan <- function(links, gr, name_column="region", split=c("ccan", "name")) {
  hits1 <- findOverlaps(links, gr, use.region="first")
  hits2 <- findOverlaps(links, gr, use.region="second")

  ccan_df <- rbind(
    data.frame(ccan=links[queryHits(hits1),]$CCAN1,
               name=mcols(gr[subjectHits(hits1),])[, name_column]),
    data.frame(ccan=links[queryHits(hits2),]$CCAN1,
               name=mcols(gr[subjectHits(hits2),])[, name_column])) |>
    unique() |>
    filter(!is.na(ccan), !is.na(name))

  if(split == "ccan") {
    final <- split(ccan_df$name, ccan_df$ccan)
  } else if (split == "name") {
    final <- split(ccan_df$ccan, ccan_df$name)
  }

  return(final)
}


get_targets_links <- function(links, re, proms, name_links="HPAP") {
  # Obtain CCANs for REs
  res <- get_ccan(links, re, name_column="region", split="name")
  res[setdiff(re$region, names(res))] <- NA
  re$CCAN <- res[re$region]

  # Obtain CCANs for gene promoters
  genes <- get_ccan(links, proms, name_column="symbol", split="ccan")
  genes[setdiff(unlist(res[!is.na(res)]), names(genes))] <- NA

  # Add information into RE mcols
  mcols(re)[[paste0("genes_", name_links)]] <- lapply(re$CCAN, function(x) unlist(genes[as.character(unlist(x))]))
  return(re)
}

re <- get_targets_links(links_hpap, retsi,
                        proms = get_promoters_protein_coding(TxDb.Hsapiens.UCSC.hg38.knownGene),
                        name_links="HPAP")
re <- tidyr::unnest(data.frame(re), cols = c(genes_HPAP))
re <- re %>% as.data.frame() %>% regioneR::toGRanges()


# ADD GETSI
getsi2add <- as.data.frame(getsi) |>
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) |>
  dplyr::rename(GETSI_entropy = norm_entropy)

result <- re |>
  data.frame() |>
  dplyr::rename(RETSI_entropy = norm_entropy) |>
  dplyr::right_join(getsi2add, by = c("genes_HPAP" = "symbol", "cell_type")) |>
  dplyr::select(seqnames, start, end, cell_type, annotation, distanceToTSS,
                genes_HPAP, region, RETSI, RETSI_entropy,
                GETSI, GETSI_entropy) %>%
  left_join(
    getsi %>%
      data.frame() %>%
      dplyr::select(symbol, seqnames, start, end) %>% distinct(),
    by = c("genes_HPAP" = "symbol")
  )  %>%
  mutate(
    seqnames = coalesce(seqnames.x, seqnames.y),
    start = coalesce(start.x, start.y),
    end = coalesce(end.x, end.y)
  ) %>%
  dplyr::select(-seqnames.x, -seqnames.y, -start.x, -start.y, -end.x, -end.y) %>%
  dplyr::select(seqnames, start, end, everything())

SPICEY_ANNOTATED_COACC <- result %>% data.frame() %>% regioneR::toGRanges() # 552230
# saveRDS(SPICEY_ANNOTATED_COACC, "/homes/users/gfuentes/scratch/projects/spicey_old/CYT/data/CYT_SPICEY_COACC.rds")
