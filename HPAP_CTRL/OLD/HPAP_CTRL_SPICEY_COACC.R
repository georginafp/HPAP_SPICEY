library(GenomicRanges)
library(InteractionSet)
library(dplyr)

links <- readRDS("/homes/users/gfuentes/shared_data/projects-data/t1d-snps/_targets/objects/LINKS_HPAP")
links <- links[links$coacc >= 0.25,]

gr_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_GETSI.rds")
gr_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL_RETSI.rds")



# 1. Function to extract CCAN groups
get_ccan <- function(links, gr, name_column, split = c("ccan", "name")) {
  split <- match.arg(split)

  # Overlaps with both anchors
  hits1 <- findOverlaps(anchors(links, "first"), gr)
  hits2 <- findOverlaps(anchors(links, "second"), gr)

  ccan_df <- rbind(
    data.frame(ccan = mcols(links)$CCAN1[queryHits(hits1)],
               name = mcols(gr)[[name_column]][subjectHits(hits1)]),
    data.frame(ccan = mcols(links)$CCAN2[queryHits(hits2)],
               name = mcols(gr)[[name_column]][subjectHits(hits2)])
  ) |>
    unique() |>
    filter(!is.na(ccan), !is.na(name))

  if (split == "ccan") {
    split(ccan_df$name, ccan_df$ccan)
  } else {
    split(ccan_df$ccan, ccan_df$name)
  }
}


# 2. Main function to annotate REs with linked genes
get_targets_links <- function(links, retsi, getsi,
                              name_column_peaks = "region",
                              name_column_genes = "symbol",
                              name_links = "HPAP") {

  # Step 1: Link REs to CCAN
  re_ccan <- get_ccan(links, retsi,
                      name_column = name_column_peaks, split = "name")
  missing <- setdiff(mcols(retsi)[[name_column_peaks]], names(re_ccan))
  re_ccan[missing] <- NA
  retsi$CCAN <- re_ccan[mcols(retsi)[[name_column_peaks]]]

  # Step 2: Link Genes to CCAN
  gene_ccan <- get_ccan(links, getsi,
                        name_column = name_column_genes, split = "ccan")
  all_re_ccans <- unique(unlist(re_ccan[!is.na(re_ccan)]))
  gene_ccan[setdiff(all_re_ccans, names(gene_ccan))] <- NA

  # Step 3: Add gene list to retsi
  retsi$genes_multiome <- lapply(retsi$CCAN, function(ccans) {
    if (is.null(ccans)) return(NA)
    unique(unlist(gene_ccan[as.character(ccans)]))
  })

  return(retsi)
}




mcols(gr_atac) <- mcols(gr_atac) %>%
  data.frame() %>%
  dplyr::select(-c("width", "strand"))

SPICEY_ANNOTATED_COACC <- get_targets_links(links, gr_atac, gr_rna)
saveRDS(SPICEY_ANNOTATED_COACC, file="/homes/users/gfuentes/scratch/projects/spicey/SPICEY_ANNOTATED_COACC.rds")



## ----
spicey_anntoated_coacc <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/SPICEY_ANNOTATED_COACC.rds")
