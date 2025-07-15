##..............##
## Paper SPICEY ##
##..............##

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
library(InteractionSet)
library(stringr)

## ATAC -------------------------------------------------------------------------
merged_df <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA_ATAC_HPAP-CTRL.rds")
merged_df <- unlist(GRangesList(merged_df))

retsi_gr <- merged_df %>%
  data.frame() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand)) %>%
  regioneR::toGRanges()

anno <- ChIPseeker::annotatePeak(retsi_gr,
                                 TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                 verbose=FALSE)

retsi_gr$distanceToTSS <- data.frame(anno)$distanceToTSS
retsi_gr$annotation <- "Distal"
retsi_gr$annotation[abs(retsi_gr$distanceToTSS)<=2e3] <- "Promoter"

## annotate peaks
genes <- readRDS("/homes/users/gfuentes/scratch/projects/coculture/data/_targets/objects/GENES")$gr
retsi_gr <- retsi_gr %>%
  plyranges::mutate(nearestGeneSymbol = genes$external_gene_name[nearest(retsi_gr, promoters(genes, 1, 0))])
names(retsi_gr) <- NULL


### RETSI computation ----
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = 2^(avg_FC),
         p_val_adj = p.adjust(p_val, method = "fdr")  # Adjust p-value using FDR
  ) %>%            # remove negative sign
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


### Entropy computation ----
# epsilon <- 1e-2
# low_thresh <- 0.4
#
# retsi_entropy <- retsi_gr %>%
#   as.data.frame() %>%
#   group_by(region) %>%
#   mutate(
#     # Replace values < 0.5 or NA with epsilon
#     RETSI_clean = ifelse(is.na(RETSI) | RETSI < low_thresh, epsilon, RETSI),
#     # RETSI_clean = RETSI + epsilon,  # smooth RETSI to avoid zeros
#
#     # Normalize RETSI values to probabilities
#     RETSI_sum = sum(RETSI, na.rm = TRUE),
#     RETSI_prob = RETSI / RETSI_sum,
#
#     # Final renormalization
#     RETSI_prob = RETSI_prob / sum(RETSI_prob, na.rm = TRUE),
#
#     # Entropy component
#     entropy_component = -RETSI_prob * log2(RETSI_prob)
#   ) %>%
#   summarise(
#     entropy = sum(entropy_component, na.rm = TRUE),
#     n_cell_types = n(),
#     norm_entropy = entropy / log2(n_cell_types),
#     RETSI_prob_check = sum(RETSI_prob, na.rm = TRUE)
#   ) %>%
#   ungroup()


retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean RETSI values
    # RETSI_clean = ifelse(is.na(RETSI) | RETSI < low_thresh, epsilon, RETSI),
    RETSI_sum = sum(RETSI, na.rm = TRUE),
    RETSI_prob = RETSI / RETSI_sum,

    # Final renormalization
    RETSI_prob = RETSI_prob / sum(RETSI_prob, na.rm = TRUE),

    # Entropy component
    entropy_component = -RETSI_prob * log2(RETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE)
  ) %>%
  # mutate(
  #   H_min = min(entropy, na.rm = TRUE),
  #   H_max = max(entropy, na.rm = TRUE),
  #   norm_entropy = (entropy - H_min) / (H_max - H_min)
  # ) %>%
  # dplyr::select(-c(H_min, H_max)) %>%
  # as.data.frame()
  mutate(
    norm_entropy = 1 - exp(-entropy)  # Exponential-based normalization
  ) %>%
  as.data.frame()


retsi_atac <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "region") %>%
  mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()
# saveRDS(retsi_atac, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds")




## RNA --------------------------------------------------------------------------
gr_list <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA_RNA_HPAP-CTRL.rds")


gr_df_list <- map(gr_list, ~ as.data.frame(.x))
source("/homes/users/gfuentes/scratch/projects/spicey/R/utils.R")
gr_list <- lapply(names(gr_df_list), function(name) {

  df <- gr_df_list[[name]]
  df$symbol <- rownames(df)
  df$ensembl_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys=df$symbol,
                                         column="ENSEMBL",
                                         keytype="SYMBOL",
                                         multiVals="first")

  df_annot <- dplyr::left_join(df,
                               biomart_genes()$df %>% dplyr::select(chromosome_name,
                                                                    start_position,,
                                                                    end_position,
                                                                    strand,
                                                                    ensembl_gene_id,
                                                                    gene_biotype),
                               by=c(ensembl_id="ensembl_gene_id"))

  if (nrow(df_annot) == 0) {
    message("Missing entry annotation for ", name)
    return(NULL)
  }

  df_annot %>%
    dplyr::filter(!is.na(chromosome_name),
                  !is.na(start_position),
                  !is.na(end_position),
                  gene_biotype == "protein_coding") %>%
    dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
    dplyr::select(chromosome_name,
                  start_position,
                  end_position,
                  strand,
                  everything()) %>%
    dplyr::mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
    makeGRangesFromDataFrame(seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field=c("start", "start_position"),
                             end.field=c("end", "stop", "end_position"),
                             strand.field="strand",
                             keep.extra.columns = T) %>%
    sort()
})


# Rename list with cell types
names(gr_list) <- names(gr_df_list)
gr_list <- Filter(Negate(is.null), gr_list) # Remove NULLs


# Bind together in one GRanges
retsi_gr <- imap(gr_list, ~ {
  mcols(.x)$cell_type <- .y  # .x = GRanges, .y = name
  .x
}) %>%
  purrr::reduce(c) #214896


### GETSI computation -----------------------------------------------------------
retsi_gr <- retsi_gr %>%
  as.data.frame() %>%
  mutate(region = paste0(seqnames, ":", start, "-", end)) %>%
  mutate(avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
         avg_FC = 2^(avg_FC),
         # avg_FC = ifelse(avg_FC < 0, 0, avg_FC)  # Clip negative values
         p_val_adj = p.adjust(p_val, method = "fdr")  # Adjust p-value using FDR
  ) %>%            # remove negative sign
  group_by(symbol) %>%
  mutate(
    N = n(),
    max_log2FC = max(avg_FC, na.rm = TRUE),
    max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)  # Avoid zero division
  ) %>%
  ungroup() %>%
  mutate(p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
         weight = -log10(p_val_adj)) %>%
  group_by(cell_type, symbol) %>%
  mutate(
    norm_log2FC = (avg_FC / max_log2FC),
    GETSI = (norm_log2FC * weight) / (N - 1)
    # GETSI = log1p(GETSI)
    # GETSI = (norm_log2FC * ifelse(p_val_adj <= 0.05, 1, 0)) / (N - 1)  # Set weight to 1 or 0 based on p_val_adj
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


### Entropy computation ----
# epsilon <- 1e-2
# low_thresh <- 0.4
#
# retsi_entropy <- retsi_gr %>%
#   as.data.frame() %>%
#   group_by(symbol) %>%
#   mutate(
#     # Replace values < 0.5 or NA with epsilon
#     GETSI_clean = ifelse(is.na(GETSI) | GETSI < low_thresh, epsilon, GETSI),
#     # GETSI_clean = GETSI + epsilon,  # smooth RETSI to avoid zeros
#
#     # Normalize GETSI values to probabilities
#     GETSI_sum = sum(GETSI, na.rm = TRUE),
#     GETSI_prob = GETSI / GETSI_sum,
#
#     # Final renormalization
#     GETSI_prob = GETSI_prob / sum(GETSI_prob, na.rm = TRUE),
#
#     # Entropy component
#     entropy_component = -GETSI_prob * log2(GETSI_prob)
#   ) %>%
#   summarise(
#     entropy = sum(entropy_component, na.rm = TRUE),
#     n_cell_types = n(),
#     norm_entropy = entropy / log2(n_cell_types),
#     RETSI_prob_check = sum(GETSI_prob, na.rm = TRUE)
#   ) %>%
#   ungroup()


retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(symbol) %>%
  mutate(
    # Clean RETSI values
    # RETSI_clean = ifelse(is.na(RETSI) | RETSI < low_thresh, epsilon, RETSI),
    GETSI_sum = sum(GETSI, na.rm = TRUE),
    GETSI_prob = GETSI / GETSI_sum,

    # Final renormalization
    GETSI_prob = GETSI_prob / sum(GETSI_prob, na.rm = TRUE),

    # Entropy component
    entropy_component = -GETSI_prob * log2(GETSI_prob)
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE)
  ) %>%
  # mutate(
  #   H_min = min(entropy, na.rm = TRUE),
  #   H_max = max(entropy, na.rm = TRUE),
  #   norm_entropy = (entropy - H_min) / (H_max - H_min)
  # ) %>%
  # dplyr::select(-c(H_min, H_max)) %>%
  # as.data.frame()
  mutate(
    norm_entropy = 1 - exp(-entropy)  # Exponential-based normalization
  ) %>%
  as.data.frame()


retsi_rna <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "symbol") %>%
  mutate(GETSI = log1p(GETSI)) %>%
  regioneR::toGRanges()
# saveRDS(retsi_rna, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds")



## SPICEY bind -----------------------------------------------------------------
### Nearest gene ---------------------------------------------------------------
# retsi_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds")
# retsi_rna <- retsi_rna %>%
#   as.data.frame() %>%
#   dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
#   dplyr::rename(GETSI_entropy = norm_entropy) %>%
#   distinct()
#
# # retsi_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds")
# retsi_atac <- retsi_atac %>%
#   as.data.frame() %>%
#   dplyr::select(seqnames, start, end, annotation, nearestGeneSymbol, cell_type, RETSI, norm_entropy) %>%
#   dplyr::rename(symbol = nearestGeneSymbol,
#                 RETSI_entropy = norm_entropy)
#
#
# retsi_gr_final <- retsi_atac %>%
#   left_join(retsi_rna %>%
#               dplyr::select(symbol, GETSI, cell_type, GETSI_entropy),  # Select relevant columns
#             by = c("symbol", "cell_type")) %>%
#   regioneR::toGRanges()


### Coaccessible links ---------------------------------------------------------
links_df <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/data/CTRL_LINKS.rds")
links_df <- links_df[links_df$coaccess > 0.6,]

# Filter out rows with valid peaks
links_df <- links_df %>%
  filter(!is.na(Peak1) & !is.na(Peak2))

# Function to convert "chrX-start-end" to GRanges
parse_peak <- function(peak_str) {
  parts <- str_split(peak_str, "-", simplify = TRUE)
  GRanges(seqnames = parts[, 1],
          ranges = IRanges(as.numeric(parts[, 2]),
                           as.numeric(parts[, 3])))
}

anchor1 <- parse_peak(links_df$Peak1)
anchor2 <- parse_peak(links_df$Peak2)

# Build GInteractions object
links <- GInteractions(anchor1, anchor2)
mcols(links)$coaccess <- links_df$coaccess
mcols(links)$CCAN1 <- links_df$CCAN1
mcols(links)$CCAN2 <- links_df$CCAN2


gr_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds") %>%
  data.frame() %>%
  dplyr::select(-c("width", "strand")) %>%
  regioneR::toGRanges()

gr_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds") %>%
  data.frame() %>%
  dplyr::select(-c("width", "strand")) %>%
  regioneR::toGRanges()


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


# get_targets_links <- function(links, retsi, getsi,
#                               name_column_peaks = "region",
#                               name_column_genes = "symbol") {
#
#   # Step 1: Link REs to CCAN
#   re_ccan <- get_ccan(links, retsi,
#                       name_column = name_column_peaks, split = "name")
#
#   missing <- setdiff(mcols(retsi)[[name_column_peaks]], names(re_ccan))
#   re_ccan[missing] <- NA
#   retsi$CCAN <- re_ccan[mcols(retsi)[[name_column_peaks]]]
#
#   # Step 2: Link Genes to CCAN
#   gene_ccan <- get_ccan(links, getsi,
#                         name_column = name_column_genes, split = "ccan")
#   all_re_ccans <- unique(unlist(re_ccan[!is.na(re_ccan)]))
#   gene_ccan[setdiff(all_re_ccans, names(gene_ccan))] <- NA
#
#   # Step 3: Add gene list to retsi
#   retsi$genes_coacc <- lapply(retsi$CCAN, function(ccans) {
#     if (is.null(ccans)) return(NA)
#     unique(unlist(gene_ccan[as.character(ccans)]))
#   })
#
#   # Step 4: Unnest genes_coacc so each gene gets its own row
#   retsi_df <- as_tibble(retsi)
#   retsi_df_unnested <- retsi_df %>%
#     unnest(genes_coacc, keep_empty = TRUE)
#
#   # Step 5: Convert back to GRanges
#   retsi_gr_unnested <- GRanges(
#     seqnames = retsi_df_unnested$seqnames,
#     ranges = IRanges(start = retsi_df_unnested$start, end = retsi_df_unnested$end),
#     strand = retsi_df_unnested$strand
#   )
#   mcols(retsi_gr_unnested) <- retsi_df_unnested %>%
#     select(-seqnames, -start, -end, -strand)
#
#   return(retsi_gr_unnested)
# }


get_targets_links <- function(links, retsi, getsi,
                              name_column_peaks = "region",
                              name_column_genes = "symbol") {

  # Step 1: Link REs to CCAN
  re_names <- mcols(retsi)[[name_column_peaks]]
  re_ccan <- get_ccan(links, retsi, name_column = name_column_peaks, split = "name")
  re_ccan_vec <- re_ccan[re_names]
  re_ccan_vec[is.na(re_ccan_vec)] <- NA
  retsi$CCAN <- re_ccan_vec

  # Step 2: Link Genes to CCAN
  gene_ccan <- get_ccan(links, getsi, name_column = name_column_genes, split = "ccan")
  used_ccans <- unique(unlist(re_ccan_vec, use.names = FALSE))
  gene_ccan <- gene_ccan[names(gene_ccan) %in% used_ccans]

  # Step 3: Add gene list to retsi
  retsi$genes_coacc <- lapply(retsi$CCAN, function(ccans) {
    if (length(ccans) == 0 || all(is.na(ccans))) return(NA_character_)
    unique(unlist(gene_ccan[as.character(ccans)], use.names = FALSE))
  })

  # Step 4: Unnest genes_coacc to one gene per row, preserving GRanges structure
  keep <- lengths(retsi$genes_coacc) > 0
  gr_list <- rep(retsi[keep], lengths(retsi$genes_coacc[keep]))
  mcols(gr_list)$genes_coacc <- unlist(retsi$genes_coacc[keep], use.names = FALSE)

  return(gr_list)
}


SPICEY_ANNOTATED_COACC <- get_targets_links(links, gr_atac, gr_rna)


gr_rna <- gr_rna %>%
  data.frame() %>%
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
  dplyr::rename(GETSI_entropy = norm_entropy) %>%
  distinct()

SPICEY_ANNOTATED_COACC <- SPICEY_ANNOTATED_COACC %>%
  data.frame() %>%
  dplyr::select(seqnames, start, end, annotation, genes_coacc, cell_type, RETSI, norm_entropy) %>%
  dplyr::rename(RETSI_entropy = norm_entropy)



# Convert GRanges to tibbles
SPICEY_ANNOTATED_COACC <- as_tibble(SPICEY_ANNOTATED_COACC)
gr_rna <- as_tibble(gr_rna)


# Add GETSI
SPICEY_ANNOTATED_COACC <- SPICEY_ANNOTATED_COACC %>%
  dplyr::left_join(gr_rna, by = c("genes_coacc" = "symbol", "cell_type")) %>%
  dplyr::rename(
    GETSI_coacc = GETSI)

# saveRDS(SPICEY_ANNOTATED_COACC, file="/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/CTRL_HPAP_SPICEY_ANNOTATED_COACC.rds")

