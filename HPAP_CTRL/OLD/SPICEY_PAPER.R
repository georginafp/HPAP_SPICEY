##################
## Paper SPICEY ##
##################

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
merged_df <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA/ATAC_CTRL_SEURAT_WL.rds")
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
epsilon <- 1e-6
near_zero_thresh <- 1e-2  # Lower the near_zero_thresh to treat more values as near-zero

retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean RETSI values
    RETSI_clean = ifelse(is.na(RETSI) | RETSI < near_zero_thresh, epsilon, RETSI),

    # Check if all RETSI values are effectively 0
    all_near_zero = all(RETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all values are near-zero, use uniform; otherwise, normalize RETSI_clean
    RETSI_adj = ifelse(all_near_zero, 1 / n(), RETSI_clean),
    RETSI_sum = sum(RETSI_adj, na.rm = TRUE),
    RETSI_prob = RETSI_adj / RETSI_sum,

    # Final renormalization in case of numeric drift
    RETSI_prob = RETSI_prob / sum(RETSI_prob, na.rm = TRUE),

    # Entropy computation
    entropy_component = -RETSI_prob * log2(RETSI_prob),
    RETSI_prob_check = sum(RETSI_prob, na.rm = TRUE)  # Sanity check
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types),
    RETSI_prob_check = first(RETSI_prob_check)
  ) %>%
  ungroup()


retsi_atac <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "region") %>%
  mutate(RETSI = log1p(RETSI)) %>%
  regioneR::toGRanges()
saveRDS(retsi_atac, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds")




## RNA --------------------------------------------------------------------------
gr_list <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/DA/RNA_CTRL_SEURAT_WL.rds")


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
    # GETSI = (norm_log2FC * ifelse(p_val_adj <= 0.05, 1, 0)) / (N - 1)  # Set weight to 1 or 0 based on p_val_adj
  ) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand, N)) %>%
  as.data.frame() %>%
  regioneR::toGRanges()


### Entropy computation ----
epsilon <- 1e-6
near_zero_thresh <- 1e-2  # Lower the near_zero_thresh to treat more values as near-zero

retsi_entropy <- retsi_gr %>%
  as.data.frame() %>%
  group_by(region) %>%
  mutate(
    # Clean GETSI values
    GETSI_clean = ifelse(is.na(GETSI) | GETSI < near_zero_thresh, epsilon, GETSI),

    # Check if all GETSI values are effectively 0
    all_near_zero = all(GETSI_clean == epsilon)
  ) %>%
  mutate(
    # If all values are near-zero, use uniform; otherwise, normalize GETSI_clean
    GETSI_adj = ifelse(all_near_zero, 1 / n(), GETSI_clean),
    GETSI_sum = sum(GETSI_adj, na.rm = TRUE),
    GETSI_prob = GETSI_adj / GETSI_sum,

    # Final renormalization in case of numeric drift
    GETSI_prob = GETSI_prob / sum(GETSI_prob, na.rm = TRUE),

    # Entropy computation
    entropy_component = -GETSI_prob * log2(GETSI_prob),
    GETSI_prob_check = sum(GETSI_prob, na.rm = TRUE)  # Sanity check
  ) %>%
  summarise(
    entropy = sum(entropy_component, na.rm = TRUE),
    n_cell_types = n(),
    norm_entropy = entropy / log2(n_cell_types),
    GETSI_prob_check = first(GETSI_prob_check)
  ) %>%
  ungroup()

retsi_rna <- retsi_gr %>%
  data.frame() %>%
  left_join(retsi_entropy, by = "region") %>%
  mutate(GETSI = log1p(GETSI)) %>%
  regioneR::toGRanges()
saveRDS(retsi_rna, "/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds")



## SPICEY bind -----------------------------------------------------------------
### Nearest gene ---------------------------------------------------------------
retsi_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds")
retsi_rna <- retsi_rna %>%
  as.data.frame() %>%
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
  dplyr::rename(GETSI_entropy = norm_entropy) %>%
  distinct()

retsi_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds")
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


### Coaccessible links ---------------------------------------------------------
links_df <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/CTRL_LINKS.rds")
links_df <- links_df[links_df$coaccess >= 0.25,]

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


gr_rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_GETSI.rds") %>%
  data.frame() %>%
  dplyr::select(-c("width", "strand")) %>%
  regioneR::toGRanges()

gr_atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/HPAP_CTRL_RETSI.rds") %>%
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


get_targets_links <- function(links, retsi, getsi,
                              name_column_peaks = "region",
                              name_column_genes = "symbol") {

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
  retsi$genes_coacc <- lapply(retsi$CCAN, function(ccans) {
    if (is.null(ccans)) return(NA)
    unique(unlist(gene_ccan[as.character(ccans)]))
  })

  # Step 4: Unnest genes_coacc so each gene gets its own row
  retsi_df <- as_tibble(retsi)
  retsi_df_unnested <- retsi_df %>%
    unnest(genes_coacc, keep_empty = TRUE)

  # Step 5: Convert back to GRanges
  retsi_gr_unnested <- GRanges(
    seqnames = retsi_df_unnested$seqnames,
    ranges = IRanges(start = retsi_df_unnested$start, end = retsi_df_unnested$end),
    strand = retsi_df_unnested$strand
  )
  mcols(retsi_gr_unnested) <- retsi_df_unnested %>%
    select(-seqnames, -start, -end, -strand)

  return(retsi_gr_unnested)
}

SPICEY_ANNOTATED_COACC <- get_targets_links(links, gr_atac, gr_rna)


gr_rna <- gr_rna %>%
  data.frame() %>%
  dplyr::select(symbol, GETSI, cell_type, norm_entropy) %>%
  dplyr::rename(GETSI_entropy = norm_entropy) %>%
  distinct()

SPICEY_ANNOTATED_COACC <- SPICEY_ANNOTATED_COACC %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, annotation, nearestGeneSymbol, genes_coacc, cell_type, RETSI, norm_entropy) %>%
  dplyr::rename(RETSI_entropy = norm_entropy)



# Convert GRanges to tibbles
SPICEY_ANNOTATED_COACC <- as_tibble(SPICEY_ANNOTATED_COACC)
gr_rna <- as_tibble(gr_rna)


# Add GETSI
SPICEY_ANNOTATED_COACC <- SPICEY_ANNOTATED_COACC %>%
  left_join(gr_rna, by = c("nearestGeneSymbol" = "symbol", "cell_type")) %>%
  rename(
    GETSI_nearest = GETSI,
    GETSI_entropy_nearest = GETSI_entropy
  ) %>%
  left_join(gr_rna, by = c("genes_coacc" = "symbol", "cell_type")) %>%
  rename(
    GETSI_coacc = GETSI,
    GETSI_entropy_coacc = GETSI_entropy
  )

saveRDS(SPICEY_ANNOTATED_COACC, file="/homes/users/gfuentes/scratch/projects/spicey/HPAP_CTRL/data/CTRL_HPAP_SPICEY_ANNOTATED_COACC.rds")

