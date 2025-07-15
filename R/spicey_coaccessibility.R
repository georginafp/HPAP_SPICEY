#' Get Candidate Cis-Regulatory Elements (CCAN) for a cell type
#'
#' @param retsi_gr GRanges with RETSI scores and cell_type metadata
#' @param cell_type Character specifying cell type to filter
#' @param ret_si_threshold Numeric RETSI score threshold
#' @return GRanges of CCAN peaks for the cell type above threshold
#' @export
get_ccan <- function(retsi_gr, cell_type, ret_si_threshold = 0.25) {
  subset <- retsi_gr[retsi_gr$cell_type == cell_type & retsi_gr$RETSI > ret_si_threshold]
  return(subset)
}

#' Link peaks to target genes via proximity or other criteria
#'
#' @param peaks_gr GRanges of candidate regulatory elements
#' @param genes_gr GRanges of annotated genes
#' @param max_dist Maximum distance for linking (default 1e6)
#' @return Data frame of links with peak and gene info
#' @importFrom GenomicRanges distanceToNearest
#' @export
get_targets_links <- function(peaks_gr, genes_gr, max_dist = 1e6) {
  nearest_hits <- GenomicRanges::distanceToNearest(peaks_gr, genes_gr)
  links_df <- data.frame(
    peak_idx = queryHits(nearest_hits),
    gene_idx = subjectHits(nearest_hits),
    distance = mcols(nearest_hits)$distance
  )
  links_df <- links_df[links_df$distance <= max_dist, ]

  peak_df <- as.data.frame(peaks_gr)[links_df$peak_idx, ]
  gene_df <- as.data.frame(genes_gr)[links_df$gene_idx, ]

  links_df <- cbind(peak_df, gene_df)
  return(links_df)
}
