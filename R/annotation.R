#' Annotate peaks with nearest gene and distance to TSS
#'
#' @param peaks_gr GRanges object of peaks
#' @param genes_gr GRanges object of genes with gene symbols (external_gene_name)
#' @return GRanges annotated with distanceToTSS, annotation, nearestGeneSymbol
#' @importFrom ChIPseeker annotatePeak
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom plyranges mutate nearest promoters
#' @import regioneR
#' @export
annotate_peaks <- function(peaks_gr, genes_gr) {
  anno <- ChIPseeker::annotatePeak(peaks_gr,
                                   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   verbose = FALSE)

  peaks_gr$distanceToTSS <- data.frame(anno)$distanceToTSS
  peaks_gr$annotation <- ifelse(abs(peaks_gr$distanceToTSS) <= 2000, "Promoter", "Distal")

  peaks_gr <- peaks_gr %>%
    plyranges::mutate(nearestGeneSymbol = genes_gr$external_gene_name[nearest(peaks_gr, promoters(genes_gr, 1, 0))])

  names(peaks_gr) <- NULL
  return(peaks_gr)
}

#' Parse genomic region string "chr:start-end" to GRanges
#'
#' @param region_str Character string in format "chr:start-end"
#' @return GRanges object
#' @importFrom regioneR toGRanges
#' @export
parse_peak_string <- function(region_str) {
  regioneR::toGRanges(region_str)
}
