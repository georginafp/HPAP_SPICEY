##..............................................................................
# BiomaRt gene data                                                         ####
##..............................................................................
#'
#' @return List with two elements:
#'         - `df`. Data.frame with gene information to use for the RNA-seq
#'         pipeline.
#'         - `gr`. GRanges with protein coding genes to use for the epigenome
#'         pipeline and analyses.
biomart_genes <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "https://www.ensembl.org",
                           path = "/biomart/martservice",
                           dataset = "hsapiens_gene_ensembl")
  genes <- biomaRt::getBM(attributes = c("chromosome_name",
                                         "start_position", "end_position", "strand",
                                         "ensembl_gene_id", "external_gene_name", "gene_biotype",
                                         "percentage_gene_gc_content"),
                          useCache = TRUE, mart = mart)
  message("Got it!")

  ## Make GRanges
  genes$strand[genes$strand == -1] <- "-"
  genes$strand[genes$strand == 1] <- "+"

  genes_gr <- regioneR::toGRanges(genes)
  strand(genes_gr) <- genes_gr$strand
  mcols(genes_gr) <- mcols(genes_gr)[, -1]
  genes_gr <- GenomeInfoDb::keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
  genes_gr <- GenomeInfoDb::dropSeqlevels(genes_gr, "MT", pruning.mode = "coarse")
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"

  genes_gr <- genes_gr %>%
    plyranges::filter(gene_biotype == "protein_coding")

  genes_gr <- genes_gr[,c("ensembl_gene_id", "external_gene_name")]

  ## Make list with both
  genes_list <- list(df=genes,
                     gr=genes_gr)
  return(genes_list)

}



# biomart_genes <- function() {
#   mart <- biomaRt::useEnsembl(biomart = "ensembl",
#                               dataset = "hsapiens_gene_ensembl",
#                               mirror = "useast")
#
#   genes <- biomaRt::getBM(attributes = c("chromosome_name",
#                                          "start_position", "end_position", "strand",
#                                          "ensembl_gene_id", "external_gene_name", "gene_biotype",
#                                          "percentage_gene_gc_content"),
#                           useCache = TRUE, mart = mart)
#   message("Got it!")
#
#   ## Make GRanges
#   genes$strand[genes$strand == -1] <- "-"
#   genes$strand[genes$strand == 1] <- "+"
#
#   genes_gr <- regioneR::toGRanges(genes)
#   strand(genes_gr) <- genes_gr$strand
#   mcols(genes_gr) <- mcols(genes_gr)[, -1]
#   genes_gr <- GenomeInfoDb::keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
#   genes_gr <- GenomeInfoDb::dropSeqlevels(genes_gr, "MT", pruning.mode = "coarse")
#   GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
#
#   genes_gr <- genes_gr %>%
#     plyranges::filter(gene_biotype == "protein_coding")
#
#   genes_gr <- genes_gr[,c("ensembl_gene_id", "external_gene_name")]
#
#   ## Make list with both
#   genes_list <- list(df=genes,
#                      gr=genes_gr)
#   return(genes_list)
# }


#' LiftOver coordinates from hg19 to hg38
#'
#' This function takes a `GRanges` object with coordinates based on the hg19
#' genome assembly and lifts them over to hg38 using a chain file.
#' Only uniquely mapped regions (regions mapping to exactly one location) are retained.
#'
#' @param x A `GRanges` object with hg19 coordinates.
#'
#' @return A `GRanges` object with lifted coordinates in the hg38 assembly,
#'         keeping only regions that map uniquely.
#'
#' @details
#' The function uses the `rtracklayer::liftOver` function and the hg19 to hg38
#' UCSC chain file. Non-uniquely mapping regions are discarded.
#'
#' @examples
#' \dontrun{
#'   library(GenomicRanges)
#'   gr_hg19 <- GRanges(seqnames = "chr1", ranges = IRanges(1000000, 1000100))
#'   gr_hg38 <- liftOver_hg19_hg38(gr_hg19)
#' }
#'
#' @export
liftOver_hg19_hg38 <- function(x) {
  hg19 <- rtracklayer::liftOver(x, chain=rtracklayer::import.chain("/homes/users/gfuentes/shared_data/refs/liftOver/hg19ToHg38.over.chain"))
  lens <- sapply(hg19, length)
  hg38_filt <- unlist(hg19[lens == 1])
  hg38 <- base::unique(hg38_filt)
  hg38
}




