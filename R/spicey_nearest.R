#' Bind two GRanges by nearest neighbor, merging metadata columns
#'
#' @param query_gr GRanges to be joined (e.g. peaks)
#' @param subject_gr GRanges to join with (e.g. genes)
#' @param max_dist Maximum distance to consider a match (default = 1e6)
#' @param ignore_strand Logical, whether to ignore strand when finding nearest (default TRUE)
#' @return GRanges with metadata columns combined from both inputs,
#'         with subject metadata suffixed by ".subject" to avoid conflicts
#' @importFrom GenomicRanges nearest mcols
#' @export
bind_nearest <- function(query_gr, subject_gr, max_dist = 1e6, ignore_strand = TRUE) {

  # Find nearest index
  nearest_idx <- GenomicRanges::nearest(query_gr, subject_gr, ignore.strand = ignore_strand)

  # Calculate distances to nearest
  dist_vec <- GenomicRanges::distance(query_gr, subject_gr[nearest_idx])

  # Filter by max_dist threshold
  valid_hits <- which(dist_vec <= max_dist)

  # Subset query and subject GRanges by valid hits
  query_sub <- query_gr[valid_hits]
  subject_sub <- subject_gr[nearest_idx[valid_hits]]

  # Prepare metadata for merging
  mcols_query <- as.data.frame(mcols(query_sub))
  mcols_subject <- as.data.frame(mcols(subject_sub))

  # Add suffix to subject metadata colnames to avoid clashes
  colnames(mcols_subject) <- paste0(colnames(mcols_subject), ".subject")

  # Combine metadata
  combined_mcols <- cbind(mcols_query, mcols_subject)

  # Assign combined metadata to query_gr subset
  mcols(query_sub) <- combined_mcols

  # Return combined GRanges with metadata from both sets
  return(query_sub)
}
