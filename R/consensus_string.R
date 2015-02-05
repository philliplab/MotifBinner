#' Constructs a consensus string using the specified technique
#' @param seqs The aligned sequences to construct a consensus of
#' @param technique The technique to use
#' @param params Extra parameters passed to the consensus string construction
#' technique
#' @export

construct_consensus <- function(seqs, technique = 'Biostrings::consensusString', params = list()){
  seq_name <- names(seqs)[1]
  seq_name <- paste(seq_name, length(seqs), sep = '_')
  x <- list(seq_name = DNAString(consensusString(seqs, ambig = 'N')))
  names(x) <- seq_name
  return(x)
}
