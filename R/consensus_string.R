#' Constructs a consensus string using the specified technique
#' 
#' This function is specifically designed to work in the situation of binning
#' sequences based on PIDS, so it handles sequence names in a specific way. It
#' assumes that the name of the first sequence is the name that needs to be
#' given to the consensus sequence. Furthermore, it will append and underscore
#' and the number of sequences that were used to construct the consensus to the
#' end of this name.
#'
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
