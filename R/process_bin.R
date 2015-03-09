#' Given a groups of sequences that were binned together, throw out outliers,
#' align them and construct a consensus sequence
#'
#' The output will be a list with a DNAString. The name of the DNAString in the
#' list will be the bin name and will contain information about how many
#' sequences were used to construct the consensus. See construct_consensus.
#' @param seqs The sequences that were binned together
#' @param classification_technique The technique to use to search for
#' mislabelled sequences
#' @param classification_params The parameters for the classification technique
#' @param alignment_technique The alignment technique to use
#' @param alignment_params The alignment parameters to use
#' @param consensus_technique The technique to use when generating the
#' consensus sequence
#' @param consensus_params The parameters for the consensus generator
#' @param remove_gaps If set to TRUE (the default, then gaps will be removed
#' from the consensus sequences)
#' @export

process_bin <- function(seqs, classification_technique = 'infovar_balance',
                        classification_params = list(threshold = 1, 
                                                     start_threshold = 0.02, 
                                                     max_sequences = 100),
                        alignment_technique = 'muscle',
                        alignment_params = list(),
                        consensus_technique = 'Biostrings::consensusString',
                        consensus_params = list(),
                        remove_gaps = TRUE){
  if (is.list(seqs)){
    if (all(sort(c('src', 'out')) == sort(names(seqs)))){
      seqs <- c(seqs$src, seqs$out)
    }
  }
  x <- classify_bin(seqs, technique = classification_technique, 
                    params = classification_params)
  x <- x$src
  if (length(x) <= 1){
    return(DNAStringSet(NULL))
  } else {
    x <- align_sequences(x, technique = alignment_technique,
                         params = alignment_params)
    x <- construct_consensus(x, technique = consensus_technique,
                             params = consensus_params)
    if (remove_gaps){
      x[[1]] <- gsub('-', '', x[[1]])
    }
    return(x)
  }
}
