# Goes from sequences binned together to a consensus sequence

#' Given a groups of sequences that were binned together, throw out outliers,
#' align them and construct a consensus sequence
#' @param seqs The sequences that were binned together
#' @param classification_technique The technique to use to search for
#' mislabelled sequences
#' @param classification_params The parameters for the classification technique
#' @param alignment_technique The alignment technique to use
#' @param alignment_params The alignment parameters to use
#' @param consensus_technique The technique to use when generating the
#' consensus sequence
#' @param consensus_params The parameters for the consensus generator
#' @export

process_bin <- function(seqs, classification_technique = 'infovar_balance',
                        classification_params = list(threshold = 1, start_threshold = 0.02),
                        alignment_technique = 'muscle',
                        alignment_params = list(),
                        consensus_technique = 'Biostrings::consensusString',
                        consensus_params = list()){
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
    return(x)
  }
}
