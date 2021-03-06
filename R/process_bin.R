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
#' @param strip_uids Remove the unique identifiers from the sequence. It is not
#' intelligent. The names will be split on '_' and the first and last pieces
#' will be kept.
#' @export

process_bin <- function(seqs, classification_technique = 'infovar_balance',
                        classification_params = list(threshold = 1, 
                                                     start_threshold = 0.02, 
                                                     max_sequences = 100),
                        alignment_technique = 'muscle',
                        alignment_params = list(),
                        consensus_technique = 'Biostrings::consensusString',
                        consensus_params = list(),
                        remove_gaps = TRUE,
                        strip_uids = FALSE){
  if (is.list(seqs)){
    if (all(sort(c('src', 'out')) == sort(names(seqs)))){
      seqs <- c(seqs$src, seqs$out)
    }
  }
  result <- list()
  if (length(names(seqs)) != length(unique(names(seqs)))){
    warning('sequence names must be unique for process_bin - adding unique id')
    names(seqs) <- paste(names(seqs), 1:length(seqs), sep = '_')
  }
  x <- classify_bin(seqs, technique = classification_technique, 
                    params = classification_params)
  result$src <- x$src
  result$out <- x$out
  result$dmat <- x$dmat
  x <- x$src
  if (length(x) <= 1){
    result$consensus <- DNAStringSet(NULL)
  } else {
    x <- align_sequences(x, technique = alignment_technique,
                         params = alignment_params)
    result$alignment <- x
    x <- construct_consensus(x, technique = consensus_technique,
                             params = consensus_params)
    if (remove_gaps){
      x[[1]] <- gsub('-', '', x[[1]])
    }
    if (strip_uids){
      y <- strsplit(names(x), '_')[[1]]
      names(x) <- paste(y[1], y[length(y)], sep = '_')
    }
    result$consensus <- x
  }
  return(result)
}
