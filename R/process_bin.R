# Goes from sequences binned together to a consensus sequence

#' Given a groups of sequences that were binned together, throw out oitliers,
#' align them and construct a consensus sequence
#' @param seqs The sequences that were binned together
#' @export

process_bin <- function(seqs){
  if (is.list(seqs)){
    if (all(sort(c('src', 'out')) == sort(names(seqs)))){
      seqs <- c(seqs$src, seqs$out)
    }
  }

  x <- classify_bin(seqs, technique = 'infovar_balance', params = list(threshold = 1))
  x <- x$src

  x <- align_sequences(x)

  x <- construct_consensus(x)
  return(x)
}
