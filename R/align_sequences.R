# Aligns sequences using external tools
# A convenient wrapper is provided for benchmarking purposes.

#' A wrapper for a number of alignment tools
#'
#' @param seqs The sequences to align
#' @param technique Which alignment program to use
#' @param params A list of parameters to pass to the alignment program
#' @export

align_sequences <- function(seqs, technique = 'muscle', params = list()){
  if (is.null(names(seqs))){
    names(seqs) <- paste('seq', 1:length(seqs), sep='_')
  }
  seqs <- as.DNAbin(seqs)
  aligned_seqs <- muscle(seqs)
  seqs <- DNAStringSet(toupper(apply(as.character(aligned_seqs), 1, paste0, collapse="")))
  return(seqs)
}
