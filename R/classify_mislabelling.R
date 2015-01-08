#' Classify a bin's sequences into 'src' and 'out'
#'
#' Using a specified method, split a bin into a list of two DNAStringSets. One
#' containing the sequences that we believe represents reads of the molecule of
#' interest and the other representing reads that were from a different
#' molecule.
#'
#' @return
#' A list with two elements 'src' and 'out'
#' @param bin The input bin as a single DNAStringSet.
#' @param technique A string selecting which technique to use for the
#' classification
#' @param params A list of parameters used by the specific
#' classification technique
#' @export

classify_bin <- function(bin, technique = 'random', params = list(n=0.2)){
}

#' Randomly classify a number of unique sequences to the 'out' group
#'
#' A very poor strategy included as a baseline. 
#'
#' First all the unique sequences in a bin is found. Then a percentage of then
#' is assigned to the unique-out group. Then all sequences in the original bin
#' that match sequences in the unique-out group is moved to the out bin.
#'
#' @param bin The input bin as a single DNAStringSet.
#' @param n The proportion of the unique sequences in the bin to remove.
#' @export

classify_bin_random <- function(bin, n){
  ubin <- unique(bin)
  n_elements <- length(ubin)
  n_to_remove <- trunc(n * n_elements)
  to_remove <- sample(1:n_elements, n_to_remove)
  usrc <- ubin[-to_remove]
  uout <- ubin[to_remove]

  return(list('src' = bin[bin %in% usrc],
              'out' = bin[bin %in% uout]))
}
