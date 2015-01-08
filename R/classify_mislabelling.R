#' Classify a bin's sequences into 'in' and 'out'
#'
#' Using a specified method, split a bin into a list of two DNAStringSets. One
#' containing the sequences that we believe represents reads of the molecule of
#' interest and the other representing reads that were from a different
#' molecule.
#'
#' @value
#' A list with two elements 'in' and 'out'
#' @param bin The input bin as a single DNAStringSet.
#' @param technique A string selecting which technique to use for the
#' classification
#' @param params A list of parameters used by the specific
#' classification technique
#' @export

classify_bin <- function(bin, technique = 'random', params = list(n=0.2)){
}

#' Randomly classify a number of sequences to the 'out' group
#'
#' A very poor strategy included as a baseline
#' @inheritParams classify_bin
#' @export

classify_bin_random <- function(bin, params){
  return(1)
}
