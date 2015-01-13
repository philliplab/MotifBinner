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
#' classification techniques
#' @export

classify_bin <- function(bin, technique = 'random', params = list(n=0.2)){
  if (is.list(bin)){
    if (all(sort(c('src', 'out')) == sort(names(bin)))){
      bin <- c(bin$src, bin$out)
    }
  }
  params[['bin']] <- bin
  classified <- FALSE
  if (technique == 'random'){
    classified <- do.call(classify_bin_random, params)
  }
  if (technique == 'infovar_balance'){
    classified <- do.call(classify_bin_infovar_balance, params)
  }
  return(classified)
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
  if (n_to_remove > 0) {
    to_remove <- sample(1:n_elements, n_to_remove)
    usrc <- ubin[-to_remove]
    uout <- ubin[to_remove]
  } else {
    usrc <- ubin
    uout <- DNAStringSet(NULL)
  }

  return(list('src' = bin[bin %in% usrc],
              'out' = bin[bin %in% uout]))
}

#' Finds the outliers in a bin by balancing the amount of information that is
#' lost with the variance reduction that is realized.
#'
#' This approach removes outliers from a bin by finding the sequence with the
#' highest average distance to all other sequences in the bin. It is then
#' removed. The ratio of the reduction in the average distance between the
#' sequences and the reduction in the information available is then computed.
#' The process will continue until either all data has been removed or the
#' ratio drops below some threshold
#'
#' @param bin The input bin as a single DNAStringSet.
#' @param threshold Outlier sequences are removed from the bin until the
#' distanct / information ratio drops below this threshold.
#' @export

classify_bin_infovar_balance <- function(bin, threshold){
  dists <- NULL
  bin_dists <- stringDist(bin)
  dmat <- as.matrix(bin_dists)
  row.names(dmat) <- 1:nrow(dmat)
  removed_sequences <- NULL
  orig_dmat <- dmat
  pdr_psr_rat <- 1000
  while((nrow(unique(dmat)) > 1) & (pdr_psr_rat > threshold)){
    dists <- c(dists, mean(dmat))
    dvec <- apply(dmat, 1, sum)
    max_indx <- which(dvec == max(dvec))
    new_dmat <- dmat[-max_indx, -max_indx]
    psr <- length(max_indx) / nrow(dmat)
    if (nrow(new_dmat) == 0){
      pdr <- pdr # Since all data was removed, set this to the amount of data
                 # that was left in the prev iteration
    } else {
      pdr <- (mean(dmat) - mean(new_dmat)) / mean(orig_dmat)
    }
    pdr_psr_rat <- pdr / psr

    if (pdr_psr_rat > threshold){
      removed_sequences <- c(removed_sequences, row.names(dmat)[max_indx])
      dmat <- new_dmat
    }
  }
  src_seq <- bin[as.integer(row.names(dmat))]
  out_seq <- bin[as.integer(removed_sequences)]
  return(list(src = src_seq,
              out = out_seq))
}
