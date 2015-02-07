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

classify_bin <- function(bin, technique = 'random', params = list()){
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
  if (technique == 'most_frequent'){
    classified <- do.call(classify_bin_most_frequent, params)
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
#' sequences and the reduction in the information available (where the number
#' of DNA sequences in the bin is taken as the measure of the amount of
#' information) is then computed.  The process will continue until either all
#' data has been removed or the ratio drops below some threshold
#'
#' @param bin The input bin as a single DNAStringSet.
#' @param threshold Outlier sequences are removed from the bin until the
#' distanct / information ratio drops below this threshold.
#' @param  start_threshold Only start the classification if maximum distance in
#' the sample is greater then this number of bases normalized to the length of
#' the sequences. A value of 0.01 means that the procedure will only start if
#' the maximum distance between two sequences is greater then 1 if the
#' sequences is exactly 100 bases long. TODO further normalize for number of
#' sequences
#' @param max_sequences The maximum number of sequences to use for the
#' computation of the distance matrix. If more sequences than this is present,
#' then randomly select this many sequences and run the classification
#' algorithm on them.
#' @export

classify_bin_infovar_balance <- function(bin, threshold, start_threshold = 0, 
                                         max_sequences = 100){
  discarded <- DNAStringSet(NULL)
  if (length(bin) > max_sequences){
    picks <- sample(1:length(bin), max_sequences, replace = FALSE)
    discarded <- bin[-picks]
    bin <- bin[picks]
  }
  seq_length <- min(nchar(bin))
  dists <- NULL
  bin_dists <- stringDist(bin)
  if (max(bin_dists)/seq_length < start_threshold){
    return(list(src = bin,
                out = DNAStringSet(NULL)))
  }
  dmat <- as.matrix(bin_dists)
  row.names(dmat) <- 1:nrow(dmat)
  removed_sequences <- NULL
  orig_dmat <- dmat
  pdr_psr_rat <- 1000
  pdr <- 1000*1000
  nrow_uniq_dmat_gt_1 <- nrow(unique(dmat)) > 1
  while((nrow_uniq_dmat_gt_1) & (pdr_psr_rat > threshold)){
    dists <- c(dists, mean(dmat))
    dvec <- apply(dmat, 1, sum)
    max_indx <- which(dvec == max(dvec))
    new_dmat <- dmat[-max_indx, -max_indx]
    psr <- length(max_indx) / nrow(dmat)
    if (is.null(nrow(new_dmat))){
      pdr <- pdr
      nrow_uniq_dmat_gt_1 <- FALSE
    } else {
      if (nrow(new_dmat) == 0){
        pdr <- pdr # Since all data was removed, set this to the amount of data
                   # that was left in the prev iteration
      } else {
        pdr <- (mean(dmat) - mean(new_dmat)) / mean(orig_dmat)
      }
      nrow_uniq_dmat_gt_1 <- nrow(unique(new_dmat)) > 1
    }
    pdr_psr_rat <- pdr / psr

    if (pdr_psr_rat > threshold){
      removed_sequences <- c(removed_sequences, row.names(dmat)[max_indx])
      dmat <- new_dmat
    }
  }
  src_seq <- bin[as.integer(row.names(dmat))]
  out_seq <- c(bin[as.integer(removed_sequences)], discarded)
  return(list(src = src_seq,
              out = out_seq))
}

#' Classifies the sequences in a set by only choosing the most frequent
#' sequence.
#'
#' Given a set of sequences, pick the most frequently occurring sequence as the
#' only sequence that belongs to that bin and label all other sequences as
#' outliers.
#' @param bin The input bin as a single DNAStringSet
#' @export

classify_bin_most_frequent <- function(bin){
  tab <- table(bin)
  max_occ <- max(tab)
  src_seq <- names(tab)[which(tab == max_occ)[1]]
  src <- bin[bin == src_seq]
  out <- bin[bin != src_seq]
  return(list(src = src,
              out = out))
}
