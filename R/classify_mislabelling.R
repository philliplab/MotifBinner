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
  if (technique == 'absolute'){
    classified <- do.call(classify_absolute, params)
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
#' data has been removed or the ratio drops below some threshold. If the ratio
#' will drop be low the threshold if the next sequence(s) is removed, then the
#' process stops, so that the process will stop before the ratio goes under the
#' threshold.
#'
#' Both the reduction in the average distance to all other sequences and the
#' reduction in the amount of information available is computed as a percentage
#' relative to the original input data. The formula for the average reduction
#' in distances is (mean(new_dmat) - mean(prev_dmat))/mean(orig_dmat) where
#' new_dmat is the new distance matrix constructed from removing the next
#' sequence(s), pre_dmat is the distance matrix constructed in the previous
#' step and orig_dmat is the distance matrix constructed on the original bin
#' passed to the function. The percentage reduction in information available is
#' the number of sequences that will be removed in this step over the number of
#' sequences in the input data set.
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
  # NOTE: It will be incredably hard to refactor this function
  # unless you move to oo so that you can assign multiple variables in one
  # method call
  discarded <- DNAStringSet(NULL)
  if (length(bin) > max_sequences){
    picks <- sample(1:length(bin), max_sequences, replace = FALSE)
    discarded <- bin[-picks]
    bin <- bin[picks]
  }
  seq_length <- min(nchar(bin))
  bin_dists <- stringDist(bin)
  if (max(bin_dists)/seq_length < start_threshold){
    return(list(src = bin,
                out = DNAStringSet(NULL)))
  }
  dmat <- as.matrix(bin_dists)
  row.names(dmat) <- 1:nrow(dmat)
  removed_sequences <- NULL
  orig_dmat <- dmat
  pdr_psr_rat <- 1000 # percentage distance reduction percentage sequence reduction ratio
  pdr <- 1000*1000 # percentage distance reduction
  nrow_uniq_dmat_gt_1 <- nrow(unique(dmat)) > 1
  while((nrow_uniq_dmat_gt_1) & (pdr_psr_rat > threshold)){
    dvec <- apply(dmat, 1, sum)
    max_indx <- which(dvec == max(dvec))
    new_dmat <- dmat[-max_indx, -max_indx]
    psr <- length(max_indx) / nrow(dmat) # percentage sequence reduction
    if (is.null(nrow(new_dmat)) | (nrow(new_dmat) == 0)){
      pdr <- pdr
      nrow_uniq_dmat_gt_1 <- FALSE
    } else {
      pdr <- (mean(dmat) - mean(new_dmat)) / mean(orig_dmat)
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
              out = out_seq,
              dmat = orig_dmat))
}

#' Classifies the sequences in a set by only choosing the most frequent
#' sequence.
#'
#' Given a set of sequences, pick the most frequently occurring sequence as the
#' only sequence that belongs to that bin and label all other sequences as
#' outliers.
#' 
#' If there is more than one sequence that occurs the maximum number of times,
#' it will return all the sequences that occur that number of times.
#'
#' @param bin The input bin as a single DNAStringSet
#' @export

classify_bin_most_frequent <- function(bin){
  tab <- table(bin)
  max_occ <- max(tab)
  src_seq <- names(tab)[which(tab == max_occ)]
  src <- bin[bin %in% src_seq]
  out <- bin[!(bin %in% src_seq)]
  return(list(src = src,
              out = out))
}

#' Removes the most outlying sequences in a bin until the maximum distance
#' between any two sequences reaches a threshold
#'
#' Thresholds are expressed as the probability that any given letter is an
#' error. Thus, before applying the thresholds to the distances between the
#' sequences, they are doubled.
#'
#' The rationale for doubling the threshold is that a sequence has a read
#' error if there is an error in 1 of the bases the sequence sequence while the
#' distance between 2 sequences is 1 if there is a read error in any of the
#' bases of the two sequences under consideration. Obviously its more subtle
#' than this, see the benchmarking document for more details.
#'
#' Note that as the thresholds get bigger, the behaviour will get strange.
#' This is because it assumes that the likihood of the same mutation occuring
#' in the two sequences being very low. However, as the error rate gets much
#' higher, that assumption becomes invalid and strange things happen. This will
#' not be fixed, since this library is not designed to work in an environment
#' where the error rates are high.
#'
#' @param bin The input bin as a single DNAStringSet.
#' @param threshold Outlier sequences are removed from the bin until the
#' maximum distance between any two sequences drops below this threshold.
#' @param  start_threshold Only start the classification if the maximum
#' distance between and two sequences in the bin is greater than this.
#' @param max_sequences The maximum number of sequences to use for the
#' computation of the distance matrix. If more sequences than this is present,
#' then randomly select this many sequences and run the classification
#' algorithm on them.
#' @export

classify_absolute <- function(bin, threshold=0.01, start_threshold = 0.02, 
                              max_sequences = 100, max_iterations = 1000){
  threshold <- threshold*2
  start_threshold <- start_threshold*2
  discarded <- DNAStringSet(NULL)
  if (length(bin) > max_sequences){
    picks <- sample(1:length(bin), max_sequences, replace = FALSE)
    discarded <- bin[-picks]
    bin <- bin[picks]
  }
  seq_length <- min(nchar(bin))
  bin_dists <- stringDist(bin)
  dmat <- as.matrix(bin_dists)
  row.names(dmat) <- 1:nrow(dmat)
  removed_sequences <- NULL
  orig_dmat <- dmat
  if (max(bin_dists)/seq_length < start_threshold){
    return(list(src = bin,
                out = DNAStringSet(NULL),
                dmat = orig_dmat))
  }
  max_dist_below_threshold <- max(dmat)/seq_length < threshold
  counter <- 0
  while(!max_dist_below_threshold){
    counter <- counter + 1
    if (counter > max_iterations){
      stop('Max iterations reached')
    }
    dvec <- apply(dmat, 1, sum)
    max_indx <- which(dvec == max(dvec))
    new_dmat <- dmat[-max_indx, -max_indx]
    if (is.null(nrow(new_dmat)) ){
      max_dist_below_threshold <- TRUE
    } else if (nrow(new_dmat) == 0){
      max_dist_below_threshold <- TRUE
    } else {
      max_dist_below_threshold <- max(new_dmat)/seq_length < threshold
    }

    removed_sequences <- c(removed_sequences, row.names(dmat)[max_indx])
    dmat <- new_dmat
  }
  src_seq <- bin[as.integer(row.names(dmat))]
  out_seq <- c(bin[as.integer(removed_sequences)], discarded)
  return(list(src = src_seq,
              out = out_seq,
              dmat = orig_dmat))
}

