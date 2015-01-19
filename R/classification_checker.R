#' Checks that the output from a classifier satisfies basic constraints
#'
#' Makes sure that the following holds:
#' \itemize{
#'  \item{Output from the classifier is a list with 'in' and 'out' elements}
#'  \item{Each of the in and out elements are DNAStringSets}
#'  \item{None of the sequences in $in is in $out}
#'  \item{None of the sequences in $out is in $in}
#'  \item{Check that all the input sequences are present in the output data at
#'  the right frequencies}
#'  \item{Check that no new sequences were introduced during the classification
#'  step}
#' }
#' 
#' @param bin The input bin that was classified
#' @param classied The result from a classifier
#' @export

check_classification <- function(bin, classified){
  if (!(class(classified) == 'list')){
    stop("classified must be a list")
  }
  if (!all(sort(names(classified)) == c("out", "src"))){
    stop("classified must contain and only contain the 'src' and 'out' elements")
  }
  if (class(classified$src)!='DNAStringSet'){
    stop("'src' element of classified must be a DNSStringSet")
  }
  if (class(classified$out)!='DNAStringSet'){
    stop("'out' element of classified must be a DNSStringSet")
  }
  if (any(classified$src %in% classified$out)){
    stop("None of the 'src' sequences may be in the 'out' sequences")
  }
  tbin <- table(bin)
  tclass <- table(c(classified$src, classified$out))
  for (seq_name in names(tbin)){
    if (tbin[[seq_name]] != tclass[[seq_name]]){
      stop('All input sequences must be present at the same level in the output')
    }
  }
  if (any(!(unique(names(tclass)) %in% unique(names(tbin))))){
    stop('Classification process may not introduce new sequences')
  }
  return(TRUE)
}

#' Given a test dataset and a classification strategy, apply the strategy and compute metrics
#'
#' @param test_bin The input test bin as a list of two DNAStringSets.
#' @param technique A string selecting which technique to use for the
#' classification
#' @param params A list of parameters used by the specific
#' classification techniques
#' @export

score_classification <- function(test_bin, technique, params){
  bin <- c(test_bin$src, test_bin$out)
  ptm <- proc.time()
  classified <- classify_bin(bin, technique, params)
  time_taken <- proc.time() - ptm
  time_taken <- time_taken['user.self']
  names(time_taken) <- NULL
  # Compute Sensitivity
  tp <- sum(classified$src %in% test_bin$src)
  total_src <- length(test_bin$src)
  if ((total_src == 0) & (tp == 0)){
    sn <- 1
  } else {
    sn <- tp/total_src
  }
  # Compute Specificity
  tn <- sum(classified$out %in% test_bin$out)
  total_out <- length(test_bin$out)
  if ((total_out == 0) & (tn == 0)){
    sp <- 1
  } else {
    sp <- tn / total_out
  }
  # Compute max distance
  if (length(unique(classified$src)) < 2){
    max_dist <- 0
  } else {
    max_dist <- max(stringDist(unique(classified$src)))
    seq_length <- max(width(classified$src))
    max_dist <- (max_dist / seq_length) * 100
  }
  input_n <- length(bin)
  output_n <- length(classified$src)
  true_n <- length(test_bin$src)
  return(list(sn = sn, 
              sp = sp,
              combined = ((sn - 1)^2 + (sp -1)^2)^(1/2),
              max_dist = max_dist,
              time_taken = time_taken,
              input_n = input_n,
              output_n = output_n,
              true_n = true_n))
}

#' Given a list of datasets and a classification strategy, apply the strategy
#' to all datasets and compute the metrics
#' @param test_bins The list of input test bins.
#' @param technique A string selecting which technique to use for the
#' classification
#' @param params A list of parameters used by the specific
#' classification techniques
#' @export

score_all_classifications <- function(test_bins, technique, params){
  metrics <- data.frame(name = character(0),
                        sn = numeric(0),
                        sp = numeric(0),
                        combined = numeric(0),
                        max_dist = numeric(0),
                        time_taken = numeric(0),
                        input_n = numeric(0),
                        output_n = numeric(0),
                        true_n = numeric(0))
  for (i in 1:length(test_bins)){
    test_bin <- test_bins[[i]]
    bin_name <- names(test_bins)[[i]]
    x <- score_classification(test_bin, technique, params)
    metrics <- rbind(metrics,
                     data.frame(name = bin_name,
                                sn = x$sn,
                                sp = x$sp,
                                combined = x$combined,
                                max_dist = x$max_dist,
                                time_taken = x$time_taken,
                                input_n = x$input_n,
                                output_n = x$output_n,
                                true_n = x$true_n))
  }
  metrics <- rbind(metrics,
                   data.frame(name = 'summary',
                              sn = mean(metrics$sn),
                              sp = mean(metrics$sp),
                              combined = mean(metrics$combined),
                              max_dist = NA_real_,
                              time_taken = mean(metrics$time_taken),
                              input_n = NA_real_,
                              output_n = NA_real_,
                              true_n = NA_real_))
  return(metrics)
}
