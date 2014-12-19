#' Groups sequences together into bins based on their names
#' @export

bin_by_name <- function(seq_dat){
  bins <- list()
  for (motif in unique(names(seq_dat))){
    bins[[motif]] <- seq_dat[names(seq_dat) == motif]
  }
  return(bins)
}

#' Checks a bin for contamination
#' @export

bin_contaminated <- function(bin){
  # Draw kinky diagrams
  # Put avg dist on y
  # number of sequences on x
  # and systematically throw out only the most outlierish sequence
  if (length(bin) < 2){
    return(FALSE)
  } else {
    str_dists <- stringDist(bin)
    return(list(max = max(str_dists),
                avg = mean(str_dists),
                min = min(str_dists)))
  }
}
