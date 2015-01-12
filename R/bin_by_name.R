#' Groups sequences together into bins based on their names
#' @param seq_dat The sequence data whose names indicate what bins they belong
#' to
#' @param add_uniq_id If True, an integer will be appended to
#' the end of each sequence's name so that all identical sequences in a bin
#' gets the same number and sequences who are not identical will get different
#' numbers. Another integer will be added to ensure that each sequence name is
#' unique
#' @export

bin_by_name <- function(seq_dat, add_uniq_id = FALSE){
  bins <- list()
  for (motif in unique(names(seq_dat))){
    cbin <- seq_dat[names(seq_dat) == motif]
    if (add_uniq_id){
      names(cbin) <- paste(names(cbin), match(cbin, unique(cbin)),
                           sep = "_")
      names(cbin) <- paste(names(cbin), 1:length(cbin), sep = "_")
      cbin <- cbin[order(names(cbin))]
    }
    bins[[motif]] <- cbin
  }
  return(bins)
}
