#' Groups sequences together into bins based on their names
#' @export

bin_by_name <- function(seq_dat){
  bins <- list()
  for (motif in unique(names(seq_dat))){
    bins[[motif]] <- seq_dat[names(seq_dat) == motif]
  }
  return(bins)
}
