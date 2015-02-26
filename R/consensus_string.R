#' Constructs a consensus string using the specified technique
#' 
#' This function is specifically designed to work in the situation of binning
#' sequences based on PIDS, so it handles sequence names in a specific way. It
#' assumes that the name of the first sequence is the name that needs to be
#' given to the consensus sequence. Furthermore, it will append and underscore
#' and the number of sequences that were used to construct the consensus to the
#' end of this name.
#'
#' @param seqs The aligned sequences to construct a consensus of
#' @param technique The technique to use
#' @param params Extra parameters passed to the consensus string construction
#' technique
#' @export

construct_consensus <- function(seqs, technique = 'Biostrings::consensusString', params = list()){
  seq_name <- names(seqs)[1]
  seq_name <- paste(seq_name, length(seqs), sep = '_')
  if (technique == 'Biostrings::consensusString'){
    params$x <- seqs
    con_str <- do.call(consensusString, params)
    x <- list(seq_name = DNAString(con_str))
  } else if (technique == 'mostConsensusString'){
    params$seqs <- seqs
    con_str <- do.call(mostConsensusString, params)
    x <- list(seq_name = DNAString(con_str))
  } else if (technique == 'twoStepConsensusString'){
    params$seqs <- seqs
    con_str <- do.call(twoStepConsensusString, params)
    x <- list(seq_name = DNAString(con_str))
  }
  names(x) <- seq_name
  return(x)
}

#' A custom consensus string constructor that uses a two step threshold method
#' to construct consensus sequences.
#'
#' First, threshold_1 is applied to search for all positions where a letter
#' occur in that many or more sequences. (threshold_1 must be bigger than 0.5).
#' The letter that occurs at least that frequently will then get put in the
#' consensus sequence.
#'
#' For those positions where there is not a letter that meets threshold_1,
#' threshold_2 will then be used to remove all letters that occur in less than
#' threshold_2 sequences. Using the letters that remain, an ambiguity character
#' will be constructed.
#' 
#' @param seqs The aligned sequences to construct a consensus of
#' @param threshold The percentage of sequences that must have a letter in a
#' certain positon for it to be included in the consensus string
#' @export

twoStepConsensusString <- function(seqs, threshold_1 = 0.501, threshold_2 = 0.25){
  stopifnot(threshold_1 > 0.5)
  conm <- consensusMatrix(seqs)/length(seqs)
  src_mat <- matrix(row.names(conm), ncol = ncol(conm), nrow = nrow(conm))
  not_thres <- which(apply(conm>0.5, 2, max)==0)
  if (length(not_thres) > 0){
    old_school_cons <- consensusString(seqs, threshold = threshold_2)
    for (i in seq_along(not_thres)){
      let <- substr(old_school_cons, not_thres[i], not_thres[i])
      conm[which(row.names(conm) == let), not_thres[i]] <- 1
    }
  }
  cons <- paste(src_mat[conm>0.5], sep="", collapse="")
  return(DNAString(cons))
}

#' A custom consensus string constructor that just uses the letters that occurs
#' most frequently at each position to construct the consensus sequence
#'
#' In the case that two or more letter occurs the same number of times, an
#' ambiguity character is constructed from them
#' 
#' @param seqs The aligned sequences to construct a consensus of
#' @export

mostConsensusString <- function(seqs){
  conm <- consensusMatrix(seqs)/length(seqs)
  row_maxes <- apply(conm, 2, max)
  new_seq <- rep('+', ncol(conm))
  for (i in 1:ncol(conm)){
    lets <- sort(row.names(conm)[conm[,i] == row_maxes[i]])
    lets <- paste(lets, sep="", collapse="")
    amb_char <- names(IUPAC_CODE_MAP)[IUPAC_CODE_MAP == lets]
    new_seq[i] <- amb_char
  }
  return(DNAString(paste(new_seq, sep="", collapse="")))
}
