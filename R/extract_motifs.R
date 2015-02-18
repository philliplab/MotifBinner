#' Extracts motifs from a set of reads
#' @param seq_data The sequences whose motifs must be extracted
#' @param prefix The prefix that is used to identify the motif
#' @param suffix The suffix that is used to identify the motif
#' @param max.mismatch See ?vmatchPattern
#' @param fixed See ?vmatchPattern
#' @export

extract_motifs <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
                          fixed = FALSE){
  if (is.null(names(seq_data))){
    names(seq_data) <- paste('seq', 1:length(seq_data), sep="_")
  }
  seq_data <- DNAStringSet(gsub("[^ACGT]", "+", seq_data))
  motif_n <- paste(rep("N", motif_length), collapse="")
  padded_motif <- DNAString(paste0(prefix, motif_n, suffix))
  matches <- vmatchPattern(padded_motif, 
                           seq_data, 
                           max.mismatch = max.mismatch, 
                           with.indels=FALSE, 
                           fixed = fixed)

  matching_seq <- sapply(matches, length)

  if (all(matching_seq == rep(0, length(matching_seq)))){
    return(list(matched_seq = DNAStringSet(NULL),
                unmatched_seq = seq_data))
  }

  #stopifnot(all(matching_seq < 2))
  matching_seq <- matching_seq != 0
  matched_seq <- seq_data[matching_seq]
  unmatched_seq <- seq_data[!matching_seq]
  
  last_match <- list()
  for (i in seq_along(matches)){
    nr <- length(matches[[i]])
    if (nr > 0){
      last_match[[names(matches)[i]]] <- matches[[i]][nr]
    }
  }

  matches <- IRanges(start=unlist(lapply(last_match, start)),
                     end=unlist(lapply(last_match, end)),
                     names=names(last_match))

  shifted_matches <- matches
  start(shifted_matches) <- start(shifted_matches) + nchar(prefix)
  end(shifted_matches) <- end(shifted_matches) - nchar(suffix)
  motifs <- padAndClip(matched_seq, shifted_matches, Lpadding.letter="+", 
                       Rpadding.letter="+")
  end(matches) <- start(matches) - 1
  start(matches) <- 1
  motif_free_matched_seq <- padAndClip(matched_seq, matches, Lpadding.letter="+", 
                       Rpadding.letter="+")
  names(motif_free_matched_seq) <- motifs
  return(list(matched_seq = motif_free_matched_seq,
              unmatched_seq = unmatched_seq))
}
