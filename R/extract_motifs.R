#' Extracts motifs from a set of reads
#' @export

extract_motifs <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
                          fixed = FALSE){
  seq_data <- DNAStringSet(gsub("[^ACGT]", "+", seq_data))
  motif_n <- paste(rep("N", motif_length), collapse="")
  padded_motif <- DNAString(paste0(prefix, motif_n, suffix))
  matches <- vmatchPattern(padded_motif, 
                           seq_data, 
                           max.mismatch = max.mismatch, 
                           with.indels=FALSE, 
                           fixed = fixed)
  matching_seq <- sapply(matches, length)
  stopifnot(all(matching_seq < 2))
  matching_seq <- matching_seq != 0
  matched_seq <- seq_data[matching_seq]
  unmatched_seq <- seq_data[!matching_seq]
  matches <- unlist(matches)
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
