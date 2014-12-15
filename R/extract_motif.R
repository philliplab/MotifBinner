library(Biostrings)
library(ShortRead)

x <- readFastq("~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq")
x <- x@sread
seq_data <- x
prefix <- "CCAGCTGGTTATGCGATTCTMARGTG"
suffix <- "CTGAGCGTGTGGCAAGGCCC"
motif_length <- 9
max.mismatch <- 5
fixed = FALSE

extract_motif <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
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
  names(matched_seq) <- motifs
  return(list(matched_seq = matched_seq,
              unmatched_seq = unmatched_seq))
}
