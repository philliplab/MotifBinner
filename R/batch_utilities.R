# Useful functions for applying functions to files or folders

#' Bins a given FASTA file and outputs each bin as a seperate file
#'
#' @param file_name The file name
#' @param add_sequence_uniqueness_id If True, an integer will be appended to
#' the end of each sequence's name so that all identical sequences in a bin
#' gets the same number and sequences who are not identical will get different
#' numbers.
#' @param number_of_front_bases_to_discard The number of bases to remove from
#' the front of the sequence
#' @param prefix See ?extract_motifs
#' @param suffix See ?extract_motifs
#' @param motif_length See ?extract_motifs
#' @param max.mismatch See ?extract_motifs
#' @param fixed See ?extract_motifs
#' @param write_files If this is a directory, the bins will be written to that
#' folder as FASTA files.
#' @export

bin_file <- function(file_name = "~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq", 
                     add_sequence_uniqueness_id = TRUE,
                     number_of_front_bases_to_discard = 28,
                     prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
                     suffix = "CTGAGCGTGTGGCAAGGCCC",
                     motif_length = 9,
                     max.mismatch = 5,
                     fixed = FALSE,
                     write_files = FALSE
                     ){

  x <- readFastq(file_name)
  x <- x@sread
  x <- padAndClip(x, IRanges(start = number_of_front_bases_to_discard, 
                             end=width(x)), 
                  Lpadding.letter="+", Rpadding.letter="+")
  seq_data <- x
  y <- extract_motifs(seq_data, prefix, suffix, motif_length, max_mismatch, fixed)

  bin_seqs <- bin_by_name(y$matched_seq, add_sequence_uniqueness_id)

  if (write_files != FALSE){
    dir.create(write_files, showWarnings=FALSE)
    print('writing files')
    for (i in names(bin_seqs)){
      n_seqs <- sprintf("%05.0f", length(bin_seqs[[i]]))
      writeXStringSet(bin_seqs[[i]],
                      file.path(write_files, paste0("Bin_", n_seqs, 
                                                    "_", i, ".FASTA")))
    }
  }
  return(bin_seqs)
}

