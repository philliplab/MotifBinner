# Useful functions for applying functions to files or folders

#' @include align_sequences.R
#' @include bin_by_name.R
#' @include classification_checker.R
#' @include classify_mislabelling.R
#' @include consensus_string.R
#' @include datasets.R
#' @include extract_motifs.R
#' @include imports.R
#' @include process_bin.R
NULL

paramz <- list(file_name = '~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq',
               prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
                         suffix = "CTGAGCGTGTGGCAAGGCCC",
                         motif_length = 9,
                         max.mismatch = 5,
                         fixed = FALSE,
                         add_uniq_id = T,
                         classification_technique = 'absolute',
                         classification_params = list(threshold = 8/600, 
                                                      start_threshold = 8/600, 
                                                      max_sequences = 100),
                         alignment_technique = 'muscle',
                         alignment_params = list(),
                         consensus_technique = 'mostConsensusString',
                         consensus_params = list(),
                         remove_gaps = TRUE)

#' Processes a file into consensus bins
#' @param file_name The file name
#' @export

process_file <- function(file_name,
                         prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
                         suffix = "CTGAGCGTGTGGCAAGGCCC",
                         motif_length = 9,
                         max.mismatch = 5,
                         fixed = FALSE,
                         add_uniq_id = T,
                         classification_technique = 'absolute',
                         classification_params = list(threshold = 8/600, 
                                                      start_threshold = 8/600, 
                                                      max_sequences = 100),
                         alignment_technique = 'muscle',
                         alignment_params = list(),
                         consensus_technique = 'mostConsensusString',
                         consensus_params = list(),
                         remove_gaps = TRUE){
  report_dat <- list()
  rsf_dat <- list(file_name = file_name)
  seq_dat <- read_sequence_file(file_name)
  rsf_dat$seq_dat <- seq_dat
  report_dat$rsf_dat <- rsf_dat

  em_dat <- list(seq_dat = seq_dat,
                 prefix = prefix,
                 suffix = suffix,
                 motif_length = motif_length,
                 max.mismatch = max.mismatch,
                 fixed = fixed)
  motif_dat <- do.call(extract_motifs, em_dat)
  em_dat$motif_dat <- motif_dat
  report_dat$em_dat <- em_dat

  bbn_dat <- list(seq_dat = motif_dat$matched_seq,
                  add_uniq_id = add_uniq_id)
  bin_seqs <- do.call(bin_by_name, bbn_dat)
  bin_seqs <- randomize_list(bin_seqs)
  bbn_dat$bin_seqs <- bin_seqs
  report_dat$bbn_dat <- bbn_dat

  pb_dat <- list()
  pb_dat$classification_technique <- classification_technique
  pb_dat$classification_params <- classification_params 
  pb_dat$alignment_technique <- alignment_technique 
  pb_dat$alignment_params <- alignment_params 
  pb_dat$consensus_technique <- consensus_technique
  pb_dat$consensus_params <- consensus_params
  pb_dat$remove_gaps <- remove_gaps 
  pb_out <- list()
  pb_dat$pb_out <- NULL
  for (bin_name in seq_along(bin_seqs)){
    print(bin_name)
    pb_dat$seqs <- bin_seqs[[bin_name]]
    pb_out[[bin_name]] <- do.call(process_bin, pb_dat)
  } 
  #  Just process the ouput into a friendlier data structure
  consensuses <- DNAStringSet()
  for (i in seq_along(pb_out)){
    if (length(pb_out[[i]]$consensus) > 0){
      dss <- DNAStringSet(pb_out[[i]]$consensus[[1]])
      names(dss) <- names(pb_out[[i]]$consensus)
      consensuses <- c(consensuses, dss)
    }
  }
  pb_dat$consensuses <- consensuses
  pb_dat$pb_out <- pb_out
  report_dat$pb_dat <- pb_dat

  return(report_dat)
}

#' Provides some settings useful for testing
get_test_settings <- function(){
  return(list(file_name = '~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq',
              prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
              suffix = "CTGAGCGTGTGGCAAGGCCC",
              motif_length = 9,
              max.mismatch = 5,
              fixed = FALSE,
              add_uniq_id = T))
}

#' Reads in a sequence file given a file name.
#'
#' Detects if it is a FASTA or FASTQ file and reads it appropriately. Quality
#' data will be dropped.
#'
#' @param file_name The file name
#' @export

read_sequence_file <- function(file_name){
  if (grepl('fastq', file_name, ignore.case = TRUE)){
    x <- readFastq(file_name)
    x <- x@sread
  } else {
    x <- readDNAStringSet(file_name)
  }
  x <- gsub('\\+', 'N', x)
  x <- DNAStringSet(x)
  return(x)
}

process_all_bins <- function(bins){
  y <- randomize_list(bins)
#  sfInit(parallel = TRUE, cpus = 6)
#  sfLibrary(MotifBinner)
#  pbins <- sfLapply(y, process_bin, consensus_technique = 'mostConsensusString')
#  sfStop()
  pbins <- list()
  for (bin_name in seq_along(y)){
    print(bin_name)
    pbins[[bin_name]] <- process_bin(y[[bin_name]], consensus_technique = 'mostConsensusString')
  } 
  #  Just process the ouput into a friendlier data structure
  consensuses <- DNAStringSet()
  for (i in seq_along(pbins)){
    if (length(pbins[[i]]) > 0){
      dss <- DNAStringSet(pbins[[i]][[1]])
      names(dss) <- names(pbins[[i]])
      consensuses <- c(consensuses, dss)
    }
  }
  if (is.null(names(consensuses))){
    names(consensuses) <- paste('seq', 1:length(consensuses), sep = '_')
  }
  return(consensuses)
}

#' Bins a given FASTA file and outputs each bin as a seperate file
#'
#' @param file_name The file name
#' @param add_uniq_id If True, an integer will be appended to
#' the end of each sequence's name so that all identical sequences in a bin
#' gets the same number and sequences who are not identical will get different
#' numbers.
#' @param number_of_front_bases_to_discard The number of bases to remove from
#' the front of the sequence. The first few bases are part of the primer
#' sequence and needs to be trimmed off since they do not contain any extra
#' information. This varies between sequencing approaches, so the parameter
#' should be set with knowledge of the process. The primer can contain
#' degeneracies, so its better to chop it off earlier rather than later.
#' @param prefix See ?extract_motifs
#' @param suffix See ?extract_motifs
#' @param motif_length See ?extract_motifs
#' @param max.mismatch See ?extract_motifs
#' @param fixed See ?extract_motifs
#' @param write_files If this is a directory, the bins will be written to that
#' folder as FASTA files.
#' @export

bin_file <- function(file_name = "~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq", 
                     add_uniq_id = TRUE,
                     number_of_front_bases_to_discard = 28,
                     prefix = "CCAGCTGGTTATGCGATTCTMARGTG",
                     suffix = "CTGAGCGTGTGGCAAGGCCC",
                     motif_length = 9,
                     max.mismatch = 5,
                     fixed = FALSE,
                     write_files = FALSE
                     ){

  if (grepl('fastq', file_name)){
    x <- readFastq(file_name)
    x <- x@sread
  } else {
    x <- readDNAStringSet(file_name)
  }
  x <- gsub('+', 'N', x)

  x <- padAndClip(x, IRanges(start = number_of_front_bases_to_discard, 
                             end=width(x)), 
                  Lpadding.letter="+", Rpadding.letter="+")
  seq_data <- x
  y <- extract_motifs(seq_data, prefix, suffix, motif_length, max.mismatch, fixed)
  y <- y$matched_seq
  attr(y@ranges@NAMES, "names") <- NULL

  bin_seqs <- bin_by_name(y, add_uniq_id)

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


#' Randomizes the order of the items of a list
#'
#' @param x The list to randomize
#' @export

randomize_list <- function(x){
  if (is.null(names(x))){
    names(x) <- 1:length(x)
  }
  y <- list()
  for (i in sample(names(x), length(names(x)), replace = FALSE)){
    y <- c(y, x[[i]])
  }
  return(y)
}

#' Bins and constructs consensus sequences for an entire fastq file
#'
#' Running process_bin in parallel yields 3x improvement in execution speed on
#' my laptop (213 sec -> 64 sec)
#'
#' @param file_name The file name
#' @export

file_to_consensus <- function(file_name = "~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq",
                              ...){
  x <- bin_file(file_name, write_files = FALSE, ...)
  y <- randomize_list(x)
  sfInit(parallel = TRUE, cpus = 6)
  sfLibrary(MotifBinner)
  z <- sfLapply(y, process_bin)
  sfStop()
  # Just process the ouput into a friendlier data structure
  consensuses <- DNAStringSet()
  for (i in seq_along(z)){
    if (length(z[[i]]) > 0){
      dss <- DNAStringSet(z[[i]][[1]])
      names(dss) <- names(z[[i]])
      consensuses <- c(consensuses, dss)
    }
  }
  if (is.null(names(consensuses))){
    names(consensuses) <- paste('seq', 1:length(consensuses), sep = '_')
  }
  return(consensuses)
}

#' Reads a classified binned file and splits it into bins
#' @param file_name Name of the file to process
#' @param prefixes The prefixes that indicate the different bins
#' @export

read_classified_binned_file <- function(file_name, prefixes = c('src', 'out')){
  # file_name <- "/home/phillipl/projects/MotifBinner/data/sample_classified_bins/classified/Bin_00030_CTGGAACCT.FASTA"
  seq_dat <- readDNAStringSet(file_name)
  bins <- list()
  for (prefix in prefixes){
    bins[[prefix]] <- character(0)
    for (i in 1:length(seq_dat)){
      seq_name <- names(seq_dat)[i]
      if (length(grep(paste0('^', prefix), seq_name)) == 1){
        bins[[prefix]] <- c(bins[[prefix]], as.character(seq_dat[i]))
      }
    }
  }
  stopifnot(length(seq_dat) == sum(unlist(lapply(bins, length))))
  return(bins)
}

#' Reads in all files from a folder, assuming that they are classified binned
#' files
#' @param binned_folder The folder containing the binned classified files
#' @export

dput_classified_binned_folder <- function(binned_folder){
  # binned_folder <- "/home/phillipl/projects/MotifBinner/data/sample_classified_bins/binfiasco"
  files <- list.files(binned_folder)
  all_binned_files <- list()
  for (file_name in files){
    file_path <- file.path(binned_folder, file_name)
    bins <- read_classified_binned_file(file_path)
    all_binned_files[[file_name]] <- bins
  }
  return(all_binned_files)
}
