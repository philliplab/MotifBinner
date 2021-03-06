#' A wrapper for extract motifs parallel that will iteratively allow for more
#' mismatches in the prefix and suffix surrounding the pids.
#' @param seq_data The sequences whose motifs must be extracted
#' @param prefix The prefix that is used to identify the motif
#' @param suffix The suffix that is used to identify the motif
#' @param motif_length The length of the motif that forms the pid.
#' @param max.mismatch See ?vmatchPattern
#' @param fixed See ?vmatchPattern
#' @param ncpu The number of cores to use
#' @param job_size The number of sequences to group into a single job
#' @param max.suffix.chop The maximum amount whereby to shorten the suffix
#' defaults to 0.5*length(suffix)
#' @param max.mismatch_start What is the minimum number of mismatches to search
#' for?
#' @param verbose If set to true a progress will be printed after each loop
#' iteration. Also passed through to the extract_motifs_par
#' function.
#' @export

extract_motifs_iterative <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
                          fixed = FALSE, ncpu = 6, job_size = NULL, max.suffix.chop = NULL,
                          max.mismatch_start = 0, verbose = FALSE){
  if (is.null(max.suffix.chop)){
    max.suffix.chop <- trunc(nchar(suffix)/3)
  }
  print(max.suffix.chop)
  matched_seq <- DNAStringSet(NULL)
  unmatched_seq <- seq_data
  params <- list(prefix = prefix,
                 suffix = suffix,
                 motif_length = motif_length,
                 fixed = fixed,
                 ncpu = ncpu,
                 job_size = job_size,
                 verbose = verbose)
  start_time <- Sys.time()
  for (i in max.mismatch_start:max.mismatch){
    params$max.mismatch <- i
    params$seq_data <- unmatched_seq
    if (verbose) {print(params)}
    result <- do.call(extract_motifs_par, params)
    matched_seq <- c(matched_seq, result$matched_seq)
    unmatched_seq <- result$unmatched_seq
    if (verbose){
      print(c(length(seq_data), length(matched_seq), 
              length(unmatched_seq), round(Sys.time()-start_time, 2)))
    }
  }
  if (max.suffix.chop > 0){
    for (i in 0:max.suffix.chop){
      params$max.mismatch <- max.mismatch
      params$suffix <- substr(suffix, 1, nchar(suffix) - i)
      params$seq_data <- unmatched_seq
      if (verbose) {print(params)}
      result <- do.call(extract_motifs_par, params)
      matched_seq <- c(matched_seq, result$matched_seq)
      unmatched_seq <- result$unmatched_seq
      if (verbose){
        print(c(length(seq_data), length(matched_seq), 
                length(unmatched_seq), round(Sys.time()-start_time, 2)))
      }
    }
  }
  return(list(matched_seq = matched_seq,
              unmatched_seq = unmatched_seq))
}

#' A wrapper for extract motifs that will execute it in a parallel nature using
#' the specified number of cpus
#' @param seq_data The sequences whose motifs must be extracted
#' @param prefix The prefix that is used to identify the motif
#' @param suffix The suffix that is used to identify the motif
#' @param motif_length The length of the motif that forms the pid.
#' @param max.mismatch See ?vmatchPattern
#' @param fixed See ?vmatchPattern
#' @param ncpu The number of cores to use
#' @param job_size The number of sequences to group into a single job
#' @param verbose If set to true a progress report will be given via file
#' creation in the sessions temp dir
#' @export

extract_motifs_par <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
                          fixed = FALSE, ncpu = 6, job_size = NULL, verbose=FALSE){
  if (is.null(job_size)){
    job_size <- ceiling(sqrt(length(seq_data)))
  }
  if (job_size < 10){
    job_size <- 10
  }
  seq_sets <- list()
  for (i in 1:(ceiling(length(seq_data)/job_size))){
    seq_sets[[i]] <- seq_data[(((i-1)*job_size)+1):min(i*job_size, length(seq_data))]
  }

  params <- list(prefix = prefix,
                 suffix = suffix,
                 motif_length = motif_length,
                 max.mismatch = max.mismatch,
                 fixed = fixed)
  all_params <- list()
  registerDoMC(cores=ncpu)
  list_results <- foreach(i=seq_along(seq_sets)) %dopar% {
    if (verbose){
      random_file_name <- paste('motsearch_', i, '_', length(seq_sets), '.fasta', sep = '')
      tmpfile_name <- file.path(tempdir(), random_file_name)
      file.create(tmpfile_name)
    }
    params <- list(prefix = prefix,
                   suffix = suffix,
                   motif_length = motif_length,
                   max.mismatch = max.mismatch,
                   fixed = fixed,
                   seq_data = seq_sets[[i]])
    y <- do.call(extract_motifs, params)
    y
  }
  matched_seq <- DNAStringSet(NULL)
  unmatched_seq <- DNAStringSet(NULL)
  for (i in seq_along(list_results)){
    matched_seq <- c(matched_seq, list_results[[i]]$matched_seq)
    unmatched_seq <- c(unmatched_seq, list_results[[i]]$unmatched_seq)
  }
  return(list(matched_seq = matched_seq,
              unmatched_seq = unmatched_seq))
}

#' Extracts motifs from a set of reads
#' @param seq_data The sequences whose motifs must be extracted
#' @param prefix The prefix that is used to identify the motif
#' @param suffix The suffix that is used to identify the motif
#' @param motif_length The length of the motif that forms the pid.
#' @param max.mismatch See ?vmatchPattern
#' @param fixed See ?vmatchPattern
#' @export

extract_motifs <- function(seq_data, prefix, suffix, motif_length, max.mismatch = 5,
                          fixed = FALSE){
  seq_data <- clean_seq_data(seq_data)

  motif_n <- paste(rep("N", motif_length), collapse="")
  padded_motif <- DNAString(paste0(prefix, motif_n, suffix))
  pattern_length <- nchar(padded_motif)-4
  matches <- vmatchPattern(padded_motif, 
                           seq_data, 
                           max.mismatch = max.mismatch, 
                           with.indels=FALSE, 
                           fixed = fixed)

  last_match <- clean_matches(matches, pattern_length)
  matching_seq <- names(seq_data) %in% names(last_match)

  if (sum(matching_seq) == 0){
    return(list(matched_seq = DNAStringSet(NULL),
                unmatched_seq = seq_data))
  }

  matched_seq <- seq_data[matching_seq]
  unmatched_seq <- seq_data[!matching_seq]
  
  matches <- IRanges(start=unlist(lapply(last_match, start)),
                     end=unlist(lapply(last_match, end)),
                     names=names(last_match))

  motif_free_matched_seq <- remove_motifs(matches, prefix, suffix, matched_seq)
  
  return(list(matched_seq = motif_free_matched_seq,
              unmatched_seq = unmatched_seq))
}

#' Internal Function: Removes motifs from sequences given match information
#' 
#' @param matched A set of IRanges describing the matched motif locations
#' @param prefix The prefix that is used to identify the motif
#' @param suffix The suffix that is used to identify the motif
#' @param matched_seq The sequences in which there were motifs.

remove_motifs <- function(matches, prefix, suffix, matched_seq){
  shifted_matches <- matches
  start(shifted_matches) <- start(shifted_matches) + nchar(prefix)
  end(shifted_matches) <- end(shifted_matches) - nchar(suffix)
  motifs <- padAndClip(matched_seq, shifted_matches, Lpadding.letter="+", 
                       Rpadding.letter="+")
# old code used to debug invalid motifs that were picked up due to error in
# iterative search function
#  invalid_motifs <- grep("\\+", motifs)
#  if (length(invalid_motifs)>0){
#    random_file_name <- paste('badmot_', paste(sample(c(LETTERS, letters), 20), collapse=""), '.fasta', sep = '')
#    tmpfile_name <- file.path(tempdir(), random_file_name)
#    writeXStringSet(DNAStringSet(matched_seq[invalid_motifs]),
#                    tmpfile_name)
#  }
  end(matches) <- start(matches) - 1
  start(matches) <- 1
  motif_free_matched_seq <- padAndClip(matched_seq, matches, Lpadding.letter="+", 
                       Rpadding.letter="+")
  names(motif_free_matched_seq) <- motifs
  return(motif_free_matched_seq)
}

#' Internal Function: Clean sequence data for motif extraction
#' @param seq_data The sequences whose motifs must be extracted

clean_seq_data <- function(seq_data){
  if (is.null(names(seq_data))){
    names(seq_data) <- paste('seq', 1:length(seq_data), sep="_")
  }
  if (length(unique(names(seq_data))) != length(names(seq_data))){
    names(seq_data) <- paste(names(seq_data), 1:length(seq_data), sep = '_')
  }
  seq_data <- DNAStringSet(gsub("[^ACGT]", "+", seq_data))
  return(seq_data)
}

#' Internal Function: Clean matches to the motifs
#' @param matches The matches to the motifs which potentially contains multiple
#' @param pattern_length The minimum length that the match must be in order to
#' be kept

clean_matches <- function(matches, pattern_length){
  last_match <- list()
  for (i in seq_along(matches)){
    c_matches <- matches[[i]][width(matches[[i]]) >= pattern_length]
    nr <- length(c_matches)
    if (nr > 0){
      last_match[[names(matches)[i]]] <- c_matches[nr]
    }
  }
  return(last_match)
}
