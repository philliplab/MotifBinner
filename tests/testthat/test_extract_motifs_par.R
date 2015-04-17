context('Extract Motifs Parallel - Simple')

test_that('The motif extractor works when motifs identifiers match perfectly', {
  seq_data <- DNAStringSet(c(
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'AACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = ""),
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'AACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = "")
                        ))
  extracted <- extract_motifs_par(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=5)
  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$matched_seq)[2]=='ACGTACGT', is_true())
  expect_that(width(extracted$matched_seq)[1] == 20, is_true())

  extracted <- extract_motifs_par(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=0)
  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$matched_seq)[2]=='ACGTACGT', is_true())
})

test_that('The motif extractor works when motif identifiers have errors', {

  seq_data <- DNAStringSet(c(
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'CACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = ""),
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'CACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = "")
                        ))
  extracted <- extract_motifs_par(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=0)

  expect_that(length(extracted$matched_seq) == 0, is_true())
  expect_that(length(extracted$unmatched_seq) == 2, is_true())
  expect_that(names(extracted$unmatched_seq)[1]=='seq_1', is_true())
  expect_that(names(extracted$unmatched_seq)[2]=='seq_2', is_true())

  extracted <- extract_motifs_par(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=1)

  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$matched_seq)[2]=='ACGTACGT', is_true())
  expect_that(length(extracted$unmatched_seq) == 0, is_true())
  
  seq_data <- DNAStringSet(c(
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'CACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = ""),
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'CTCGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = "")
                        ))
  names(seq_data) <- c('s1', 's2')
  extracted <- extract_motifs_par(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=1)

  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$unmatched_seq)[1]=='s2', is_true())
  expect_that(length(extracted$unmatched_seq) == 1, is_true())
  expect_that(width(extracted$matched_seq) < width(extracted$unmatched_seq), is_true())
})

context('Extract Motifs Parallel - datasets')

test_that('results between single core and parallel versions of motif finder matches', {
  prefix <- "CCAGCTGGTTATGCGATTCTCAGGTG"
  suffix <- "CTGAGCGTGTGGCAAGGCCC"
  for (i in 1:4){
    in_data <- DNAStringSet(NULL)
    params <- list(seq_len = 500,
                   pid_len = 9,
                   prefix_len = prefix,
                   suffix_len = suffix,
                   prefix_snps = 1,
                   suffix_snps = 0,
                   suffix_chop = 1)
    for (j in 1:90){
      pid_search <- do.call(gen_pid_search_scenario, params)
      in_data <- c(in_data, DNAStringSet(pid_search$seq_dat))
    } 
    params$pid_len <- 7
    for (j in 1:10){
      pid_search <- do.call(gen_pid_search_scenario, params)
      in_data <- c(in_data, DNAStringSet(pid_search$seq_dat))
    } 
    em <- extract_motifs(in_data, 
                         prefix = prefix, 
                         motif_length = 9,
                         suffix = suffix, 
                         max.mismatch = 5)
    em_par <- extract_motifs_par(in_data, 
                                 prefix = prefix, 
                                 motif_length = 9,
                                 suffix = suffix, 
                                 max.mismatch = 5)
    expect_that(all(names(em$matched) %in% names(em_par$matched)), is_true())
    expect_that(all(names(em_par$matched) %in% names(em$matched)), is_true())

    expect_that(all(em$matched %in% em_par$matched), is_true())
    expect_that(all(em_par$matched %in% em$matched), is_true())
  }
})

