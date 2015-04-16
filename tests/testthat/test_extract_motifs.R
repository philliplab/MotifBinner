context('Extract Motifs - simple sequences')

test_that('The motif extractor works when motifs identifiers match perfectly', {
  seq_data <- DNAStringSet(c(
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'AACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = ""),
          paste(c(paste0(sample(c('A', 'C', 'G', 'T'), 20, replace=T), collapse=""),
                  'AACGAATTAA', 'ACGTACGT', 'CCAACCGCTC'), collapse = "")
                        ))
  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=5)
  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$matched_seq)[2]=='ACGTACGT', is_true())
  expect_that(width(extracted$matched_seq)[1] == 20, is_true())

  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
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

  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=0)
  expect_that(length(extracted$matched_seq) == 0, is_true())
  expect_that(length(extracted$unmatched_seq) == 2, is_true())
  expect_that(names(extracted$unmatched_seq)[1]=='seq_1', is_true())
  expect_that(names(extracted$unmatched_seq)[2]=='seq_2', is_true())

  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
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

  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=1)
  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$unmatched_seq)[1]=='s2', is_true())
  expect_that(length(extracted$unmatched_seq) == 1, is_true())
  expect_that(width(extracted$matched_seq) < width(extracted$unmatched_seq), is_true())
})

context('extract motifs - MiSeq like')

test_that('The motif extractor works with different sequence lengths', {
  params <- list(pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:5){ # do some repeats of each test
    for (len in c(100, 300, 400, 500, 600, 700, 1000)){
      c_params <- params
      c_params$seq_len <- len
      pid_search <- do.call(gen_pid_search_scenario, c_params)
      em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                           prefix = pid_search$prefix, 
                           motif_length = c_params$pid_len,
                           suffix = pid_search$suffix, 
                           max.mismatch = 4)
      expect_that(names(em$matched_seq) == pid_search$pid, is_true())
    }
  }
})

