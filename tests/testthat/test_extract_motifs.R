set.seed(1)

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
  for (i in 1:3){ # do some repeats of each test
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

test_that('The motif extractor works with different pid lengths', {
  params <- list(seq_len = 500,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:3){ # do some repeats of each test
    for (len in 5:15){
      c_params <- params
      c_params$pid_len <- len
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

test_that('The motif extractor works with different prefix lengths', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:3){ # do some repeats of each test
    for (len in seq(19, 31, by=2)){
      c_params <- params
      c_params$prefix_len <- len
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

test_that('The motif extractor works with different suffix lengths', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:3){ # do some repeats of each test
    for (len in 13:23){
      c_params <- params
      c_params$suffix_len <- len
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

test_that('The motif extractor works with different number of snps in the prefix', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:3){ # do some repeats of each test
    for (len in 0:4){
      c_params <- params
      c_params$prefix_snps <- len
      pid_search <- do.call(gen_pid_search_scenario, c_params)
      em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                           prefix = pid_search$prefix, 
                           motif_length = c_params$pid_len,
                           suffix = pid_search$suffix, 
                           max.mismatch = 5)
      expect_that(names(em$matched_seq) == pid_search$pid, is_true())
    }
  }
})

test_that('The motif extractor works with different number of snps in the suffix', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_chop = 1)
  for (i in 1:3){ # do some repeats of each test
    for (len in 0:3){
      c_params <- params
      c_params$suffix_snps <- len
      pid_search <- do.call(gen_pid_search_scenario, c_params)
      em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                           prefix = pid_search$prefix, 
                           motif_length = c_params$pid_len,
                           suffix = pid_search$suffix, 
                           max.mismatch = 5)
      expect_that(names(em$matched_seq) == pid_search$pid, is_true())
    }
  }
})

test_that('The motif extractor works with different number of chops from the suffix', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0)
  for (i in 1:3){ # do some repeats of each test
    for (len in 0:4){
      c_params <- params
      c_params$suffix_chop <- len
      pid_search <- do.call(gen_pid_search_scenario, c_params)
      em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                           prefix = pid_search$prefix, 
                           motif_length = c_params$pid_len,
                           suffix = pid_search$suffix, 
                           max.mismatch = 5)
      expect_that(names(em$matched_seq) == pid_search$pid, is_true())
    }
  }
})

context('Extract Motifs - indels in prefix and pid')

test_that('what happens with indels in pid', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:5){ # do some repeats of each test
    c_params <- params
    wrong_pid_len <- 8
    pid_search <- do.call(gen_pid_search_scenario, c_params)
    em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                         prefix = pid_search$prefix, 
                         motif_length = wrong_pid_len,
                         suffix = pid_search$suffix, 
                         max.mismatch = 5)
    expect_that(length(em$matched_seq) == 0, is_true())
  }

  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:5){ # do some repeats of each test
    c_params <- params
    wrong_pid_len <- 10
    pid_search <- do.call(gen_pid_search_scenario, c_params)
    em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                         prefix = pid_search$prefix, 
                         motif_length = wrong_pid_len,
                         suffix = pid_search$suffix, 
                         max.mismatch = 5)
    expect_that(length(em$matched_seq) == 0, is_true())
  }
})

test_that('what happens with indels in pid', {
  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:5){ # do some repeats of each test
    c_params <- params
    pid_search <- do.call(gen_pid_search_scenario, c_params)
    p <- pid_search$prefix
    wrong_prefix <- DNAStringSet(paste0(substr(p, 1, 11), 'A', substr(p, 12, 100)))
    em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                         prefix = wrong_prefix, 
                         motif_length = c_params$pid_len,
                         suffix = pid_search$suffix, 
                         max.mismatch = 5)
    expect_that(length(em$matched_seq) == 0, is_true())
  }

  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 1,
                 suffix_snps = 0,
                 suffix_chop = 1)
  for (i in 1:5){ # do some repeats of each test
    c_params <- params
    pid_search <- do.call(gen_pid_search_scenario, c_params)
    p <- pid_search$prefix
    wrong_prefix <- DNAStringSet(paste0(substr(p, 1, 10), substr(p, 12, 100)))
    em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                         prefix = wrong_prefix, 
                         motif_length = c_params$pid_len,
                         suffix = pid_search$suffix, 
                         max.mismatch = 5)
    expect_that(length(em$matched_seq) == 0, is_true())
  }

  params <- list(seq_len = 500,
                 pid_len = 9,
                 prefix_len = 27,
                 suffix_len = 18,
                 prefix_snps = 0,
                 suffix_snps = 0,
                 suffix_chop = 0)
  for (i in 1:5){ # do some repeats of each test
    c_params <- params
    pid_search <- do.call(gen_pid_search_scenario, c_params)
    p <- pid_search$prefix
    wrong_prefix <- DNAStringSet(paste0(substr(p, 1, 2), substr(p, 4, 100)))
    em <- extract_motifs(DNAStringSet(pid_search$seq_dat), 
                         prefix = wrong_prefix, 
                         motif_length = c_params$pid_len,
                         suffix = pid_search$suffix, 
                         max.mismatch = 5)
    expect_that(names(em$matched_seq) == pid_search$pid, is_true())
  }
})

