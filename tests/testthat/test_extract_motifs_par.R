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
#  for (i in 1:3){
#    in_data <- DNAStringSet(NULL)
#    params <- list(seq_len = 500,
#                   pid_len = 9,
#                   prefix_len = 27,
#                   suffix_len = 18,
#                   prefix_snps = 1,
#                   suffix_snps = 0,
#                   suffix_chop = 1)
#    for (j in 1:90){
#      pid_search <- do.call(gen_pid_search_scenario, params)
#      in_data <- c(in_data, DNAStringSet(pid_search$seq_dat))
#    } 
#    params <- list(seq_len = 500,
#                   pid_len = 7,
#                   prefix_len = 27,
#                   suffix_len = 18,
#                   prefix_snps = 1,
#                   suffix_snps = 0,
#                   suffix_chop = 1)
#    for (j in 1:10){
#      pid_search <- do.call(gen_pid_search_scenario, params)
#      in_data <- c(in_data, DNAStringSet(pid_search$seq_dat))
#    } 
#    em <- extract_motifs(in_data, 
#                         prefix = wrong_prefix, 
#                         motif_length = 9,
#                         suffix = pid_search$suffix, 
#                         max.mismatch = 5)
#  }
})

