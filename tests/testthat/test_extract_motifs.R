context('Extract Motifs')

test_that('The motif extractor is working', {
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

  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=0)
  expect_that(names(extracted$matched_seq)[1]=='ACGTACGT', is_true())
  expect_that(names(extracted$matched_seq)[2]=='ACGTACGT', is_true())

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
  extracted <- extract_motifs(seq_data, prefix = 'AACGAATTAA', motif_length=8, 
                              suffix = 'CCAACCGCTC', max.mismatch=1)

})

