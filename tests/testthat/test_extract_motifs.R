context('Extract Motifs')

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

test_that('The motif extractor works with realistic prefixes and suffixes', {
  prefix <- "CAGYACAGTACAATGTACACATGGAAT" 
  suffix <- "CTGAGCGTGTGGCAAGGC" 
  pid <- 'AAAGGCAAA'
  motif_length <- nchar(pid)
  max.mismatch <- 0
  seq_data <- DNAStringSet(paste0(gen_seq(50), prefix, pid, suffix))
  seq_data <- randomize_ambig(seq_data[1])
  
  extracted <- extract_motifs(seq_data, prefix = prefix, 
                              motif_length = motif_length, 
                              suffix = suffix, max.mismatch = max.mismatch)

  
})

