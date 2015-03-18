context('Process Bin')

test_that('The bin processor works', {
  bin <- DNAStringSet(c(rep('AA', 5)))
  pbin <- process_bin(bin)

  out_names <- c('src', 'out', 'alignment', 'consensus')
  expect_that(all(out_names %in% names(pbin)), is_true())
  expect_that(as.character(pbin$consensus[[1]]) == 'AA', is_true())

  bin <- DNAStringSet(c(rep('AAAAAAAAAAAAAAAAAAA', 5),
                        'CCCCCCCCCCCCCCCCCCC'))
  pbin <- process_bin(bin)
  expect_that(as.character(pbin$consensus[[1]]) == 'AAAAAAAAAAAAAAAAAAA', is_true())

  bin <- DNAStringSet(c(rep('AAAAAAAAAAAAAAAAAAA', 5),
                        'CCCCCCCCCCCCCCCCCCC',
                        'AAAAAAAAAAAAAAAAAAG'))
  pbin <- process_bin(bin)
  expect_that(as.character(pbin$consensus[[1]]) == 'AAAAAAAAAAAAAAAAAAA', is_true())
})

test_that('gaps are inserted and removed appropriately by process_bin', {
  bin <- DNAStringSet(c(rep('AAAAAAAAAAAACCCCCCAAAAAAAAAAAAA', 8),
                        rep('AAAAAAAAAAAAAACCCCCCAAAAAAAAAAAAA', 4)))
  classification_params <- list(threshold = 1, 
                                start_threshold = 1, 
                                max_sequences = 100)
  consensus_technique = 'mostConsensusString'
  pbin <- process_bin(bin, classification_params = classification_params,
                      consensus_technique = consensus_technique)
  expect_that(nchar(pbin$consensus[[1]]), equals(31))
  pbin <- process_bin(bin, classification_params = classification_params,
                      consensus_technique = consensus_technique,
                      remove_gaps = FALSE)
  expect_that(nchar(pbin$consensus[[1]]), equals(33))
  gaps_inserted <- sum(strsplit(as.character(pbin$consensus[[1]]), split='')[[1]] == '-')
  expect_that(gaps_inserted, equals(2))
})
