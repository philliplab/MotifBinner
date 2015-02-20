context('Process Bin')

test_that('The bin processor works', {
  bin <- DNAStringSet(c(rep('AA', 5)))
  pbin <- process_bin(bin)
  expect_that(as.character(pbin[[1]]) == 'AA', is_true())

  bin <- DNAStringSet(c(rep('AAAAAAAAAAAAAAAAAAA', 5),
                        'CCCCCCCCCCCCCCCCCCC'))
  pbin <- process_bin(bin)
  expect_that(as.character(pbin[[1]]) == 'AAAAAAAAAAAAAAAAAAA', is_true())

  bin <- DNAStringSet(c(rep('AAAAAAAAAAAAAAAAAAA', 5),
                        'CCCCCCCCCCCCCCCCCCC',
                        'AAAAAAAAAAAAAAAAAAG'))
  pbin <- process_bin(bin)
  expect_that(as.character(pbin[[1]]) == 'AAAAAAAAAAAAAAAAAAA', is_true())
})

