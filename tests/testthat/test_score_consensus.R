context('consensus checker')

test_that('The consensus checker works', {
  bin <- list(src = DNAStringSet(c(rep('AAAAA', 5))),
              out = DNAStringSet(NULL),
              true_consensus = DNAString('AAAAA'))
  pbin <- score_consensus(bin)
  expect_that(pbin$edit_dist, equals(0))
  expect_that(as.character(pbin$result[[1]]), equals('AAAAA'))

  bin <- list(src = DNAStringSet(c(rep('AAAAA', 5))),
              out = DNAStringSet(NULL),
              true_consensus = DNAString('AACAA'))
  pbin <- score_consensus(bin)
  expect_that(pbin$edit_dist, equals(1))
})


