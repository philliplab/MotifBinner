context('Alignment')

test_that('The aligner works', {
  bin <- DNAStringSet(c(rep('AA', 5), rep('CC', 10), 
                        rep('GG', 2), rep('TT', 1)))
  classified <- classify_bin_random(bin, 0.25)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) > 0, is_true())
  expect_that(length(classified$out) > 0, is_true())

  classified <- classify_bin(bin, technique='random', params=list(n= 0.25))
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) > 0, is_true())
  expect_that(length(classified$out) > 0, is_true())

  classified <- classify_bin_random(bin, 0)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$out) == 0, is_true())
})

