context('Most Frequent Classifier')

test_that('The most frequent classifier works', {
  bin <- DNAStringSet(c(rep('AA', 5), rep('CC', 10), 
                        rep('GG', 2), rep('TT', 1)))
  classified <- classify_bin_most_frequent(bin)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 10, is_true())
  expect_that(length(classified$out) == 8, is_true())

  classified <- classify_bin(bin, technique='most_frequent', params=list())
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 10, is_true())
  expect_that(length(classified$out) == 8, is_true())
})
