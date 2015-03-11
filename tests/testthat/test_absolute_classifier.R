context('Absolute threshold based classifier')

test_that('The abs thres based classifier works', {
  bin <- DNAStringSet(c(rep('AAAAA', 5), rep('AAAAC', 10), 
                        rep('AAAAG', 2), rep('TTTTT', 1)))
  classified <- classify_absolute(bin, threshold=1)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 17, is_true())
  expect_that(length(classified$out) == 1, is_true())

})

test_that('it stops when no decision can be made about which sequence to remove', {
  bin <- DNAStringSet(c('AAAAA', 'CCCCC', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_absolute(bin, threshold=1)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 0, is_true())
  expect_that(length(classified$out) == 4, is_true())
})

test_that('all equally distant sequences are always removed simultaneously', {
  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_absolute(bin, threshold=1.001)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 4, is_true())
  expect_that(length(classified$out) == 0, is_true())

  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_bin_infovar_balance(bin, threshold=1)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 2, is_true())
  expect_that(length(classified$out) == 2, is_true())
})
