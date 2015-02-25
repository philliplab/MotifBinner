context('Infomation Variance Balance Classifier')

test_that('The infovar balance classifier works', {
  bin <- DNAStringSet(c(rep('AAAAA', 5), rep('AAAAC', 10), 
                        rep('AAAAG', 2), rep('TTTTT', 1)))
  classified <- classify_bin_infovar_balance(bin, threshold=1)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 17, is_true())
  expect_that(length(classified$out) == 1, is_true())

  classified <- classify_bin(bin, technique='infovar_balance', params=list(threshold = 1))
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 17, is_true())
  expect_that(length(classified$out) == 1, is_true())

  classified <- classify_bin_infovar_balance(bin, threshold = 1, start_threshold = 0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 17, is_true())
  expect_that(length(classified$out) == 1, is_true())

  classified <- classify_bin_infovar_balance(bin, threshold = .1, start_threshold = 0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 10, is_true())
  expect_that(length(classified$out) == 8, is_true())

  classified <- classify_bin_infovar_balance(bin, threshold = 1, start_threshold = 5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == length(bin), is_true())
  expect_that(length(classified$out) == 0, is_true())

  bin <- DNAStringSet(c(rep('AAAAA', 5), rep('AAAAC', 10), 
                        rep('AAAAG', 2), rep('AATTT', 1)))
  classified <- classify_bin_infovar_balance(bin, threshold = 1, start_threshold = 0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 10, is_true())
  expect_that(length(classified$out) == 8, is_true())

  bin <- DNAStringSet(c(rep('AAAAA', 5), rep('AAAAC', 10), 
                        rep('AAAAG', 2), rep('AAATT', 1)))
  classified <- classify_bin_infovar_balance(bin, threshold = 1, start_threshold = 0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == length(bin), is_true())
  expect_that(length(classified$out) == 0, is_true())
})

test_that('it stops when no decision can be made about which sequence to remove', {
  bin <- DNAStringSet(c('AAAAA', 'CCCCC', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_bin_infovar_balance(bin, threshold=1)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 0, is_true())
  expect_that(length(classified$out) == 4, is_true())
})

test_that('all equally distant sequences are always removed simultaneously', {
  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_bin_infovar_balance(bin, threshold=2.00)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 4, is_true())
  expect_that(length(classified$out) == 0, is_true())

  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_bin_infovar_balance(bin, threshold=1.9999999)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 2, is_true())
  expect_that(length(classified$out) == 2, is_true())
})
