context('Absolute threshold based classifier')

test_that('The abs thres based classifier works', {
  bin <- DNAStringSet(c(rep('AAAAA', 5), rep('AAAAC', 10), 
                        rep('AAAAG', 2), rep('TTTTT', 1)))
  classified <- classify_absolute(bin, threshold=0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 17, is_true())
  expect_that(length(classified$out) == 1, is_true())
  expect_that(all(dim(classified$dmat) == c(18, 18)), is_true())
  expect_that(max(classified$dmat), equals(5))

  bin <- DNAStringSet(rep('AAAAA', 5))
  classified <- classify_absolute(bin, threshold=0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 5, is_true())
  expect_that(length(classified$out) == 0, is_true())

  bin <- DNAStringSet(c(rep('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 5),
                      rep('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC', 2),
                      rep('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT', 2),
                      rep('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACG', 2)
                      ))
  # stop when less than 2 errors in 50
  classified <- classify_absolute(bin, threshold=0.999/50, start_threshold = 0.999/50)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 9, is_true())
  expect_that(length(classified$out) == 2, is_true())
  expect_that(sort(unique(as.numeric(classified$dmat))), equals(0:2))

  # start when more than 3 errors in 50
  classified <- classify_absolute(bin, threshold=0.999/50, start_threshold = 1.5/50)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 11, is_true())
  expect_that(length(classified$out) == 0, is_true())

  # stop when less than 1 errors in 50
  classified <- classify_absolute(bin, threshold=0.499/50, start_threshold = 0.499/50)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 5, is_true())
  expect_that(length(classified$out) == 6, is_true())

  classified <- classify_bin(bin, technique='absolute', 
                             params = list(threshold=0.499/50, 
                                           start_threshold = 0.499/50))
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 5, is_true())
  expect_that(length(classified$out) == 6, is_true())
})

test_that('all equally distant sequences are always removed simultaneously', {
  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_absolute(bin, threshold=0.501)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 4, is_true())
  expect_that(length(classified$out) == 0, is_true())

  bin <- DNAStringSet(c('AAAAA', 'AAAAA', 
                        'GGGGG', 'TTTTT'))
  classified <- classify_bin_infovar_balance(bin, threshold=0.5)
  expect_that(check_classification(bin, classified), is_true())
  expect_that(length(classified$src) == 2, is_true())
  expect_that(length(classified$out) == 2, is_true())
})
