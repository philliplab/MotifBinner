context('Classification Checker')

test_that('The classification check detects various issues', {
  bin <- DNAStringSet('AA')
  classified <- list('src' = bin,
                     'out' = DNAStringSet(NULL))
  expect_that(check_classification(bin, classified), is_true())
})

