context('Classification Checker')

test_that('The classification check allows valid classifications through', {
  bin <- DNAStringSet('AA')
  classified <- list('src' = bin,
                     'out' = DNAStringSet(NULL))
  expect_that(check_classification(bin, classified), is_true())

})

test_that('The classification check detects various issues', {

  bin <- DNAStringSet('AA')
  classified <- list('src' = bin,
                     'out' = bin)

  expect_error(check_classification(bin, 1),
               "classified must be a list")
  
  expect_error(check_classification(bin, list('hello' = 1)),
               "classified must contain and only contain the 'src' and 'out' elements")

  expect_error(check_classification(bin, list('src' = 1,
                                              'out' = 2)),
               "'src' element of classified must be a DNSStringSet")

  expect_error(check_classification(bin, classified), 
               "None of the 'src' sequences may be in the 'out' sequences")

  bin <- DNAStringSet(c('AA', 'AA', 'AB'))
  classified <- list('src' = bin[1:2],
                     'out' = c(bin[3], bin[3]))
  expect_error(check_classification(bin, classified), 
               "All input sequences must be present at the same level in the output")

  bin <- DNAStringSet(c('AA', 'AA', 'AB'))
  classified <- list('src' = bin[1:2],
                     'out' = c(bin[3], DNAStringSet('AC')))
  expect_error(check_classification(bin, classified), 
               "Classification process may not introduce new sequences")
})


