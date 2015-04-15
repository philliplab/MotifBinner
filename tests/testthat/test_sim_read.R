context('sim_read')

test_that('gen_seq works as expected', {
  seq_dat <- gen_seq(400)
  expect_that(nchar(seq_dat), equals(400))
  uniq_lets <- sort(unique(strsplit(seq_dat, '')[[1]]))
  expect_that(uniq_lets, equals(c('A', 'C', 'G', 'T')))
})
