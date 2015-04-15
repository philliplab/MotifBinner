context('sim_read')

test_that('gen_seq works as expected', {
  seq_dat <- gen_seq(400)
  expect_that(nchar(seq_dat), equals(400))
  uniq_lets <- sort(unique(strsplit(seq_dat, '')[[1]]))
  expect_that(uniq_lets, equals(c('A', 'C', 'G', 'T')))
})

test_that('randomize_ambig preserves input class reasonably', {
  seq_dat <- 'AA'
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(class(seq_dat), equals(class(proc_seq)))

  seq_dat <- DNAString('AA')
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(class(proc_seq) == 'DNAStringSet', is_true())

  seq_dat <- DNAStringSet('AA')
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(class(seq_dat), equals(class(proc_seq)))
})

