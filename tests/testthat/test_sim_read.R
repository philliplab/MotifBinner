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

test_that('degeneracy is removed', {
  seq_dat <- 'AY'
  targets <- c('AC', 'AT')
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(proc_seq %in% targets, is_true())

  seq_dat <- 'AK'
  targets <- c('AG', 'AT')
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(proc_seq %in% targets, is_true())

  seq_dat <- 'MK'
  targets <- c('AG', 'AT', 'CG', 'CT')
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(proc_seq %in% targets, is_true())
})

test_that('multiple input sequences are handled correctly', {
  seq_dat <- DNAStringSet(c('AA', 'AY'))
  targets <- list('1' = DNAStringSet(c('AA', 'AC')),
                  '2' = DNAStringSet(c('AA', 'AT')))
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(all(proc_seq == targets[[1]]) | 
              all(proc_seq == targets[[2]]), 
              is_true())

  seq_dat <- DNAStringSet(c('WA', 'CG'))
  targets <- list('1' = DNAStringSet(c('AA', 'CG')),
                  '2' = DNAStringSet(c('TA', 'CG')))
  proc_seq <- randomize_ambig(seq_dat)
  expect_that(all(proc_seq == targets[[1]]) | 
              all(proc_seq == targets[[2]]), 
              is_true())

})
