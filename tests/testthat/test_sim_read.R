context('gen_seq')

test_that('gen_seq works as expected', {
  seq_dat <- gen_seq(400)
  expect_that(nchar(seq_dat), equals(400))
  uniq_lets <- sort(unique(strsplit(as.character(seq_dat), '')[[1]]))
  expect_that(uniq_lets, equals(c('A', 'C', 'G', 'T')))
})

context('randomize_ambig')

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

context('mutate_base')

test_that('mutate_base works', {
  x <- 'A'
  y <- mutate_base(x)
  expect_that(y != x, is_true())
  expect_that(y %in% c('A', 'C', 'G', 'T'), is_true())
  expect_that(length(y), equals(1))
  expect_that(nchar(y), equals(1))

  x <- 'C'
  y <- mutate_base(x)
  expect_that(y != x, is_true())
  expect_that(y %in% c('A', 'C', 'G', 'T'), is_true())
  expect_that(length(y), equals(1))
  expect_that(nchar(y), equals(1))

  x <- 'CC'
  expect_that(mutate_base(x), throws_error())
  x <- c('T','A')
  expect_that(mutate_base(x), throws_error())
})

context('add_snps')

test_that('length of sequence is longer than number of desired snps', {
  seq_dat <- DNAStringSet(c('AAAA', 'CCCC'))
  expect_that(class(add_snps(seq_dat, 1)) == 'DNAStringSet', is_true())
  expect_that(class(add_snps(seq_dat, 4)) == 'DNAStringSet', is_true())

  expect_that(add_snps(seq_dat, 5), throws_error())
  expect_that(add_snps(seq_dat, 100), throws_error())
})

test_that('the correct number of snps are inserted', {
  seq_dat <- DNAStringSet('ACGT')
  mut_seq <- add_snps(seq_dat, 1)
  seq_dist <- stringDist(DNAStringSet(c(as.character(seq_dat), 
                                        as.character(mut_seq))))
  expect_that(seq_dist == 1, is_true())


  seq_dat <- DNAStringSet('ACGT')
  mut_seq <- add_snps(seq_dat, 3)
  seq_dist <- stringDist(DNAStringSet(c(as.character(seq_dat), 
                                        as.character(mut_seq))))
  expect_that(seq_dist == 3, is_true())

  seq_dat <- DNAStringSet(c('ACGTACGT', 'GGTTCCAA'))
  mut_seq <- add_snps(seq_dat, 2)

  seq_dist1 <- stringDist(DNAStringSet(c(as.character(seq_dat[1]), 
                                         as.character(mut_seq[1]))))
  seq_dist2 <- stringDist(DNAStringSet(c(as.character(seq_dat[2]), 
                                         as.character(mut_seq[2]))))
  expect_that(seq_dist1 == 2, is_true())
  expect_that(seq_dist2 == 2, is_true())
})

test_that('zero snps can be handled', {
  seq_dat <- DNAStringSet('ACGT')
  mut_seq <- add_snps(seq_dat, 0)
  seq_dist <- stringDist(DNAStringSet(c(as.character(seq_dat), 
                                        as.character(mut_seq))))
  expect_that(seq_dist == 0, is_true())
})

context('gen_pid_search_scenario ')

test_that('the output is in the correct format', {
  pid_search <- gen_pid_search_scenario(seq_len = 500, prefix_len = 25, 
                                        pid_len = 9, suffix_len = 15, 
                                        prefix_snps = 0, suffix_snps = 0, 
                                        suffix_chop = 0)
  expect_that('seq_dat' %in% names(pid_search), is_true())
  expect_that('prefix' %in% names(pid_search), is_true())
  expect_that('suffix' %in% names(pid_search), is_true())
  expect_that('pid' %in% names(pid_search), is_true())
  expect_that(nchar(pid_search$seq_dat), equals(500 + 25 + 9 + 15))
  expect_that(nchar(pid_search$prefix), equals(25))
  expect_that(nchar(pid_search$suffix), equals(15))
  expect_that(nchar(pid_search$pid), equals(9))
})

test_that('prefix and suffix mutation works', {
})
