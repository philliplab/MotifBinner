context('Alignment')

test_that('The aligner works', {
  seqs <- DNAStringSet(rep(paste(rep('A', 50), collapse=""), 5))
  aligned <- align_sequences(seqs)
  expect_that(all(seqs == unname(aligned)), is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                         'AAAAAAAAAAATGCGTTAGCGCGTTTTTTTTTT'))
  aligned <- align_sequences(seqs)
  target <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                           'AAAAAAAAAAATG---CGTTAGCGCGTTTTTTTTTT'))
  expect_that(all(target == aligned), is_true())
})

