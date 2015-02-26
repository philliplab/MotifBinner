context('Consensus String Constructor')

test_that('The consensus string constructor works', {
  seqs <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                         'AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT'))
  consen <- construct_consensus(seqs)
  target <- DNAString('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT')
  expect_that(consen[[1]] == target, is_true())

  names(seqs) <- c('s1', 's2')
  consen <- construct_consensus(seqs)
  expect_that(names(consen), equals('s1_2'))

  seqs <- DNAStringSet('AAAAAAAAA')
  consen <- construct_consensus(seqs)
  target <- DNAString('AAAAAAAAA')
  expect_that(consen[[1]] == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- construct_consensus(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- construct_consensus(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- construct_consensus(seqs)
  target <- DNAString('AAAAAAAAAAA')
  expect_that(consen[[1]] == target, is_true())

  seqs <- DNAStringSet(c('AAAA-AAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- construct_consensus(seqs, params=list(ambiguityMap='N'))
  target <- DNAString('AAAAAAAAAAA')
  expect_that(consen[[1]] == target, is_true())

  seqs <- DNAStringSet(c('AAAA-AAAAAA',
                         'AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- construct_consensus(seqs, params=list(ambiguityMap='N'))
  target <- DNAString('AAAANAAAAAN')
  expect_that(consen[[1]] == target, is_true())
})
