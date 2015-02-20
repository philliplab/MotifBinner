context('Most Consensus String Constructor')

test_that('The most consensus string constructor works', {
  seqs <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                         'AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT'))
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet('AAAAAAAAA')
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAA')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT'))
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())


  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAC','AAAAAAAAAAG','AAAAAAAAAAT'))
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAA')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT',
                         'AAAAAAAAAAC'))
  consen <- mostConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())

})
