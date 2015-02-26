context('Two Step Consensus String Constructor')

test_that('The Two Step consensus string constructor works', {
  seqs <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                         'AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT'))
  consen <- twoStepConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet('AAAAAAAAA')
  consen <- twoStepConsensusString(seqs)
  target <- DNAString('AAAAAAAAA')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- twoStepConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT'))
  consen <- twoStepConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())


  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAC','AAAAAAAAAAG','AAAAAAAAAAT'))
  consen <- twoStepConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAA')
  expect_that(consen == target, is_true())
})

test_that('two step consensus string can be called via construct_consensus', {
  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT',
                         'AAAAAAAAAAC'))
  consen <- construct_consensus(seqs, technique = 'twoStepConsensusString', params = list())
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_7', is_true())
})
