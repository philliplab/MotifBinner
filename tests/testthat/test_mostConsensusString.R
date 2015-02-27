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

test_that('most consensus string can be called via construct_consensus', {
  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT',
                         'AAAAAAAAAAC'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_7', is_true())
})

test_that('most consensus string plays nice with gaps',{
  seqs <- DNAStringSet(c('-AAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         'AAAAAAAAAAT','AAAAAAAAAAT','AAAAAAAAAAT',
                         'AAAAAAAAAAC'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_7', is_true())

  seqs <- DNAStringSet(c('-AAAAAAAAAA','-AAAAAAAAAA','AAAAAAAAAAA',
                         '-AAAAAAAAAT','-AAAAAAAAAT','AAAAAAAAAAT',
                         '-AAAAAAAAAC'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_7', is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA','AAAAAAAAAAA','AAAAAAAAAAA',
                         '-AAAAAAAAAT','-AAAAAAAAAT','-AAAAAAAAAT'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_7', is_true())

})
