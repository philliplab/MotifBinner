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
  seqs <- DNAStringSet(c('AAA','AAA','AAA', 'AAA'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAA')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_4', is_true())

  seqs <- DNAStringSet(c('AAA','AAA','AAA', 'AA-'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AAA')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_4', is_true())

  seqs <- DNAStringSet(c('AAA','AAA','AA-', 'AA-'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AA+')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_4', is_true())

  seqs <- DNAStringSet(c('AAA','AA-','AA-', 'AA-'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AA-')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_4', is_true())

  seqs <- DNAStringSet(c('AA-','AA-','AA-', 'AA-'))
  consen <- construct_consensus(seqs, technique = 'mostConsensusString', params = list())
  target <- DNAString('AA-')
  expect_that(consen[[1]] == target, is_true())
  expect_that(names(consen)[1] == '_4', is_true())
})
