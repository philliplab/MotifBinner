context('Easy Consensus String Constructor')

test_that('The easy consensus string constructor works', {
  seqs <- DNAStringSet(c('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT',
                         'AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT'))
  consen <- easyConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAATGAAACGTTAGCGCGTTTTTTTTTT')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet('AAAAAAAAA')
  consen <- easyConsensusString(seqs)
  target <- DNAString('AAAAAAAAA')
  expect_that(consen == target, is_true())

  seqs <- DNAStringSet(c('AAAAAAAAAAA',
                         'AAAAAAAAAAT'))
  consen <- easyConsensusString(seqs)
  target <- DNAString('AAAAAAAAAAW')
  expect_that(consen == target, is_true())
})
