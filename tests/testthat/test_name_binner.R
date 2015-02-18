context('Name Binner')

test_that('The name binner works with 2 bins', {
  motif_extracted <- DNAStringSet(c('AAA', 'AAA', 'AAC', 'CCC'))
  names(motif_extracted) <- c('a', 'a', 'a', 'b')
  binned <- bin_by_name(motif_extracted)
  expect_that(all(sort(names(binned)) == c('a', 'b')), is_true())
  expect_that(length(binned$a), equals(3))
  expect_that(length(binned$b), equals(1))
  expect_that(unname(as.character(binned$b)), equals('CCC'))

})

test_that('The name binner works with 1 bins', {
  motif_extracted <- DNAStringSet(c('AAA', 'AAA', 'AAC', 'CCC'))
  names(motif_extracted) <- c('a', 'a', 'a', 'a')
  binned <- bin_by_name(motif_extracted)
  expect_that(all(sort(names(binned)) == c('a')), is_true())
  expect_that(length(binned$a), equals(4))
})

test_that('The name binner works with size 1 bins', {
  motif_extracted <- DNAStringSet(c('AAA', 'TAA', 'AAC', 'CCC'))
  names(motif_extracted) <- c('a', 'b', 'c', 'd')
  binned <- bin_by_name(motif_extracted)
  expect_that(all(sort(names(binned)) == c('a', 'b', 'c', 'd')), is_true())
  expect_that(length(binned$a), equals(1))
  expect_that(length(binned$b), equals(1))
  expect_that(unname(as.character(binned$b)), equals('TAA'))
})

