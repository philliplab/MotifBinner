library(Biostrings)
library(ShortRead)

x <- DNAStringSet('ACGTGTGCGCGATGCGAATTTACG')
x1 <- DNAStringSet('ACGTGTGCGCGATCGATTACG')

qual <- as.integer(rep(38, nchar(x)))
qual1 <- as.integer(rep(38, nchar(x1)))

ill_qual <- IlluminaQuality(qual)
ill_qual1 <- IlluminaQuality(qual1)
ugh <- BStringSet(ill_qual)
sr_q <- IntegerQuality(qual)

y <- QualityScaledDNAStringSet(x, ill_qual)

writeXStringSet(y, '/tmp/x.fastq')

z <- grep('A', '', y)

yay <- ShortReadQ(sread = x, quality = ugh, id = BStringSet('ogggieboogie'))
yay <- ShortReadQ(sread = x, quality = ill_qual, id = BStringSet('ogggieboogie'))
yay1 <- ShortReadQ(sread = c(x, x1), quality = c(ill_qual, ill_qual1), 
                   id = BStringSet(c('ogggieboogie', 'x1')))

sread(yay1)

writeFastq(yay, '/tmp/yay.fastq', compress=FALSE)
