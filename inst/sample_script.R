library(MotifBinner)

x <- readFastq("~/projects/MotifBinner/data/CAP177_2040_v1merged.fastq")
x <- x@sread
seq_data <- x
prefix <- "CCAGCTGGTTATGCGATTCTMARGTG"
suffix <- "CTGAGCGTGTGGCAAGGCCC"
motif_length <- 9
max.mismatch <- 5
fixed = FALSE

y <- extract_motifs(seq_data, prefix, suffix, motif_length, max.mismatch, fixed)

print(y)

print(table(table(names(y$matched_seq))))
