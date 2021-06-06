# generate splits for array slurm for annotating SNPeff resutls (script -4b)

contigs <- 30271
cores <- 10
spacing <- floor(30271/cores)

starts <- 0:(cores-1)*spacing + 1

ends <- c(1:(cores-1)*spacing,contigs)

cbind(starts, ends)

write.table(cbind(starts, ends), file = "../input/contig_splits.txt", row.names = FALSE, col.names = FALSE)