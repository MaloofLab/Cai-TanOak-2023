# generate splits for array slurm for annotating SNPeff results (script -4b)
# note contigs are in order of size, so we can't just split them and give the first 1/10 to 1 array, etc.

library(tidyverse)

cores <- 10

contig.lengths <-
  system("bcftools view -h '/Volumes/GoogleDrive/Shared drives/TanOak/new_12_2021_assembly/cohort.all.genotyped.snpEff.vcf.gz' | grep contig", inter=TRUE) %>%
  as_tibble() %>%
  mutate(value=str_remove_all(value, "(##contig=<ID=)|(length=)|>")) %>%
  separate(value, into = c("name", "length"), sep=",", convert = TRUE ) %>%
  mutate(split=rep(1:cores, length.out=nrow(.)))

write_csv(contig.lengths, "../input/contig_splits.csv")
