# This is a modified version of script 04_ParseSnpeEfSNPs and 06_Analyze_SnpEffSNPs.  
# While the other version filters SNPs just to retain those of moderate and high effect, 
# here we will keep all (After filtering for quality, etc).  Will then run fisher test and anova tests
# annotate

# This is modified to be run as a script on Barbera

# setup

.libPaths("/share/malooflab/Julin/R/x86_64-pc-linux-gnu-library/4.1")

library(tidyverse)
library(VariantAnnotation)

datadir <- "../input"
outdir <- "../output"

taskID <- commandArgs(trailingOnly = TRUE)[1]

cat("task: ", taskID, "\n")

tab <- TabixFile(file.path(datadir, "cohort.all.genotyped.snpEff.vcf.gz"))

tab_seqnames <- seqnamesTabix(tab)

## load contig splits file 
contig.splits <- read_csv("../input/contig_splits.csv")

contig.splits_tab <- contig.splits %>%
  filter(name %in% tab_seqnames)

if(last > nrow(contig.lengths_tab)) last <- nrow(contig.lengths_tab)

outfile <- str_c("Annotated_SNPeff_ALL_quant", first, last,".csv.gz", sep="_")

if(file.exists(file.path(outdir, outfile))) file.remove(file=file.path(outdir,outfile))

## get blast info
header <- c("query",
            "subject",
            "percentID",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore"
)
blast <- read_csv("../output/TanOakVsA.t.blast_out.csv", col_names = header)
blast.best <- blast %>% group_by(query) %>%
  filter(row_number(dplyr::desc(bitscore))==1)


# get annotation.  Download from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSubscriber_Data_Releases%2FTAIR_Data_20201231

atDesc <- read_tsv("../input/Araport11_functional_descriptions_20201231.txt.gz") %>%
  mutate(name = str_remove(name, "\\..*$")) %>%
  rename_all(funs(str_c("At_", .))) %>%
  filter(!duplicated(At_name)) %>% #not ideal
  dplyr::select(-At_gene_model_type)

atSymbol <- read_tsv("../input/gene_aliases_20201231.txt.gz") %>%
  rename_all(funs(str_c("At_", .))) %>%
  filter(!duplicated(At_name)) #not ideal

blast.best <- blast.best %>%
  mutate(AGI = str_remove(subject, "\\..*$")) %>%
  left_join(atSymbol, by = c("AGI" = "At_name")) %>%
  left_join(atDesc, by = c("AGI" = "At_name")) %>%
    dplyr::select(query, subject, percentID, length, starts_with("At_")) %>%
    mutate(query=str_remove(query, "-mRNA-.*"))

## get pheno data
pheno <- read_csv("../input/TanOakResistance.csv") %>%
  mutate(nameMatch=str_replace_all(Tanoak, "-","."),
         nameMatch=str_replace(nameMatch, "SM\\.(52|54)\\.(42|81|97|37)", "\\1.\\2")) %>%
  arrange(Rank_Sorting)

## fisher p function
fisherp <- function(data){ #get p value from fisher test
  if(length(unique(data$GT))==1) 
    return(NA) 
  tb <- table(data$GT, data$Tolerance_Resistance)
  fisher.test(tb) %>%
    magrittr::extract("p.value") %>% 
    unlist() %>%
    return()
}

## kendall p function

kendallp <- function(data) { # get p value from kendall ranked correlation test
  if(length(unique(data$GT))==1) 
    return(NA) 
  cor.test(as.numeric(as.factor(data$GT)), 
           data$Rank_Sorting, 
           method = "kendall",
           exact=FALSE)$p.value
  }

## anova p function
anovap <- function(data){ #get p value from anova test
  if(length(unique(data$GT))==1) 
    return(NA) 
  aov(Rank_Sorting ~ GT, data=data) %>% 
    summary()%>%
    magrittr::extract2(1) %>%
    magrittr::extract(1, "Pr(>F)") %>%
    return()
}

## loop through the contigs

for(i in first:last) {

  # read it
  vcf <-  try(
    readVcf(tab, 
                  param=GRanges(
                    seqnames = contig.lengths$name[i],
                    ranges = IRanges(
                      start = 1,
                      end=contig.lengths$length[i])))
  )
  
  if(class(vcf)=="try-error") next() #some small contigs aren't in vcf and throw an error.
  
  # filter it
  vcf <- vcf[rowRanges(vcf)$QUAL > 50 &
               sapply(info(vcf)$AC, length) == 1 & #stick with biallelic SNPs for now
               info(vcf)$AN == 28] # require complete genotype info for now 
  
  if(nrow(vcf)==0) next()
  
  # convert to tibble
  vcf.filter.VR <- as(vcf, "VRanges") 
  VRtibble <- tibble(
    seqname=as.character(seqnames(vcf.filter.VR)),
    start=start(vcf.filter.VR),
    end=end(vcf.filter.VR),
    sample=as.character(sampleNames(vcf.filter.VR)),
    ref=ref(vcf.filter.VR),
    alt=alt(vcf.filter.VR),
    RD=refDepth(vcf.filter.VR),
    AD=altDepth(vcf.filter.VR),
    as.data.frame(mcols(vcf.filter.VR)[c("QUAL", "GT", "GQ", "PL", "ANN" )])
  ) %>%
    mutate(snpID=str_c(seqname,start,sep=":"),
           GT=str_replace(GT, stringr::fixed("|"), "/")) 
  
  rm(vcf, vcf.filter.VR)
  gc()
  
  
  #join with pheno data
  VRtibble <- VRtibble %>% inner_join(pheno, by=c("sample"="nameMatch"))
  VRtibble <- VRtibble %>% dplyr::select(seqname, start, end, snpID, sample, ref, alt, GT, Tolerance_Resistance, Rank_Sorting, ANN) %>%
    group_by(snpID) %>%
    nest(data=c(sample, GT, Tolerance_Resistance, Rank_Sorting))
  
  # calculate fisher p value
  VRtibble <- VRtibble %>%
    mutate(fisher.p=map_dbl(data, fisherp),
          # anova.p=map_dbl(data, anovap),
           kendall.p=map_dbl(data, kendallp)) %>% 
    dplyr::select(-data)
  
   # reformat snpEff annotation
  annnames <- c("allele", "effect", "impact", "GeneName", "GeneID", "FeatureType", "FeatureID", "TranscriptType", "Rank_Total", "HGVS.c", "HGVS.p", "cDNApos_len", "CDSpos_len", "protpos_len", "distance", "warn")
  
  VRtibble <- VRtibble %>% ungroup() %>% 
    hoist(ANN, ANN1=1, ANN2=2, ANN3=3, ANN4=4, .remove = FALSE) %>% # ANN not always removed even if TRUE, perhaps if it has more than 5 components
    dplyr::select(-ANN) 

  VRtibble <- VRtibble %>%
    pivot_longer(cols = starts_with("ANN"), names_to = "ANN_ID", values_to = "ANN", values_drop_na = TRUE) 
  
  VRtibble <- VRtibble %>%  separate(ANN, into=annnames, sep = "\\|") 
  
  # filter to remove redundant entries
  VRtibble <- VRtibble %>% # only keep one entry per SNP.  keep the one with the most impact.
    group_by(snpID) %>%
    mutate(impact=factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
    arrange(impact) %>%
    filter(!duplicated(snpID)) %>%
    arrange(start) %>% dplyr::select(seqname:kendall.p, effect, GeneName, cDNApos_len, distance )
  
  # add blast annotation
  VRtibble <- VRtibble %>% 
    left_join(blast.best, by= c("GeneName" = "query")) %>%
    dplyr::select(seqname, start, end, ref, alt, fisher.p:kendall.p, GeneName, everything())
  
  # write it
  write_csv(VRtibble, file=file.path(outdir, outfile), append = TRUE)
}


