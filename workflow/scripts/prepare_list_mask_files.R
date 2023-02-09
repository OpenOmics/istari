#!/usr/bin/env Rscript
library(dplyr)
library(magrittr)
library(data.table)
# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
FILE2 <- args[3]
FILE3 <- args[4]
FILE4 <- args[5]
FILE5 <- args[6]
FILE6 <- args[7]
FILE7 <- args[8]
FILE8 <- args[9]

#### first set of files
# read in annotation file
setwd(DIR)
annot_file1 = fread(FILE1, header = F)
annot_file2=subset(annot_file1, annot_file1$V2 !=".")
annot_file2$V12 <- gsub(":", "_", annot_file2$V1)
fwrite(annot_file2[,c('V12', 'V2', 'V3')], FILE2, row.names = F, col.names = F, sep = '\t', quote = F)

# create list file
set_file1 = annot_file2[,c('V1', 'V2')]
set_file2 = set_file1[!duplicated(set_file1$V2),]
set_file2$CHR = stringr::str_split(set_file2$V1, ':', simplify = T)[,1]
set_file2$BP = stringr::str_split(set_file2$V1, ':', simplify = T)[,2]
set_file2$SNPS = NA
for(i in 1:nrow(set_file2)){
  set_file2$SNPS[i] = paste(annot_file2$V1[annot_file2$V2==set_file2$V2[i]], collapse = ',')
}

set_file2$SNPS2 <- gsub(":", "_", set_file2$SNPS)
fwrite(set_file2[,c('V2', 'CHR', 'BP','SNPS2')], FILE3, row.names = F, col.names = F, sep = '\t', quote = F)

#create mask file
masks1=data.table(unique(annot_file2$V3))
masks1
all1=data.table(Col1=rep("All", nrow(masks1)), Col2=masks1$V1)
All_df1 <- all1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
missenses1=data.table(masks1[V1 %like% "missense_variant"])
missense1=data.table(Col1=rep("Missense", nrow(missenses1)), Col2=missenses1$V1)
Missense_df1 <- missense1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
badness1=masks1[masks1$V1 != "missense_variant", ]
Badness1=data.table(Col1=rep("Badness", nrow(badness1)), Col2=badness1$V1)
Badness_df1 <- Badness1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
mask_df1=rbind(Missense_df1, Badness_df1, All_df1)
write.table(mask_df1, FILE4, row.names=F, col.names=F, quote=F, sep=' ')


## second set of files
annot_file1 = fread(FILE5, header = F)
annot_file2=subset(annot_file1, annot_file1$V2 !=".")
annot_file2$V12 <- gsub(":", "_", annot_file2$V1)
fwrite(annot_file2[,c('V12', 'V2', 'V3')], FILE6, row.names = F, col.names = F, sep = '\t', quote = F)

# create list file
set_file1 = annot_file2[,c('V1', 'V2')]
set_file2 = set_file1[!duplicated(set_file1$V2),]
set_file2$CHR = stringr::str_split(set_file2$V1, ':', simplify = T)[,1]
set_file2$BP = stringr::str_split(set_file2$V1, ':', simplify = T)[,2]
set_file2$SNPS = NA
for(i in 1:nrow(set_file2)){
  set_file2$SNPS[i] = paste(annot_file2$V1[annot_file2$V2==set_file2$V2[i]], collapse = ',')
}

set_file2$SNPS2 <- gsub(":", "_", set_file2$SNPS)
fwrite(set_file2[,c('V2', 'CHR', 'BP','SNPS2')], FILE7, row.names = F, col.names = F, sep = '\t', quote = F)

#create mask file
masks1=data.table(unique(annot_file2$V3))
masks1
all1=data.table(Col1=rep("All", nrow(masks1)), Col2=masks1$V1)
All_df1 <- all1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
missenses1=data.table(masks1[V1 %like% "missense_variant"])
missense1=data.table(Col1=rep("Missense", nrow(missenses1)), Col2=missenses1$V1)
Missense_df1 <- missense1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
badness1=masks1[masks1$V1 != "missense_variant", ]
Badness1=data.table(Col1=rep("Badness", nrow(badness1)), Col2=badness1$V1)
Badness_df1 <- Badness1 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
mask_df1=rbind(Missense_df1, Badness_df1, All_df1)
write.table(mask_df1, FILE8, row.names=F, col.names=F, quote=F, sep=' ')
