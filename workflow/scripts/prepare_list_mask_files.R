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
#### first set of files 
# read in annotation file
setwd(DIR)
annot_file1 = fread(FILE1, header = F)

# create list file
set_file1 = annot_file1[,c('V1', 'V2')]
set_file1 = set_file1[!duplicated(set_file1$V2),]
set_file1$CHR = stringr::str_split(set_file1$V1, ':', simplify = T)[,1]
set_file1$BP = stringr::str_split(set_file1$V1, ':', simplify = T)[,2]
set_file1$SNPS = NA
for(i in 1:nrow(set_file1)){
  set_file1$SNPS[i] = paste(annot_file1$V1[annot_file1$V2==annot_file1$V2[i]], collapse = ',')
}
fwrite(set_file1[,c('V2', 'CHR', 'BP','SNPS')], FILE2, row.names = F, col.names = F, sep = '\t', quote = F)

#create mask file
masks1=data.table(unique(annot_file1$V3))
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
write.table(mask_df1, FILE3, row.names=F, col.names=F, quote=F, sep=' ')


## second set of files
library(data.table)
# read in annotation file
annot_file2 = fread(FILE4, header = F)

# create list file
set_file2 = annot_file2[,c('V1', 'V2')]
set_file2 = set_file2[!duplicated(set_file2$V2),]
set_file2$CHR = stringr::str_split(set_file2$V1, ':', simplify = T)[,1]
set_file2$BP = stringr::str_split(set_file2$V1, ':', simplify = T)[,2]
set_file2$SNPS = NA
for(i in 1:nrow(set_file2)){
  set_file2$SNPS[i] = paste(annot_file2$V1[annot_file2$V2==annot_file2$V2[i]], collapse = ',')
}
fwrite(set_file2[,c('V2', 'CHR', 'BP','SNPS')], FILE5, row.names = F, col.names = F, sep = '\t', quote = F)

#create mask file
masks2=data.table(unique(annot_file2$V3))
masks2
all2=data.table(Col1=rep("All", nrow(masks2)), Col2=masks2$V1)
All_df2 <- all2 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
missenses2=data.table(masks2[V1 %like% "missense_variant"])
missense2=data.table(Col1=rep("Missense", nrow(missenses2)), Col2=missenses2$V1)
Missense_df2 <- missense2 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
badness2=masks2[masks2$V1 != "missense_variant", ]    
badness2=data.table(Col1=rep("Badness", nrow(badness2)), Col2=badness2$V1)
Badness_df2 <- badness2 %>%
  dplyr::group_by(Col1) %>%
  dplyr::summarise(Col2 = paste(Col2, collapse = ","))
mask_df2=rbind(Missense_df2, Badness_df2, All_df2)
write.table(mask_df2, FILE6, row.names=F, col.names=F, quote=F, sep=' ')
