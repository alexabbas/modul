# Alex Abbas
# 201114 ... Saturday :(
# Make SummarizedExperiment of MODUL RNAseq and clinical data

library(tidyverse)
library(limma)
library(SummarizedExperiment)

source("~/Google Drive/R/aaStuff.R")

#countsTable = read.csv("/Volumes/OBD/avastin/mo29112_modul/ngs/outdata/ngs_pipeline/ea_obd_mo29112_ngs_rna_targrna_rnaaccess_20200814/GeneExonic_StrandedCounts.csv")
countsTable = read.csv("~/Desktop/GeneExonic_StrandedCounts.csv")
clinDat = read.csv("~/Google Drive/crc/MODUL/data/adk200413.csv")

clinPatids = paste0("PAT",clinDat$SUBJID)
countsPatids = gsub("count.","PAT",colnames(countsTable)[4:dim(countsTable)[2]])
#setTally(clinPatids,countsPatids)

# hgxDat = read.csv("~/Google Drive/crc/MODUL/data/MO29112_HGX_Results_Prelim_20200702135253_MP 6.csv")
# hgxPatids = paste0("PAT",hgxDat$Patient.ID)
# setTally(hgxPatids,countsPatids)

colData = clinDat %>% mutate(
  SUBJID = paste0("PAT",SUBJID),
  TRT = factor(gsub("+","",TRT,fixed=T)),
  RESINDI = factor(gsub("/","",RESINDI,fixed=T)),
  KRAS = factor(KRAS),
  SEX = factor(SEX),
  RNR = factor(RNR),
  PFS.event = 1-PFSCNSR,
  PFSmod = ifelse(PFSCNSR==0,PFSTTE,100),
  progressEarly = ifelse(PFSmod<median(PFSmod),"early","notEarly")
) %>% select(-BEP)
colData = colData[match(countsPatids,colData$SUBJID),]

# filter RNAseq features down to named genes. For duplicate symbols, keep highest expressed.
countsWithSymbol = countsTable %>% filter(HGNC_Symbol!="")
countsUniqSymbol = countsWithSymbol %>%
  mutate(sum = rowSums(countsWithSymbol[,4:dim(countsWithSymbol)[2]])) %>%
  arrange(desc(sum)) %>% group_by(HGNC_Symbol) %>% dplyr::slice(n=1) %>% ungroup() %>% select(-sum)

counts = as.matrix(countsUniqSymbol[,4:dim(countsUniqSymbol)[2]])
colnames(counts) = countsPatids
rownames(counts) = countsUniqSymbol$HGNC_Symbol
v = voom(counts)$E
rownames(v) = rownames(counts)

rowData = data.frame(countsUniqSymbol %>% select(HGNC_Symbol,ENSEMBL_ID))

modulSE = SummarizedExperiment(assays=list(counts=counts,voom=v),
                     rowData=rowData, colData=colData)

save(modulSE,file="~/Google Drive/crc/MODUL/data/modulSE.RData")


