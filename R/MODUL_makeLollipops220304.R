# Alex Abbas
# 220304
# make lollipop plots for MODUL

library(tidyverse)
library(readxl)
library(g3viz)

setwd("/Volumes/GoogleDrive/Shared drives/GI Biomarkers [gRED OBD]/CRC/MODUL Cohort 1/figures/")
source("../R/loadData220218.R")

# maf.file <- system.file("extdata", "TCGA.BRCA.varscan.somatic.maf.gz", package = "g3viz")
# mutation.dat <- readMAF(maf.file)

coh99toexp = function(x) { x[x=="Coh99"] = "Experimental"; x }

oaplMafGain = oapl %>%
  mutate(Arm = coh99toexp(Arm)) %>% 
  filter(SUBJECT.ID %in% C1PdSubjs) %>%
  mutate(variant = gsub("-","",paste0(SV.PROTEIN.CHANGE,CNA.TYPE,REARR.DESCRIPTION))) %>%
  select(SUBJECT.ID,GENE,variant,baseEnd,Arm,SV.PROTEIN.CHANGE,SOMATIC.STATUS.FUNCTIONAL.IMPACT) %>% unique() %>% 
  group_by(SUBJECT.ID,GENE,variant,Arm,SV.PROTEIN.CHANGE,SOMATIC.STATUS.FUNCTIONAL.IMPACT) %>% summarise(baseEnd = toString(baseEnd)) %>% ungroup() %>%
  filter(!grepl(",",baseEnd)) %>% 
  mutate(gainloss = ifelse(baseEnd == "base","loss","gain")) %>% unique() %>% 
  filter(gainloss == "gain" & !is.na(GENE)) %>%
  group_by(SUBJECT.ID,variant) %>% slice(n=1) %>% ungroup() %>%
  #  select(GENE,SV.PROTEIN.CHANGE,SOMATIC.STATUS.FUNCTIONAL.IMPACT) %>%
  mutate(
    AA_Position = gsub("\\D","",SV.PROTEIN.CHANGE),
    HGVSp_Short = paste0("p.",SV.PROTEIN.CHANGE),
    Hugo_Symbol = GENE,
    Impact = SOMATIC.STATUS.FUNCTIONAL.IMPACT
  )

oneLollipop = function(dat,agene) {
  dat = dat %>% filter(AA_Position != "")
  g3Lollipop(dat,
             gene.symbol = agene,
             plot.options = g3Lollipop.theme( theme.name = "nature2" ), # default, cbioportal, nature, nature2, dark, blue, ggplot2, and simple
             output.filename = paste0(agene,"_lollipop"),
             factor.col = "Impact"
  )
}

oneLollipop(oaplMafGain,"BRAF")
oneLollipop(oaplMafGain,"MAP2K1")
oneLollipop(oaplMafGain,"KRAS")
oneLollipop(oaplMafGain,"NRAS")

#maf.file <- system.file("extdata", "TCGA.COAD.varscan.somatic.maf.gz", package = "g3viz")
mutation.dat <- readMAF("~/Downloads/gdc_download_20220310_204838.443922.tar")
