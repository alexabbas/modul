# Alex Abbas
# 210608
# Update sample matching table for MODUL FMI with 210309 manifest info 

library(tidyverse)

ctdna_dir="/Volumes/OBD/avastin/mo29112_modul/ngs/results/ctdna/data"
sampleMatch = read.csv(paste(ctdna_dir,"modul_FMIsamples.csv",sep="/"))
newManifest = read.csv("~/Google Drive/crc/MODUL/data/210309_MO29112_BPA-RET-20-2047_BRMRI-5905_20xPooledPlasmaSamples_FMI.csv")

newPats = setdiff(newManifest$Subject.ID,sampleMatch$subject_id_patient_id)
newSampleMatch = sampleMatch %>% dplyr::rename(Subject.ID=subject_id_patient_id) %>% bind_rows(
  data.frame(Subject.ID=newPats)
)

newPD = newManifest %>% filter(Visit=="PD")
newSampleMatch %>% filter(Subject.ID %in% newPD$Subject.ID)
newSampleMatch[match(newPD$Subject.ID,newSampleMatch$Subject.ID),"plasma_PD"] = newPD$Partner.Sample.ID
newSampleMatch %>% filter(Subject.ID %in% newPD$Subject.ID)

newSampleMatch = newSampleMatch %>% dplyr::rename(SUBJECT.ID = Subject.ID)
write.csv(newSampleMatch,paste(ctdna_dir,"modul_FMIsamples_210608.csv",sep="/"),na="")
