# Alex Abbas
# 200324
# Make analysis-ready MODUL data table from SAS files and HGX mutation table

library(tidyverse)
library(xlsx)
library(haven)

wb = loadWorkbook("~/Google Drive/crc/MODUL/data/MO29112_IxRs_ResultsTracker_Ch2_02-JAN-2017.xlsx")
sheetnames = names(getSheets(wb))
hgx = lapply(sheetnames,function(i){
  read.xlsx("~/Google Drive/crc/MODUL/data/MO29112_IxRs_ResultsTracker_Ch2_02-JAN-2017.xlsx",i)
})
names(hgx) = sheetnames

sasfiles = list.files("~/Google Drive/crc/MODUL/data/ADaM Jul 2019 (9)/",full.names = T)
sasshorts = gsub(".sas7bdat","",list.files("~/Google Drive/crc/MODUL/data/ADaM Jul 2019 (9)/"))
sasdats = lapply(sasfiles,read_sas)
names(sasdats) = sasshorts

sapply(sasdats,colnames)
colnames(sasdats[["adtte"]])

# RESINDC Tumor Response at the end of Induction treatment phase (ITP)  based on investigator assessment (CRF)
# RESIGR1 Tum Resp end of induction treatment phase (ITP)-  based on investigator assessment (CRF) Goup 1 (this variable is pooling PR and CR together) 
# RESINDI Tumor Response at the end of ITP  induction treatment phase  based on ixRS . There is a difference between IxRS and CRF because investigators made mistakes during randomization
# CNSR (0 = event)
# PARAMCD
#   OS
#   PFSM

krasHash = c(
  "Mutation Detected" = "MUT",
  "Mutation Not Detected" = "WT",
  "Invalid" = NA
)
kras = hgx[["KRAS"]] %>%
  mutate(KRAS = krasHash[as.character(Result)]) %>%
  select(subject.ID, KRAS) %>%
  rename(SUBJID = subject.ID)
kras$KRAS[kras$KRAS=="NA"] = NA

adtte = sasdats[["adtte"]] %>% mutate(
  AMO = ADY/365*12,   # convert days to months
  TRT = gsub("FP+","",TRT02P,fixed=T)
)
adtteWide = full_join(
  adtte %>%
    filter(PARAMCD == "OS" & COHORT == "Cohort 2") %>% 
    select(SUBJID,AGE,SEX,ARM,TRT,RESINDI,CNSR,AMO) %>%
    rename(OSTTE = AMO, OSCNSR = CNSR),
  adtte %>%
    filter(PARAMCD == "PFSM" & COHORT == "Cohort 2") %>% 
    select(SUBJID,CNSR,AMO) %>%
    rename(PFSTTE = AMO, PFSCNSR = CNSR)
) %>%
  mutate(
    RNR = ifelse(grepl("R",RESINDI),"R","NR")
  )

adk = full_join(kras,adtteWide) %>%
  mutate(BEP = grepl("\\w",KRAS)) #%>% 
  # filter(!is.na(KRAS))
write.csv(adk,"~/Google Drive/crc/MODUL/data/adk200413.csv",na="",row.names=F)
  