# Alex Abbas
# 220218
# R chunk for loading MODUL cohort 1 data

# 220523
#  now points to May 19 FMI data, which includes 13 new on-tx samples
#
# 220218
#  excludes subjects c("2804", "4850", "7842", "4441", "9565", "6403", "7203", "8401")
#  pseudonymizes subjects to C/E numbers


library(tidyverse)
library(readxl)

#setwd("/Volumes/GoogleDrive/Shared drives/GI Biomarkers [gRED OBD]/CRC/MODUL Cohort 1/R/")

### Lesion data
# Read lesion data file Cohort_1_RECIST.csv to get ADT (date), ADY (relative date), and LESDIAM (lesion diameter).

lesionData = read.csv("../data/Cohort_1_RECIST.csv")
lesionData = lesionData %>% mutate(
  ADT = as.Date(ADT, format = "%d%b%Y"),
  ADY = as.numeric(ADY),
  LESDIAM = as.numeric(LESDIAM)
)

#  data_dir="/gne/data/obdroot/avastin/mo29112_modul/ngs/results/ctdna/data"
ctDnaDir <- "../data/ctDNA FMI/"
tissueDir <- "../data/FMI/"
Aug30fmiDir <- "/Volumes/OBD/avastin/mo29112_modul/ngs/rawdata/fmi_obd_mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_std_20210830/"

patsToExclude = c("2804", "4850", "7842", "4441", "9565", "6403", "7203", "8401")
ArmToTx = c( "Experimental" = "Experimental", "Control" = "Control", "Coh99" = "Experimental")

### Patient Data
# Load patient data from MO29112_Ch1 clinical-biomarker dashboard (2).xlsx to get study arm and best overall RECIST response. Exclude patients 2804, 4850, 7842, 4441, 9565, 6403, and 8401 because of insufficient samples.
patientData <- read_excel("../data/MO29112_Ch1 clinical-biomarker dashboard (2).xlsx")
patientData <- patientData %>%
  dplyr::rename(Arm = ARM, Response = BOR, SUBJECT.ID = "Subject ID") %>%
  select(SUBJECT.ID, Arm, Response) %>% 
  mutate(SUBJECT.ID = as.character(SUBJECT.ID), Tx = ArmToTx[Arm])
patientData = patientData %>% filter( !(SUBJECT.ID %in% patsToExclude) )

# Load more patient data from Clinical data cohort1.xlsx to get RECIST status by date.
clinicalData <- read_excel("../data/Clinical data cohort1.xlsx") %>% 
  dplyr::rename(
    SUBJECT.ID = SUBJID
  )
clinicalData = clinicalData %>% filter( !(SUBJECT.ID %in% patsToExclude) )
clinicalData = clinicalData %>% mutate(
  ADT = as.Date(ADT, format = "%d%b%Y"),
  ADY = as.numeric(ADY)
)
clinicalData$ARM[is.na(clinicalData$ARM)] = "X"

### ctDNA data
# Load ctDNA data, both OAPL and OSPL.
# OAPL: mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_n241_one-alteration-per-line-plus_20220125.txt
# OSPL: mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_n241_one-sample-per-line_20220125.xlsx
oapl_ctdna_File <- paste0(
  ifelse ( F & dir.exists(Aug30fmiDir), Aug30fmiDir, ctDnaDir ),
  "mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_n254_one-alteration-per-line-plus_20220519.txt"
)
liqOapl <- read.delim(oapl_ctdna_File,sep="\t") %>% select(-SUBJECT.ID)
ospl_ctdna_file = paste0(ctDnaDir,"mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_n254_one-sample-per-line_20220519.xlsx")
liqOspl <- read_excel(ospl_ctdna_file) %>%
  dplyr::rename(FMI.SAMPLE.ID = "FMI SAMPLE ID", SAMPLE.ID = "SAMPLE ID", SUBJECT.ID = "SUBJECT ID") %>% 
  mutate(SAMPLE.ID = as.character(SAMPLE.ID)) %>% select(-SUBJECT.ID)

# Add cycle and visit data to ctDNA by joining to sample annotation in mo29112_F1L_sample_data_220217.xlsx by internal sample ID (e.g. 6201895392).
sampleMeta_file = paste0(ctDnaDir,"mo29112_F1L_sample_data_220523.xlsx")
sampleMeta = read_excel(sampleMeta_file)
sampleMeta = sampleMeta %>% mutate(SUBJECT.ID = gsub("\\.0","",as.character(SUBJECT.ID)))
#sampleMeta[sampleMeta$Visit=="induction","Visit"] = "C1D1IP"

osplOnlySamples = setdiff(liqOspl$SAMPLE.ID,sampleMeta$SAMPLE.ID)
if (length(osplOnlySamples)>0) {
  stop(paste("Samples not found in mo29112_F1L_sample_data_220523.xlsx:\n",paste(osplOnlySamples,collapse="\n")))
}

liqOspl = liqOspl %>%
  left_join(sampleMeta %>% select(-FMI.SAMPLE.ID))

liqOapl <- liqOapl %>%
  mutate(SAMPLE.ID = as.character(SAMPLE.ID)) %>% 
  left_join(liqOspl %>% select(FMI.SAMPLE.ID, SUBJECT.ID, Cycle, Visit))

### Tissue data
# Similarly load tissue FMI data, both OAPL and OSPL.
# OAPL: mo29112_ngs_dna_targdna_foundationone_20-2041_n65_one-alteration-per-line_20200615.txt
# OSPL: mo29112_ngs_dna_targdna_foundationone_20-2041_n65_one-sample-per-line_20200615.xlsx
oapl_tissue_file <- paste0(
  ifelse ( F & dir.exists(Aug30fmiDir), Aug30fmiDir, tissueDir ),
  "mo29112_ngs_dna_targdna_foundationone_20-2041_n65_one-alteration-per-line_20200615.txt"
)
tissueOapl = read.delim(oapl_tissue_file,sep="\t") %>% 
  mutate(SUBJECT.ID = as.character(SUBJECT.ID))
ospl_tissue_file <- paste0(tissueDir,"mo29112_ngs_dna_targdna_foundationone_20-2041_n65_one-sample-per-line_20200615.xlsx")
tissueOspl <- read_excel(ospl_tissue_file) %>%
  dplyr::rename(SAMPLE.ID = "SAMPLE ID", SUBJECT.ID = "SUBJECT ID") %>% 
  mutate(SUBJECT.ID = as.character(SUBJECT.ID))

# Stack all tissue and plasma FMI data
oapl <- bind_rows(
  tissueOapl %>% mutate(Visit = "Tissue"),
  liqOapl %>% mutate(SAMPLE.ID=as.character(SAMPLE.ID))
) %>%
  mutate(
    variant = gsub("-","",paste0(SV.PROTEIN.CHANGE,CNA.TYPE,REARR.DESCRIPTION)),
    baseEnd = ifelse(grepl("C1|Tissue|induction|IP",Visit),"base",
                     ifelse(grepl("ontx|PD|EOT",Visit),"end",NA)),
    #    baseEnd = ifelse(grepl("C1",Visit),"C1D1",
    #                     ifelse(grepl("PD",Visit),"PD",NA)),
    SUBJECT.ID = gsub("\\.0","",as.character(SUBJECT.ID)),
    TMB.SCORE = as.numeric(TMB.SCORE),
    SV.PERCENT.READS = as.numeric(SV.PERCENT.READS)
  ) %>%
  left_join(
    patientData
  ) %>%
  filter( !(SUBJECT.ID %in% patsToExclude) ) %>% 
  mutate(Arm = factor(Arm,levels = c("Control","Experimental","Coh99")))

ospl = bind_rows(
  liqOspl %>% mutate( SAMPLE.ID = as.character(SAMPLE.ID) ),
  tissueOspl %>% mutate( Visit = "Tissue", SUBJECT.ID = as.character(SUBJECT.ID) )
)

# Join OSPL data to patient data by SUBJECT.ID.
ospl = patientData %>% mutate( SUBJECT.ID = as.character(SUBJECT.ID) ) %>%
  left_join(ospl) %>% 
  mutate(SUBJECT.ID = gsub("\\.0","",SUBJECT.ID),
         baseEnd = ifelse(grepl("C1|Tissue|induction|IP",Visit),"base",
                          ifelse(grepl("ontx|PD|EOT",Visit),"end",NA)))
ospl = ospl %>% filter( !(SUBJECT.ID %in% patsToExclude) & !(is.na(SAMPLE.ID)))

# Joining is complete, so create shorter patient IDs here
revNum = function(x) { as.integer(vapply(lapply(strsplit(as.character(x), "", fixed = TRUE), rev), paste, collapse = "", FUN.VALUE = character(1L))) }

shortArm = c(
  Coh99 = "EPD",
  Control = "C",
  Experimental = "E"
)
subjToNewIDdf = ospl %>%
  select(SUBJECT.ID,Arm) %>% unique() %>% 
  mutate(revid = revNum(SUBJECT.ID)) %>%
  arrange(revid) %>%
  group_by(Arm) %>%
  mutate(n = 1:n()) %>%
  ungroup() %>% 
  mutate(newID = paste0(shortArm[Arm],n))

newIDout = subjToNewIDdf %>% select(SUBJECT.ID,Arm,newID)
#write.csv(newIDout,"/Volumes/GoogleDrive/Shared drives/GI Biomarkers [gRED OBD]/CRC/MODUL Cohort 1/data/oldIDtoNewID.csv")

subjToNewID = setNames(subjToNewIDdf$newID,subjToNewIDdf$SUBJECT.ID)
ospl = ospl %>% mutate(SUBJECT.ID = subjToNewID[SUBJECT.ID])
oapl = oapl %>% mutate(SUBJECT.ID = subjToNewID[SUBJECT.ID])
clinicalData = clinicalData %>% mutate(SUBJECT.ID = subjToNewID[SUBJECT.ID])
patientData = patientData %>% mutate(SUBJECT.ID = subjToNewID[SUBJECT.ID])
lesionData = lesionData %>% mutate(SUBJECT.ID = subjToNewID[as.character(SUBJID)])

# update BOR for two patients
patientData[patientData$SUBJECT.ID %in% c("EPD3","EPD8"),"Response"] = "NE"

# Calculate event times in weeks by subtracting one from the cycle number, multiplying it by 2, and subtracting 16 if it's an IP cycle.
cycleToTime = function (Cycle) {
  #  Cycle[Cycle == "C1D1IP"] = -16
  Cycle[grepl("IP",Cycle)] = ( as.numeric(gsub("C(\\d{1,2})D1IP","\\1",Cycle[grepl("IP",Cycle)])) - 1) * 2 - 16
  Cycle[Cycle == "Tissue"] = -22    # this is a hack to place the tissue "timepoint" in the right place on the plots
  Cycle[grepl("MP",Cycle)] = ( as.numeric(gsub("C(\\d{1,2})D1MP","\\1",Cycle[grepl("MP",Cycle)])) - 1) * 2
  numericCycle = as.numeric(Cycle)
  if(any(is.na(numericCycle))) {
    cycleNA = Cycle[is.na(numericCycle)]
    stop(paste(cycleNA,collapse=" "))
  }
  return(numericCycle)
}
cleanCycle = function(cycleDf) {
  cycleDf[cycleDf$Visit=="Tissue","Cycle"] = "Tissue"
  cycleDf = cycleDf %>% mutate(
    Cycle = gsub("^(C\\d+)1MP$","\\1D1MP",Cycle),
    Cycle = gsub("^(C\\d+D\\d+)$","\\1MP",Cycle)
  )
}
oapl = cleanCycle(oapl)
oapl <- oapl %>% mutate( timepoint = cycleToTime(Cycle) )
ospl = cleanCycle(ospl)
ospl <- ospl %>% mutate( timepoint = cycleToTime(Cycle) )

# complete rows missing an SNV because it wasn't detected at that timepoint
# oapl = oapl %>%
#   group_by(Arm, SUBJECT.ID) %>% complete(nesting(Cycle, Visit), nesting(GENE, variant, VARIANT.TYPE), fill = list(
#   SV.PERCENT.READS = 0)) %>% ungroup()
# #a = oapl2 %>% filter(SUBJECT.ID=="3323")
# oapl[oapl$Visit=="Tissue","Cycle"] = "Tissue"
# oapl <- oapl %>% mutate( timepoint = cycleToTime(Cycle) )

# split MSI and TMB data from variants
msiScore <- oapl %>% filter(VARIANT.TYPE=="BIOMARKER-MSI SCORE")
msiStatus <- oapl %>% filter(VARIANT.TYPE=="BIOMARKER-MSI STATUS")
tmb <- oapl %>% filter(VARIANT.TYPE=="BIOMARKER-TMB")
oaplv <- oapl %>% filter(!grepl("BIOMARKER",VARIANT.TYPE))

# re-subset tissue and liq
tissueOapl <- oaplv %>% filter(Visit == "Tissue")
tissueOaplAll <- tissueOapl
tissueOapl <- tissueOapl %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "unknown")

liqOapl <- oaplv %>% filter(Visit != "Tissue")
liqOaplAll <- liqOapl
liqOapl <- liqOapl %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "unknown")

# Remove unknown impact variants.
#  (now that tissueOaplAll and liqOaplAll are made)
oaplv <- oaplv %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "unknown")

pathwayGenes <- c("EGFR", "KRAS", "NRAS", "BRAF", "MAP2K1", "NF1")
