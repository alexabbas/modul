# Alex Abbas
# 220304
# get data from manifests for calculating TM/ml

library(tidyverse)
library(readxl)

plasmaDetail = read_excel("/Volumes/GoogleDrive/Shared drives/GI Biomarkers [gRED OBD]/CRC/MODUL Cohort 1/data/manifests/MO29112_FMI PLASMA Tracker Detail_02-09-22_JC.xlsx")
plasmaDetail = plasmaDetail %>%
  rename(
    RocheID = `Roche Specimen Id`,
    Yield = `DNA Yield (ng)`,
    Concentration = `DNA Concentration (ng/uL)`
  ) %>% 
  select(RocheID, Yield, Concentration) %>% 
  mutate()