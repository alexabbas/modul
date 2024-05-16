# Alex Abbas
# 220506
# Make Table 1 for MODUL paper

library(tidyverse)
library(table1)

# setwd("/Volumes/GoogleDrive/Shared drives/GI Biomarkers [gRED OBD]/CRC/MODUL Cohort 1/results/")

source("../R/loadData220218.R")

clinicalData = clinicalData %>%
  select(SUBJECT.ID,AGE,SEX,ARM) %>% 
  unique() %>% 
  mutate(
    BEP = factor(ifelse(SUBJECT.ID %in% ospl$SUBJECT.ID,"BEP","non-BEP")),
    AGE = as.numeric(AGE)
  ) %>%
  filter(
    ARM != "X"
  ) %>% 
    rename(
      Sex = SEX,
      Age = AGE
    )

units(clinicalData$Age) = "years"
table1(~ Sex + Age | BEP + ARM, data = clinicalData, topclass="Rtable1-zebra")

# strata = c(list(Total=clinicalData), split(clinicalData, clinicalData$BEP))
# labels = list(variables=list(Sex="Sex",Age="Age (years)"))
# table1(strata, labels)
