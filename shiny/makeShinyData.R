# Alex Abbas
# 240428
# Make data for MODUL shiny app

library(tidyverse)

setwd("/Volumes/Macintosh HD/Users/alex/Google Drive (alexabbas@gmail.com)/GNE/MODUL/shiny/")

source("../R/loadData220218.R")

oapl = oaplv1
write.csv(oapl,"MODUL/data/oapl.csv")

ospl = ospl %>% rename(tmbscore = `TMB SCORE`) %>% mutate(tmbscore = as.numeric(tmbscore))
write.csv(ospl,"MODUL/data/ospl.csv")

clinicalData = clinicalData %>% select(SUBJECT.ID,ADY,AVALC,APHASE)
write.csv(clinicalData,"MODUL/data/clin.csv")

write.csv(patientData,"MODUL/data/pat.csv")
