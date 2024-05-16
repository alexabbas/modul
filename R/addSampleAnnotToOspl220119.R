# Alex Abbas
# 220119
# add our sample annotation to OSPL

library(tidyverse)
library(readxl)

sample_annot = read_excel("~/Desktop/mo29112_ctDNA_sample_annotation_220119.xlsx")
liqOspl = read_excel("~/Desktop/mo29112_ngs_dna_targdna_foundationoneliquidcdx_20-2047_n214_one-sample-per-line_20211214.xlsx")

out = liqOspl %>% full_join(sample_annot %>% select(-'FMI SAMPLE ID'),by="SAMPLE ID")
write.csv(out,"~/Desktop/mo29112_F1Lcdx_n214_ospl_20211214.csv",row.names = F, na="")
