# Alex Abbas
# 240428
# Develop oncoprint for shiny

library(tidyverse)
library(ComplexHeatmap)

setwd("~/Google Drive (alexabbas@gmail.com)/GNE/MODUL/shiny/MODUL/")

oapl = read.csv("data/oapl.csv")
ospl = read.csv("data/ospl.csv")
clin = read.csv("data/clin.csv")
pat = read.csv("data/pat.csv")

genes = c("EGFR","KRAS","NRAS","BRAF","MAP2K1","NF1")

oapl = oapl %>% 
  mutate(
    mapkMut = GENE %in% genes & variant != "V600E"
  )
osplPath = oapl %>%
  select(SUBJECT.ID,SAMPLE.ID,mapkMut,Cycle,Visit,timepoint,Tx) %>% 
  group_by(SUBJECT.ID,Cycle,Visit,timepoint,Tx,SAMPLE.ID) %>%
  summarise(state = ifelse(any(mapkMut),"Positive","Negative")) %>% 
  mutate(
    timepointMonth = timepoint / 52 * 12,
    timepointDay = timepoint / 52 * 365,
    Event = ifelse(Visit=="Tissue","Archival Biopsy","ctDNA")
  ) %>% 
  ungroup()

# fill in events with no mutations in oapl
toAdd = ospl %>%
  filter( !(SAMPLE.ID %in% osplPath$SAMPLE.ID) ) %>%
  select( SUBJECT.ID, Cycle, Visit, Tx, timepoint) %>% 
  dplyr::mutate(
    timepointMonth = timepoint / 52 * 12,
    timepointDay = timepoint / 52 * 365,
    Event = "ctDNA",
    state = "Negative"
  )
osplPath = osplPath %>% bind_rows(toAdd)

# add clinical timecourse response data
clinToBind = clin %>% 
  dplyr::rename(
    state = AVALC
  ) %>% 
  mutate(
    timepointMonth = as.numeric(ADY) / 365 * 12,
    timepointDay = as.numeric(ADY),
    Event = "Scan"
  ) %>% 
  filter(
    SUBJECT.ID %in% osplPath$SUBJECT.ID
  ) %>% 
  left_join(osplPath %>% select(SUBJECT.ID,Tx) %>% unique()) %>% 
  mutate(
    PostTx = ifelse(APHASE == "Post-trt" | Event == "ctDNA", "yes", "no")
  )
osplClin = bind_rows(clinToBind,osplPath)

# make a swimmer plot with it
strokeColors = c(
  "Negative" = "gray20",
  "Positive" = "gray20",
  #  "Submitted" = "#4488FF",
  "yes" = "gray20",
  "no" = "#00000000",
  "CR" = "#22CC00",
  "PR" = "#DDDD00",
  "SD" = "#DD8811",
  "PD" = "#CC3366"
)

tdf = "triangle down filled"
tf = "triangle filled"
sq = "square filled"
cf = "circle filled"

eventPch = c(
  "Archival Biopsy" = cf,
  "ctDNA" = tf,
  "Scan" = sq
)

stateFills = c(
  "Negative" = "white",
  "Positive" = "gray20",
  "Submitted" = "#4488FF",
  "Complete Response (CR)" = "#22CC00",
  "Partial Response (PR)" = "#DDDD00",
  "Stable Disease (SD)" = "#DD8811",
  "Progressive Disease (PD)" = "#CC3366",
  "Not Evaluable (NE)" = "gray50"
)

osplClin$state = factor(osplClin$state, levels = names(stateFills))
osplClin = osplClin %>% left_join(patientData)

osplClin = osplClin %>% mutate(timepointDay2 = timepointDay / 15.5)
osplClin = osplClin %>% mutate(Response = factor(Response, levels = names(strokeColors)))
osplClin = osplClin %>% mutate(SUBJECT.ID = factor(SUBJECT.ID, levels = unique(osplClin$SUBJECT.ID[order(osplClin$Response)])))

gp <- ggplot(osplClin, aes(x = timepointDay2, y = SUBJECT.ID, fill = state, alpha = state, col = PostTx, shape = Event, group = SUBJECT.ID)) +
  geom_line(aes(color = Response), size = 1.5) +
  geom_point(data = osplClin %>% filter(Event == "ctDNA"), size = 1.9, position = position_nudge( y = .2 ) ) +
  geom_point(data = osplClin %>% filter(Event == "Archival Biopsy"), size = 2.5 ) +
  geom_point(data = osplClin %>% filter(Event == "Scan"), size = 2.3, position = position_nudge( y = -.18 ) ) +
  theme_classic() +
  facet_grid(rows = "Tx", space = "free", scales = "free" ) +
  xlab("Time relative to start of maintenance phase (days)") +
  scale_color_manual(
    values = strokeColors,
    breaks = c("yes","no"),
    guide = guide_legend(override.aes = list(shape = sq, linetype = 0, fill = "gray80"), title = "Post-tx")) +
  scale_fill_manual(
    values = stateFills,
    name = "RECIST",
    breaks = c(  "Complete Response (CR)", "Partial Response (PR)", "Stable Disease (SD)", "Progressive Disease (PD)", "Not Evaluable (NE)"),
    guide = guide_legend(override.aes = list(color = "gray50", linetype = 1, shape = sq) ) ) +
  scale_shape_manual(values = eventPch) +
  scale_alpha_manual(
    name = "Mutation status",
    values = c(1, 1, 1, 1, 1, 1, 1, 1),
    breaks = c("Positive", "Negative"),#, "Submitted"),
    guide = guide_legend(
      override.aes = list(
        linetype = 0,
        shape = tf,
        fill = c("gray20","white"),#,"#4488FF"),
        color = c("gray20","gray20")#,"#4488FF")
      ))) #+
gp
