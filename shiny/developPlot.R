# Alex Abbas
# 240428
# Develop oncoprint for shiny

library(tidyverse)
library(ComplexHeatmap)

setwd("~/Google Drive (alexabbas@gmail.com)/GNE/MODUL/shiny/MODUL/")

oapl = read.csv("data/oapl.csv")
ospl = read.csv("data/ospl.csv")

mapkGenes = c("EGFR","KRAS","NRAS","BRAF","MAP2K1","NF1")

modulPrint = function(oapl, ospl, marks = T, subjGroups = c("C","E"), doTmb = T, doBor = T, genes = mapkGenes) {
  subjs = unique(ospl$SUBJECT.ID)
  subjs = grep(paste(subjGroups,collapse="|"),subjs,value=T)
  subjs = subjs[order(as.numeric(gsub("\\D","",subjs)))]
  
  # remove some variants?
  #  like: oapl = oapl %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "unknown")
  
  # initialize the matrix with samples tested. Call them all wildtype first.
  makeVisitMat = function(avisit) {
    # Make a full data matrix with everything initialized as wildtype
    mat <- matrix("", nrow=length(genes), ncol=length(subjs))
    dimnames(mat) <- list(genes,subjs)
    for( subjectid in unique( ospl %>% filter( Visit %in% avisit ) %>% select( SUBJECT.ID )) ) {
      mat[genes,subjectid] = "wildtype"
    }
    
    # change from wildtype to V600E based on oapl
    pathwayOaplv600e <- oapl %>% filter(GENE %in% genes & Visit %in% avisit & variant == "V600E")
    x <- apply(pathwayOaplv600e,1,function(oa){
      mat[oa["GENE"],oa["SUBJECT.ID"]] <<- "V600E"
    })
    
    # change from wildtype to other variant based on oapl
    pathwayOapl <- oapl %>%
      filter(GENE %in% genes & Visit %in% avisit & variant != "V600E") %>% 
      arrange(desc(VARIANT.TYPE))
    x <- apply(pathwayOapl,1,function(oa){
      mat[oa["GENE"],oa["SUBJECT.ID"]] <<- oa["VARIANT.TYPE"]
    })
    
    # handle multiple variants
    pathway_multiSamples <- unlist(oapl %>%
                                     filter(GENE %in% genes & Visit %in% avisit & variant != "V600E") %>% 
                                     select(GENE,SUBJECT.ID,Visit,VARIANT.TYPE,SAMPLE.ID) %>%
                                     unique() %>% 
                                     group_by(GENE,SUBJECT.ID,Visit,SAMPLE.ID) %>% 
                                     summarise(n = n()) %>% filter(n == 2) %>%
                                     ungroup() %>% 
                                     select(SAMPLE.ID))
    pathwayOapl_multi = oapl %>% 
      filter(GENE %in% genes & Visit %in% avisit & variant != "V600E") %>% 
      filter(SAMPLE.ID %in% pathway_multiSamples)
    x <- apply(pathwayOapl_multi,1,function(oa){
      mat[oa["GENE"],oa["SUBJECT.ID"]] <<- "multiple"
    })
    
    mat
  }
  matTissue = makeVisitMat("Tissue")
  matC1 = makeVisitMat(c("C1D1IP","induction"))
  matPD = makeVisitMat(c("C1D1MP","ontx","EOT","PD"))
  mat = rbind(matTissue,matC1,matPD)
  c1PdRowSplit = factor(
    c(
      rep("Archival\nTissue",times = dim(matTissue)[1]),
      rep("Baseline\nPlasma",times = dim(matC1)[1]),
      rep("Final\nPlasma",times = dim(matPD)[1])
    ),
    levels = c(
      "Archival\nTissue",
      "Baseline\nPlasma",
      "Final\nPlasma"
    )
  )
  
  responseCols <- c(
    "CR" = "#22CC00",
    "PR" = "#DDDD00",
    "SD" = "#DD8811",
    "PD" = "#CC3366"
  )
  
  tissueTmbDat = tmb %>% filter(Visit == "Tissue")
  tissueTmbDat = tissueTmbDat[match(colnames(mat),tissueTmbDat$SUBJECT.ID),]
  patientHeatDat <- patientData[match(colnames(mat),patientData$SUBJECT.ID),] %>% 
    mutate(
      tissueTMB = as.numeric(tissueTmbDat$TMB.SCORE),
      Response = factor(Response,levels = names(responseCols))
    )
  subjSplit = patientHeatDat %>% select(SUBJECT.ID,"Tx") %>% unique()
  splitFactor = factor(
    unlist(subjSplit[match(subjs,subjSplit$SUBJECT.ID),"Tx"]),
    levels = c("Control","Experimental")
  )
  
  oncoCols <- c(
    "short-variant" = "forestgreen",
    "V600E" = "#BBCCBB",
    "copy-number-alteration" = "darkorange",
    rearrangement = "purple",
    multiple = "magenta",
    wildtype = "gray90",
    background = "white"
  )
  
  alterFun <- lapply(names(oncoCols),function(vartype) { function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = oncoCols[vartype], col = NA)) })
  names(alterFun) = names(oncoCols)
  
  annotCols = list(BOR = responseCols)
  
  osplForOncoprint = ospl %>% filter(SUBJECT.ID %in% colnames(mat)) %>% group_by(SUBJECT.ID) %>% slice(1)
  
  leftAnnotArgs = list()
  if (doTmb) { leftAnnotArgs[["TMB\n(Mut/MB)"]] = row_anno_barplot(osplForOncoprint$tmbscore, border = F, axis_param = c(direction="reverse")) }
  if (doBor) {
    leftAnnotArgs[["BOR"]] = osplForOncoprint$Response
    leftAnnotArgs[["col"]] = annotCols
  }
  leftAnnotArgs[["which"]] = "row"
  leftAnnot = do.call(HeatmapAnnotation,leftAnnotArgs)

  op = oncoPrint(t(mat),
                  alter_fun = alterFun,
                  col = oncoCols,
                  show_pct = F,
                  row_split = splitFactor,
                  column_split = c1PdRowSplit,
                  left_annotation = leftAnnot,
                  right_annotation = NULL,
                  top_annotation = NULL,
                  row_order = order(patientHeatDat$Response),
                  show_row_names = T,
                  show_column_names = T,
                  row_names_gp = gpar(fontsize = 9)
  )
  
  if (marks) {
    preTimes = c("Tissue","C1D1IP","induction")
    postTimes = c("C1D1MP","ontx","EOT","PD")
    
    oapl = oapl %>% mutate(subjectGeneVariant = paste(SUBJECT.ID,GENE,variant))
    
    preSGVs = (oapl %>% filter(Visit %in% preTimes))$subjectGeneVariant
    acqPathVars = oapl %>%
      filter(
        GENE %in% genes &
          variant != "V600E" &
          Visit %in% postTimes &
          !( subjectGeneVariant %in% preSGVs)
      ) %>% 
      group_by(SUBJECT.ID,GENE) %>% summarise(variants = paste(unique(variant), collapse = " ")) %>%
      ungroup() %>% 
      mutate(geneVariants = paste(GENE,variants)) %>% 
      group_by(SUBJECT.ID) %>% summarise(geneVariants = paste(geneVariants, collapse = ", ")) %>%
      ungroup() %>% 
      select(SUBJECT.ID,geneVariants) %>% unique()
    
    labels_list = unlist(lapply(acqPathVars$geneVariants, function(x) 
      paste0(strwrap(paste0(x, collapse = " "), width = 45), collapse = "\n")))
    
    variantMarks = anno_mark(
      at = match(acqPathVars$SUBJECT.ID, colnames(mat)),
      labels = labels_list,
      which = "row",
      padding = unit(3, "mm"),
      labels_gp = gpar(fontsize = 10)
    )
    op = op + rowAnnotation(mark = variantMarks)
  }
  
  op
}

modulPrint(oapl, ospl)

###

anOncoprint <- function(dat,splitVar=NA) {
  # genes <- names(sort(table(dat$Gene),decreasing = T))
  #allgenes <- unique(oapl$Gene)
  subjs <- unique(dat$SUBJECT.ID)
  mat <- matrix("",nrow=length(topGenes),ncol=length(subjs))
  dimnames(mat) <- list(topGenes,subjs)
  topOapl <- dat %>% filter(GENE %in% topGenes)
  x <- apply(topOapl,1,function(oa){
    mat[oa["GENE"],oa["SUBJECT.ID"]] <<- oa["VARIANT.TYPE"]
  })
  
  splitFactor <- NA
  if (!is.na(splitVar)) {
    subjSplit = dat %>% select(SUBJECT.ID,!!splitVar) %>% unique()
    splitFactor = subjSplit[match(subjs,subjSplit$SUBJECT.ID),splitVar]
  }
  oncoCols <- c(
    "short-variant" = "forestgreen",
    "copy-number-alteration" = "darkorange",
    rearrangement = "purple",
    # Insertion = "gray",
    # Complex = "firebrick3",
    # "Splice site" = "hotpink",
    # "Stop codon" = "gold",
    background = "white"
  )
  
  alter_fun <- lapply(names(oncoCols),function(vartype) { function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = oncoCols[vartype], col = NA)) })
  names(alter_fun) = names(oncoCols)
  
  oncoPrint(mat,
            alter_fun = alter_fun,
            col = oncoCols,
            show_pct = T,
            column_split = splitFactor,
            right_annotation = rowAnnotation(
              rbar = anno_oncoprint_barplot(
                width = unit(3, "cm"),
                axis_param = list(at = c(0, .25, .5, .75) * dim(mat)[2],
                                  labels = c("0", "25%", "50%", "75%"),
                                  side = "bottom",
                                  labels_rot = 0)))
  )
}
topGenes <- names(sort(table(tissueOapl$GENE),decreasing = T)[1:30])

tissueOapl = oapl %>% filter(Cycle=="Tissue") %>% mutate(Arm = factor(Arm,levels=c("Control","Experimental","Coh99")))
onc1 = anOncoprint(tissueOapl %>% filter(GENE %in% topGenes),splitVar="Arm")
onc1
