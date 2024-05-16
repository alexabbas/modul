library(shiny)
library(tidyverse)
library(ComplexHeatmap)

oapl = read.csv("data/oapl.csv")
ospl = read.csv("data/ospl.csv")

function(input, output, session) {

  modulPrint = function(oapl, ospl, marks = T, subjGroups = c("C","E"), doTmb = T, doBor = T, selectedGenes = mapkGenes, impact = 2) {
    subjs = unique(ospl$SUBJECT.ID)
    subjs = grep(paste(subjGroups,collapse="|"),subjs,value=T)
    subjs = subjs[order(as.numeric(gsub("\\D","",subjs)))]
    
    # filter variants by impact
    if (impact > 1) {
      oapl = oapl %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "unknown") 
    }
    if (impact > 2) {
      oapl = oapl %>% filter(SOMATIC.STATUS.FUNCTIONAL.IMPACT != "likely") 
    }
    
    # initialize the matrix with samples tested. Call them all wildtype first.
    makeVisitMat = function(avisit) {
      # Make a full data matrix with everything initialized as wildtype
      mat <- matrix("", nrow=length(selectedGenes), ncol=length(subjs))
      dimnames(mat) <- list(selectedGenes,subjs)
      for( subjectid in unique( ospl %>% filter( Visit %in% avisit ) %>% select( SUBJECT.ID )) ) {
        mat[selectedGenes,subjectid] = "wildtype"
      }
      
      # change from wildtype to V600E based on oapl
      pathwayOaplv600e <- oapl %>% filter(GENE %in% selectedGenes & Visit %in% avisit & variant == "V600E")
      x <- apply(pathwayOaplv600e,1,function(oa){
        mat[oa["GENE"],oa["SUBJECT.ID"]] <<- "V600E"
      })
      
      # change from wildtype to other variant based on oapl
      pathwayOapl <- oapl %>%
        filter(GENE %in% selectedGenes & Visit %in% avisit & variant != "V600E") %>% 
        arrange(desc(VARIANT.TYPE))
      x <- apply(pathwayOapl,1,function(oa){
        mat[oa["GENE"],oa["SUBJECT.ID"]] <<- oa["VARIANT.TYPE"]
      })
      
      # handle multiple variants
      pathway_multiSamples <- unlist(oapl %>%
                                       filter(GENE %in% selectedGenes & Visit %in% avisit & variant != "V600E") %>% 
                                       select(GENE,SUBJECT.ID,Visit,VARIANT.TYPE,SAMPLE.ID) %>%
                                       unique() %>% 
                                       group_by(GENE,SUBJECT.ID,Visit,SAMPLE.ID) %>% 
                                       summarise(n = n()) %>% filter(n == 2) %>%
                                       ungroup() %>% 
                                       select(SAMPLE.ID))
      pathwayOapl_multi = oapl %>% 
        filter(GENE %in% selectedGenes & Visit %in% avisit & variant != "V600E") %>% 
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
    
    ospl = ospl %>% mutate(Response = factor(Response,levels = names(responseCols)))
    osplForOncoprint = ospl %>% filter(SUBJECT.ID %in% colnames(mat)) %>% group_by(SUBJECT.ID) %>% slice(1)

    osplForOncoprint = osplForOncoprint[match(colnames(mat),osplForOncoprint$SUBJECT.ID),]
    subjSplit = osplForOncoprint %>% select(SUBJECT.ID,"Tx") %>% unique()
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
                   row_order = order(osplForOncoprint$Response),
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
          GENE %in% selectedGenes &
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
  
  output$oncoprint <- renderPlot({
      modulPrint(oapl, ospl,
               doBor = input$doBor, doTmb = input$doTmb, selectedGenes = input$genes,
               marks = input$marks, impact = input$impact)
  }, height = 700)
  
  
}
