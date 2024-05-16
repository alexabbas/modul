# Alex Abbas
# 240428
# Exploration of tumor mutation data from the MODUL trial

library(shiny)

topgenes = c("APC","ATM","BRAF","BRCA2","MAP2K1","MAP3K1","MLL2","MSH3","MTOR","NF1","NOTCH1","NTRK1","NTRK3","PIK3CA","SMAD4","TP53")
mapkGenes = c("EGFR","KRAS","NRAS","BRAF","MAP2K1","NF1")
genes = unique(c(topgenes,mapkGenes))

fluidPage(

    h1(id="big-heading", "MODUL Cohort 1: Tumor Mutation Exploration"),
    tags$style(HTML("#big-heading{color: white; background-color:#4d3a7d; padding:30px}")),

        p(id="sub-heading", "MODUL (MO29112, NCT02291289) was an adaptable, phase 2, signal-seeking trial testing novel agents as first-line therapy for predefined subgroups of patients with metastatic colorectal cancer."),
    p(id="sh2", "Here we examine cohort 1, subjects with BRAF V600E mutant tumors, randomized to either the experimental arm of 5-FU/LV + cetuximab + vemurafenib or a control arm of FP + bevacizumab. Subjects' tumor tissue was tested with the FMI Foundation One CDx assay, and plasma samples from C1D1 and/or PD/EOT were tested with the FMI Foundation One Liquid CDx assay. Cetuximab treatment put a negative selective pressure on V600E mutations and a positive selective pressure on compensatory mutations elsewhere in the EGFR pathway."),
    p(id="sh3", "This tool's output is similar to figure 3a that I made for the study publication."),
    
    tags$h4("Publication:"),
    tags$div(
      tags$a(id="link",href="https://pubmed.ncbi.nlm.nih.gov/36921494/", 
             "Clinical and exploratory biomarker findings from the MODUL trial (Cohorts 1, 3 and 4) of biomarker-driven maintenance therapy for metastatic colorectal cancer. Eur J Cancer. 2023 May;184:137-150. doi: 10.1016/j.ejca.2023.01.023. Epub 2023 Feb 4. PMID: 36921494.")
    ),p(),
    
    # Sidebar
    sidebarLayout(
        sidebarPanel(fluid = FALSE,
                      checkboxInput("doBor", "Best Overall Response (BOR)", value = TRUE),
                      checkboxInput("doTmb", "Tumor Mutational Burden (TMB)", value = TRUE),
                      checkboxInput("marks", "Mark Acquired Variants", value = TRUE),
                      selectInput("impact", "Filter genes by impact:", choices = c("All"=1,"Known/Likely"=2,"Known"=3), selected=2),
                      checkboxGroupInput("genes", "Genes", genes, mapkGenes)
        ),
        # Oncoprint
        mainPanel(
            plotOutput("oncoprint")
        )
    )
)
