


####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","BiocManager",
              "dplyr","DT","ggplot2","ggpubr","reshape2","tibble",
              "plotly","readr","enrichR","ggrepel","tidyr","tools","shinycssloaders")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
}
#bioconductor packages
bioCpacks <- c("clusterProfiler","GSVA","limma","enrichplot")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
    BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(packages, library, character.only = TRUE))
invisible(lapply(bioCpacks, library, character.only = TRUE))




####----User Data Input----####

#Input desired project name for webpage - will be followed by 'Expression Analysis'
ProjectName <- "CCLE Cell Line Analysis"


##--User Input File Names--##

#expression data
expr_file <- "htseq_gene_level_fpkm_T_geneName_max_1cutoff_v2.txt"

#meta data
meta_file <- "USP7_meta.tsv"
#Is there a header?
header <- TRUE

#Enriched Signatures data table
ES_tables <- c()

#If human: set TRUE
#If mouse: set FALSE
human <- TRUE

##--User Gene Set Input--##

#write in the name of your gene set list for shiny UI
userGSlist_name <- 'CellMarker Gene Sets'

#path to your gene set file .gmt or .txt/.tsv
userGS_file <- 'CellMarker_gsNsym_HS.tsv'
#Does gene set file have header?
header.gs <- TRUE

#path to your R data list object for ssGSEA
userRData_file <- 'CellMarker_GS_HS.RData'
#userRData_file <- 'CellMarker_GS_HS.RData'


## ADD CCLE input
#CCLE Master Meta
ccle_meta_file <- 'MasterMeta_v2.tsv'
#CCLE expression
ccle_expr_file <- 'CCLE_expression.tsv'

####----Read Files----####

#read ccle_meta
ccle_meta <- as.data.frame(read_tsv(ccle_meta_file))
#read ccle_expression
ccle_expr <- as.data.frame(read_csv(ccle_expr_file))
rownames(ccle_expr) <- ccle_expr[,1]
ccle_expr <- ccle_expr[,-1]

####----Selection Generation----####

#lineage selection
#ccle_lineage_choice <- "NoneSelected"
#ccle_lineage_choice <- c(ccle_lineage_choice, unique(ccle_meta$lineage))
ccle_lineage_choice <- unique(ccle_meta$lineage)
#ccle_lineage_choice <- c(ccle_lineage_choice,"All")

ccle_genecols <- colnames(ccle_meta)[c(15:ncol(ccle_meta))]
ccle_genecols <- "NoneSelected"
ccle_genecols <- c(ccle_genecols,gsub("_.*","",ccle_genecols))

#ccle_disease_choice <- "NoneSelected"
#ccle_disease_choice <- c(ccle_disease_choice, unique(ccle_meta$primary_disease))
ccle_disease_choice <- unique(ccle_meta$primary_disease)




####----Backend Data Input----####


if (human == TRUE) {
    #MSigDB gene set
    msigdb <- 'msigdb_gsNsym_HS.tsv'
    #MSigDB gene set FOR UI
    msigdb2 <- 'msigdb_gsNcat_HS.tsv'
    #gene set list for ssGSEA
    load('gs_list_HS.RData')
    #Cytokine genes for human
    CTKgenes <- c("IL2","IL12A","IL12B","IL17A","IFNA1","IFNB1","IFNG","IFNGR","CD11b",
                  "ITGAM","CD33","ENTPD1","ICOSLG","CD275","CD278","TNFSF9","TNFRSF9",
                  "CD40","CD40LG","CD70","CD27","TNFSF18","TNFRSF18","TNFSF14","TNFRSF14",
                  "TNFSF4","TNFRSF4","HLA-A","CD3","CEACAM1","CD80","CD86","CTLA4","CD276",
                  "VTCN1","PVR","CD226","TIGIT","CD96","LGALS3","LGALS3BP","LGALS9","LGALS9C",
                  "HAVCR2","HHLA2","TMIGD2","CD274","PDCD1LG2","PDCD1","VSIR")
}
if (human == FALSE) {
    #MSigDB gene set
    msigdb <- 'msigdb_gsNsym_MM.tsv'
    #MSigDB gene set FOR UI
    msigdb2 <- 'msigdb_gsNcat_MM.tsv'
    #gene set list for ssGSEA 
    load('gs_list_MM.RData')
    #Cytokiny genes for mouse
    CTKgenes <- c("Il2","Il12a","Il12b","Il17a","Ifna13","Ifnb1","Ifng","Ifngr1","Cd11b","Itgam",
                  "Cd33","Entpd1","Icosl","Icos","Tnfsf9","Tnfrsf9","Cd40","Cd40lg","Cd70","Cd27",
                  "Tnfsf18","Tnfrsf18","Tnfsf14","Tnfrsf14","Tnfsf4","Tnfrsf4","H2-K1","CD3G",
                  "Ceacam1","Cd80","Cd86","Ctla4","Cd276","Vtcn1","Pvr","Cd226","Tigit","Cd96","Lgals3",
                  "Lgals3bp","Lgals9","Lgals9c","Havcr2","Hhla2","Cd274","Pdcd1lg2","Pdcd1","Vsir")
}




####----Read and Manipulate Files----####


##read files

##reading expression data
expr <- read.delim(expr_file, sep = '\t', header = T, strip.white = T)
colnames(expr)[1] <- "Gene"
expr <- expr %>%
    drop_na()
row.names(expr) <- make.names(expr[,1], unique = T)
expr <- expr[,-1]
colnames(expr) <- gsub("[_.-]", "_", colnames(expr))
A <- as.matrix(expr)

#gene list file from expression data
Gene <- rownames(expr)
geneList <- as.data.frame(Gene)


#meta
meta <- read.delim(meta_file, sep = '\t', header = header, strip.white = T)
meta[,1] <- gsub("[_.-]", "_", meta[,1])
colnames(meta) <- c("SampleName","Group")
metagroups <- as.vector(levels(factor(meta[,2])))



#for heatmap sample selection
sampsames <- intersect(colnames(expr),meta[,1])
#ensure expression samples and meta are exact
expr <- expr[,sampsames]
meta <- meta[which(meta[,1] %in% sampsames),]


#boxplot choices based on meta groups
if (length(metagroups) == 2) {
    boxopt <- c("wilcox.test", "t.test", "none")
}
if (length(metagroups) >= 3) {
    boxopt <- c("kruskal.test", "anova", "none")
}




#Enriched Signatures
if (!is.null(ES_tables)) {
    ldf <- list()
    SigNames <- c()
    ES_Tab_List <- c()
    n=1
    for (k in 1:length(ES_tables)){
        file <- basename(ES_tables[k])
        file2 <- gsub("\\..*$","",file)
        SigNames <- c(SigNames, file2)
        ldf[[k]] <- as.data.frame(read_delim(ES_tables[k], delim = '\t'))
        ES_Tab_List <- c(ES_Tab_List,paste("ES_table",n,sep = ""))
        n=n+1
    }
    j <- 1
    for (i in 1:length(ldf)){
        names(ldf)[i] <- paste("ES_table",j,sep = "")
        j=j+1
    }
    list2env(ldf,globalenv())
}

#load blank enriched signatures table if none given by user
if (is.null(ES_tables)) {
    SigNames <- "None Loaded"
    ES_Tab_List <- c()
    n=1
    ES_table1 <- data.frame("ID" = "No Enriched Signature Tables loaded","Description" = "NA","setSize" = "NA",
                            "enrichmentScore" = "NA","NES" = "NA","pvalue" = "NA","p.adjust" = "NA","qvalues" = "NA",
                            "rank" = "NA","leading_edge" = "NA","core_enrichment" = "NA")
    ES_Tab_List <- c(ES_Tab_List,paste("ES_table",n,sep = ""))
    n=n+1
}


#MSigDB gene sets
msigdb.gsea <- read.delim(msigdb, header = T, sep = '\t')
gmt <- msigdb.gsea
#MSigDB gene sets FOR UI
msigdb.gsea2 <- read.delim(msigdb2, header = T, sep = '\t')


#tab2 User gene set
if (file_ext(userGS_file) == "gmt") {
    tab2 <- read.gmt(userGS_file)
}
if (file_ext(userGS_file) == "tsv" || file_ext(userGS_file) == "txt") {
    tab2 <- read.delim(userGS_file, header = header.gs, sep = '\t')
}
#tab2 back end
GeneSet2 <- as.data.frame(unique(tab2[,1]))
rownames(GeneSet2) <- 1:nrow(GeneSet2)
colnames(GeneSet2)[1] <- "Gene_Set"
#tab2 R Data list
loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}
gs2 <- loadRData(userRData_file)


#CV function for variance
cv <- function(x){
    (sd(x)/mean(x))*100
}

#enrichR pathway choices
enrichRtab <- listEnrichrDbs()
enrichRchoice <- enrichRtab[,"libraryName"]


####----Shiny UI----####


shinytheme("sandstone")

ui <-
    
    navbarPage(paste("{",ProjectName,"Expression Analysis }", sep=" "),
               
               ####----Intro Tab----####

               tabPanel("Intro/Methodology",
                        fluidPage(
                            mainPanel(
                                h3("Introduction"),
                                #p("With the rapid generation of large datasets as a result of advancing next-generation sequecing (NGS) technologies, the need for quick and reproducible data analysis tools has become paramount. The ability to analyze both genomic and proteomic data sets can allow for the detection of compelling genes and features that may be of interest for further analysis. The General Expression Analysis App hosted in R Shiny does not require an extensive computation background to run these analyses. With user input of expression and meta data, along with gene set we have sourced and provided, the user may produce useful anaylsis and visualizations on a web-based interface within minutes. This method of reproducible data analysis generates a number of visualization to view Gene Set Enrichment (GSEA) and differential gene expression analysis, as well as the ability to download various tables and gene sets that are produced by the App for further use. Additional apps are being developed for the EASY family. DRPPM-EASY-Integraction allows users to further analyze data obtained from the main DRPPM-EASY app and compare expression data between two matrices. DRPPM-EASY-CCLE integrates a sample selection tab which allows users to select expression and meta data from the Cancer Cell Line Encyclopedia (CCLE) based on cancer type or lineage and sample type for analysis within the DRPPM-EASY app. The flow chart below gives a layout of the EASY app's family infrastructure which we will describe in further detail."),
                                p("This is an extention of the DRPPM Expression Analysis ShinY (EASY) App which integrates data from the Cancer Cell Line Encyclopedia (CCLE) for sample selection perform differential gene expression and gene set enrichment analyses. This R Shiny app is very similar in features to the main app, expect for the addition of the first tab which allows for CCLE sample selection. Based on this selection, the expression and meta data will be subset and imported in the back end to the app for further analysis and visualization. Below is a schematic for the DRPPM-EASY-CCLE pipeline."),
                                #img(src=flowchart_file, align = "left"),
                                #imageOutput("flowchartimage"),
                                h3("Unsupervised Clustering Methods"),
                                p("Unsupervised clustering is performed by calculating the top variable genes through MAD, CV, or VAR as the variance measure. The choice of variance measure and number to top genes/probes can be determined by the user, as well as the number of clusters with k-means. The identified most variable genes can be downloaded as a table for the user. The “complete” method was chosen as the default algorithm for clustering, but other methods such as Ward, average, and centroid could be chosen based on the ‘hclust()’ function in R."),
                                h3("Differential Gene Expression Methods"),
                                p("Differential gene expression analysis is performed on the expression data between two groups defined by the provided metadata file. The samples chosen are log-transformed (log2 + 1). A model matrix is then designed for LIMMA linear modeling followed by empirical Bayes statistics to moderate gene variance and modeling of the global characteristic of the data. Of note, the current pipeline expects the matrix input quantification are pre-normalized, such as in CPM, TMM, or FPKM."),
                                h3("Gene Set Enrichment Analysis (GSEA) Methods"),
                                p("Gene Set Enrichment Analysis (GSEA) is performed through signal-to-noise calculations on the expression matrix. This calculation requires at least two types of sample groups and at least three samples for each grouping. The signal-to-noise calculation scales the difference of means by standard deviation and generates a ranked list of genes, effectively capturing the differences between phenotypes. The GSEA identifies the enrichment score (ES) for a gene set, which indicates the extent that the gene set is overrepresented at the beginning or end of the list of ranked genes. The enrichment score is calculated by walking down the ranked gene list and increasing the running-sum statistic when a gene is in the gene set and decreasing when it is not. The maximum deviation from zero that is encountered through walking the list is the ES, where a positive ES indicates the gene set is enriched at the top of the ranked list and a negative ES indicates the gene set is enriched at the bottom of the ranked list. A leading-edge gene list is provided displaying a subset of the genes in the gene set that contribute the most to the ES."),
                                h4("Gene Set Sources"),
                                p("The Molecular Signatures Database (MSigDB) is a collection of over 32,000 annotated gene sets divided into 9 major collections as well as various sub-collections [PMID: 16199517, PMID: 12808457]. The MSigDB  gene sets were downloaded with the msigdbr package and processed through R using Tidyverse packages to combine the gene sets into a data frame."),
                                p("The Cell Marker gene sets were derived from a comprehensive list of cell markers from cell types in various tissues obtained through manually curating over 100,000 published papers [PMID: 30289549]. This data was obtained as an extensive data frame and parsed to make gene sets with the tidyverse packages."),
                                p("The Library of Integrated Network-based Cellular Signatures (LINCS) L1000 gene sets were derived from the Connectivity Map (CMAP) projects data which treated 4 human cell lines with ~1300 drugs followed by profiling genome-wide mRNA expression following small molecule perturbation [PMID: 24906883].")
                                #p("Mapped reads were derived from BBSR and likely normalized to CPM. Pathway enrichment analysis was performed by GSEA [PMID: 16199517] and enrichR[PMID: 23586463]. Single-sample gene set enrichment analysis (ssGSEA) [PMID: 19847166, PMID: 30595505] was used to quantify the expression signatures. LIMMA was used to define differentially expressed genes [PMID: PMID: 25605792]"))
                                )
                            )
                        ),


	       tabPanel("Sample Selection",
                   fluidPage(
                       mainPanel(
                           fluidRow(
                               column(6,
                                      h4("Subset CCLE Samples by Lineage or Disease type:"),
                                      radioButtons("linordischeck","",
                                                         choices = list("Lineage" = 1,
                                                                        "Disease Type" = 2)),
                                      uiOutput("subsetselec"),
                                      uiOutput("sublinselec")
                                      ),
                               column(6,
                                      h4("Condition Selection:"),
                                      radioButtons("condseleccheck","",
                                                         choices = list("Sex" = 1,
                                                                        #"Sample Source" = 2,
                                                                        "Sample Collection Site" = 3,
                                                                        "Primary or Metastasis" = 4,
                                                                        "Choose Gene of Interest" = 5)),
                                      uiOutput("condselec")
                                      )
                           ),
                           fluidRow(
                               #column(6,
                                      #textInput("ccle_exprFileName","Expression File Download Name:",
                                       #         value = "ExpressionData")),
                                      #downloadButton("downloadEXPR","Download Expression Matrix"),
				      #actionButton("updateEXPR", "UpdateEXPR")),
                               column(12,
                                      #textInput("ccle_metaFileName","Meta File Download Name:",
                                      #          value = "MetaData"),
                                      #downloadButton("downloadMETA","Download Meta Data"),
				      actionButton("updateMETA", "Load/Update Selected Data")
                               )

                           ),
			   fluidRow(

                               div(DT::dataTableOutput("testtable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                               div(DT::dataTableOutput("testexprtable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
			   )
                       )
                     )
               ),
		

               #tabPanel("Intro/Methodology",
               #         fluidPage(
               #             mainPanel(
               #                 h2("GSEA and RNAseq Analysis Method"),
               #                 p("Mapped reads were derived from BBSR and likely normalized to CPM. Pathway enrichment analysis was performed by GSEA [PMID: 16199517] and enrichR[PMID: 23586463]. Single-sample gene set enrichment analysis (ssGSEA) [PMID: 19847166, PMID: 30595505] was used to quantify the expression signatures. LIMMA was used to define differentially expressed genes [PMID: PMID: 25605792]"))
               #         )
               #),
               
               ####----RNAseq Tab----####
               
               tabPanel("Standard Analysis",
                        fluidPage(
                            title = "Standard Analysis",
                            sidebarLayout(
                                sidebarPanel(
                                    
                                    #heatmap side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '22'",
                                                     tabsetPanel(id = "customs",
                                                                 tabPanel("Parameters",
                                                                          p(),
                                                                          numericInput("NumFeatures", step = 1, label = "Number of Genes", value = 100),
                                                                          h4("Clustering Parameters"),
                                                                          selectInput("VarianceMeasure", "Select Variance Measure",
                                                                                      choices = c("MAD","CV","VAR")),
                                                                          selectInput("ClusteringMethod",
                                                                                      "Select clustering Method",
                                                                                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                                                          numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 2),
                                                                          h4("Download Cluster Result:"),
                                                                          downloadButton("downloadClusters", "Download .tsv"),
                                                                          h4("Download Most Variable Genes List:"),
                                                                          downloadButton("MVGdownload", "Download MVG .tsv"),
                                                                          downloadButton("MVGdownloadgmt", "Download MVG .gmt"),
                                                                          p(),
                                                                          sliderInput("heatmapFont2.c", "Heatmap Font Size (columns):",
                                                                                      min = 5, max = 75,
                                                                                      value = 12, step = 1),
                                                                          sliderInput("heatmapFont2.r", "Heatmap Font Size (rows):",
                                                                                      min = 5, max = 75,
                                                                                      value = 9, step = 1),
                                                                          value = 222
                                                                 ),
                                                                 tabPanel("Custom Heatmap",
                                                                          p(),
                                                                          selectizeInput("heatmapGeneSelec","Gene Selection:",
                                                                                         choices = sort(as.vector(geneList[,1])),
                                                                                         multiple = T, selected = CTKgenes),
                                                                          textInput("userheatgenes", "Text Input of Gene List (space or tab delimited):", value = ""),
                                                                          selectInput("ClusteringMethod2",
                                                                                      "Select clustering Method",
                                                                                      choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                                                          selectizeInput("userheatsamp2", "Samples Selection:",
                                                                                         choices = sampsames, multiple = T, selected = sampsames),
                                                                          sliderInput("heatmapFont3.c", "Heatmap Font Size (columns):",
                                                                                      min = 5, max = 75,
                                                                                      value = 12, step = 1),
                                                                          sliderInput("heatmapFont3.r", "Heatmap Font Size (rows):",
                                                                                      min = 5, max = 75,
                                                                                      value = 12, step = 1),
                                                                          value = 444
                                                                 )
                                                     )
                                    ),
                                    
                                    #Volcano and MA plot side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '44' || input.dataset == '66'",
                                                     p(),
                                                     selectInput("comparisonA2", "Comparison: GroupA",
                                                                 choices = metagroups, selected = metagroups[1]),
                                                     selectInput("comparisonB2", "Comparison: GroupB",
                                                                 choices = metagroups, selected = metagroups[2]),
                                                     numericInput("fc_cutoff", "LogFC Threshold (Absolute Value)",
                                                                  min = 0, max = 5, step = 0.1, value = 1),
                                                     numericInput("p_cutoff", "Significance Threshold (-log10(P.Value)):",
                                                                  min = 0, max = 10, step = 0.1, value = 0.05),
                                                     numericInput("top_x", "Number of Top Hits:", value = 10,
                                                                  min = 0),
                                                     selectizeInput("userGeneSelec", "User Selected Hits:",
                                                                    choices = sort(as.vector(geneList[,1])), multiple = T, selected = "-"),
                                                     textInput("userGeneSelec2", "Text Input of Gene List (space or tab delimited):", value = "")
                                    ),
                                    
                                    #Boxplot side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '88'",
                                                     p(),
                                                     h3("Select A Gene:"),
                                                     div(DT::dataTableOutput("GeneListTable"), style = "font-size:10px"),
                                                     p(),
                                                     sliderInput("boxplot.font", "Font Size:",
                                                                 min = 5, max = 75,
                                                                 value = 15, step = 1),
                                                     sliderInput("boxplotDot", "Box Plot Dot Size:",
                                                                 min = 0, max = 5,
                                                                 value = 0.75, step = .25)
                                    ),
                                    
                                    #Gene scatter plot side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '26'",
                                                     p(),
                                                     selectInput("scatterG1","Select Gene 1",
                                                                 choices = Gene),
                                                     selectInput("scatterG2","Select Gene 2",
                                                                 choices = Gene, selected = Gene[2]),
                                                     checkboxInput("logask","Log2 Transform Expression Data")
                                    ),
                                    
                                    #EnrichR pathway side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '28'",
                                                     selectInput("comparisonA2.path", "Comparison: GroupA",
                                                                 choices = metagroups, selected = metagroups[1]),
                                                     selectInput("comparisonB2.path", "Comparison: GroupB",
                                                                 choices = metagroups, selected = metagroups[2]),
                                                     numericInput("pathpval", "Adjusted P-Value Cutoff:",
                                                                  min = 0, value = 0.05),
                                                     numericInput("pathFC", "Log2 Fold Change Cutoff:",
                                                                  min = 0, value = 1),
                                                     selectInput("SelectedPathway", "Select Pathway",
                                                                 choices = enrichRchoice)
                                    ),
                                    
                                    #DEG table side panel
                                    
                                    conditionalPanel(condition = "input.dataset == '24'",
                                                     p(),
                                                     selectInput("comparisonA2.DEG", "Comparison: GroupA",
                                                                 choices = metagroups, selected = metagroups[1]),
                                                     selectInput("comparisonB2.DEG", "Comparison: GroupB",
                                                                 choices = metagroups, selected = metagroups[2]),
                                                     h4("Download DEG Table as GMT File"),
                                                     textInput("DEGfileName", "File Name for Download:",value = "DEGgeneSet"),
                                                     numericInput("fc_cutoff2", "LogFC Threshold",
                                                                  min = 0, max = 5, step = 0.1, value = 1),
                                                     selectInput("UpDnChoice","Up-regulated or Down-regulated:",
                                                                 choices = c("UpAndDown_Regulated","Up_Regulated","Down_Regulated")),
                                                     numericInput("p_cutoff2", "Adj.P.Value Cutoff:",
                                                                  min = 0, max = 10, step = 0.1, value = 0.05),
                                                     numericInput("top_x2", "Number of Top Hits:", value = 100),
                                                     downloadButton("DEGgmtDownload", "Download DEG .gmt")
                                    ),
                                ),
                                mainPanel(
                                    tabsetPanel(
                                        id = "dataset",
                                        tabPanel("Heatmap",
                                                 withSpinner(jqui_resizable(plotOutput("heatmap1", width = "100%", height = "1000px")), type = 6),
                                                 value = 22),
                                        tabPanel("Volcano Plot",
                                                 p(),
                                                 verbatimTextOutput("VolGroupsText"),
                                                 jqui_resizable(plotOutput('Volcano3', width = "800px", height = "550px",
                                                                           hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))),
                                                 uiOutput("hover_info"),
                                                 value = 44),
                                        tabPanel("MA Plot",
                                                 p(),
                                                 verbatimTextOutput("MAGroupsText"),
                                                 jqui_resizable(plotOutput('MAPlot1', width = "800px", height = "550px",
                                                                           hover = hoverOpts("plot_hover2", delay = 10, delayType = "debounce"))),
                                                 uiOutput("hover_info2"),
                                                 value = 66),
                                        tabPanel("Box Plots",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotOutput('boxplot1', width = "100%", height = "600px")), type = 6),
                                                 value = 88),
                                        tabPanel("Gene Scatter Plot",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotlyOutput("geneScatter0", height = "500px")), type = 6),
                                                 DT::dataTableOutput("geneScatterTable"),
                                                 downloadButton("geneScatterDownload", "Download Non-log2 Transformed .tsv"),
                                                 value = 26),
                                        tabPanel("Pathway Analysis",
                                                 #verbatimTextOutput("upregpathtitle"),
                                                 #h3("Up-regulated pathway"),
                                                 uiOutput("UpRegPathLabel"),
                                                 verbatimTextOutput("upregpath_text"),
                                                 withSpinner(plotOutput('UpRegPathway1'), type = 6),
                                                 div(DT::dataTableOutput("UpRegPathwayTable1"), style = "font-size:10px"),
                                                 downloadButton("UpRegPathDownload", "Download Table .tsv"),
                                                 downloadButton("UpRegPathDownloadgmt", "Download .gmt"),
                                                 #verbatimTextOutput("dnregpathtitle"),
                                                 #h3("Down-regulated pathway"),
                                                 uiOutput("DnRegPathLabel"),
                                                 verbatimTextOutput("downregpath_text"),
                                                 withSpinner(plotOutput('DnRegPathway1'), type = 6),
                                                 div(DT::dataTableOutput("DnRegPathwayTable1"), style = "font-size:10px"),
                                                 downloadButton("DnRegPathDownload", "Download Table .tsv"),
                                                 downloadButton("DnRegPathDownloadgmt", "Download .gmt"),
                                                 p(),
                                                 value = 28),
                                        tabPanel("DEG Table",
                                                 p(),
                                                 verbatimTextOutput("degtext"),
                                                 div(DT::dataTableOutput("DEGtable1"), style = "font-size:12px"),
                                                 downloadButton("DEGtableDownload", "Download DEG table .tsv"),
                                                 value = 24)
                                    )
                                )
                            )
                        )
               ),
               
               ####----GSEA Tab----####
               
               tabPanel("GSEA Analysis",
                        fluidPage(
                            title = "GSEA Analysis",
                            sidebarLayout(
                                sidebarPanel(
                                    tabsetPanel(id = "GSEA",
                                                tabPanel("GSEA Parameters",
                                                         p(),
                                                         selectInput("comparisonA", "Comparison: GroupA",
                                                                     choices = metagroups, selected = metagroups[1]),
                                                         selectInput("comparisonB", "Comparison: GroupB",
                                                                     choices = metagroups, selected = metagroups[2]),
                                                         numericInput("userPval", "Pvalue Cutoff", value = 1.0, width = '100%'),
                                                         h4("ssGSEA Parameters"),
                                                         selectInput("ssGSEAtype","Choose ssGSEA Method",
                                                                     choices = c("ssgsea","gsva","zscore","plage")),
                                                         selectInput("boxplotcompare", "Boxplot Stat Compare Method:",
                                                                     choices = boxopt),
                                                         tabsetPanel(
                                                             id = "tables",
                                                             tabPanel("MSigDB Gene Sets",
                                                                      div(DT::dataTableOutput("msigdbTable"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                                      value = 1),
                                                             tabPanel(paste(userGSlist_name),
                                                                      div(DT::dataTableOutput("tab2table"), style = "font-size:10px; height:500px; overflow-X: scroll"),
                                                                      value = 3),
                                                             tabPanel("Use your own gene set",
                                                                      p(),
                                                                      uiOutput("user.gmt"),
                                                                      uiOutput("user.RDataButton"),
                                                                      uiOutput("RDataMessage"),
                                                                      uiOutput("user.GStable"),
                                                                      value = 5)
                                                         )
                                                ),
                                                tabPanel("Figure Parameters",
                                                         p(),
                                                         h4("Heatmap Parameters"),
                                                         sliderInput("heatmapFont1.c", "Heatmap Font Size (columns):",
                                                                     min = 5, max = 75,
                                                                     value = 12, step = 1),
                                                         sliderInput("heatmapFont1.r", "Heatmap Font Size (rows):",
                                                                     min = 5, max = 75,
                                                                     value = 10, step = 1),
                                                         h4("Box Plot Parameters"),
                                                         sliderInput("boxplotFontss", "Box Plot Font Size:",
                                                                     min = 5, max = 75,
                                                                     value = 15, step = 1),
                                                         sliderInput("boxplotDotss", "Box Plot Dot Size:",
                                                                     min = 0, max = 5,
                                                                     value = 0.75, step = .25)
                                                )
                                    )
                                ),
                                mainPanel(
                                    tabsetPanel(
                                        id = "datasettwo",
                                        tabPanel("Enrichment Plot",
                                                 h3("GSEA Enrichment Plot"),
                                                 verbatimTextOutput("NESandPval"),
                                                 withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
                                                 h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                                 downloadButton("LEGdownload", "Download .tsv"),
                                                 div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; overflow-y: scroll"),
                                                 value = 2),
                                        tabPanel("Gene Expression Heatmap",
                                                 withSpinner(jqui_resizable(plotOutput("heatmap0", width = "100%", height = "1000px")), type = 6),
                                                 value = 4),
                                        #tabPanel("GSEA Summary Table",
                                        #         p(),
                                        #         selectInput("SigTableChoice", "Select Enriched Signatures Table:",
                                        #                     choices = SigNames),
                                        #         h3("Enriched Signatures Table"),
                                        #         DT::dataTableOutput("enrich_sig_table", width = "100%"),
                                        #         downloadButton("enrich_sig_download","Download .tsv"),
                                        #         value = 6),
                                        #tabPanel("Generate Summary Table",
                                        #         p("Please note this may take several minutes depending on size and quantity of gene sets in GMT file."),
                                        #         withSpinner(DT::dataTableOutput("enrich_sig_table_gen"), type = 6),
                                        #         downloadButton("enrich_sig_download.u","Download .tsv"),
                                        #         value = 8),
                                        tabPanel("ssGSEA Boxplots",
                                                 p(),
                                                 withSpinner(jqui_resizable(plotOutput('boxplot2', width = "100%", height = "500px")), type = 6),
                                                 DT::dataTableOutput("ssGSEAtable"),
                                                 downloadButton("ssGSEAdownload", "Download .tsv"),
                                                 value = 10)
                                    )
                                )
                            )
                        )
               )
    )

