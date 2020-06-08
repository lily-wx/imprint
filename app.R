#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinythemes)
library(DT)
library(ggpubr)
library(data.table)
# setRepositories(addURLs =c( BioCsoft ="https://bioconductor.org/packages/3.8/bioc",
#                             BioCann="https://bioconductor.org/packages/3.8/data/annotation",
#                             BioCexp="https://bioconductor.org/packages/3.8/data/experiment",
#                             BioCworkflows="https://bioconductor.org/packages/3.8/workflows"))
#options("repos"=c( BioCsoft ="https://bioconductor.org/packages/3.8/bioc",
#                   BioCann="https://bioconductor.org/packages/3.8/data/annotation",
#                   BioCexp="https://bioconductor.org/packages/3.8/data/experiment",
#                   BioCworkflows="https://bioconductor.org/packages/3.8/workflows"))
library(BiocManager)
library(regioneR)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(karyoploteR)
library(shiny)
library(shinydashboard)
options(repos = BiocManager::repositories())


#读取数据，并处理
imprint<-fread("hg38.imprint_gene_all.sort.bed")
hg19.imprint<-fread("hg19.imprint_gene_all.sort.bed")
imprint_genes<-toGRanges(imprint,format="BED")
genelist<-fread("hg38.imprintgenes.signif_variant_gene_pairs.txt",header=F)
state.file <- fread("state_info.txt",header = T)
#k562hmm.file<-read.table("hg38.wgEncodeBroadHmmK562HMM.bed",header=F)
#K562.hmm <- toGRanges(k562hmm.file)
hmm.info <- fread("Roadmap.metadata.qc.jul2013 - Consolidated_EpigenomeIDs_summary_Table.tsv",header=T)
hmm.info.dat <- hmm.info[-c(1:2),c(2,15,18)]
hmm.info.tissue <- hmm.info.dat[TYPE=="PrimaryTissue"]
clinvar <- fread("clinvar.imprint_promoter.gDMR.simp.txt",header=F)
dbvar <- fread("clinical_variants_for_nstd102.promoter.txt",header=F)
gDMR <- fread("hg19.gDMR.bed")
ICGC <- fread("ICGC.promoter.gDMR.tumor.uniq.txt",header=F,quote="")
cluster <- fread("cluster_h3k27me3_result.csv",header=T)
great.m <- fread("great_maternal.tsv",header = T,sep = "\t")
great.p <- fread("great_paternal.tsv",header = T,sep="\t")
#GREAT
input.great<-rbindlist(list(great.m[,c(3,5)],great.p[,c(1,3)]),use.names = F)
input.great[,logP:=-log10(BinomP)]
input.great[, Class := c(rep("M", nrow(great.m)),rep("P", nrow(great.p)))]
#gDMR
colnames(gDMR) <- c("chr","start","end","name")
#clinvar
clinvar.info <- clinvar[,c(1,2,3,5,12,4)]
colnames(clinvar.info) <- c("chr","start","end","gene","gDMR","description")
#dbvar
dbvar.info <- dbvar[,c(1,2,3,4,6,10)]
colnames(dbvar.info) <- c("chr","start","end","mut","category","gene")
#处理snp
snps <- na.omit(unique(data.table(genelist$V3, genelist$V4, genelist$V11)))
tissue <- c("all",as.character(unique(snps$V1)))
snp <- fread("hg38.imprintgenes.signif_variant_gene_pairs.txt",header=F)
colnames(snp) <- c("gencode","gene_symbol","tissue","loc","name","tss_distance","ma_samples","ma_count","maf",
                   "pval_nominal","slope","slope_se","pval_nominal_threshold","min_pval_nominal","pval_beta")
snp.dat <- snp[,c(5,2,4,3,10,11,15)]
imprint_gene.dat <- as.data.table(matrix(unlist(strsplit(imprint_genes$V4,split = "\\|")),ncol = 4,byrow = T)[,1])
imprint_gene.dat[,class:="imprint"]
snp.dat[,class:=imprint_gene.dat[match(snp.dat$gene_symbol,V1)]$class]
#处理imprint
hg19.imprint[,class:=matrix(unlist(strsplit(as.character(hg19.imprint$V4),split = "\\|")),ncol = 4,byrow = T)[,2]]
hg19.imprint[,fun:=matrix(unlist(strsplit(as.character(hg19.imprint$V4),split = "\\|")),ncol = 4,byrow = T)[,4]]
maternal <- hg19.imprint[grep("M",hg19.imprint$class)]
paternal <- hg19.imprint[grep("P",hg19.imprint$class)]
coding <- hg19.imprint[fun=="protein_coding"]
noncoding <- hg19.imprint[fun!="protein_coding"]
#write.table(maternal[-1,1:6],"~/Desktop/毕业设计/数据/imprint_gene_beddata/hg19.maternal.bed",quote = F,row.names = F,col.names = F)
#write.table(paternal[-1,1:6],"~/Desktop/毕业设计/数据/imprint_gene_beddata/hg19.paternal.bed",quote = F,row.names = F,col.names = F)
#ICGC
icgc.mut <- ICGC[,c(14,15,16,17,18,19,23,2)]
colnames(icgc.mut) <- c("chr","start","end","Mu_Id","ref","alt","gene","sample")
#tissue spec
phaser_all <- fread("phASER.imprintgenes.combine_2.txt",header=F)
phaser_all[,V4:=paste(V2,"_",V1,"_",V3,"_",V4,sep="")]
phaser <- fread("phaser.process.combine.txt", sep=" ", header = T)
tmp <- phaser_all[,-c(1:3)]
colnames(tmp) <- colnames(phaser) 
phaser_change.dt <- as.data.table(melt(tmp, id.vars = "#contig_name_start_stop"))
phaser_change.dt <- phaser_change.dt[!is.na(value)]
sample_tissue_age <- fread("sample_tissue_age.process.txt", sep = "\t", header = F)
phaser_change.dt[,Tissue:=sample_tissue_age[match(phaser_change.dt$variable,V1)]$V2]
phaser_change.dt[,Age:=sample_tissue_age[match(phaser_change.dt$variable,V1)]$V3]
colnames(phaser_change.dt)[1] <- "contig_name_start_stop"
phaser_change.dt[, meanvalue := mean(value), by=c("contig_name_start_stop","Tissue")]
phaser_change.dt[, medianvalue := median(value), by=c("contig_name_start_stop","Tissue")]
phaser_change.dt[, max_medianvalue := max(medianvalue), by="contig_name_start_stop"]
phaser_change.dt.median <- unique(phaser_change.dt[,medianvalue,by=c("contig_name_start_stop","Tissue")])
phaser_change.dt.median[, max_medianvalue := max(medianvalue), by="contig_name_start_stop"]                                               
phaser_change.dt.median[,Tissue_count:=.(.N),by="contig_name_start_stop"]
phaser_change.dt.median[,Tissue_spec:=sum(1-medianvalue/max_medianvalue)/(Tissue_count-1),by="contig_name_start_stop"]
phaser_change.dt.median[Tissue_spec>=0.76,Flag:="TissueSpec"]
phaser_change.dt.median[Tissue_spec<0.19,Flag:="nonTissueSpec"]
phaser_change.dt.median[,Genename:=matrix(unlist(strsplit(phaser_change.dt.median$contig_name_start_stop, split = "\\_")), ncol = 4, byrow = T)[,2]]
phaser_change.dt.median[Flag=="nonTissueSpec",Median_medianvalue:=median(medianvalue),by="Genename"]
phaser_change.dt.median[Flag=="TissueSpec",Median_medianvalue:=median(medianvalue),by="Genename"]
phaser_change.dt.median[Median_medianvalue<0.25&Flag=="nonTissueSpec",Class:="nonImprint"]
phaser_change.dt.median[Median_medianvalue>=0.75&Flag=="nonTissueSpec",Class:="Imprint"]
phaser_change.dt.median[Median_medianvalue>=0.25&Median_medianvalue<0.75&Flag=="nonTissueSpec",Class:="other"]
phaser.dt <- phaser_change.dt.median[,c(1,8,2,3,6,7,10)]

# Define UI for application that draws a histogram
header=dashboardHeader(title = "RESULT")
sidebar=dashboardSidebar(
    sidebarMenu(
        menuItem("DISTRIBUTION", tabName = "distribution", icon = icon("dashboard")),
        menuItem("MECHANISM", tabName = "mechanism", icon = icon("angellist")),
        menuItem("TISSUE-SPEC", tabName = "spec", icon = icon("android")),
        menuItem("EQTL", icon = icon("th"), tabName = "eqtl",
                 badgeLabel = "new", badgeColor = "green"),
        menuItem("DISEASE", tabName = "disease", icon = icon("cog", lib = "glyphicon"))
    ),
    textInput("gene_spec","Gene Tissue Spec",value = "PEG10"),
    textInput("loc","GRCh38",value = "chr11:1700000-2151603"),
    sliderInput("zoom","zoom", min =0,max= 2,value = 1,width = "100%",step=0.25),
    tags$div(
        tags$a(href="https://genome.ucsc.edu/cgi-bin/hgLiftOver","LiftOver!"),
        style="color:#FFFFCC;text-align:center;"
             ),
    selectInput("hmm_tissue","hmm.tissue",choices=hmm.info.tissue[,2],selected = "Skeletal Muscle Male"),
    selectInput("tissue","snp.tissue",choices=tissue,selected = "all"),
    actionButton("button","GO!"),
    helpText("GRCh38 is only the genome coordinate of EQTL part!")
)
body=dashboardBody(
    tabItems(
        tabItem(tabName = "distribution",
                h2("DISTRIBUTION"),
                tabBox(width = 12,
                       tabPanel("IMPRINT",DT::dataTableOutput("imprint"),icon = icon("tree")),
                       tabPanel("gDMR",DT::dataTableOutput("gdmr")),
                       tabPanel("M/P", plotOutput("MP")),
                       tabPanel("PROTEIN",plotOutput("CD")),
                       tabPanel("FUNCTION",plotOutput("func"))
                )
        ),
         tabItem(tabName = "spec",
                 h2("TISSUE-SPEC"),
                 tabBox(width = 12,
                        tabPanel("Tissue Spec Plot",
                                 plotOutput("spec_plot")),
                        tabPanel("Tissue Spec Data",
                                 DT::dataTableOutput("spec")),
                        tabPanel("Methy",icon = icon("paw"),
                                 tags$img(src="methy.png",width="1000px",height="1000px")),
                        tabPanel("Imprint-Index", 
                                 fluidRow(
                                     column(width=12,
                                     box(title = "SPEC", width=12,collapsible=TRUE,collapsed=TRUE,
                                         tags$img(src="spec.png",width="800px",height="900px")),
                                     box(title= "Not-SPEC", width=12,collapsible=TRUE,collapsed = TRUE,
                                         tags$img(src="nonspec.ordr.png",width="600px",height="600px"))
                                     )
                                 )
                                 )
                 )
         ),
        tabItem(tabName = "mechanism",
                h2("MECHANISM"),
                tabBox(width = 12,
                       tabPanel("Cluster_Info", icon = icon("bell"),
                                DT::dataTableOutput("cluster_info")),
                       tabPanel("Cluster_Heatmap", 
                                fluidRow(
                                    column(width = 12,
                                    box(title = "cluster 1344", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="1344_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 1382", width=12,collapsible=TRUE,collapsed = TRUE,
                                        tags$img(src="1382_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 1489", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="1489_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 1956", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="1956_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 2421", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="2421_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 2431", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="2431_promoter_h3k27me3.png",width="800px",height="900px")),
                                    box(title = "cluster 2857", width=12,collapsible=TRUE,collapsed=TRUE,
                                        tags$img(src="2857_promoter_h3k27me3.png",width="800px",height="900px"))
                                    )
                                )
                       ),
                       tabPanel("LncRNA",tags$img(src="lncRNA.lr100k.genes.great.png",width="400px",height="400px"))
                )
        ),
        
        tabItem(tabName = "eqtl",
                h2("Gtex_Eqtl"),
                tabBox(
                    width = 12,
                    tabPanel("Eqtl_Plot",icon=icon("birthday-cake"),
                             fluidRow(
                                 column(width = 12,
                                 box(title = "PLOT", width=12,collapsible=TRUE,
                                     plotOutput("Plot")),
                                 box(title = "STATE", width=12,collapsible=TRUE,collapsed=TRUE,
                                     DT::dataTableOutput("state"))
                                 )
                             )
                    ),
                    tabPanel("SNP",DT::dataTableOutput("mytable")),
                    tabPanel("Sankey",list(
                                            #imageOutput("sankey"),
                                            tags$img(src="muscle.png",width="1000px",height="1000px"),
                                            tags$a(em("open link"),href="muscle_snp.html"))
                                            )
                )
        ),
        tabItem(tabName = "disease",
                h2("DISEASE"),
                tabBox(width = 12,
                       tabPanel("CLINVAR",icon=icon("cookie-bite"),
                                DT::dataTableOutput("clinvar")),
                       tabPanel("dbVAR",DT::dataTableOutput("dbvar")),
                       tabPanel("ICGC",DT::dataTableOutput("icgc")),
                       tabPanel("CTCF",
                                fluidRow(
                                box(title = "Distri.Gain", width=6,collapsible=TRUE,
                                    tags$img(src="count.gain.png",width="300px",height="300px")),
                                box(title= "Distri.Loss", width=6,collapsible=TRUE,
                                    tags$img(src="count.loss.png",width="300px",height="300px")
                                    )
                                ),
                                fluidRow(
                                    box(title = "Motif.Gain.ref", width=6,collapsible=TRUE,
                                        tags$img(src="ctcf.gain.png",width="300px",height="300px")),
                                    box(title= "Motif.Loss.ref", width=6,collapsible=TRUE,
                                        tags$img(src="ctcf.loss.png",width="300px",height="300px"))
                                ),
                                fluidRow(
                                    box(title = "Motif.Gain.mut", width=6,collapsible=TRUE,
                                        tags$img(src="ctcf.gain.alt.png",width="300px",height="300px")),
                                    box(title= "Motif.Loss.mut", width=6,collapsible=TRUE,
                                        tags$img(src="ctcf.loss.alt.png",width="300px",height="300px"))
                                )
                                )
                )
        )
    )
)
ui <- dashboardPage(header,sidebar,body)
# ui <- fluidPage(theme = shinytheme("flatly"),
# 
#     # Application title
#     titlePanel("Gtex_Eqtl"),
# 
#     # Sidebar with a slider input for number of bins 
#  
#         sidebarPanel(
#             textInput("loc","GRCh38",value = "chr11:1700000-2151603"),
#             selectInput("tissue","tissue",choices=tissue,selected = "all"),
#             actionButton("button","GO!"),
#             helpText("Example:chr11:1700000-2151603")
#         ),
# 
#         # Show a plot of the generated distribution
#         mainPanel(
#            tabsetPanel(
#                tabPanel("Eqtl_Plot",plotOutput("Plot")),
#                #tabPanel("SNP",tableOutput("SNP")),
#                tabPanel("SNP",DT::dataTableOutput("mytable")),
#                tabPanel("Sankey",list(
#                    imageOutput("sankey"),
#                    tags$a(em("open link"),href="muscle_snp.html"))
#                    )
#            )
#         )
#     
# )

server <- function(input, output,session) {

    #EQTL Plot
    output$Plot <- renderPlot({
        # generate bins based on input$bins from ui.R
        input$button
        #制作坐标
        loc <- isolate(input$loc)
        coor=unlist(strsplit(loc,split = ":"))[2]
        start=as.numeric(unlist(strsplit(coor,split = "-"))[1])
        end=as.numeric(unlist(strsplit(coor,split = "-"))[2])
        zoom <- isolate(input$zoom)
        chrom <- matrix(unlist(strsplit(loc,split = ":")),ncol=2,byrow = T)[,1]
        start.z <-start-floor((zoom-1)/2*(end-start))
        end.z <- end+floor((zoom-1)/2*(end-start))
        loc.z<- paste(chrom,":",start.z,"-",end.z,sep="")
        
        tss <- isolate(input$tissue)
        hmm.tss <- isolate(input$hmm_tissue)
        #确定展示的snp来自的组织
        if(tss=="all")
        {
            snp.se=snps
        }
        else{
            snp.se=snps[V1 %in% tss]
        }
        snp.pos <- as.data.table(matrix(unlist(strsplit(as.character(snp.se$V2), split = "_")), nrow = nrow(snp.se), byrow = T))
        colnames(snp.pos) <- c("chr", "end", "ref", "alt", "genome")
        snp.pos.df <- data.frame(chr = snp.pos$chr, 
                                 start = as.numeric(snp.pos$end) - 1, 
                                 end = as.numeric(snp.pos$end))
        snp.gr <- makeGRangesFromDataFrame(snp.pos.df)
        snp.gr$y <- snp.se$V3
        
        #确定展示的chromhmm来自的组织
        EID <- hmm.info.dat[`Standardized Epigenome name`==hmm.tss,1]
        hmm.file <- paste(EID,"_15_coreMarks_hg38lift_dense.bed",sep="")
        hmm.show.file <-fread(hmm.file,header = F,quote = F)
        hmm.show.file[,RGB:=state.file[match(hmm.show.file$V9,`COLOR CODE`)]$`RGB CODE`]
        hmm.show <- toGRanges(hmm.show.file)
        # draw the histogram with the specified number of bins
        at <- autotrack(current.track = 1:2, total.tracks = 5)
        kp<-plotKaryotype(genome="hg38",chromosomes=chrom,zoom=loc.z,main=loc.z,plot.type=4)
        kpAddLabels(kp, labels = "eqtl", srt=90, pos=3, cex=1.8, label.margin = 0.025)
        at <- autotrack(current.track = 1:2, total.tracks = 5)
        #kpPlotGenes(kp, data=imprint_genes,add.transcript.names = FALSE, r0=at$r0,r1=at$r1, cex=0.8,gene.name.position = "left")
        kpPlotMarkers(kp, data=imprint_genes,labels=gsub("\\|\\w+","",imprint_genes$V4),r0=at$r0,r1=at$r1, cex=0.8,
                      text.orientation="horizontal")
        at <- autotrack(current.track = 3, total.tracks = 5)
        #kpPlotRegions(kp, K562.hmm, col=K562.hmm$V9, r0=at$r0, r1=at$r1)
        kpPlotRegions(kp, hmm.show, col=hmm.show$RGB, r0=at$r0, r1=at$r1)
        at <- autotrack(current.track =4:5 , total.tracks = 5)
        #kpAxis(kp, ymin=min(snp.gr$y), ymax=max(snp.gr$y),r0=at$r0,r1=at$r1,side=2)
        kpAxis(kp, ymin=-3, ymax=3,r0=at$r0,r1=at$r1,side=2)
        kpPoints(kp, data=snp.gr,y=(snp.gr$y*0.4/6+at$r0),col = "blue",cex=0.5,r0=at$r0,r1=at$r1)
    })
    
    #Gtex_snp
    output$mytable = DT::renderDataTable({
        snp.dat
    })
    
    #Sankey Plot
    # output$sankey <- renderImage(
    #     {   
    #         list(src = "/Users/apple/Desktop/毕业设计/结果/eqtl/try/muscle.png",
    #          contentType = 'image/png',
    #          width = 600,
    #          height = 800,
    #          alt = "This is alternate text")}
    #     , deleteFile = FALSE)
    
    #MP distribution
    output$MP <- renderPlot(
        {
            paternal_genes<-toGRanges(paternal,format="BED")
            maternal_genes<-toGRanges(maternal,format="BED")
            kp <- plotKaryotype()
            kpPlotRegions(kp,data=maternal_genes,col="red")
            kpPlotRegions(kp,data=paternal_genes,col="blue")
        }
    )
    
    #Protein coding distribution
    output$CD <- renderPlot(
        {
            coding_genes<-toGRanges(coding,format="BED")
            noncoding_genes<-toGRanges(noncoding,format="BED")
            kp <- plotKaryotype()
            kpPlotRegions(kp,data=coding_genes,col="green")
            kpPlotRegions(kp,data=noncoding_genes,col="red")
        }
    )
    
    #Imprint Data Table
    output$imprint = DT::renderDataTable({
        hg19.imprint
    })
    
    #clinvar DATA Table
    output$clinvar = DT::renderDataTable({
        clinvar.info
    })
    #dbvar DATA Table
    output$dbvar = DT::renderDataTable({
        dbvar.info
    })
    #icgc DATA Table
    output$icgc = DT::renderDataTable({
        icgc.mut
    })
    #gDMR DATA Table
    output$gdmr = DT::renderDataTable({
        gDMR
    })
    #cluster DATA Table
    output$cluster_info = DT::renderDataTable({
        cluster
    })
    #chrom state DATA Table
    output$state = DT::renderDataTable({
        state.file
    })
    #function plot
    output$func = renderPlot(
        {
            p1 <- ggbarplot(data = input.great[logP>5], x = "Desc", y = "logP",
                      color = "Class", fill = "Class", palette = "npg",
                      xlab = "", ylab = "-10logP",sort.by.groups = F, sort.val = "desc")
            ggpar(p1, rotate = T)
        }
    )
    #gene spec
    output$spec = DT::renderDataTable({
        phaser.dt
    })
    output$spec_plot=renderPlot(
        {
            gene=input$gene_spec
            p1 <- ggbarplot(phaser.dt[Genename==gene],x="Tissue",
                            y="medianvalue",fill  = "Tissue",
                            palette = get_palette("npg",30),sort.val = "desc",
                            sort.by.groups = F)
            ggpar(p1,x.text.angle = 60,title = gene)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
