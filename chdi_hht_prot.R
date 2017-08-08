# This script needs to be excecuted from working directory
# The working directory needs to have the provided Ensembl information file
# in it, and the uncompressed directories that contain the cortex, liver,
# and striatum experimental conditions and results

# Turn off automatic factorization of character vectors. With all the sample
# and gene IDs, it is pointless, and it gets in the way of subsets based
# on vector.

stime <- proc.time()
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
options(stringsAsFactors = FALSE)

# Definitions ands parameters
genesrc <- 'All'   # Test, Cytoskeletal, Synaptic, or All
plotoption <- 'none'       # none, nopdf, makepdf, pdfonly
outliers <- c('F10Q92C_1','F2Q20C_4','F10Q7H_2','M10Q20H_2','F2Q92B_2')

rcut <- 0.5
pvmax <- 0.05
dfpcut <- 0.2
siglim <- 0
agetest <- c(2,6,10)
qltest <- c(7,175)

tissue <- c("Cerebellum","Cortex","Hippocampus","Striatum")
tcodes <- c("B","C","H","S")
lctissue <- tolower(tissue)

# Read sample info and experimental data files
if (! exists("alldata")) {
  sample_info <- read.table("Sample_name_info3.txt",header=TRUE)
  sample_info$Length <- as.integer(sub("^Q","",sample_info$Genotype))
  for (ttype in lctissue) {
    tabname <- paste0(ttype,"_data")
    fname <- paste0(ttype,"_data_samples_names2.txt")
    assign(tabname,read.table(fname,sep="\t", header=TRUE, fill=TRUE))
  }
  rm(ttype,tabname,fname)
  
  alldata <- left_join(cerebellum_data, cortex_data, by='GeneNames') %>%
    left_join(., hippocampus_data, by='GeneNames') %>%
    left_join(., striatum_data, by='GeneNames')

    rm(cerebellum_data,cortex_data,hippocampus_data,striatum_data)
  # Remove sample info for outliers that were removed
  # Likely will be a list of outliers to remove rather than missing data entries
  sample_info <- sample_info[sample_info$SampleName %in% colnames(alldata),]
  sample_info <- sample_info[! sample_info$SampleName %in% outliers,]
  
  SampleGene <- data.frame(GeneNames = alldata$GeneNames)
  SampleGene$Symbol <- gsub("_.*$","",SampleGene$GeneNames)
  SampleGene <- separate_rows(SampleGene, Symbol, sep = ";")
}

if ((! exists("lastsrc")) || (lastsrc != genesrc)) {
  if (exists("sylist")) {
    rm(sylist)
  }
  if (exists("listdata")) {
    rm(listdata)
  }
  if (exists("expinf")) {
    rm(expinf)
  }
  if (exists("qlresults")) {
    rm(qlresults)
  }
  if (exists("ageresults")) {
    rm(ageresults)
  }
  
  if (genesrc == 'Cytoskeletal') {
    sylist <- read.table("Cytoskeletal_mouse_symbol.txt",skip=1)[,1]
  } else if (genesrc == 'Synaptic') {
    sylist <- read.table("Synaptic_mouse_symbol.txt",skip=1)[,1]
  } else if (genesrc == 'Test') {
    sylist <- sample(SampleGene$Symbol,3)
  }
}

if (exists("sylist")) {
  listdata <- alldata[alldata$GeneNames %in%
    SampleGene[SampleGene$Symbol %in% sylist,1],]
} else {
  listdata <- alldata
}

# Linear regressions correlations for 2, 6, and 10 months
if (! exists("expinf")) {
  for (tct in 1:length(tissue)) {
    ttype <- tissue[tct]
    tcode <- tcodes[tct]
    for (age in agetest) {
      qlength <- subset(sample_info,Tissue == ttype & Age == age,
                        select=c(SampleName,Length))
      expr <- listdata[,c("GeneNames",qlength$SampleName)]
      expr <- expr[rowSums(!is.na(expr)) > 3,]
  
      alab <- age
      if (age < 10) alab <- paste0(0,alab)
      expcode <- paste0(tcode,"RL","A",alab)
      
      exptab <- data.frame("Exp"=expcode,"Comparison"="Regression Q-length",
        "Tissue"=ttype,"Constant"=paste0(alab," mo"),"Value"="Q-length")
      if (exists("expinf")) {
        expinf <- bind_rows(expinf,exptab)
      } else {
        expinf <- exptab
      }
      
      restab <- data.frame("GeneNames"=expr$GeneNames,
        "Symbol"=gsub("_.*$","",expr$GeneNames),"Exp"=expcode,
        "Tissue"=rep(ttype,nrow(expr)),"Age"=rep(age,nrow(expr)))
      expr <- data.matrix(expr[-1])
      restab$PearsonR <- apply(expr,1,function(y) cor(qlength$Length,y,
                               use='na.or.complete'))
      restab[,c("Intercept","Slope","Pvalue")] <-
        t(apply(expr,1,function(y) coef(summary(lm(y~qlength$Length,
                                        na.action=na.omit)))[c(1,2,8)]))
      if (exists("qlresults")) {
        qlresults <- bind_rows(qlresults,restab)
      } else {
        qlresults <- restab
      }
    }
  }
  qlresults$Pvalue[is.na(qlresults$Pvalue)] <- 1
  qlresults$PearsonR[is.na(qlresults$PearsonR)] <- 0
  qlresults$absPR <- abs(qlresults$PearsonR)
  qlresults$absRise <- abs(qlresults$Slope*168)
  
  # Linear regressions and correlations for q-lengths of 175, 111, 80, and 7 (WT)
  for (tct in 1:length(tissue)) {
    ttype <- tissue[tct]
    tcode <- tcodes[tct]
    for (length in qltest) {
      ages <- subset(sample_info,Tissue == ttype & Length == length,
                     select=c(SampleName,Age))
      expr <- listdata[,c("GeneNames",ages$SampleName)]
      expr <- expr[rowSums(!is.na(expr)) > 3,]
      
      qlab <- "Q"
      if (length < 100) qlab <- paste0(qlab,0)
      if (length < 10) qlab <- paste0(qlab,0)
      qlab <- paste0(qlab,length)
      expcode <- paste0(tcode,"RA",qlab)
      
      exptab <- data.frame("Exp"=expcode,"Comparison"="Regression Age",
        "Tissue"=ttype,"Constant"=qlab,"Value"="Age")
      expinf <- bind_rows(expinf,exptab)
      
      restab <- data.frame("GeneNames"=expr$GeneNames,
        "Symbol"=gsub("_.*$","",expr$GeneNames),"Exp"=expcode,
        "Tissue"=rep(ttype,nrow(expr)),"Length"=rep(length,nrow(expr)))
      expr <- data.matrix(expr[-1])
      
      restab$PearsonR <- apply(expr,1,function(y) cor(ages$Age,y,
                               use='na.or.complete'))
      restab[,c("Intercept","Slope","Pvalue")] <-
        t(apply(expr,1,function(y) coef(summary(lm(y~ages$Age,
                                        na.action=na.omit)))[c(1,2,8)]))
      if (exists("ageresults")) {
        ageresults <- bind_rows(ageresults,restab)
      } else {
        ageresults <- restab
      }
    }
  }
  ageresults$Pvalue[is.na(ageresults$Pvalue)] <- 1
  ageresults$PearsonR[is.na(ageresults$PearsonR)] <- 0
  ageresults$absPR <- abs(ageresults$PearsonR)
  ageresults$absRise <- abs(ageresults$Slope*8)
}

# Assign significant correlations

qlresults$SigGene <- 0
qlresults[qlresults$Pvalue < pvmax & qlresults$absPR > rcut &
            qlresults$absRise > dfpcut,"SigGene"] <- 1
qlsiglm <- qlresults %>% filter(SigGene == 1)
write.table(qlresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/prot_",genesrc,"_qlength_lm.txt"))

ageresults$SigGene <- 0
ageresults[ageresults$Pvalue < pvmax & ageresults$absPR > rcut &
            ageresults$absRise > dfpcut,"SigGene"] <- 1
agesiglm <- ageresults %>% filter(SigGene == 1)
write.table(ageresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/prot_",genesrc,"_age_lm.txt"))

ftime <- proc.time()
print("Data input plus inear regression time:")
print(ftime - stime)

# Form summary table of significant expression changes
sigtable <- bind_rows(qlsiglm %>% select(GeneNames,Exp) %>% mutate(Match = 1),
                 agesiglm %>% select(GeneNames,Exp) %>% mutate(Match = 1))

# Need to add in zero match line for all tests to account for tests that have
# zero significant matches
tmplg <- agesiglm[1,"GeneNames"]
sigtable <- bind_rows(sigtable,data.frame(GeneNames=tmplg,Exp=expinf$Exp,Match=0))
sigtable <- sigtable %>% group_by(GeneNames,Exp) %>% summarize(Match = sum(Match))

dsummary <- dcast(sigtable,GeneNames ~ Exp,value.var="Match")
dsummary[is.na(dsummary)] <- 0
dsummary$Total <- dsummary %>% select(-GeneNames) %>% rowSums()
dsummary$Symbol <- gsub("_.*$","",dsummary$GeneNames)
dsummary <- dsummary[c(1,ncol(dsummary),ncol(dsummary)-1,2:(ncol(dsummary)-2))]
dsummary <- dsummary %>% arrange(desc(Total),GeneNames)

write.table(dsummary,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/prot_",genesrc,"_summary_sig.txt"))

expinf <- expinf %>% arrange(Exp)
write.table(t(expinf),sep="\t",col.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/prot_",genesrc,"_exp_info.txt"))
hpgenes <- dsummary %>% filter(Total > siglim) %>%
  select(GeneNames) %>% unlist(use.names = F)

# Plot Intensity vs. Q Length and Age for select genes
plotprot <- function(sterm,goption="") {
  pg <- list()
  plc <- 1
  for (ttype in tissue) {
    qldata <- subset(sample_info,Tissue == ttype & Age %in% agetest,
                     select=c(SampleName,Age,Length))
    qldata$Expr <- t(listdata[listdata$GeneNames==sterm,qldata$SampleName])
    qldata <- na.omit(qldata)
    qltitle <- paste0(gsub("_"," ",sterm)," ", ttype)
    qldata$Age <- as.factor(qldata$Age)
    pg[[plc]] <- ggplot(qldata, aes(x=Length, y=Expr, color=Age, shape=Age)) +
      geom_point() + geom_smooth(method=lm,se=FALSE) + ggtitle(qltitle) +
      labs(x="Q Length",y="Intensity") + xlim(0,200) +
      theme(plot.title = element_text(hjust = 0.5))
    if (NROW(qldata) == 0) {
      pg[[plc]] <- pg[[plc]] + annotate("text",x=100,y=0.5,label = "No data")
    }
    plc <- plc + 1
  }

  plgt=textGrob(paste0(sterm, " Length Effect"),
                gp=gpar(fontsize=18,fontface="bold"))
  mgplot <- arrangeGrob(pg[[1]],pg[[2]],pg[[3]],
        pg[[4]],ncol=2,top=plgt)
  if (goption != 'pdfonly') {
    grid.newpage()
    grid.draw(mgplot)
  }
  if (goption %in% c('pdfonly','makepdf')) {
    plname <- paste0(genesrc, "/", sterm, "_qlprot.pdf")
    ggsave(plname,mgplot,width=8,height=10,units="in")
  }

  pg <- list()
  plc <- 1
  for (ttype in tissue) {
    agedata <- subset(sample_info,Tissue == ttype & Length %in% qltest,
                      select=c(SampleName,Age,Length))
    agedata$Expr <- t(listdata[listdata$GeneNames==sterm,agedata$SampleName])
    agedata <- na.omit(agedata)
    agetitle <- paste0(gsub("_"," ",sterm)," ", ttype)
    agedata$Length <- as.factor(agedata$Length)
    pg[[plc]] <- ggplot(agedata, aes(x=Age, y=Expr, color=Length, shape=Length)) +
      geom_point() + geom_smooth(method=lm,se=FALSE) + ggtitle(agetitle) +
      labs(x="Age",y="Intensity") + xlim(2,10) +
      theme(plot.title = element_text(hjust = 0.5))
    if (NROW(agedata) == 0) {
      pg[[plc]] <- pg[[plc]] + annotate("text",x=6,y=0.5,label = "No data")
    }
    plc <- plc + 1
  }
  
  plgt=textGrob(paste0(sterm, " Age Effect"),
                gp=gpar(fontsize=18,fontface="bold"))
  mgplot <- arrangeGrob(pg[[1]],pg[[2]],pg[[3]],
                        pg[[4]],ncol=2,top=plgt)
  if (goption != 'pdfonly') {
    grid.newpage()
    grid.draw(mgplot)
  }
  if (goption %in% c('pdfonly','makepdf')) {
    plname <- paste0(genesrc,"/", sterm, "_ageprot.pdf")
    ggsave(plname,mgplot,width=8,height=10,units="in")
  }
}

if (plotoption != "none") {
  for (protid in hpgenes) {
    plotprot(protid,plotoption)
  }
}

lastsrc <- genesrc

etime <- proc.time()
print("Plotting time:")
print(etime - ftime)

print("Total time:")
print(etime - stime)
