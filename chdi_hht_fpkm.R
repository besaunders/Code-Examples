# This script needs to be excecuted from working directory
# The working directory needs to have the provided Ensembl information file
# in it, and the uncompressed directories that contain the cortex, liver,
# and striatum experimental conditions and results

# Turn off automatic factorization of character vectors. With all the sample
# and gene IDs, it is pointless, and it gets in the way of subsets based
# on vector.

stime <- proc.time()
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
options(stringsAsFactors = FALSE)

# Genesrc = Cytoskeletal, Synaptic, Test, or All
genesrc <- 'All'

# Plotoption = none, nopdf, makepdf, pdfonly
plotoption <- 'none'

# Set various parameters
tissue <- c("Cerebellum","Cortex","Hippocampus","Striatum")
tcodes <- c("B","C","H","S")
lctissue <- tolower(tissue)
agetest <- c(2,6,10)
qltest <- c(7,175)
rcut <- 0.7
pvmax <- 0.05
dfpcut <- 2

# Read Ensembl name file. This assumes a unique mapping, so the file should
# not have any repeats of EnsemblID
ensembl <- read.table("ensembl_mm_symbol.txt",sep="\t",quote="",header=TRUE)

# Read sample info files
sample_info <- do.call(rbind,lapply(paste0(lctissue,"_sampleinfo.txt"),read.table))
sample_info$Length <- as.integer(sub("WT","7",sub("^Q","",sample_info$Genotype)))

for (ttype in lctissue) {
  tabname <- paste0(ttype,"_fpkm")
  fname <- paste0(ttype,"_procdata.txt")
  assign(tabname,read.table(fname,sep="\t", header=TRUE))
}
rm(ttype,tabname,fname)

if (genesrc == 'Cytoskeletal') {
  sylist <- read.table("Cytoskeletal_mouse_symbol.txt",header=TRUE)
} else if (genesrc == 'Synaptic') {
  sylist <- read.table("Synaptic_mouse_symbol.txt",header=TRUE)
} else if (genesrc == 'Test') {
  # enslist <- sample(cerebellum_fpkm$EnsemblID,5)
  enslist <- c('ENSMUSG00000020115')
}
if (exists("sylist")) {
  colnames(sylist) <- c("GeneName")
  enslist <- inner_join(sylist,ensembl,by="GeneName") %>%
    select(EnsemblID) %>% unlist(use.names = F)
}

if (exists("enslist")) {
  cerebellum_fpkm <- cerebellum_fpkm[cerebellum_fpkm$EnsemblID %in% enslist,]
  cortex_fpkm <- cortex_fpkm[cortex_fpkm$EnsemblID %in% enslist,]
  hippocampus_fpkm <- hippocampus_fpkm[hippocampus_fpkm$EnsemblID %in% enslist,]
  striatum_fpkm <- striatum_fpkm[striatum_fpkm$EnsemblID %in% enslist,]
}

fpkm <- left_join(cerebellum_fpkm, cortex_fpkm, by='EnsemblID') %>%
  left_join(., hippocampus_fpkm, by='EnsemblID') %>%
  left_join(., striatum_fpkm, by='EnsemblID')
fpkm[is.na(fpkm)] <- 0

# Remove sample info for outliers that were removed
sample_info <- sample_info[sample_info$SampleName %in% colnames(fpkm),]

# Linear regressions correlations for 2, 6, and 10 months
if (exists("expinf")) {
  rm(expinf)
}
if (exists("qlresults")) {
  rm(qlresults)
}
for (tct in 1:length(tissue)) {
  ttype <- tissue[tct]
  tcode <- tcodes[tct]
  for (age in agetest) {
    qlength <- subset(sample_info,Tissue == ttype & Age == age,
                      select=c(SampleName,Length))
    expr <- fpkm[,qlength$SampleName]
    
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
    
    restab <- data.frame("EnsemblID"=fpkm$EnsemblID,"Exp"=expcode,
      "Tissue"=rep(ttype,nrow(fpkm)),"Age"=rep(age,nrow(fpkm)))
    restab$PearsonR <- apply(expr,1,function(y) cor(qlength$Length,y))
    restab[,c("Slope","Pvalue")] <-
      t(apply(expr,1,function(y) coef(summary(lm(y~qlength$Length)))[c(2,8)]))
    if (exists("qlresults")) {
      qlresults <- bind_rows(qlresults,restab)
    } else {
      qlresults <- restab
    }
  }
}
rm(age,qlength,alab)
qlresults$Pvalue[is.na(qlresults$Pvalue)] <- 1
qlresults$PearsonR[is.na(qlresults$PearsonR)] <- 0
qlresults <- left_join(qlresults,ensembl,by="EnsemblID")
qlresults <- qlresults[c("EnsemblID","GeneName","Exp","Tissue","Age",
                       "PearsonR","Slope","Pvalue")]
qlresults$absPR <- abs(qlresults$PearsonR)
qlresults$absRise <- abs(qlresults$Slope*168)
qlresults$SigGene <- 0
qlresults[qlresults$Pvalue < pvmax & qlresults$absPR > rcut &
            qlresults$absRise > dfpcut,"SigGene"] <- 1
qlsiglm <- qlresults %>% filter(SigGene == 1)
write.table(qlresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/fpkm_",genesrc,"_expr_qlength_lm.txt"))

# Linear regressions and correlations for q-lengths of 175, 111, 80, and 7 (WT)
if (exists("ageresults")) {
  rm(ageresults)
}
for (tct in 1:length(tissue)) {
  ttype <- tissue[tct]
  tcode <- tcodes[tct]
  for (length in qltest) {
    ages <- subset(sample_info,Tissue == ttype & Length == length,
                   select=c(SampleName,Age))
    expr <- fpkm[,ages$SampleName]
    
    qlab <- "Q"
    if (length < 100) qlab <- paste0(qlab,0)
    if (length < 10) qlab <- paste0(qlab,0)
    qlab <- paste0(qlab,length)
    expcode <- paste0(tcode,"RA",qlab)
    
    exptab <- data.frame("Exp"=expcode,"Comparison"="Regression Age",
      "Tissue"=ttype,"Constant"=qlab,"Value"="Age")
    expinf <- bind_rows(expinf,exptab)

    restab <- data.frame("EnsemblID"=fpkm$EnsemblID,"Exp"=expcode,
      "Tissue"=rep(ttype,nrow(fpkm)),"Length"=rep(length,nrow(fpkm)))
    restab$PearsonR <- apply(expr,1,function(y) cor(ages$Age,y))
    restab[,c("Slope","Pvalue")] <-
      t(apply(expr,1,function(y) coef(summary(lm(y~ages$Age)))[c(2,8)]))
    if (exists("ageresults")) {
      ageresults <- bind_rows(ageresults,restab)
    } else {
      ageresults <- restab
    }
  }
}
rm(length,ttype,ages,expr,restab,exptab,qlab)
ageresults$Pvalue[is.na(ageresults$Pvalue)] <- 1
ageresults$PearsonR[is.na(ageresults$PearsonR)] <- 0
ageresults <- left_join(ageresults,ensembl,by="EnsemblID")
ageresults <- ageresults[c("EnsemblID","GeneName","Exp","Tissue","Length",
                         "PearsonR","Slope","Pvalue")]
ageresults$absPR <- abs(ageresults$PearsonR)
ageresults$absRise <- abs(ageresults$Slope*8)
ageresults$SigGene <- 0
ageresults[ageresults$Pvalue < pvmax & ageresults$absPR > rcut &
            ageresults$absRise > dfpcut,"SigGene"] <- 1
agesiglm <- ageresults %>% filter(SigGene == 1)
write.table(ageresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/fpkm_",genesrc,"_expr_age_lm.txt"))

ftime <- proc.time()
print("Linear regression time:")
print(ftime - stime)

# Form summary table of significant expression changes
sigtable <- bind_rows(qlsiglm %>% select(EnsemblID,Exp) %>% mutate(Match = 1),
                 agesiglm %>% select(EnsemblID,Exp) %>% mutate(Match = 1))
# expids <- bind_rows(ageresults %>% distinct(Exp),
#                    qlresults %>% distinct(Exp)) %>% arrange(Exp)

# Need to add in zero match line for all tests to account for tests that have
# zero significant matches
tmplg <- agesiglm[1,"EnsemblID"]
sigtable <- bind_rows(sigtable,data.frame(EnsemblID=tmplg,Exp=expinf$Exp,Match=0))
sigtable <- sigtable %>% group_by(EnsemblID,Exp) %>% summarize(Match = sum(Match))

dsummary <- dcast(sigtable,EnsemblID ~ Exp,value.var="Match")
dsummary[is.na(dsummary)] <- 0
dsummary$Total <- dsummary %>% select(-EnsemblID) %>% rowSums()
dsummary <- dsummary %>% left_join(ensembl,by="EnsemblID")
dsummary <- dsummary[c(1,ncol(dsummary),ncol(dsummary)-1,2:(ncol(dsummary)-2))]
dsummary <- dsummary %>% arrange(desc(Total),GeneName,EnsemblID)

write.table(dsummary,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/fpkm_",genesrc,"_summary_sig.txt"))

expinf <- expinf %>% arrange(Exp)
write.table(t(expinf),sep="\t",col.names=FALSE,quote=FALSE,na="",
            file=paste0(genesrc,"/fpkm_",genesrc,"_exp_info.txt"))

# hpgenes <- dsummary %>% filter(Total > 2 & !is.na(GeneName)) %>%
#  select(EnsemblID) %>% unlist(use.names = F)

if (! exists("hpgenes")) {
  hpgenes <- dsummary[1:3,1]
}

# Plot FPKM vs. Q Length and Age for select genes
plotgene <- function(sterm,goption="") {
  geneg <- ensembl[ensembl$EnsemblID == sterm | ensembl$GeneName == sterm,1]
  if(length(geneg)==0) {
    return(paste0("plotgene error: ", sterm, " is not a Ensembl ID or symbol"))
  }
  gsy <- ensembl[ensembl$EnsemblID==geneg,2]

  pg <- list()
  plc <- 1
  for (ttype in tissue) {
    qldata <- subset(sample_info,Tissue == ttype & Age %in% agetest,
                     select=c(SampleName,Age,Length))
    qldata$Expr <- t(fpkm[fpkm$EnsemblID==geneg,qldata$SampleName])
    qltitle <- paste0(gsy," ", ttype)
    qldata$Age <- as.factor(qldata$Age)
    pg[[plc]] <- ggplot(qldata, aes(x=Length, y=Expr, color=Age, shape=Age)) +
      geom_point() + geom_smooth(method=lm) + ggtitle(qltitle) +
      labs(x="Q Length",y="FPKM") +
      theme(plot.title = element_text(hjust = 0.5))
    plc <- plc + 1
  }

  plgt=textGrob(paste0(gsy, " Length Effect"),
                gp=gpar(fontsize=18,fontface="bold"))
  mgplot <- arrangeGrob(pg[[1]],pg[[2]],pg[[3]],
        pg[[4]],ncol=2,top=plgt)
  if (goption != 'pdfonly') {
    grid.newpage()
    grid.draw(mgplot)
  }
  if (goption %in% c('pdfonly','makepdf')) {
    plname <- paste0(genesrc, "/", gsy, "_qlfpkm.pdf")
    ggsave(plname,mgplot,width=8,height=10,units="in")
  }

  pg <- list()
  plc <- 1
  for (ttype in tissue) {
    agedata <- subset(sample_info,Tissue == ttype & Length %in% qltest,
                      select=c(SampleName,Age,Length))
    agedata$Expr <- t(fpkm[fpkm$EnsemblID==geneg,agedata$SampleName])
    agetitle <- paste0(gsy," ", ttype)
    agedata$Length <- as.factor(agedata$Length)
    pg[[plc]] <- ggplot(agedata, aes(x=Age, y=Expr, color=Length, shape=Length)) +
      geom_point() + geom_smooth(method=lm) + ggtitle(agetitle) +
      labs(x="Age",y="FPKM") +
      theme(plot.title = element_text(hjust = 0.5))
    plc <- plc + 1
  }
  
  plgt=textGrob(paste0(gsy, " Age Effect"),
                gp=gpar(fontsize=18,fontface="bold"))
  mgplot <- arrangeGrob(pg[[1]],pg[[2]],pg[[3]],
                        pg[[4]],ncol=2,top=plgt)
  if (goption != 'pdfonly') {
    grid.newpage()
    grid.draw(mgplot)
  }
  if (goption %in% c('pdfonly','makepdf')) {
    plname <- paste0(genesrc,"/", gsy, "_agefpkm.pdf")
    ggsave(plname,mgplot,width=8,height=10,units="in")
  }
}

if (plotoption != "none") {
  for (dgene in hpgenes) {
     plotgene(dgene,plotoption)
  }
}

etime <- proc.time()
print("Plotting time:")
print(etime - ftime)

print("Total time:")
print(etime - stime)
