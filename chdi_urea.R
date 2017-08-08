# This script needs to be excecuted from working directory
# The working directory needs to have the provided Ensembl information file
# in it, and the uncompressed directories that contain the cortex, liver,
# and striatum experimental conditions and results

# Turn off automatic factorization of character vectors. With all the sample
# and gene IDs, it is pointless, and it gets in the way of subsets based
# on vector.

stime <- proc.time()
library(DESeq2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(pheatmap)

options(stringsAsFactors = FALSE)

rcut <- 0.7
fcut <- 1.5
pvmax <- 0.05
dfpcut <- 2

# Set various parameters
dprefix <- c("cortex","liver","striatum")
# files <- c("sample_info","counts","counts_avg","counts_avg_na",
#           "fpkm","fpkm_avg","fpkm_avg_na")
files <- c("sample_info","counts","fpkm")
tissue <- c("Cortex","Liver","Striatum")
agetest <- c(6,10)
qltest <- c(7,80,111,140,175)
diffql <- c(20,80,92,111,140,175)
diffaq <- c(80,111,175)

# Genes of interest: Ass1, Asl, Otc Arg1, and Cps1
ucenz <- c("ENSMUSG00000076441","ENSMUSG00000025533","ENSMUSG00000031173",
           "ENSMUSG00000019987","ENSMUSG00000025991")

# Outliers from external analysis (can incorporate into code if desired)
outliers <- c("F2Q80C_3","M2Q80C_5","M2Q92C_3","M2Q80S_5")

# Read Ensembl name file. This assumes a unique mapping, so the file should
# not have any repeats of EnsemblID
ensembl <- read.table("ensmouse_table_uq.txt",sep="\t",quote="",header=TRUE)
ensembl$Description <- sub("\\s+$", "", ensembl$Description)

# Read data files 
for (ttype in dprefix) {
  for (ftype in files) {
    tabname <- paste0(ttype,"_",ftype)
    fname <- paste0(ttype,"\\",tabname,".csv")
    assign(tabname,read.table(fname,sep="\t", header=TRUE,row.names=1))
  }
}
rm(ttype,ftype,tabname,fname)

# Combine data from tissue types into one structure (easier processing)
# Convert Q-length to digit, for linear modeling purposes
sample_info <- rbind(cortex_sample_info,liver_sample_info,striatum_sample_info)
sample_info$Length <- as.integer(sub("WT","7",sub("^Q","",sample_info$Genotype)))
rm(cortex_sample_info,liver_sample_info,striatum_sample_info)

fpkm <- cbind(cortex_fpkm,liver_fpkm,striatum_fpkm)
rm(cortex_fpkm,liver_fpkm,striatum_fpkm)

counts <- cbind(cortex_counts,liver_counts,striatum_counts)
rm(cortex_counts,liver_counts,striatum_counts)

# Remove outliers
sample_info <- sample_info[!(rownames(sample_info) %in% outliers),]
fpkm <- fpkm[,!(names(fpkm) %in% outliers)]
counts <- counts[,!(names(counts) %in% outliers)]

# Calculate differential expression with DESeq2, comparing specified Q lengths to
# the wildtype (Q length 7), for a combination of the tissue types and ages 6 and
# 10 months.
if (exists("expinf")) {
  rm(expinf)
}
if (exists("qldiffe")) {
  rm(qldiffe)
}
for (ttype in tissue) {
  for (age in agetest) {
    qinfo <- subset(sample_info,Tissue == ttype & Age == age,select=Length)
    qinfo$Length <- as.factor(qinfo$Length)
    expr <- counts[,rownames(qinfo)]

    deds <- DESeqDataSetFromMatrix(countData = expr,
            colData = qinfo, design = ~ Length) 
    deds <- DESeq(deds)
    
    for (gtype in diffql) {
      qlab <- "Q"
      if (gtype < 100) qlab <- paste0(qlab,0)
      if (gtype < 10) qlab <- paste0(qlab,0)
      qlab <- paste0(qlab,gtype)
      alab <- age
      if (age < 10) alab <- paste0(0,alab)
      
      expcode <- paste0(substr(ttype,0,1),"DLA",alab,qlab)

      exptab <- data.frame("Exp"=expcode,"Comparison"="Differential Q-length",
        "Tissue"=ttype,"Constant"=paste0(alab," mo"),
        "Value"=paste0(qlab,":WT"))
      if (exists("expinf")) {
        expinf <- bind_rows(expinf,exptab)
      } else {
        expinf <- exptab
      }
        
      dres <- results(deds, contrast=c("Length",gtype,7))
      restab <- data.frame("EnsemblID"=rownames(dres),
        "Tissue"=rep(ttype,nrow(dres)), "Age"=rep(age,nrow(dres)),
        "Diff_QLength"=rep(gtype,nrow(dres)), "log2FC"=dres$log2FoldChange,
        "padj"=dres$padj, "Exp"=expcode)
      if (exists("qldiffe")) {
        qldiffe <- bind_rows(qldiffe,restab)
      } else {
        qldiffe <- restab
      }
    }
  }
}
rm(qinfo)

# Add gene symbols to results (expired Ensembl IDs will have none and be NA)
# Set undefined adjusted P-values to 1, and undefined log2 fold changes to 0 -
# this keeps the vectors numerical
# Indicate experiments with significant result

qldiffe <- left_join(qldiffe,ensembl,by="EnsemblID") %>% select(-Description)
qldiffe <- qldiffe[c("EnsemblID","GeneName","Exp","Tissue","Age",
  "Diff_QLength","log2FC","padj")]
qldiffe$padj[is.na(qldiffe$padj)] <- 1
qldiffe$log2FC[is.na(qldiffe$log2FC)] <- 0
qldiffe$abslfc <- abs(qldiffe$log2FC)
qldiffe$SigGene <- 0
qldiffe[qldiffe$padj < pvmax & qldiffe$abslfc > log2(fcut),"SigGene"] <- 1
qlsigde <- qldiffe %>% filter(SigGene == 1)
qlsigdect <- qlsigde %>% group_by(Tissue,Age,Diff_QLength) %>%
  summarize(Sig_Genes=n())
qlsigdect <- bind_rows(qlsigdect,qldiffe %>%
                         group_by(Tissue,Age,Diff_QLength) %>% summarize(Sig_Genes=0))
qlsigdect <- qlsigdect %>% group_by(Tissue,Age,Diff_QLength) %>%
  summarize(Sig_Genes=sum(Sig_Genes))
write.table(qlsigdect,sep="\t",row.names=FALSE,quote=FALSE,
            file="output\\diff_expr_qlsigct.txt")
write.table(qldiffe,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file="output\\diff_expr_ql.txt")

# Calculate differential expression with DESeq2, comparing ages 6 and 10 months
# to the age two months, for a combination of the tissue types and specified
# Q lengths
if (exists("agediffe")) {
  rm(agediffe)
}

for (ttype in tissue) {
  for (gtype in diffaq) {
    ainfo <- subset(sample_info,Tissue == ttype & Length == gtype,select=Age)
    ainfo$Age <- as.factor(ainfo$Age)
    expr <- counts[,rownames(ainfo)]

    deds <- DESeqDataSetFromMatrix(countData = expr,
            colData = ainfo, design = ~ Age) 
    deds <- DESeq(deds)
    
    for (age in agetest) {
      qlab <- "Q"
      if (gtype < 100) qlab <- paste0(qlab,0)
      if (gtype < 10) qlab <- paste0(qlab,0)
      qlab <- paste0(qlab,gtype)
      alab <- age
      if (age < 10) paste0(0,alab)
      expcode <- paste0(substr(ttype,0,1),"DA",qlab,"A",alab)

      exptab <- data.frame("Exp"=expcode,"Comparison"="Differential Age",
         "Tissue"=ttype,"Constant"=qlab,"Value"=paste0(alab," mo: 2 mo"))
      expinf <- bind_rows(expinf,exptab)

      dres <- results(deds, contrast=c("Age",age,2))
      restab <- data.frame("EnsemblID"=rownames(dres),
        "Tissue"=rep(ttype,nrow(dres)),"Length"=rep(gtype,nrow(dres)),
        "Diff_Age"=rep(age,nrow(dres)),"log2FC"=dres$log2FoldChange,
        "padj"=dres$padj, "Exp"=expcode)
      if (exists("agediffe")) {
        agediffe <- bind_rows(agediffe,restab)
      } else {
        agediffe <- restab
      }
    }
  }
}
rm(ainfo,gtype,ttype,age,expr,deds,restab,dres,qlab,alab,exptab)

agediffe <- left_join(agediffe,ensembl,by="EnsemblID") %>% select(-Description)
agediffe <- agediffe[c("EnsemblID","GeneName","Exp","Tissue","Length",
                       "Diff_Age","log2FC","padj")]
agediffe$padj[is.na(agediffe$padj)] <- 1
agediffe$log2FC[is.na(agediffe$log2FC)] <- 0
agediffe$abslfc <- abs(agediffe$log2FC)
agediffe$SigGene <- 0
agediffe[agediffe$padj < pvmax & agediffe$abslfc > log2(fcut),"SigGene"] <- 1
agesigde <- agediffe %>% filter(SigGene == 1)
agesigdect <-agesigde %>% group_by(Tissue,Length,Diff_Age) %>%
  summarize(Sig_Genes=n())
agesigdect <- bind_rows(agesigdect,agediffe %>%
                          group_by(Tissue,Length,Diff_Age) %>% summarize(Sig_Genes=0))
agesigdect <- agesigdect %>% group_by(Tissue,Length,Diff_Age) %>%
  summarize(Sig_Genes=sum(Sig_Genes))
write.table(agesigdect,sep="\t",row.names=FALSE,quote=FALSE,
            file="output\\diff_expr_agesigct.txt")
write.table(agediffe,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file="output\\diff_expr_age.txt")

dtime <- proc.time()
print("DESeq2 time:")
print(dtime - stime)

# Linear regressions correlations for 10 months and 6 months
if (exists("qlresults")) {
  rm(qlresults)
}
for (ttype in tissue) {
  for (age in agetest) {
    qlength <- subset(sample_info,Tissue == ttype & Age == age,select=Length)
    expr <- fpkm[,rownames(qlength)]
    
    alab <- age
    if (age < 10) paste0(0,alab)
    expcode <- paste0(substr(ttype,0,1),"RL","A",alab)
    
    exptab <- data.frame("Exp"=expcode,"Comparison"="Regression Q-length",
      "Tissue"=ttype,"Constant"=paste0(alab," mo"),"Value"="Q-length")
    expinf <- bind_rows(expinf,exptab)

    restab <- data.frame("EnsemblID"=rownames(fpkm),"Exp"=expcode,
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
qlresults <- left_join(qlresults,ensembl,by="EnsemblID") %>% select(-Description)
qlresults <- qlresults[c("EnsemblID","GeneName","Exp","Tissue","Age",
                       "PearsonR","Slope","Pvalue")]
qlresults$absPR <- abs(qlresults$PearsonR)
qlresults$absRise <- abs(qlresults$Slope*168)
qlresults$SigGene <- 0
qlresults[qlresults$Pvalue < pvmax & qlresults$absPR > rcut &
            qlresults$absRise > dfpcut,"SigGene"] <- 1
qlsiglm <- qlresults %>% filter(SigGene == 1)
write.table(qlresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file="output\\expr_qlength_lm.txt")

# Linear regressions and correlations for q-lengths of 175, 111, 80, and 7 (WT)
if (exists("ageresults")) {
  rm(ageresults)
}
for (ttype in tissue) {
  for (length in qltest) {
    ages <- subset(sample_info,Tissue == ttype & Length == length,select=Age)
    expr <- fpkm[,rownames(ages)]
    
    qlab <- "Q"
    if (length < 100) qlab <- paste0(qlab,0)
    if (length < 10) qlab <- paste0(qlab,0)
    qlab <- paste0(qlab,length)
    expcode <- paste0(substr(ttype,0,1),"RA",qlab)
    
    exptab <- data.frame("Exp"=expcode,"Comparison"="Regression Age",
      "Tissue"=ttype,"Constant"=qlab,"Value"="Age")
    expinf <- bind_rows(expinf,exptab)
                         
    restab <- data.frame("EnsemblID"=rownames(fpkm),"Exp"=expcode,
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
ageresults <- left_join(ageresults,ensembl,by="EnsemblID") %>% select(-Description)
ageresults <- ageresults[c("EnsemblID","GeneName","Exp","Tissue","Length",
                         "PearsonR","Slope","Pvalue")]
ageresults$absPR <- abs(ageresults$PearsonR)
ageresults$absRise <- abs(ageresults$Slope*8)
ageresults$SigGene <- 0
ageresults[ageresults$Pvalue < pvmax & ageresults$absPR > rcut &
            ageresults$absRise > dfpcut,"SigGene"] <- 1
agesiglm <- ageresults %>% filter(SigGene == 1)
write.table(ageresults,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file="output\\expr_age_lm.txt")
ftime <- proc.time()

print("Linear regression time:")
print(ftime - dtime)

# Form summary table of significant expression changes
sigtable <- bind_rows(qlsigde %>% select(EnsemblID,Exp) %>% mutate(Match = 1),
                 agesigde %>% select(EnsemblID,Exp) %>% mutate(Match = 1),
                 qlsiglm %>% select(EnsemblID,Exp) %>% mutate(Match = 1),
                 agesiglm %>% select(EnsemblID,Exp) %>% mutate(Match = 1))
expids <- bind_rows(agediffe %>% distinct(Exp),ageresults %>% distinct(Exp),
  qldiffe %>% distinct(Exp),qlresults %>% distinct(Exp)) %>% arrange(Exp)

# Need to had in zero match line for all tests to account for tests that have
# zero significant matches
tmplg <- agesigde[1,"EnsemblID"]
sigtable <- bind_rows(sigtable,data.frame(EnsemblID=tmplg,Exp=expids,Match=0))
sigtable <- sigtable %>% group_by(EnsemblID,Exp) %>% summarize(Match = sum(Match))

dsummary <- dcast(sigtable,EnsemblID ~ Exp,value.var="Match")
dsummary[is.na(dsummary)] <- 0
dsummary$Total <- dsummary %>% select(-EnsemblID) %>% rowSums()
dsummary <- dsummary %>% left_join(ensembl,by="EnsemblID") %>% select(-Description)
dsummary <- dsummary[c(1,ncol(dsummary),ncol(dsummary)-1,2:(ncol(dsummary)-2))]
dsummary <- dsummary %>% arrange(desc(Total),GeneName,EnsemblID)

write.table(dsummary,sep="\t",row.names=FALSE,quote=FALSE,na="",
            file="output\\summary_sig.txt")

expinf <- expinf %>% arrange(Exp)
write.table(t(expinf),sep="\t",col.names=FALSE,quote=FALSE,na="",
            file="output\\exp_info.txt")

hpgenes <- dsummary %>% filter(Total > 6 & !is.na(GeneName)) %>%
  select(EnsemblID) %>% unlist(use.names = F)
hpgenes <- c(hpgenes,ucenz)

# Draw log2FC heatmaps from a select set of genes
drheatmap <- function(genelist) {
  if(length(genelist)==0) {
    return(paste0("Heatap error: gene list has no length"))
  }  
  for (ttype in tissue) {
    qlgraph <- qldiffe %>% subset(EnsemblID %in% genelist) %>%
      filter(Tissue==ttype) %>% select(GeneName,Age,Diff_QLength,log2FC)
    qlgraph$Age <- factor(qlgraph$Age,agetest)
    qlgraph$Diff_QLength <- factor(qlgraph$Diff_QLength,diffql)
    qlhmat <- acast(qlgraph,GeneName ~ Age + Diff_QLength,value.var="log2FC")
    qlgrlab <- data.frame("Age"=gsub("_.*","",colnames(qlhmat)),
                          "DiffQL"=gsub(".*_","",colnames(qlhmat)))
    qlgrlab$Age <- paste0(qlgrlab$Age," mo")
    qlgrlab$Age <- factor(qlgrlab$Age,c("6 mo","10 mo"))
    qlgrlab$DiffQL <- paste0("Q",qlgrlab$DiffQL,":WT")
    qlgrlab$DiffQL <- factor(qlgrlab$DiffQL,c("Q20:WT","Q80:WT","Q92:WT",
                                              "Q111:WT","Q140:WT","Q175:WT"))
    rownames(qlgrlab) <- colnames(qlhmat)
    qlgrtitle <- paste0(ttype," Q-Length Differential Expression")
    pheatmap(qlhmat,cluster_cols=FALSE,annotation_col=qlgrlab,show_colnames=F,
             show_rownames=F,main=qlgrtitle)
    qlpdf <- paste0("output\\",ttype,"_qlhmap.pdf")
    pheatmap(qlhmat,cluster_cols=FALSE,annotation_col=qlgrlab,show_colnames=F,
             show_rownames=F,main=qlgrtitle,filename=qlpdf)
    
    agegraph <- agediffe %>% subset(EnsemblID %in% genelist) %>%
      filter(Tissue==ttype) %>% select(GeneName,Diff_Age,Length,log2FC)
    agegraph$Length <- factor(agegraph$Length,diffaq)
    agegraph$Diff_Age <- factor(agegraph$Diff_Age,agetest)
    agehmap <- acast(agegraph,GeneName ~ Diff_Age + Length,value.var="log2FC")
    agegrlab <- data.frame("DiffAge"=gsub("_.*","",colnames(agehmap)),
                           "QLength"=gsub(".*_","",colnames(agehmap)))
    agegrlab$DiffAge <- paste0(agegrlab$DiffAge," mo:2 mo")
    agegrlab$DiffAge <- factor(agegrlab$DiffAge,c("6 mo:2 mo","10 mo:2 mo"))
    agegrlab$QLength <- factor(agegrlab$QLength,diffaq)
    rownames(agegrlab) <- colnames(agehmap)
    agegrtitle <- paste0(ttype," Age Differential Expression")
    pheatmap(agehmap,cluster_cols=FALSE,annotation_col=agegrlab,show_colnames=F,
             show_rownames=F,main=agegrtitle)
    agepdf <- paste0("output\\",ttype,"_agehmap.pdf")
    pheatmap(agehmap,cluster_cols=FALSE,annotation_col=agegrlab,show_colnames=F,
             show_rownames=F,main=agegrtitle,filename=agepdf)
  }
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
                     select=c(Age,Length))
    qldata$Expr <- t(fpkm[geneg,rownames(qldata)])
    qltitle <- paste0(gsy," ", ttype)
    qldata$Age <- as.factor(qldata$Age)
    pg[[plc]] <- ggplot(qldata, aes(x=Length, y=Expr, color=Age, shape=Age)) +
      geom_point() + geom_smooth(method=lm) + ggtitle(qltitle) +
      labs(x="Q Length",y="FPKM") +
      theme(plot.title = element_text(hjust = 0.5))
    plc <- plc + 1

    agedata <- subset(sample_info,Tissue == ttype & Length %in% qltest,
                      select=c(Age,Length))
    agedata$Expr <- t(fpkm[geneg,rownames(agedata)])
    agetitle <- paste0(gsy," ", ttype)
    agedata$Length <- as.factor(agedata$Length)
    pg[[plc]] <- ggplot(agedata, aes(x=Age, y=Expr, color=Length, shape=Length)) +
      geom_point() + geom_smooth(method=lm) + ggtitle(agetitle) +
      labs(x="Age",y="FPKM") +
      theme(plot.title = element_text(hjust = 0.5))
    plc <- plc + 1
  }

  plgt=textGrob(paste0(gsy, " Distribution"),
                gp=gpar(fontsize=18,fontface="bold"))
  mgplot <- arrangeGrob(pg[[1]],pg[[2]],pg[[3]],
        pg[[4]],pg[[5]],pg[[6]],ncol=2,top=plgt)
  if (goption != 'pdfonly') {
    grid.newpage()
    grid.draw(mgplot)
  }
  if (goption %in% c('pdfonly','makepdf')) {
    plname <- paste0("output\\geneplots\\", gsy, "_fpkm.pdf")
    ggsave(plname,mgplot,width=8,height=10,units="in")
  }
}

drheatmap(hpgenes)

for (dgene in hpgenes) {
  plotgene(dgene,"pdfonly")
}

etime <- proc.time()
print("Plotting time:")
print(etime - ftime)

print("Total time:")
print(etime - stime)
