require(RMySQL)
options("stringsAsFactors" = FALSE, warn = -1)
library(dplyr)

#wesdb <- 'WES_LAST'
#wesdb <- 'WES_NEW'

wesdb <- 'WES'

afcut <- 0.005

stime <- proc.time()

# Family IDs:  eg. 556, etc.
# Can be argument in Rscript command

if (! exists("Family")) Family <- commandArgs(trailingOnly=TRUE)[1]

if (is.na(Family)) stop("Enter Family ID on command line")

ipFile <- paste0(Family,'_input.txt')
if (! file.exists(ipFile)) {
  stop(paste0("Error: ",ipFile," not in working directory"))
}
famInfo <- read.table(ipFile,col.names=c('external_id','relation'))

con <- dbConnect(MySQL(),host="gleesondb.ucsd.edu",username='schview',
                 password='schview',dbname=wesdb)

getidstr <- paste0("select GleesonID, W.id as WES_id, sex_PDB, sex_COV, ",
  "affected, P.plate_id, L.name as Plate ",
  "from Annotations.Patients P, Annotations.Plates L, ", wesdb,
  ".Patients W ", "where W.external_id = ")
exteq <- " and P.external_id = W.external_id and L.id = P.plate_id"

for (i in 1:nrow(famInfo)) {
  extid <- famInfo[i,"external_id"]
  if (extid == 0) next
  getext <- paste0(getidstr,extid,exteq)
  dbret <- dbFetch(dbSendQuery(con,getext),n=-1)
  if (is.na(dbret[1,'WES_id'])) {
    stop(paste0("Error: external_id ",extid," is not in WES database"))
  }
  sex <- dbret[1,'sex_PDB']
  if (is.na(sex)) {
    sex <- dbret[1,'sex_COV']
  }
  if (is.na(sex)) {
    sex <- 'U'
  }
  famInfo[i,'GleesonID'] <- dbret[1,"GleesonID"]
  famInfo[i,'WES_id'] <- dbret[1,"WES_id"]
  famInfo[i,'Plate'] <- dbret[1,"Plate"]
  famInfo[i,'sex'] <- sex
  famInfo[i,'affected'] <- dbret[1,"affected"]
}

infoFile <- paste0(Family,'_info.txt')
sink(infoFile,split=TRUE)
print("Family member info:")
print(famInfo)
cat("\n")

gPrbvar <- "select C1.variant_id, V.CHROM, V.pos, V.REF as cRef, V.ALT as cAlt,
    V.clinvar_rs as rsid, V.1000Gp3_AF as freq100G, V.ExAC_AF as freqExAC,
    S.nHomSharesUnAff as homoShares, S.nHetSharesUnAff as heteroShares,
    V.MutationTaster_score as mutationTaster, V.GERP___RS as gerp,
    V.CADD_phred as cadd, C1.GT as zygProband, C1.DP as totDepthP,
    V.phyloP20way_mammalian_rankscore as phyloP,
    V.phastCons20way_mammalian as phast, C1.AD_ALT as varDepthP,
    C1.GQ as qualP, C1.denovo, S.AF as sharesAF
  from Calls C1, Variants V, Shares S
  where C1.patient_id = "

gPrbvar <- paste0(gPrbvar,famInfo[1,'WES_id'],"
  and V.id = C1.variant_id
  and S.variant_id = C1.variant_id")

Prbvar <- dbFetch(dbSendQuery(con,gPrbvar),n=-1)

print(paste0("Variants in Proband: ",nrow(Prbvar)))

# Remove Proband 0/1 variants on X chromosome (23) if proband is male
# Pseudoautosomal Regions:
# PAR1: chrX:10001-2781479 and chrY:10001-2781479
# PAR2: chrX:155701383-156030895 and chrY:56887903-57217415
if (famInfo[1,'sex'] == "M") {
  Prbvar <- Prbvar[!(Prbvar$CHROM == 23 & Prbvar$zygProband==1 &
              ((Prbvar$pos < 10001) | ((Prbvar$pos > 2781479) &
               (Prbvar$pos < 155701383)) | (Prbvar$pos > 156030895))),]
  print(paste0("Variants after X heterozygous removal: ",nrow(Prbvar)))
}
cat("\n")
Prbvar[,"HG38:chrPos"] <- paste0(Prbvar$CHROM,':',Prbvar$pos)
Prbvar$zygProband[Prbvar$zygProband == 1] <- "0/1"
Prbvar$zygProband[Prbvar$zygProband == 2] <- "1/1"

pdct <- nrow(Prbvar[Prbvar$denovo==1 & (! is.na(Prbvar$sharesAF)) &
  (Prbvar$sharesAF < afcut),c("HG38:chrPos","sharesAF")])
print(paste0("Potential denovo variants with AF < ",afcut,': ',pdct))

# Grab calls for all family members
zlab <- c("zygProband")
tdlab <- c("totDepthP")
vdlab <- c("varDepthP")
qlab <- c("qualP")
for (i in 2:nrow(famInfo)) {
  trel <- famInfo[i,'relation']
  twid <- famInfo[i,'WES_id']
  zlab[i] <- paste0('zyg',trel)
  tdlab[i] <- paste0('totDepth',trel)
  vdlab[i] <- paste0('varDepth',trel)
  qlab[i] <- paste0('qual',trel)
  if (is.na(twid)) next
  gCall <- paste0("select C1.variant_id, C1.GT as ", zlab[i], ", C1.DP as ",
    tdlab[i], ", C1.AD_ALT as ", vdlab[i], ", C1.GQ as ", qlab[i],
    " from Calls C1 WHERE C1.patient_id = ", famInfo[i,'WES_id'])
  CallDat <- dbFetch(dbSendQuery(con,gCall),n=-1)

  CallDat[,zlab[i]][CallDat[,zlab[i]] == 1] <- "0/1"
  CallDat[,zlab[i]][CallDat[,zlab[i]] == 2] <- "1/1"
  Prbvar <- left_join(Prbvar,CallDat,by="variant_id")
}

dtime <- proc.time()
cat("\nDatabase querying:\n")
dtime - stime
cat("\n")

# Get Homozygous WT coverage info
GTdf <- Prbvar[,c('CHROM','pos',zlab[2:length(zlab)])]
getGTv <- unique(GTdf[rowSums(is.na(GTdf)) > 0,c('CHROM','pos')])
getGTn <- nrow(getGTv)
print(paste0("Variants to get from VCF: ", getGTn))

opbcf <- paste0(Family,'_bcf.txt')
if ((NROW(getGTv) > 0) && (! file.exists(opbcf))) {
  print("Retrieving homozygoups WT info from VCF file")
  extstr <- paste0(famInfo[2:nrow(famInfo),"external_id"],collapse=",")
  idvcf <- paste0(Family,'_getgt.txt')
  write.table(getGTv,file=idvcf,sep="\t",quote=F,row.names=F,col.names=F)
  getgt <- paste0('bcftools query -f ',
      '\'%CHROM:%POS[\\t%GT\\t%DP\\t%AD\\t%GQ]\\n\'',
      ' -s ', extstr, ' -R ', idvcf,
      ' /raid/data/plates/current_WES_VCF.gz > ', opbcf)
#      '\'%CHROM:%POS[\\t%SAMPLE\\t%GT\\t%DP\\t%AD\\t%GQ]\\n\'',
  sysop <- system(getgt,intern = TRUE)
} else if (NROW(getGTv) > 0) {
  print("Retrieving homozygoups WT info from previous VCF scan")
} else {
  print("No VCF GT information needs retrieval")
}
if (file.exists(opbcf)) {
  GTres <- read.table(opbcf)
  for (i in 2:nrow(famInfo)) {
    xad <- 4*(i - 1)
    GTres[,xad] <- sub('^.*,([0-9]+)$','\\1',GTres[,xad])
  }
  colnames(GTres)[1] <- 'HG38:chrPos'
  mgcol <- ncol(Prbvar)
  Prbvar <- left_join(Prbvar,GTres,by="HG38:chrPos")
  for (j in 2:nrow(famInfo)) {
    xj <- 4*(j-2) + mgcol + 1
    Prbvar[is.na(Prbvar[,zlab[j]]) & !is.na(Prbvar[,xj]) &
      (Prbvar[,xj] == '0/0'),tdlab[j]] <- Prbvar[is.na(Prbvar[,zlab[j]]) &
      !is.na(Prbvar[,xj]) & (Prbvar[,xj] == '0/0'),1+xj]
    Prbvar[is.na(Prbvar[,zlab[j]]) & !is.na(Prbvar[,xj]) &
      (Prbvar[,xj] == '0/0'),vdlab[j]] <- Prbvar[is.na(Prbvar[,zlab[j]]) &
      !is.na(Prbvar[,xj]) & (Prbvar[,xj] == '0/0'),2+xj]
    Prbvar[is.na(Prbvar[,zlab[j]]) & !is.na(Prbvar[,xj]) &
      (Prbvar[,xj] == '0/0'),qlab[j]] <- Prbvar[is.na(Prbvar[,zlab[j]]) &
      !is.na(Prbvar[,xj]) & (Prbvar[,xj] == '0/0'),3+xj]
    Prbvar[is.na(Prbvar[,zlab[j]]) & !is.na(Prbvar[,xj]) &
      (Prbvar[,xj] == '0/0'),zlab[j]] <- '0/0'
  }
}
vtime <- proc.time()
cat("\nVCF retrieval and GT merge:\n")
vtime - dtime
cat("\n")

# Remove heterozygous variants in X chromosome for males
for (i in 2:nrow(famInfo)) {
  if (famInfo[i,'sex'] == 'M') {
    Prbvar[which(Prbvar$CHROM == 23 & Prbvar[,zlab[i]]=='0/1' &
      ((Prbvar$pos < 10001) | ((Prbvar$pos > 2781479) &
       (Prbvar$pos < 155701383)) | (Prbvar$pos > 156030895))),
       zlab[i]] <- NA
  }
}

gIso <- "select I.id as isoform_id, I.variant_id, I.EnsemblID,
  I.effect, I.Polyphen2_HVAR_score as polyPhen, I.SIFT_score as sift,
  I.codon_position as pPosition, I.REF_aa as pRef, I.ALT_aa as pAlt
  from Calls C, Isoforms I
  WHERE C.patient_id = "
gIso <- paste0(gIso,famInfo[1,'WES_id'],"
  and I.variant_id = C.variant_id")
Isovar <- dbFetch(dbSendQuery(con,gIso),n=-1)

getEns <- paste0("select distinct EnsemblID, name as hgncSymbol,",
  "MIM_disease as omimDiseaseNames from Annotations.Genes_hg38")
Ensembl <- dbFetch(dbSendQuery(con,getEns),n=-1)
Ensembl$omimDiseaseNames <- gsub(".$","",Ensembl$omimDiseaseNames)
Ensembl$omimDiseaseNames <- gsub("\\[MIM:.*?\\]","",Ensembl$omimDiseaseNames)
Ensembl$omimDiseaseNames <-gsub(";","; ",Ensembl$omimDiseaseNames)

Isomerge <- left_join(Isovar,Ensembl,by="EnsemblID")
Isomerge <- Isomerge[!is.na(Isomerge$hgncSymbol),]
print(paste0(nrow(Isomerge), " isoform rows"))

GranthamS <- dbFetch(dbSendQuery(con,paste0("select REF_aa as pRef, ",
  "ALT_aa as pAlt, grantham_score as grantham from ",
  "Annotations.GranthamScore")),n=-1)

Isomerge <- left_join(Isomerge,Prbvar,by="variant_id") %>%
  filter(!is.na(`HG38:chrPos`))

Isomerge <- left_join(Isomerge,GranthamS,by=c("pRef","pAlt"))
print(paste0(nrow(Isomerge), " variant source rows"))

SCfields <- c("CHROM", "pos", "hgncSymbol", "geneNameLong", "HG38:chrPos",
  "cSeqAnnotation", "cPosition", "cRef", "cAlt", "pSeqAnnotation",
  "pPosition", "pRef", "pAlt", "rsid", zlab, "effect", "freq1000G",
  "freqExAC", "homoShares", "heteroShares", "omimNumber", "omimDiseaseNames",
  "variantAccession", "variantPathogenicity", "polyPhen", "mutationTaster",
  "sift", "gerp", "grantham", "phat", "phast", "phyloP", "strandBias",
  "knownSplice")
for (i in 1:NROW(tdlab)) {
  SCfields <- c(SCfields,tdlab[i],vdlab[i],qlab[i])
}

Sifields = colnames(Isomerge)
for (scol in SCfields) {
  if (!(scol %in% Sifields)) {
    Isomerge[,scol] <- NA
  }
}

SimulOP <- Isomerge[,SCfields]
SimulOP <- SimulOP %>% unique() %>% arrange(CHROM,pos,hgncSymbol)
print(paste0(nrow(SimulOP), " output rows"))
SimulOP$CHROM <- NULL
SimulOP$pos <- NULL

opFile <- paste0(Family,'_SCtable.tsv')
cat("fileformat=gleesonVersion2\n",file=opFile)
write.table(SimulOP,file=opFile,sep='\t',quote=F,row.names=F,na="",append=TRUE)

ftime <- proc.time()

cat("\nAnnotations + Output:\n")
ftime - vtime
cat("\nTotal:\n")
ftime - stime
