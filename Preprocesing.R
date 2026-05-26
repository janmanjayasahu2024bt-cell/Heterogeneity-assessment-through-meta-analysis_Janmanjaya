setwd("C:/Users/College/Desktop/Dissertation/Restart/Predictive_Model")


getwd()
library(biomaRt)
library(limma)
library(tibble)
library(dplyr)
library(plyr)
library(HGNChelper)
library(oligo)
library(stringr)
library(affy)
BiocManager::install("lumi")
BiocManager::install("lumiHumanIDMapping")
library(lumi)
library(lumiHumanIDMapping)


#18312

# LIST .CEL FILES
cels.gse18312<-list.celfiles("RAW/", pattern="CEL")

cel_path <- "C:/Users/College/Desktop/Dissertation/Restart/18312/RAW"  

cel_files <- list.files(
  path = cel_path,
  pattern = "[cC][eE][lL](\\.gz)?$",
  full.names = TRUE
)

RAW18312 <- read.celfiles(cel_files)

# PERFORM RMA
GSE18312.rma <- oligo::rma(RAW18312)

RAW.GSE18312.DF = as.data.frame(exprs(GSE18312.rma))

#Load phenodata file
my_data_18312 <- read.csv("C:/Users/College/Desktop/Dissertation/Restart/18312/18312_Ph_data.csv")

# Editing column names of the raw expression matrix df

colnames(RAW.GSE18312.DF)=my_data_18312$GSM.ID
colnames(RAW.GSE18312.DF)

#/ GENE ANNOTATION

featureData(GSE18312.rma) = getNetAffx(GSE18312.rma, "transcript")

Ann.GSE18312.DF= as.data.frame(GSE18312.rma@featureData@data)
dim(Ann.GSE18312.DF)

# Check for unassigned ProbeIDs

table(is.na(Ann.GSE18312.DF$geneassignment))


# Keep only ProbeIDs and annotated gene info

Ann.GSE18312.DF=Ann.GSE18312.DF[,c(1,2,8)]


# Remove rows with unassigned ProbeIDs

Ann.GSE18312.omit=na.omit(Ann.GSE18312.DF)
dim(Ann.GSE18312.omit)#17638
table(is.na(Ann.GSE18312.omit$geneassignment))


# Extract gene symbols in a separate column
#Stringr package should be enabled

Ann.GSE18312.omit[,c(4:6)]=str_split_fixed(Ann.GSE18312.omit$geneassignment, "//",3)

Ann.GSE18312.omit=Ann.GSE18312.omit[,c(1,2,5)] # keep relevant columns


# Strip split may have white space.

Ann.GSE18312.omit$V5=trimws(Ann.GSE18312.omit$V5)

#/ CHECK FOR HGNC SYMBOLS

gsymbols=Ann.GSE18312.omit$V5
check.symbols=checkGeneSymbols(gsymbols, unmapped.as.na = TRUE, map= NULL, species="human")

# Remove unassigned genes

table(is.na(check.symbols$Suggested.Symbol))
check.symbols.df=na.omit(check.symbols)
table(is.na(check.symbols.df$Suggested.Symbol))

# Add approved gene symbols to the modified annotated matrix

Ann.GSE18312.omit.HGNC=merge(Ann.GSE18312.omit, check.symbols.df, by.x= "V5", by.y= "x")


# Add HGNC approved gene symbols to raw matrix

RAW.GSE18312.HGNC=rownames_to_column(RAW.GSE18312.DF, "transcriptclusterid")

RAW.GSE18312.final=merge(RAW.GSE18312.HGNC,Ann.GSE18312.omit.HGNC, by= "transcriptclusterid")

#/ LIMMA AVEREPS

GSE18312.avereps = as.data.frame(limma::avereps(RAW.GSE18312.final, RAW.GSE18312.final$Suggested.Symbol))

GSE18312.avereps.df = column_to_rownames(GSE18312.avereps, "Suggested.Symbol")
colnames(GSE18312.avereps.df)

GSE18312.avereps.final= GSE18312.avereps.df[,-c(1,23:28)]
dim(GSE18312.avereps.final)#17131

saveRDS(GSE18312.avereps.final, file = "18312_normalized.rds")

#27383


#/ LIST .CEL FILES

cels.gse27383<-list.celfiles("RAW/", pattern="CEL")

#/ READING RAW DATA

#RAW<-read.celfiles(list.files(pattern = '*CEL', full.names = TRUE))
#or
#RAW <- oligo::read.celfiles(cels.gse27383)
      # or oligo, depending on your platform

cel_path <- "C:/Users/College/Desktop/Dissertation/Restart/27383/RAW"  

cel_files <- list.files(
  path = cel_path,
  pattern = "[cC][eE][lL](\\.gz)?$",
  full.names = TRUE
)

RAW_27383 <- read.celfiles(cel_files)


#/ FETCHING PHENOTYPE DATA

##{ph.data<-read.table("", header = TRUE, sep = "\t", nrows = 21)}
#or
data_27383 <- read.csv("C:/Users/College/Desktop/Dissertation/Restart/27383/27383_Ph_Data.csv", header = TRUE, stringsAsFactors = FALSE)


#/ PERFORM RMA
#GSE27383.rma= rma(RAW, target="core")
GSE27383.rma <- oligo::rma(RAW_27383)

#/ RAW EXPRESSION MATRIX AS DATA FRAME
RAW.GSE27383.DF = as.data.frame(exprs(GSE27383.rma))

##Load phenodata file
my_data_27383 <- read.csv("C:/Users/College/Desktop/Dissertation/Restart/27383/27383_Ph_Data.csv")

# Editing column names of the raw expression matrix df
colnames(RAW.GSE27383.DF)=my_data_27383$Accession
colnames(RAW.GSE27383.DF)



#/ GENE ANNOTATION (Finished till here)

#featureData(GSE27383.rma) = getNetAffx(GSE27383.rma, "transcript")
#instead of this alternate code was used
##annotation <- AnnotationDbi::select(
#hgu133plus2.db,
#keys = rownames(expr),
#columns = c("SYMBOL", "GENENAME", "ENTREZID"),
#keytype = "PROBEID"
#)

GSE27383.rma.df<- as.data.frame(exprs(GSE27383.rma))


Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
GSE27383.pIDs=rownames(GSE27383.rma.df)
GS.pIDs.27383= getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
                     filters = "affy_hg_u133_plus_2",
                     values = GSE27383.pIDs,
                     mart = Martfunction)


# Column names edited

colnames(GSE27383.rma.df)=my_data_27383$Accession
colnames(GSE27383.rma.df)


#/ REMOVE UNASSIGNED PROBE-IDS

GS.pIDs.27383.df=GS.pIDs.27383[!(GS.pIDs.27383$hgnc_symbol==""),]

#/ VERIFY GENE SYMBOLS WITH HGNC

geneSymbols.verified<-checkGeneSymbols(GS.pIDs.27383.df$hgnc_symbol, species = "human")


# Exchange the columns with corrected values (from verified to unverified gene symbol df)

GS.pIDs.27383.df$hgnc_symbol<-geneSymbols.verified$Suggested.Symbol


# Creating blank spaces in place of NA and then removing the blank spaces 

GS.pIDs.27383.df[is.na(GS.pIDs.27383.df)] <- ""
GS.pIDs.27383.df=GS.pIDs.27383.df[!(GS.pIDs.27383.df$hgnc_symbol==""),]


#/ BRING ROW-NAMES TO COLUMN IN RAW EXPRESSSION DATA

GSE27383.rma.df=rownames_to_column(GSE27383.rma.df, "PIDs")
colnames(GS.pIDs.27383.df)[1]= "PIDs"



#/ ADD GENE SYMBOLS TO NORMALIZED DATA

GSE27383.geneSymbols= merge(GS.pIDs.27383.df,GSE27383.rma.df, by= "PIDs")

dim(GSE27383.geneSymbols)#43042


#/ AVERAGE OUT DUPLICATE GENE SYMBOLS ROW-WISE USING LIMMA AVEREPS

GSE27383.avreps=as.data.frame((limma::avereps(GSE27383.geneSymbols, GSE27383.geneSymbols$hgnc_symbol)))

dim(GSE27383.avreps)#21532



#/ HAVE UNIQUE GENE SYMBOLS AS ROW IDS

GSE27383.avreps.df=column_to_rownames(GSE27383.avreps, "hgnc_symbol")
GSE27383.avreps.final=GSE27383.avreps.df[-c(1)] # remove unwanted column
dim(GSE27383.avreps.final)#21532


abc = mutate_all(GSE27383.avreps.final, function(x) as.numeric(as.character(x)))

saveRDS(abc, file = "27383_normalized.rds")

#38481

#/ LISTING FILES IN THE DIRECTORY

txt.gse38481<-list.files("C:/Users/College/Desktop/Dissertation/Restart/38481/GSE38481/", pattern="txt")


#/ READING RAW DATA

RAW.DATA.GSE38481 <- lumiR("C:/Users/College/Desktop/Dissertation/Restart/38481/GSE38481_non_normalised.txt", convertNuID = TRUE, lib.mapping = 'lumiHumanIDMapping', QC = TRUE,  dec = '.', parseColumnName = TRUE, columnNameGrepPattern = list(exprs='AVG_SIGNAL', se.exprs='NA', detection='DETECTION'))
dim(RAW.DATA.GSE38481)
#24526

#/ FILTERING EXPRESSED PROBES by taking the average of 3 probes<0.05

EXPRESSED.PROBES<- rowSums(RAW.DATA.GSE38481@assayData$detection<0.05)>=3


# Subtract from main object or retain only expressed probes in expression matrix

GSE38481.RAW=RAW.DATA.GSE38481[EXPRESSED.PROBES,]
dim(GSE38481.RAW) #16023

## FETCHING PHENOTYPE DATA

phData_38481<- read.csv("C:/Users/College/Desktop/Dissertation/Restart/38481/38481_Ph_data.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(GSE38481.RAW)= phData_38481$X.Sample_geo_accession
colnames(GSE38481.RAW)


#/ PERFORMING BACKGROUND CORRECTION

#BCOR.GSE38481<- lumiExpresso(GSE38481.RAW, bg.correct = TRUE,
#                            variance.stabilize = FALSE,
#                           normalize =TRUE,
#                          verbose = TRUE)

BCOR.GSE38481 <- lumiN(lumiT(lumiB(
  GSE38481.RAW,method="bgAdjust"),method="log2"),method="quantile")

dim(BCOR.GSE38481)
#16023

# Extract background corrected matrix

BCOR.GSE38481.matrix= exprs(BCOR.GSE38481)
dim(BCOR.GSE38481.matrix)


# Normalization
#library(lumi)
##file.lumi = lumiR(BCOR.GSE38481.matrix)
##lumiExpr = lumiExpresso(file.lumi, bg.correct = TRUE, normalise = TRUE, verbose = TRUE)
##exprs <- exprs(lumiExpr)
##pvalue <- detection(lumiExpr)

# 1. Convert your matrix to a LumiBatch object
# Replace 'your_matrix' with the name of your variable (e.g., BCOR.GSE38481.matrix)
#lumi_batch <- new("LumiBatch", exprs = BCOR.GSE38481.matrix)

# 2. Log2 Transformation
# This ensures data is normally distributed
#lumi_log2 <- lumiT(lumi_batch, method = "log2")

# 3. Quantile Normalization
# This makes the distribution of all samples identical
#lumi_norm <- lumiN(lumi_log2, method = "quantile")

# 4. Extract the final normalized matrix
#final_matrix <- exprs(BCOR.GSE38481.matrix.2)





# Check for negative values and replace them with zero

bcd=apply(BCOR.GSE38481.matrix, 1, function(row) any(row<0))
length(which(bcd))
# 1700 rows with negative values


BCOR.GSE38481.matrix[BCOR.GSE38481.matrix<0]=0

BCOR.GSE38481.df= as.data.frame(BCOR.GSE38481.matrix)

dim(BCOR.GSE38481.df) 
#16023

## GENE ANNOTATION

BCOR.GSE38481.df= rownames_to_column(BCOR.GSE38481.df)

GENE.SYMBOLS<- nuID2targetID(BCOR.GSE38481.df$rowname, lib.mapping = "lumiHumanIDMapping")

GENE.SYMBOLS= as.data.frame(GENE.SYMBOLS)

GENE.SYMBOLS.sub = GENE.SYMBOLS$GENE.SYMBOLS

## CHECKING HGNC SYMBOLS

GENE.SYMBOLS[,2:4]<-checkGeneSymbols(GENE.SYMBOLS.sub, species = "human", unmapped.as.na = TRUE)

GENE.SYMBOL.ASSIGNED=na.omit(GENE.SYMBOLS)

table(is.na(GENE.SYMBOL.ASSIGNED$GENE.SYMBOLS))

dim(GENE.SYMBOL.ASSIGNED) 
#15799


# Remove rows with unassigned gene symbols

GENE.SYMBOL.ASSIGNED= rownames_to_column(GENE.SYMBOL.ASSIGNED, "rowname")


GENES.MERGED= merge(GENE.SYMBOL.ASSIGNED, BCOR.GSE38481.df, by="rowname")

dim(GENES.MERGED)

#/ LIMMA AVEREPS

GSE38481.avreps=as.data.frame((limma::avereps(GENES.MERGED, GENES.MERGED$Suggested.Symbol)))

GSE38481.avreps= GSE38481.avreps[,-c(1:4)]
dim(GSE38481.avreps
)#12647

GSE38481.avreps.final= column_to_rownames(GSE38481.avreps, "Suggested.Symbol")
dim(GSE38481.avreps.final)
# 12647

is.numeric(GSE38481.avreps.final)
GSE38481.avreps.final[]= lapply(GSE38481.avreps.final, function(x) as.numeric(as.character(x)))

saveRDS(GSE38481.avreps.final, file = "38481_normalized.rds")


#38484

#/ LISTING FILES IN THE DIRECTORY

txt.gse38484<-list.files("C:/Users/College/Desktop/DIssertation/Restart/38484/GSE38484/", pattern="txt")

#/ READING RAW DATA
RAW.DATA.GSE38484 <- lumiR("C:/Users/College/Desktop/DIssertation/Restart/38484/GSE38484_non_normalized.txt", convertNuID = TRUE, 
                           lib.mapping = 'lumiHumanIDMapping', QC = TRUE,  dec = '.', parseColumnName = TRUE,
                           columnNameGrepPattern = list(exprs='AVG_SIGNAL', se.exprs='NA', detection='DETECTION'))
dim(RAW.DATA.GSE38484)#48742

#/ FILTERING EXPRESSED PROBES by taking the average of 3 probes<0.05

EXPRESSED.PROBES_38484<- rowSums(RAW.DATA.GSE38484@assayData$detection<0.05)>=3

# Subtract from main object or retain only expressed probes in expression matrix

GSE38484.RAW=RAW.DATA.GSE38484[EXPRESSED.PROBES_38484,]
dim(GSE38484.RAW) #38975

#/ FETCHING PHENOTYPE DATA

phData_38484<- read.csv("C:/Users/College/Desktop/DIssertation/Restart/38484/38484_Phenodata.csv")

colnames(GSE38484.RAW)= phData_38484$Geo_accession
colnames(GSE38484.RAW)

#/ PERFORMING BACKGROUND CORRECTION

#BCOR.GSE38484<- lumiB(GSE38484.RAW, method = 'bgAdjust')


BCOR.GSE38484 <- lumiN(lumiT(lumiB(
  GSE38484.RAW,method="bgAdjust"),method="log2"),method="quantile")

dim(BCOR.GSE38484)#38975

# Extracting background corrected matrix

BCOR.GSE38484.matrix= exprs(BCOR.GSE38484)
dim(BCOR.GSE38484.matrix)

# Check for negative values and replace them with zero

cde=apply(BCOR.GSE38484.matrix, 1, function(row) any(row<0))
length(which(cde)) # 27742 rows with negative values

BCOR.GSE38484.matrix[BCOR.GSE38484.matrix<0]=0

BCOR.GSE38484.df= as.data.frame(BCOR.GSE38484.matrix)

dim(BCOR.GSE38484.df) #38975

#/ GENE ANNOTATION

BCOR.GSE38484.df= rownames_to_column(BCOR.GSE38484.df)

GENE.SYMBOLS<- nuID2targetID(BCOR.GSE38484.df$rowname, lib.mapping = "lumiHumanIDMapping")

GENE.SYMBOLS= as.data.frame(GENE.SYMBOLS)

GENE.SYMBOLS.sub = GENE.SYMBOLS$GENE.SYMBOLS

#/ CHECKING HGNC SYMBOLS

GENE.SYMBOLS[,2:4]<-checkGeneSymbols(GENE.SYMBOLS.sub, species = "human", unmapped.as.na = TRUE)

GENE.SYMBOL.ASSIGNED=na.omit(GENE.SYMBOLS)

table(is.na(GENE.SYMBOL.ASSIGNED$GENE.SYMBOLS))

dim(GENE.SYMBOL.ASSIGNED) #25373

# Remove rows with unassigned gene symbols

GENE.SYMBOL.ASSIGNED= rownames_to_column(GENE.SYMBOL.ASSIGNED, "rowname")

#BCOR.GSE38484.df= rownames_to_column(BCOR.GSE38484.df, "rowname")

GENES.MERGED= merge(GENE.SYMBOL.ASSIGNED, BCOR.GSE38484.df, by="rowname")

dim(GENES.MERGED)#25373

#/ LIMMA AVEREPS

GSE38484.avreps=as.data.frame(limma::avereps(GENES.MERGED, GENES.MERGED$Suggested.Symbol))

GSE38484.avreps= GSE38484.avreps[,-c(1:4)]
dim(GSE38484.avreps)#17233

GSE38484.avreps.final= column_to_rownames(GSE38484.avreps, "Suggested.Symbol")
dim(GSE38484.avreps.final)# 17233

is.numeric(GSE38484.avreps.final)
GSE38484.avreps.final[]= lapply(GSE38484.avreps.final, function(x) as.numeric(as.character(x)))

saveRDS(GSE38484.avreps.final, file = "38484_normalized.rds")

#48072

#/ READING RAW DATA
RAW.GSE48072 <- lumiR("C:/Users/College/Desktop/Dissertation/Restart/48072/Rawdata.txt", convertNuID = TRUE, 
                      lib.mapping = 'lumiHumanIDMapping', QC = TRUE,  dec = '.', parseColumnName = TRUE,
                      columnNameGrepPattern = list(exprs='SAMPLE', se.exprs='NA', detection='Detection'))


#/ FILTERING EXPRESSED PROBES by taking the average of 3 probes<0.05

EXP.PROBES=rowSums(RAW.GSE48072@assayData$detection<0.05)>=3

# Subtract from main object or retain only expressed probes in expression matrix

GSE48072.RAW=RAW.GSE48072[EXP.PROBES,]

dim(GSE48072.RAW)# 30355

y=as.data.frame(GSE48072.RAW)

#/ FETCHING PHENOTYPE DATA
phData<- read.csv("C:/Users/College/Desktop/Dissertation/Restart/48072/48072_Phenodata.csv")

colnames(GSE48072.RAW)=phData$Geo_accession
colnames(GSE48072.RAW)

z=as.data.frame(GSE48072.RAW)
#/ PERFORMING BACKGROUND CORRECTION

BCOR.GSE48072<- lumiB(GSE48072.RAW)
dim(BCOR.GSE48072)#30355

expr <- exprs(BCOR.GSE48072)

# Remove invalid values
expr[expr <= 0] <- NA

# Remove probes with any NA
expr <- expr[complete.cases(expr), ]

# Now log2 is safe
expr <- log2(expr)

# Now normalize
expr <- normalizeBetweenArrays(expr, method = "quantile")

# Final check
sum(is.na(expr))
sum(is.infinite(expr))



##BCOR.GSE48072 <- lumiN(lumiT(lumiB(
#  GSE48072.RAW,method="bgAdjust"),method="log2"),method="quantile")
#dim(BCOR.GSE48072)


#expr <- exprs(GSE48072.RAW)
#expr <- log2(expr + 1)
#expr.qn <- normalizeBetweenArrays(expr, method = "quantile")


# Extracting background corrected matrix

BCOR.GSE48072.matrix= expr
dim(BCOR.GSE48072.matrix)


BCOR.GSE48072.df= as.data.frame(BCOR.GSE48072.matrix)

BCOR.GSE48072.df= rownames_to_column(BCOR.GSE48072.df)

GENE.SYMBOLS<- nuID2targetID(BCOR.GSE48072.df$rowname, lib.mapping = "lumiHumanIDMapping")

GENE.SYMBOLS= as.data.frame(GENE.SYMBOLS)

GENE.SYMBOLS.sub = GENE.SYMBOLS$GENE.SYMBOLS

#/ CHECKING HGNC SYMBOLS

GENE.SYMBOLS[,2:4]<-checkGeneSymbols(GENE.SYMBOLS.sub, species = "human", unmapped.as.na = TRUE)

GENE.SYMBOL.ASSIGNED=na.omit(GENE.SYMBOLS)

table(is.na(GENE.SYMBOL.ASSIGNED$GENE.SYMBOLS))

dim(GENE.SYMBOL.ASSIGNED) #21521

# Remove rows with unassigned gene symbols

GENE.SYMBOL.ASSIGNED= rownames_to_column(GENE.SYMBOL.ASSIGNED, "rowname")

BCOR.GSE48072.df= rownames_to_column(BCOR.GSE48072.df, "rownames")

BCOR.GSE48072.df$rownames <- NULL


GENES.MERGED= merge(GENE.SYMBOL.ASSIGNED, BCOR.GSE48072.df, by="rowname")

dim(GENES.MERGED)#21521

#/ LIMMA AVEREPS

GSE48072.avreps=as.data.frame((limma::avereps(GENES.MERGED, GENES.MERGED$Suggested.Symbol)))

GSE48072.avreps= GSE48072.avreps[,-c(1:4)]
dim(GSE48072.avreps)#15155

GSE48072.avreps.final= column_to_rownames(GSE48072.avreps, "Suggested.Symbol")
dim(GSE48072.avreps.final)# 15155

is.numeric(GSE48072.avreps.final)
GSE48072.avreps.final[]= lapply(GSE48072.avreps.final, function(x) as.numeric(as.character(x)))

GSE48072.avreps.final=as.data.frame(GSE48072.avreps.final)
boxplot(GSE48072.avreps.final)

saveRDS(GSE48072.avreps.final, file = "48072_normalized.rds")


#/ LISTING FILES
TXT.GSE54913<-list.files("C:/Users/College/Desktop/Dissertation/Restart/54913/GSE54913/", pattern="txt", full.names = TRUE)

#/ READING RAW DATA

GSE54913 <- read.delim("C:/Users/College/Desktop/Dissertation/Restart/54913/GSE54913_mRNA_raw_for_GEO.txt", header = TRUE, sep ="\t")


phData<- read.csv("C:/Users/College/Desktop/Dissertation/Restart/54913/54913_Ph_data.csv")

#/ MODIFYING THE RAW DATA TABLE
GSE54913.mod1= GSE54913[,c(1,3)]
GSE54913.mod2= GSE54913[,-c(2:14)]


expr_raw <- GSE54913.mod2[, -1]
expr_raw[] <- lapply(expr_raw, function(x) as.numeric(as.character(x)))

#/ CHECKING HGNC SYMBOLS

GSE54913.mod1[,3:5]<-checkGeneSymbols(GSE54913.mod1$GeneSymbol, species = "human", unmapped.as.na = TRUE)

table(is.na(GSE54913.mod1$Suggested.Symbol))

# Remove rows with unassigned gene symbols

GSE54913.mod1.df=na.omit(GSE54913.mod1)
table(is.na(GSE54913.mod1.df$Suggested.Symbol))

# Keep only PROBE NAMES and suggested symbols

GSE54913.mod1.df1=GSE54913.mod1.df[,c(1,5)]

# Add HGNC symbols to normalized expression matrix

GSE54913.HGNC= merge(GSE54913.mod2, GSE54913.mod1.df1, by="ProbeName")

# Shifting last column to second


GSE54913.HGNC <- GSE54913.HGNC %>% relocate(Suggested.Symbol, .before = NC72)


GSE54913 <- read.delim("C:/Users/College/Documents/GSE54913_mRNA_raw_for_GEO.txt", header = TRUE, sep ="\t")


# Remove Probe Names

GSE54913.FINAL=GSE54913.HGNC[,-c(1)]


#/ BACKGROUND CORRECTION

GSE54913.FINAL.BC= limma::backgroundCorrect(GSE54913.FINAL[,c(2:31)], method = "normexp")


expr <- GSE54913.FINAL.BC
expr <- log2(expr)
expr <- normalizeBetweenArrays(expr, method = "quantile")


# Merge again with Gene symbols

GSE54913.FINAL.BC=expr

GSE54913.GS= cbind(GSE54913.FINAL.BC, GSE54913.FINAL[,1])
#GSE54913.GS <- GSE54913.GS %>% relocate(V31, .before = NC72)
GSE54913.GS.df=as.data.frame(GSE54913.GS)
GSE54913.GS.df <- GSE54913.GS.df %>% relocate(V31, .before = NC72)
colnames(GSE54913.GS.df)[colnames(GSE54913.GS.df) == "V31"] <- "Suggested.Symbol"

dim(GSE54913.GS.df)# 17038 31

#/ LIMMA AVEREPS
GSE54913.avereps= as.data.frame(limma::avereps(GSE54913.GS.df, GSE54913.GS.df$Suggested.Symbol))

dim(GSE54913.avereps)# 13023 31

#/ HAVE UNIQUE GENE SYMBOLS AS ROW IDS

GSE54913.avereps.df=column_to_rownames(GSE54913.avereps, "Suggested.Symbol")
colnames(GSE54913.avereps.df)= phData$Geo_accession
colnames(GSE54913.avereps.df)

dim(GSE54913.avereps.df)# 13023 30

#/ EXPORT DATA

class(GSE54913.avereps.df)

#is.numeric(GSE54913.avreps.final)
#GSE54913.avreps.final[]= lapply(GSE54913.avreps.df, function(x) as.numeric(as.character(x)))


df_num <- GSE54913.avereps.df
df_num[] <- lapply(df_num, function(x) as.numeric(as.character(x)))


saveRDS(df_num, file="54913_normalized.rds")

# Kum

Kum_WN <- readRDS("C:/Users/College/Desktop/Dissertation/Restart/Kum/Kum_WN.rds")

class(Kum_WN)

is.numeric(Kum_WN)

expr <- as.matrix(Kum_WN)

is.numeric(expr)

expr.log2 <- log2(expr + 1)

expr.norm <- normalizeBetweenArrays(expr.log2, method = "quantile")


#boxplot()#for com=nfirmation

expr.norm.df <- as.data.frame(expr.norm)

saveRDS(expr.norm, file = "Kum_normalized.rds")

