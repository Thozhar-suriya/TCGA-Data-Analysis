``` R

### TCGA Data Downloads ##

## Installation of BioCmanager and TCGAbiolinks ##
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("TCGAbiolinks")

## Updated Version ##
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

##if not use devtools or remotes##

### Importing Libraries ### 
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(tidyr)

## Project ##
## TCGA-STAD ##

### Query for mRNA ##

## The following options are used to search mRNA results using TCGAbiolinks:
## Harmonized database: data aligned against hg38 Reference Genome
## data.category: "Transcriptome Profiling"
## data.type: "Gene Expression Quantification"
## workflow.type = "STAR - Counts"

query.exp.hg38 <- GDCquery(
    project = "TCGA-STAD", 
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts",
    )
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
    query = query.exp.hg38,
    save = TRUE, 
    save.filename = "exp.rda"
)

### Query for miRNA ##

library(TCGAbiolinks)
query.mirna <- GDCquery(
    project = "TCGA-STAD", 
    experimental.strategy = "miRNA-Seq",
    data.category = "Transcriptome Profiling", 
    data.type = "miRNA Expression Quantification"
)
GDCdownload(query.mirna)
mirna <- GDCprepare(
    query = query.mirna,
    save = TRUE, 
    save.filename = "mirna.rda"
)

### Query for DNA methylation ##

# HumanMethylation450
query_met.hg38 <- GDCquery(
    project= "TCGA-STAD", 
    data.category = "DNA Methylation", 
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450", 
)
GDCdownload(query_met.hg38, 
                method = "api",
                directory = "./TCGA")
        
                       
 data.hg38 <- GDCprepare(query_met.hg38, 
                                save = FALSE, 
                                    summarizedExperiment = TRUE, 
                                                directory = "./TCGA",)

# Methylation Data Preparation
data.hg38
data.hg38 %>% rowRanges %>% as.data.frame %>% head
data.hg38 %>% colData %>% as.data.frame
data.hg38 %>% assay %>% head %>% as.data.frame
data.hg38 %>% rowData %>% as.data.frame %>% head

#Saving the Object
saveRDS(object = data.hg38,
        file = "data.hg38.RDS",
        compress = FALSE)

# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.hg38.RDS")

# Making assay Matrix
met <- as.data.frame(assay(data.hg38))

# Making Clincal Matrix
clinical <- data.frame(data.hg38@colData)

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
probe.na <- rowSums(is.na(met))
table(probe.na == 0)

# chose those has not NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

#####################################
## If clinical data rownmaes is duplicate ##
clin <- read.csv("Methyl_rowRanges.csv", sep = ",")
rownames(clin) <- NULL
rownames(clin) <- make.names(clin$X, unique = TRUE)
################################

## remove probes that match to chromosome  X and Y 
clin <- subset(clinical,subset = !as.character(seqnqmes(clinical)) %in% c("chrM","chrX","chrY", "*"))
Extract <- rownames(clin)

head(Extract)
Methyl_without_X_Y <- met[Extract,] %>% data.frame()
dim(Methyl_without_X_Y)

##drop NA##
Methyl_without_X_Y <- Methyl_without_X_Y %>% drop_na()
write.csv(Methyl_without_X_Y, "Methyl_without_X_Y.csv")


# mRNA Pre-processing ##  
## Data pre-processing and normalization ##
data <- read.csv("mRNA.csv", header  = TRUE, sep = ",", row.names = 1)
head(data, 5)
dim(data)

## Removing Rows containing NA
data <- data %>% drop_na()

## Removing Rows containing zeros

data <- data[apply(data, 1, function(row) all(row !=0 )), ]
write.csv(data, "mRNA_without_NA_Zero.csv")

## Converting to Gene Symbol ##

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- read.csv("mRNA_without_NA_Zero.csv")
head(df)
genes <- df$ensembl_gene_id
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=genes, mart= mart)
G_list
A <- merge(df, G_list, by = "ensembl_gene_id")
write.csv(A, "mRNA_withsymbol.csv")

# miRNA Pre-processing ##  
## Data pre-processing and normalization ##
data <- read.csv("miRNA.csv", header  = TRUE, sep = ",", row.names = 1)
head(data, 5)
dim(data)

## Removing Rows containing NA
data <- data %>% drop_na()

## Removing Rows containing all zeros ##
data <- data[rowSums(data[])>0,]
write.csv(data, "mRNA_without_NA_Zero.csv")

### Integration of the Omics based on Patient Data ##



```


