``` R

### TCGA Data Downloads ##

## Installation of BioCmanager and TCGAbiolinks ##
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

### Importing Libraries ### 

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

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
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)

# Methylation Data Preparation

## Data pre-processing and normalization ##
data <- read.csv("mRNA.csv", header  = TRUE, sep = ",", row.names = 1)
head(data, 5)
dim(data)

## Removing Rows containing NA
data <- data %>% drop_na()

## Removing Rows containing zeros
data <- data[apply(data, 1, function(row) all(row !=0 )), ]








```


