### TCGA Data Downloads ##

```R
## Installation of BioCmanager and TCGAbiolinks ##
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

### Importing Libraries ### 

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)




```
