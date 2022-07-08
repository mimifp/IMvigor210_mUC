############ Metastasic urothelial cancer RNAseq EGAS00001002556 ############
# This worfkflow was performed in R 3.3.3

setwd("/Users/mimiferreiro/Documents/GitHub/tfm_mUC/2_data/1_IMvigor210")

#### 1.  PREPARE ENVIROMENT ####
# Downloading urothelial samples package (in terminal)
#download_path <- ("../../2_data/1_IMvigor210")
#setwd(download_path)
#wget http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.0.tar.gz

# Load functions
source("http://bioconductor.org/biocLite.R")
source("../1_functions/remove_duplicated_genes.R")

install.packages("IMvigor210CoreBiologies_1.0.0.tar.gz", 
                 repos=NULL)

# Load packages
library(IMvigor210CoreBiologies)
library(dplyr)

#### 2. EXPLORE AND PREPARE DATA ####

# Load CountDataSet
data(cds)

# Check CountDataSet
head(counts(cds))
head(fData(cds))
fData <- fData(cds)
head(pData(cds))
pData <- pData(cds)

# Load NChannelSet
data(fmone)
ls.str(assayData(fmone))
head(pData(fmone))
head(assayDataElement(fmone, "known_short"))

# Convert factor variables to numeric
exmat <- counts(cds)
entrez_id <- rownames(exmat) 
exmat_trans <- as.data.frame(exmat)
exmat_trans <- lapply(exmat_trans, function(x) as.numeric(as.character(x)))
exmat <- as.data.frame(exmat_trans)
sapply(exmat, class)

# Set entez_id as rownames again
rownames(exmat) <- entrez_id
sample_id <- colnames(counts(cds))

# Make correspondence between ENTREZ ID - SYMBOL ID
fData <- fData %>% filter(!is.na(symbol) | (!is.na(entrez_id)))
counts_annot <- cbind(exmat, fData$symbol[match(rownames(exmat), fData$entrez_id)])
counts_annot <- as.data.frame(counts_annot)

# Change symbol column name and modifying genes names
counts_annot$gene_id <- counts_annot[,ncol(counts_annot)]
counts_annot$gene_id <- gsub("\\-", "_", counts_annot$gene_id)
counts_annot$gene_id <- gsub("\\.", "_", counts_annot$gene_id)
exmat <- counts_annot %>% dplyr::select(gene_id, dplyr::everything())
exmat <- exmat[, -ncol(exmat)]

# Remove empty values in symbol and gene_id column
which((exmat$gene_id == "") == TRUE)
which((fData$symbol == "") == TRUE)
exmat <- exmat %>% filter(!gene_id=="") # 31086 genes
fData <- fData %>% filter(!symbol=="")

# Checking if there are duplicates
id_dup <- which(duplicated(exmat$gene_id) == TRUE) # 1 gene (CSNK1E)
length(id_dup)

# Remove duplicated genes
exmat <- remove_duplicated_genes(exmat)
fData <- fData[-id_dup,]

# Count matrix (exmat) contains 348 samples and information for 31085 genes. Counts = number of reads
# that align with a sample for each gene. 
dim(exmat)

# Save datasets
# Change working directory
write.csv(exmat, "exmat_IMvigor210.csv", row.names = FALSE)
write.csv(fData, "fData_IMvigor210.csv", row.names = FALSE)
write.csv(pData, "pData_IMvigor210.csv", row.names = FALSE)