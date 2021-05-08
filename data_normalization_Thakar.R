# install packages needed for programmer part (Project-5 redo by Abhishek)
#install.packages("BiocManager")
#BiocManager::install("affy")
#BiocManager::install("affyPLM")
#BiocManager::install("sva")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("hgu133plus2.db")

###### set working directory to where the .CEL files are 
setwd("/projectnb2/bf528/users/dachshund/project_1/samples/all_samples")
#To see working directory
# getwd()

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)

###### set working directory to where the .CEL files are 
#setwd("/projectnb/bf528/users/dachshund/project_1/samples/all_samples")
# getwd()


####RMA expression ---

# read .CEL files in current directory and normalize them
data <- ReadAffy()
normalized_data <- rma(data)

#### fitPLM ---
#Computing relative log expression and normalized unscaled standard error for microarray samples
Pset <- fitPLM(data, normalize=TRUE, background=TRUE)
#pset)


#### RLE ---
# compute the relative log expression (RLE) of the microarray samples

# RLE for each sample
RLE_scores <- RLE(Pset, type="stats")
NUSE_scores <- NUSE(Pset,type="stats")

#Create histograms for RLE and NUSE stats
hist(RLE_scores['median',],xlab="Median RLE", main="Distribution of RLE scores")
hist(NUSE_scores['median',],xlab="Median NUSE", main="Distribution of NUSE scores")
#Using ComBat to correct for batch effects. Only when the batch covariate is known
matrix <- exprs(normalized_data)
dim(matrix)

# path to annotation file on SCC 
meta_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")

# Batch effects include both Center and RNA extraction method and have been merged
# into a single variable called *normalizationcombatbatch* in the annotation file
batch = meta_data$normalizationcombatbatch
mod = model.matrix(~normalizationcombatmod, data=meta_data)
dim(batch)

#corrected batch effects
batch_corrected <- ComBat(dat=matrix,batch=batch, mod=mod)
write.csv(batch_corrected, file='expressiondata.csv')

## Principal component analysis (PCA) for corrected batch and normalized data----
transpose_batch_corrected = t(batch_corrected)
rescaled_transposed <- scale(transpose_batch_corrected)

# return the scaled transposed data back to its original orientation
untransposed_scaled <- t(rescaled_transposed)

PCA_output <- prcomp(untransposed_scaled,scale.= FALSE, center=FALSE)
PC1 <- PCA_output$rotation[,1]
PC2 <- PCA_output$rotation[,2]

summarized <- summary(PCA_output)
summarized$importance
plot(PC2, PC1, 
     xlab = paste0("PC2 (", round(summarized$importance[2,2]*100, digits = 2), "%)"),
     ylab = paste0("PC1 (", round(summarized$importance[2,1]*100, digits = 2), "%)"))
boxplot(PC1)
which(PC1 > mean(PC1) +3*sd(PC1) | PC1 < mean(PC2)-3*sd(PC1))

boxplot(PC2)
which(PC2 > mean(PC2) +3*sd(PC2) | PC2 < mean(PC2)-3*sd(PC2))

id <- which(!(PC2 > mean(PC2) +3*sd(PC2) | PC2 < mean(PC2) - 3*sd(PC2) | PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1) - 3*sd(PC1)))

#outliers PCA
PCA_no_outliers <- prcomp(untransposed_scaled[,id],scale.=FALSE,center=FALSE)
PC1_no_outliers <- PCA_no_outliers$rotation[,1]
PC2_no_outliers <- PCA_no_outliers$rotation[,2]
summarized_no_outliers <- summary(PCA_no_outliers)
summarized_no_outliers$importance
plot(PC2_no_outliers, PC1_no_outliers, 
     xlab = paste0("PC2-no outliers (", round(summarized_no_outliers$importance[2,2]*100, digits = 2), "%)"),
     ylab = paste0("PC1-no outliers (", round(summarized_no_outliers$importance[2,1]*100, digits = 2), "%)"))

final_data <- (batch_corrected[,id])
write.csv(final_data, file="final_data.csv")


