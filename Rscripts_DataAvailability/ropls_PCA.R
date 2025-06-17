library(devtools)
library(tidyverse)
library(ropls)
data(sacurine)
names(sacurine)
attach(sacurine)
view(dataMatrix)
view(sampleMetadata)

##################### Astrocytes leva1-2 ##########################################
dataMatrix <- read.delim("C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/dataMatrix.txt")

# Transpose rows and columns--------------------------------
dataMatrix <- rownames_to_column(dataMatrix, var = "Protein.Group") %>% as_tibble()
dataMatrix2 <- data.frame(t(dataMatrix[-1]))
colnames(dataMatrix2) <- dataMatrix[, 1]
dataMatrix <- as.matrix(dataMatrix2)

# Metadata
sampleMetadata <- read.delim("C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/sampleMetadata.txt")
sampleMetadata <- sampleMetadata %>% remove_rownames %>% column_to_rownames(var="Sample")

# crear lista 
astro <- list(dataMatrix=dataMatrix, sampleMetadata=sampleMetadata)
saveRDS(astro, "C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/astro_ropls.RDS")

names(astro)
attach(astro)
view(dataMatrix)
view(sampleMetadata)

# PCA
astro.pca <- opls(dataMatrix)

group <- sampleMetadata[, "Group"]
virus <- sampleMetadata[, "Virus"]
age <- sampleMetadata[, "Age"]
batch <- sampleMetadata[, "Batch"]
rep_biol <- sampleMetadata[, "Rep_biol"]
infection <- sampleMetadata[, "Infection"]
plot(astro.pca,
     typeVc = "x-score",
     parAsColFcVn = infection)


# Partial least-squares: PLS and PLS-DA 
astro.plsda.inf <- opls(dataMatrix, infection)
# our model with 1 predictive (pre) components has significant and quite high R2Y and Q2Y values.

# Orthogonal partial least squares: OPLS and OPLS-DA
astro.oplsda.inf <- opls(dataMatrix, infection,
                        predI = 1, orthoI = NA)

# assess the predictive performance of our model. We first train the model on a subset of the samples
astro.oplsda.train <- opls(dataMatrix, infection,
                        predI = 1, orthoI = NA,
                        subset = "odd")
