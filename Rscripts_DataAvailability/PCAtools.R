library(dplyr)
library(tidyverse)
library(PCAtools)
library(cowplot)
library(ggplotify)

############ Example 1 ###########################################################
pca_quant <- read_delim("C:/Users/Laura/Downloads/astrocitos_COVID19/figs_Age_dependVirus/Diff Analysis/pca_matrix.txt", 
                        delim = "\t", trim_ws = TRUE)

# Convertir Protein.group-->rownames y el dataframe-->matrix
pca_quant <- pca_quant %>% remove_rownames %>% column_to_rownames(var="Protein.Group")
pca_quant$First.Protein.Description <- NULL
pca_quant$Genes <- NULL
pca_quant$Protein.Names <- NULL
pca_quant$Protein.Ids <- NULL
pca_quant$KEGG <- NULL
pca_quant$GOCC <- NULL
pca_quant$GOBP <- NULL
pca_quant$GOMF <- NULL
pca_quant <- as.matrix(pca_quant)

# Metadata
md <- read_delim("C:/Users/Laura/Downloads/astrocitos_COVID19/figs_Age_dependVirus/Diff Analysis/pca_metadata.txt", 
                        delim = "\t", trim_ws = TRUE)
md <- md %>% remove_rownames %>% column_to_rownames(var="...1")

# Crear el pca
pca1 <- pca(pca_quant, metadata = md, removeVar = 0.1)
getComponents(pca1)
getVars(pca1)
getLoadings(pca1)

# Visualizing results
screeplot(pca1, hline=80, axisLabSize=14, titleLabSize=14)

biplot(pca1, labSize=3, pointSize=3, sizeLoadingsNames=5, axisLabSize=14, titleLabSize=14)

biplot(pca1, 
       legendPosition='right', colby='Age', colkey = c('cent'='blue', 'adu'='red'),
       shape='Virus', shapekey=c('CTRL'=3, 'GAM'=0, 'WU'=1),
       #hline = 0, vline = c(-10, 0, 10),
       #vlineType = c('dashed', 'solid', 'dashed'),
       labSize=3, pointSize=4, sizeLoadingsNames=3, axisLabSize=14, titleLabSize=14)

biplot(pca1, showLoadings = TRUE, lab = NULL)

biplot(pca1, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
       legendPosition = 'right')

biplot(pca1, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
       shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')

pairsplot(pca1, triangle = TRUE)

plotloadings(pca1, drawConnectors=TRUE)

eigencorplot(p, components = getComponents(p, 1:10),
             metavars = c('ESR', 'CRP'))


############ Astrocytes leva 1-2 ###########################################################
dataMatrix <- read.delim("C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/dataMatrix.txt")

# Convertir Protein.group-->rownames y el dataframe-->matrix
dataMatrix <- dataMatrix %>% remove_rownames %>% column_to_rownames(var="Protein.Group")
dataMatrix <- as.matrix(dataMatrix)

# Metadata
sampleMetadata <- read.delim("C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/sampleMetadata.txt")
sampleMetadata <- sampleMetadata %>% remove_rownames %>% column_to_rownames(var="Sample")

# Crear el pca
pca1 <- pca(dataMatrix, metadata = sampleMetadata, removeVar = 0.1)
saveRDS(pca1, "C:/Users/Laura/Downloads/astrocitos_COVID19/_leva1_2_reanalysis/pca_PCAtools.RDS")
getComponents(pca1)
getVars(pca1)
getLoadings(pca1)

# Visualizing results
a <- screeplot(pca1, hline=80, axisLabSize=14, titleLabSize=14, title = 'SCREE plot')

biplot(pca1, labSize=3, pointSize=3, sizeLoadingsNames=5, axisLabSize=14, titleLabSize=14)

c <- biplot(pca1, 
       legendPosition='right', colby='Age', colkey = c('cent'='blue', 'adu'='red'),
       shape='Virus', shapekey=c('CTRL'=3, 'GAM'=0, 'WU'=1, 'pool'=8),
       #hline = 0, vline = c(-10, 0, 10),
       #vlineType = c('dashed', 'solid', 'dashed'),
       labSize=3, pointSize=3, sizeLoadingsNames=3, axisLabSize=14, titleLabSize=14, title = 'PCA bi-plot')

biplot(pca1, showLoadings = TRUE, lab = NULL)

biplot(pca1, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
       legendPosition = 'right')

biplot(pca1, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
       shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')

pairsplot(pca1, triangle = TRUE)
b <- pairsplot(pca1, components = getComponents(pca1, c(1:3)),
          triangle = TRUE, trianglelabSize = 12, titleLabSize = 12,
          hline = 0, vline = 0,
          pointSize = 3, gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Infection',
          title = 'Pairs plot by infection', plotaxes = FALSE,
          margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
          returnPlot = FALSE)

plotloadings(pca1, drawConnectors=TRUE)

eigencorplot(pca1, components = getComponents(p, 1:10),
             metavars = c('ESR', 'CRP'))

d <- plotloadings(pca1, rangeRetain = 0.01, labSize = 3,
                          title = 'Loadings plot', axisLabSize = 12,
                          subtitle = 'PC1, PC2, PC3, PC4, PC5',
                          caption = 'Top 1% variables',
                          shape = 23, shapeSizeRange = c(2, 5),
                          col = c('limegreen', 'black', 'red3'),
                          legendPosition = 'none',
                          drawConnectors = FALSE,
                          returnPlot = FALSE)

plot_grid(a, b, c, d)
