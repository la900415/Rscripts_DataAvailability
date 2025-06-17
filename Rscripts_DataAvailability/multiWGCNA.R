library(multiWGCNA)
library(multiWGCNAdata)
library(ggplot2)
library(ggalluvial)
library(ExperimentHub)
library(generics)
library(BiocGenerics)
library(AnnotationHub)
library(BiocFileCache)
library(dbplyr)
library(SummarizedExperiment)


####### 3 Load astrocyte Ribotag RNA-seq data ###################
eh = ExperimentHub()

eh_query = query(eh, c("multiWGCNAdata"))
astrocyte_se = eh_query[["EH8223"]]

# Collect the metadata in the sampleTable; the first column must be named "Sample"
sampleTable = colData(astrocyte_se)
View(sampleTable)

# Check the data
assays(astrocyte_se)[[1]][1:5, 1:5]
View(assays(astrocyte_se)[[1]])

# Define our conditions for trait 1 (disease) and 2 (brain region)
conditions1 = unique(sampleTable[,2])
conditions2 = unique(sampleTable[,3])

view(conditions1)
view(conditions2)

####### 4 Construct the combined networks and all the sub-networks (EAE, WT, and each region) ###########
# Same parameters as Tommasini and Fogel. BMC Bioinformatics
# We now perform network construction, module eigengene calculation, module-trait correlation. Let’s use power = 12 since we used this in our paper (Tommasini and Fogel. BMC Bioinformatics. 2023.) for all the networks.
astrocyte_networks = constructNetworks(astrocyte_se, sampleTable, conditions1, conditions2, 
                                       networkType = "signed", TOMType = "unsigned", 
                                       power = 12, minModuleSize = 100, maxBlockSize = 25000,
                                       reassignThreshold = 0, minKMEtoStay = 0, mergeCutHeight = 0,
                                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                                       deepSplit = 4, verbose = 3)

# Load pre-computed astrocyte networks
astrocyte_networks = eh_query[["EH8222"]] 
saveRDS(astrocyte_networks, file = "astrocyte_networks.RDS")

# Check one of the WGCNA objects
View(astrocyte_networks[["combined"]])
astrocyte_networks[["combined"]]

####### 5 Compare modules by overlap ########
# Next, we compare modules (by hypergeometric overlap) across conditions. We’ll save the results in a list.
# Save results to a list
results = list()
results$overlaps = iterate(astrocyte_networks, overlapComparisons, plot=F) #plot=T for visualizing
saveRDS(results, file = "results.RDS")

# Check the overlaps, ie between the EAE and wildtype networks
head(results$overlaps$EAE_vs_WT$overlap)
View(results$overlaps$EAE_vs_WT$overlap)

######## 6 Identify a module of interest ###########
# Then, we perform differential module expression analysis to detect modules with disease-associated expression patterns. This incorporates the linear model described in the paper and tests for significance using ANOVA.

# Run differential module expression analysis (DME) on combined networks
results$diffModExp = runDME(astrocyte_networks[["combined"]], 
                            sampleTable,
                            p.adjust = "fdr", 
                            refCondition = "Region", 
                            testCondition = "Disease") 
# plot=TRUE, 
# out="ANOVA_DME.pdf")
saveRDS(results, file = "results.RDS")

# Check results sorted by disease association FDR
results$diffModExp[order(results$diffModExp$Disease),]

# You can check the expression of module M13 from Tommasini and Fogel. BMC Bioinformatics. 2023 like this. Note that the values reported in the bottom panel title are p-values and not adjusted for multiple comparisons like in results$diffModExp
diffModuleExpression(astrocyte_networks[["combined"]], 
                     geneList = topNGenes(astrocyte_networks[[1]], "combined_013"), 
                     design = sampleTable,
                     test = "ANOVA",
                     plotTitle = "Module 013",
                     plot = TRUE)


###### 7 Draw the multiWGCNA network ##########
# We can now check to see if M13 is present in any of the sub-networks. An easy way to do this is using the network-network correspondences 
# from hypergeometric overlap. These are stored in results$overlaps. We can plot these in a convenient visualization scheme that also organizes
# the three levels of the multiWGCNA analysis: 1) combined network, 2) EAE and wildtype networks, and 3) the four regional networks.
drawMultiWGCNAnetwork(astrocyte_networks, 
                      results$overlaps, 
                      "combined_013", 
                      design = sampleTable, 
                      overlapCutoff = 0, 
                      padjCutoff = 1, 
                      removeOutliers = TRUE, 
                      alpha = 1e-50, 
                      layout = NULL, 
                      hjust = 0.4, 
                      vjust = 0.3, 
                      width = 0.5)

# We can identify the EAE module that corresponds to M13 using the overlap analysis:
bidirectionalBestMatches(results$overlaps$combined_vs_EAE)

###### 8 Observe differential co-expression of top module genes ########
# We can visually check that combined_013/EAE_015 genes are co-expressed in EAE and not co-expressed in WT samples.

# Get expression data for top 20 genes in EAE_015 module
datExpr = GetDatExpr(astrocyte_networks[[1]], 
                     genes = topNGenes(astrocyte_networks$EAE, "EAE_015", 20))

# Plot
coexpressionLineGraph(datExpr, splitBy = 1.5, fontSize = 2.5) + 
  geom_vline(xintercept = 20.5, linetype='dashed')

######### 9 Follow up with a preservation analysis ##########
# To enable multi-threading
# library(doParallel)
# library(WGCNA)
# nCores = 2
# registerDoParallel(cores = nCores)
# enableWGCNAThreads(nThreads = nCores)

# Turn off multi-threading
# registerDoSEQ()
# disableWGCNAThreads()

# Calculate preservation statistics
results$preservation=iterate(astrocyte_networks[c("EAE", "WT")], 
                             preservationComparisons, 
                             write=FALSE, 
                             plot=TRUE, 
                             nPermutations=2)


######## 10 Determining if preservation value is significant ########
# Then, one can perform a permutation procedure that estimates the probability of observing a 
# disease (or wildtype) module with this preservation score in the wildtype (or disease) 
# setting (PreservationPermutationTest). The test is designed to control for the other condition in 
# the sampleTable. In this case, it will equally distribute the samples belonging to each anatomical 
# region when testing preservation of this disease module in the wildtype samples. 
# This is the slowest step! We recommend to let this run on a computing cluster overnight.
options(paged.print = FALSE)
options(paged.print = FALSE)
results$permutation.test = PreservationPermutationTest(astrocyte_networks$combined@datExpr[sample(17000,3000),], 
                                                       sampleTable, 
                                                       constructNetworksIn = "EAE", # Construct networks using EAE samples
                                                       testPreservationIn = "WT", # Test preservation of disease samples in WT samples
                                                       nPermutations = 10, # Number of permutations for permutation test
                                                       nPresPermutations = 10, # Number of permutations for modulePreservation function
                                                       
                                                       # WGCNA parameters for re-sampled networks (should be the same as used for network construction)
                                                       networkType = "signed", TOMType = "unsigned", 
                                                       power = 12, minModuleSize = 100, maxBlockSize = 25000,
                                                       reassignThreshold = 0, minKMEtoStay = 0, mergeCutHeight = 0,
                                                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                                                       deepSplit = 4, verbose = 3
)

# Load pre-computed results
data(permutationTestResults) 

# Remove outlier modules
permutationTestResultsFiltered = lapply(permutationTestResults, function(x) x[!x$is.outlier.module,])

# Extract the preservation score distribution
results$scores.summary = PreservationScoreDistribution(permutationTestResultsFiltered, 
                                                       moduleOfInterestSize = 303 # The size of the module of interest (dM15)
)

# Observed preservation score of dM15
observed.score = 9.16261490617938

# How many times did we observe a score lower than or equal to this observed score?
z.summary.dist = results$scores.summary$z.summary
below=length(z.summary.dist[z.summary.dist <= observed.score])
probability= below/100
message("Probability of observing a score of ", round(observed.score, 2), " is ", probability)
#> Probability of observing a score of 9.16 is 0.01
saveRDS(results, file = "results.RDS")

# Plot on a histogram
ggplot(results$scores.summary, aes(x=z.summary)) + 
  geom_histogram(color="black", fill="white", bins = 15)+
  xlab("Preservation score")+
  ylab("Frequency")+
  geom_vline(xintercept=observed.score, color="red3", linetype="solid")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))







