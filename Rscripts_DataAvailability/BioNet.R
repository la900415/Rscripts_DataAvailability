library(dplyr)
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)

############ Bionet example ######################################################
# load the BioNet package and the required data sets, containing a human protein-protein interaction network and p-values derived from DEA
# aggregate these two p-values into one p-value.
pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)

# A subnetwork of the complete network is derived, containing all the proteins which are represented by probesets on the microarray. And self-loops are removed.
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet

# score each node of the network and FDR is set to 0.001
fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)

# use a fast heuristic approach to calculate an approximation to the optimal scoring subnetwork
module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label

# Both 2D and 3D module visualization procedures are available in BioNet. For a 3D visualization, see section 3.4. 
# Alternatively, the network could be easily exported in Cytoscape format
plotModule(module, scores = scores, diff.expr = logFC)
saveNetwork(module, file = "C:/Users/Laura/Downloads/ALL_module", type = "XGMML")

library("rgl")
plot3dModule(module, scores = scores, diff.expr = logFC)


####### Astrocytes ###########################
dea2 <- readRDS("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/clusterprofiler/dea_signif_standarderror.RDS")
GcCc <- dea2 %>% filter(contrast == "Gc/Cc" & signif=="TRUE") %>% arrange(desc(log2fc))
keytypes(org.Hs.eg.db)
ids2 <- bitr(GcCc$protein_id, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ids <- bitr(GcCc$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ids$SYMBOL2 <- ids$SYMBOL
GcCc2 <- merge(GcCc, ids, by.x="symbol", by.y="SYMBOL2", all.x=TRUE, all.y=FALSE) %>% arrange(desc(log2fc))
write_csv(GcCc2, "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/clusterprofiler/GcCc2.csv")
GcCc2 <- read_excel("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/clusterprofiler/GcCc2.xlsx", 
                    col_types = c("text", "text", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text", "text", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text"))
GcCc2$symbol_entrezid <-  paste0(GcCc2$SYMBOL, "(", GcCc2$ENTREZID, ")")
saveRDS(GcCc2, "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/clusterprofiler/GcCc2.RDS")

GcCc2_pval <- GcCc2 %>% select(SYMBOL, ENTREZID, pvalue_adjust)
GcCc2_pval$symbol_entrezid <-  paste0(GcCc2_pval$SYMBOL, "(", GcCc2_pval$ENTREZID, ")")
GcCc2_pval <- GcCc2_pval %>% select(symbol_entrezid, pvalue_adjust)
GcCc2_pval <- GcCc2_pval %>% column_to_rownames(var="symbol_entrezid")
pval2 <- cbind(x=GcCc2_pval$pvalue_adjust)
row.names(pval2) <- GcCc2$symbol_entrezid
pval <- aggrPvals(pval2, order = 1, plot = FALSE)

subnet <- subNetwork(GcCc2$symbol_entrezid, interactome)
subnet <- rmSelfLoops(subnet)
subnet

fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)

module <- runFastHeinz(subnet, scores)
logFC <- GcCc2$log2fc
names(logFC) <- GcCc2$symbol_entrezid

plotModule(module, scores = scores, diff.expr = logFC)
saveNetwork(module, file = "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/cytoscape/GcCc_bionet", type = "XGMML")

