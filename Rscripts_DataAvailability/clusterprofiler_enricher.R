library(org.Hs.eg.db)
library(clusterProfiler) 
library(GOplot)
library(DOSE) 
library(enrichplot)
library(PerseusR)
library(msdap)
library(goat)
library(AnnotationDbi)
library(GO.db)
library(msigdbr)
#BiocManager::install("pathview")
library(pathview)
library(ReactomePA)

library(readxl)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(eulerr)
library(readxl) 
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)
library(ggplot2) 
library(ComplexUpset)
library(ggrepel)
library(ggfortify)

############# MayanLab COVID19 GeneSets gmt file ##############################
url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=COVID-19_Related_Gene_Sets_2021"
download.file(url, destfile = "COVID19_GeneSets_2021.gmt")
covid19_gs2021 <- read.gmt("COVID19_GeneSets_2021.gmt")

############# MayanLab ClinVar GeneSets gmt file ##############################
url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ClinVar_2019"
download.file(url, destfile = "ClinVar_2019.gmt")
clinvar_2019 <- read.gmt("ClinVar_2019.gmt")

################ WikiPathways ###################################################
## downloaded from https://wikipathways-data.wmcloud.org/current/gmt/
url <- "https://data.wikipathways.org/current/gmt/"
download.file(url, destfile = "wikipathways-20240910-gmt-Homo_sapiens.gmt")
wp <- read.gmt.wp("wikipathways-20240910-gmt-Homo_sapiens.gmt")
ewp <- GSEA(geneList, TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")])


############# Cell marker file ##############################
cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download/Cell_marker_Human.xlsx')
## instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
  dplyr::select(cellName, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
  tidyr::unnest()

############# Molecular Signatures Database (MSigDB) ##############################
# H: hallmark gene sets, C1: positional gene sets, C2: curated gene sets, 
# C3: motif gene sets, C4: computational gene sets, C5: GO gene sets, C6: oncogenic signatures, 
# C7: immunologic signatures
library(msigdbr)
msigdbr_show_species()
m_t2g <- msigdbr(species = "Homo sapiens", category = "C8" ) %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g) 
#ORA
em <- enricher(gene, TERM2GENE=m_t2g)
head(em)

#GSEA
em2 <- GSEA(geneList, TERM2GENE = C3_t2g)
head(em2)

#########################################################
GcCc2 <- GcCc2 %>% arrange(desc(log2fc))
GcCc2_list <- cbind(GcCc2$log2fc)
rownames(GcCc2_list) <- GcCc2$symbol

cov_gsea_GcCc <- GSEA(GcCc2_list, TERM2GENE = covid19_gs2021)
cov_ora_GcCc <- enricher(GcCc2$symbol, TERM2GENE=covid19_gs2021)
heatplot(cov_ora_GcCc, showCategory=5)
cnetplot(cov_ora_GcCc, foldChange=GcCc2_list)

library(enrichplot)
GcCc2_list <- cbind(GcCc2$log2fc)
rownames(GcCc2_list) <- GcCc2$ENTREZID
de <-  GcCc2$ENTREZID
edo_GcCc <- enrichDGN(de)      #ORA
edo2_GcCc <- gseDO(GcCc2_list, organism = "hsa") #GSEA

dotplot(edo_GcCc, showCategory=30)
barplot(edo_GcCc, showCategory=20) 

edox_GcCc <- setReadable(edo_GcCc, 'org.Hs.eg.db', 'ENTREZID')

cnetplot(edox_GcCc, foldChange=GcCc2_list)
cnetplot(edox_GcCc, categorySize="pvalue", foldChange=GcCc2_list)
cnetplot(edox_GcCc, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 

heatplot(edox_GcCc, showCategory=5)
heatplot(edox_GcCc, foldChange=geneList, showCategory=5)


# KEGG
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
kk <- enrichKEGG(gene         = GcCc2$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk2a <- setReadable(kk2a, 'org.Hs.eg.db', 'ENTREZID')

kk2a <- gseKEGG(geneList     = genelist_entrezid,
               organism     = 'hsa',
               keyType = "kegg",
               exponent = 1,
               minGSSize = 3,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = TRUE,
               use_internal_data = FALSE,
               seed = FALSE)
browseKEGG(kk2, 'hsa05171')
browseKEGG(kk, 'hsa03050')

library("pathview")
hsa05171 <- pathview(gene.data  = genelist_entrezid,
                     pathway.id = "hsa05171",
                     species    = "hsa",
                     limit      = list(gene=max(abs(genelist_entrezid)), cpd=1))
# ReactomePA
rct <- enrichPathway(gene=GcCc2$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
viewPathway("Regulation of expression of SLITs and ROBOs", 
            readable = TRUE, 
            foldChange = geneList)
heatplot(rct, showCategory=5)

######### test of geneList for GSEA ##############################
genelist <- GcCc2$log2fc[GcCc2$log2fc > 0]
names(genelist) <- GcCc2$ENTREZID [GcCc2$log2fc > 0]
genelist <- na.omit(genelist)
genelist <- sort(genelist, decreasing = TRUE)
view(genelist)

GcCc2$log2fc <- sort(GcCc2$log2fc, decreasing = TRUE)
genelist <- GcCc2$log2fc
names(genelist) <- GcCc2$symbol
genelist <- sort(genelist, decreasing = TRUE)
view(genelist)

GcCc2$log2fc <- sort(GcCc2$log2fc, decreasing = TRUE)
genelist_entrezid <- GcCc2$log2fc
names(genelist_entrezid) <- GcCc2$ENTREZID
genelist_entrezid <- sort(genelist_entrezid, decreasing = TRUE)
view(genelist_entrezid)

keytypes(org.Hs.eg.db)

# GO
ego <- enrichGO(gene          = GcCc2$symbol,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

gse <- gseGO(geneList= genelist, 
             OrgDb= org.Hs.eg.db,
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE)
gse <- pairwise_termsim(gse)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) #Dotplot
emapplot(gse, showCategory = 10) #Enrichment Map
cnetplot(gse, categorySize="pvalue", foldChange=genelist, showCategory = 15) #Category Netplot
ridgeplot(gse) + labs(x = "enrichment distribution") #Ridgeplot


# KEGG
kk2 <- gseKEGG(geneList= genelist_entrezid,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2b <- gseKEGG(geneList     = genelist,
                    organism     = 'hsa',
                    keyType = "kegg",
                    exponent = 1,
                    minGSSize = 3,
                    maxGSSize = 500,
                    eps = 1e-10,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    verbose = TRUE,
                    use_internal_data = FALSE,
                    seed = FALSE)
kk2 <- pairwise_termsim(kk2)
dotplot(kk2, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign) #Dotplot
emapplot(kk2) #Enrichment Map
cnetplot(kk2, categorySize="pvalue", foldChange=genelist_entrezid, showCategory = 10) #Category Netplot
ridgeplot(kk2) + labs(x = "enrichment distribution") #Ridgeplot
dme <- pathview(gene.data=genelist_entrezid, pathway.id="hsa05171", species = "hsa")

GaCa <- dea_signif %>% filter(contrast == "Ga/Ca" & signif=="TRUE")
GaCa$log2fc <- sort(GaCa$log2fc, decreasing = TRUE)
genelistGaCa_entrezid <- GaCa$log2fc
x <- GaCa$symbol
ids <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
GaCa2 <- GaCa %>% left_join(ids, by=c("symbol"="SYMBOL"))
write_csv(GaCa2, "GaCa2.csv")
saveRDS(GaCa2, "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/clusterprofiler/GaCa2.RDS")

GaCa2$log2fc <- sort(GaCa2$log2fc, decreasing = TRUE)
genelistGaCa_entrezid <- GaCa2$log2fc
names(genelistGaCa_entrezid) <- GaCa2$ENTREZID
genelistGaCa_entrezid <- sort(genelistGaCa_entrezid, decreasing = TRUE)
view(genelistGaCa_entrezid)

kk2_GaCa <- gseKEGG(geneList     = genelistGaCa_entrezid,
                organism     = 'hsa',
                keyType = "kegg",
                exponent = 1,
                minGSSize = 3,
                maxGSSize = 500,
                eps = 1e-10,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = TRUE,
                use_internal_data = FALSE,
                seed = FALSE)
kk2_GaCa <- pairwise_termsim(kk2_GaCa)

browseKEGG(kk2_GaCa, 'hsa05171') #COVID19
browseKEGG(kk2_GaCa, 'hsa01100') #Metabolic pathway


hsa05171_GaCa <- pathview(gene.data  = genelistGaCa_entrezid,
                     pathway.id = "hsa05171",
                     species    = "hsa",
                     limit      = list(gene=max(abs(genelistGaCa_entrezid)), cpd=1))

hsa01100_GaCa <- pathview(gene.data  = genelistGaCa_entrezid,
                          pathway.id = "hsa01100",
                          species    = "hsa",
                          limit      = list(gene=max(abs(genelistGaCa_entrezid)), cpd=1))
dotplot(kk2_GaCa, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign) #Dotplot
emapplot(kk2_GaCa) #Enrichment Map
cnetplot(kk2_GaCa, categorySize="pvalue", foldChange=genelistGaCa_entrezid, showCategory = 5) #Category Netplot
ridgeplot(kk2_GaCa, showCategory = 15) + labs(x = "enrichment distribution") #Ridgeplot
