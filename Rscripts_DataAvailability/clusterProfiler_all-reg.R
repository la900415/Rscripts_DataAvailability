library(org.Mm.eg.db)
library(clusterProfiler) 

DAPs_allregions <- read_excel("C:/Users/Luis Ariel/Downloads/DAPs_allregions.xlsx")

ck_ego_all <- compareCluster(UniProt_code~Regions, data=DAPs_allregions, fun="enrichGO", ont="BP",
                            pvalueCutoff=0.05, keyType="UNIPROT", OrgDb=org.Mm.eg.db, pAdjustMethod="BH", 
                            readable=F)
ck_ego_all2 <- setReadable(ck_ego_all, OrgDb = org.Mm.eg.db, keyType="UNIPROT")
saveRDS(object=ck_ego_all, file = "C:/Users/Luis Ariel/Downloads/ck_ego_all.RDS")

dotplot(ck_ego_all, showCategory=10,includeAll=F,label_format=70) #+
  theme(axis.text.x = element_text(size = 14, angle=90, hjust=1.1, vjust=1.2),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16) )

#Remove redundancy
ck_ego_all3 <- simplify(ck_ego_all, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL)
ck_ego_all_GOlevel4 <- gofilter(ck_ego_all, level = 4)

dotplot(ck_ego_all_GOlevel4, showCategory=10,includeAll=F,label_format=70)
