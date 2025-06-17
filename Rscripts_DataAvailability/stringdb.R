library(STRINGdb)
string_db <- STRINGdb$new( version="12.0", species=9606, score_threshold=200, network_type="full", input_directory="C:/Users/Laura/Downloads")
data(diff_exp_example1)
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )
view(example1_mapped)

hits <- example1_mapped$STRING_id[1:200]
view(hits)
string_db$plot_network( hits )

#2 Payload mechanism
 # filter by p-value and add a color column
   # (i.e. green down-regulated gened and red for up-regulated genes)
   example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05), logFcColStr="logFC" )
 # post payload information to the STRING server
   payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color )
 # display a STRING network png with the "halo"
   string_db$plot_network( hits, payload_id=payload_id )

#4 CLUSTERING   
 # get clusters
clustersList <- string_db$get_clusters(example1_mapped$STRING_id[1:600])
# plot first 4 clusters
library(igraph)
par(mfrow=c(2,2)) 
for(i in seq(1:4)){
    string_db$plot_network(clustersList[[i]]) 
  }

STRINGdb$help("get_graph")
STRINGdb$help("export")
STRINGdb$help("plot_network")


######################### stringdb for Cc/Ca ###############################################
string_db <- STRINGdb$new( version="12.0", species=9606, score_threshold=200, network_type="full", input_directory="C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-08-15_00-45-32_9excluded_msqrob/stringdbR")
CcCa <- dea3 %>% filter(contrast=="Cc/Ca" & signif=="TRUE") %>% dplyr::select(pvalue, foldchange.log2, ENTREZID)
CcCa_mapped <- string_db$map (CcCa, "ENTREZID", removeUnmappedRows=TRUE)
view(CcCa_mapped)

CcCa_ppi <- CcCa_mapped$STRING_id
view(CcCa_ppi)
string_db$plot_network (CcCa_ppi)