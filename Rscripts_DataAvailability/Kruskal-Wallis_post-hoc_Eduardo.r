#### Kruskal-Wallis com post-hoc test ####
library(dunn.test)

inp_ptn$proc$kruskal$kw_table <- cbind.data.frame(inp_ptn$rawdata$groups, t(inp_ptn$proc$norm$transf_zscore))
inp_ptn$proc$kruskal$kw_results_posthoc <- lapply(inp_ptn$proc$kruskal$kw_table[-1], function(x) kruskal.test(x ~ inp_ptn$proc$kruskal$kw_table[,1]))

inp_ptn$proc$kruskal$pvalues <- sapply(inp_ptn$proc$kruskal$kw_results, function(model) model$p.value)
inp_ptn$proc$kruskal$qvalues <- p.adjust(inp_ptn$proc$kruskal$pvalues, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

##################  Dunn post-hoc test  ###################################################################
inp_ptn$proc$kruskal$DunnT <- lapply(inp_ptn$proc$kruskal$kw_table[-1], function(x) {
  dunn.test(x,inp_ptn$proc$kruskal$kw_table[,1],method="bh",kw=TRUE,alpha=0.05)})

##################  Extracting the signif pairs from Dunn poshoc  ###################################################################
inp_ptn$proc$kruskal$DunnT_comparisons <- sapply(inp_ptn$proc$kruskal$DunnT, function(x) x$comparisons[x$P < 0.05 | x$P.adjusted < 0.05 ] )
inp_ptn$proc$kruskal$DunnT_comparisons <- sapply(inp_ptn$proc$kruskal$DunnT_comparisons, function(x) paste(x, collapse = ";"))

################# Convert the results into a data frame for easier viewing ##################
inp_ptn$proc$kruskal$kw_table_output <- cbind.data.frame(inp_ptn$rawdata$abs, 
                                                   inp_ptn$rawdata$proteinid_all,
                                                   inp_ptn$proc$kruskal$pvalues,
                                                   inp_ptn$proc$kruskal$qvalues,
                                                   inp_ptn$proc$kruskal$DunnT_comparisons)
colnames(inp_ptn$proc$kruskal$kw_table_output)[21] <- "protein_id"
colnames(inp_ptn$proc$kruskal$kw_table_output)[22] <- "KW_pvalue_interaction"
colnames(inp_ptn$proc$kruskal$kw_table_output)[23] <- "KW_qvalue_interaction"
colnames(inp_ptn$proc$kruskal$kw_table_output)[24] <- "DunnT_poshoc_signif_pairs_interaction"

  



