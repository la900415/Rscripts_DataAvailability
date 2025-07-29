#### Kruskal-Wallis com post-hoc test ####

install.packages(dunn.test)
library(dunn.test)

pv_lim <- 0.1 # significance of post-hoc test/alpha

inp_ptn$proc$kruskal$kw_table <- cbind.data.frame(inp_ptn$rawdata$groups ,t(inp_ptn$proc$norm$transf_zscore))
inp_ptn$proc$kruskal$kw_results_posthoc <- lapply(inp_ptn$proc$kruskal$kw_table[-1], function(x) kruskal.test(x ~ inp_ptn$proc$kruskal$kw_table[,1]))

inp_ptn$proc$kruskal$pvalues <- sapply(inp_ptn$proc$kruskal$kw_results, function(model) model$p.value)
inp_ptn$proc$kruskal$qvalues <- p.adjust(inp_ptn$proc$kruskal$pvalues, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

inp_ptn$proc$kruskal$DunnT <- lapply(inp_ptn$proc$kruskal$kw_table[-1], function(x) {
  dunn.test(x,inp_ptn$proc$kruskal$kw_table[,1],method="bh",kw=TRUE,alpha=pv_lim)})

