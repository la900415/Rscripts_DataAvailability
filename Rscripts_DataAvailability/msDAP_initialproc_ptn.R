library(tidyverse)
library(dunn.test)



setwd("E:/RESEARCH/ygor/DN/2_tecido_FPspeclib_BEST/msdap_results/2025-03-16_17-11-57_vwmb_BEST")

# Script to open filtered protein abundances in log2 from msDAP,
# makes some value imputations and calculate ANOVA, 
# generating table for further analyses


#### Select protein_abundance__global data filter.tsv file ####

inp_file <- file.choose() # protein_abundance__global data filter.tsv
resfolder <- inp_file %>% dirname() %>% normalizePath( winslash = "\\") 

output_files <- paste0(resfolder,"\\Normalization_Eduard\\")

setwd(resfolder)
if (!dir.exists(output_files)){dir.create(paste0(resfolder,"\\Normalization_Eduard\\"))}

#### Input data ####
inp_ptn <- list()
inp_table <- read.csv2(inp_file, sep = "\t", na.strings = "",blank.lines.skip=TRUE)
inp_samplefile <- read.csv2(file.choose(), sep = "\t",header = TRUE, na.strings = "",blank.lines.skip=TRUE) # samples.tsv


inp_ptn$rawdata$proteinid_all <- inp_table[,grep('^protein_id*',colnames(inp_table),invert=FALSE)]
inp_ptn$rawdata$abs <- inp_table[,grep('^QE*',colnames(inp_table),invert=FALSE)]
inp_ptn$rawdata$abs <- inp_ptn$rawdata$abs %>% sapply(function (x) as.numeric(x))
inp_ptn$rawdata$samplenames <- inp_table[grep('^QE*',colnames(inp_table),invert=FALSE)] %>% names()
inp_ptn$rawdata$sn_simpl <- inp_ptn$rawdata$samplenames %>% sapply(function (x) sub("^.*_(\\d+)$", "\\1", x) )
inp_ptn$rawdata$interaction <- inp_ptn$rawdata$samplenames %>% sapply(function (x) {
  inp_samplefile$group[inp_samplefile$sample_id %in% x]
} )
inp_ptn$rawdata$maneuver <- inp_ptn$rawdata$samplenames %>% sapply(function (x) {
  inp_samplefile$maneuver[inp_samplefile$sample_id %in% x]
} )

#### Kruskal-Wallis test ####

inp_ptn$proc$kruskal$kw_tab_interact <- cbind.data.frame(inp_ptn$rawdata$interaction, t(inp_ptn$rawdata$abs))
inp_ptn$proc$kruskal$kw_tab_maneuver <- cbind.data.frame(inp_ptn$rawdata$maneuver, t(inp_ptn$rawdata$abs))

inp_ptn$proc$kruskal$kw_result_interact <- lapply(inp_ptn$proc$kruskal$kw_tab_interact[-1], 
                                                         function(x) kruskal.test(x ~ inp_ptn$proc$kruskal$kw_tab_interact[,1]))
inp_ptn$proc$kruskal$kw_result_maneuver <- lapply(inp_ptn$proc$kruskal$kw_tab_maneuver[-1], 
                                                         function(x) kruskal.test(x ~ inp_ptn$proc$kruskal$kw_tab_maneuver[,1]))


inp_ptn$proc$kruskal$kw_pvalues_interact <- sapply(inp_ptn$proc$kruskal$kw_result_interact, function(model) model$p.value)
inp_ptn$proc$kruskal$kw_qvalues_interact <- p.adjust(inp_ptn$proc$kruskal$kw_pvalues_interact, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

inp_ptn$proc$kruskal$kw_pvalues_maneuver <- sapply(inp_ptn$proc$kruskal$kw_result_maneuver, function(model) model$p.value)
inp_ptn$proc$kruskal$kw_qvalues_maneuver <- p.adjust(inp_ptn$proc$kruskal$kw_pvalues_maneuver, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

##################  Dunn post-hoc test  ###################################################################
inp_ptn$proc$kruskal$DunnT_interact <- lapply(inp_ptn$proc$kruskal$kw_tab_interact[-1], function(x) {
  dunn.test(x,inp_ptn$proc$kruskal$kw_tab_interact[,1],method="bh",kw=TRUE,alpha=0.05)})

inp_ptn$proc$kruskal$DunnT_maneuver <- lapply(inp_ptn$proc$kruskal$kw_tab_maneuver[-1], function(x) {
  dunn.test(x,inp_ptn$proc$kruskal$kw_tab_maneuver[,1],method="bh",kw=TRUE,alpha=0.05)})

##################  Extracting the signif pairs from Dunn poshoc  ###################################################################
inp_ptn$proc$kruskal$DunnT_signif_pairs_interact <- sapply(inp_ptn$proc$kruskal$DunnT_interact, function(x) x$comparisons[x$P < 0.05 | x$P.adjusted < 0.05 ] )
inp_ptn$proc$kruskal$DunnT_signif_pairs_interact <- sapply(inp_ptn$proc$kruskal$DunnT_signif_pairs_interact, function(x) paste(x, collapse = ";"))

inp_ptn$proc$kruskal$DunnT_signif_pairs_maneuver <- sapply(inp_ptn$proc$kruskal$DunnT_maneuver, function(x) x$comparisons[x$P < 0.05 | x$P.adjusted < 0.05 ] )
inp_ptn$proc$kruskal$DunnT_signif_pairs_maneuver <- sapply(inp_ptn$proc$kruskal$DunnT_signif_pairs_maneuver, function(x) paste(x, collapse = ";"))

################# Convert the results into a data frame for easier viewing ##################
inp_ptn$proc$kw_dunn_tab_output <- cbind.data.frame(inp_ptn$rawdata$abs, 
                                                         inp_ptn$rawdata$proteinid_all,
                                                         inp_ptn$proc$kruskal$kw_pvalues_interact,
                                                         inp_ptn$proc$kruskal$kw_qvalues_interact,
                                                         inp_ptn$proc$kruskal$DunnT_signif_pairs_interact,
                                                         inp_ptn$proc$kruskal$kw_pvalues_maneuver,
                                                         inp_ptn$proc$kruskal$kw_qvalues_maneuver,
                                                         inp_ptn$proc$kruskal$DunnT_signif_pairs_maneuver)
colnames(inp_ptn$proc$kw_dunn_tab_output)[21] <- "protein_id"
colnames(inp_ptn$proc$kw_dunn_tab_output)[22] <- "KW_pvalue_interaction"
colnames(inp_ptn$proc$kw_dunn_tab_output)[23] <- "KW_qvalue_interaction"
colnames(inp_ptn$proc$kw_dunn_tab_output)[24] <- "DunnT_signif_pairs_interaction"
colnames(inp_ptn$proc$kw_dunn_tab_output)[25] <- "KW_pvalue_maneuver"
colnames(inp_ptn$proc$kw_dunn_tab_output)[26] <- "KW_qvalue_maneuver"
colnames(inp_ptn$proc$kw_dunn_tab_output)[27] <- "DunnT_signif_pairs_maneuver"

saveRDS(inp_ptn, "inp_ptn.RDS")

###############################################################################################
######### Using centered data from perseus using global_filt                          #########
###############################################################################################
pdata <- readRDS("E:/RESEARCH/ygor/DN/2_tecido_FPspeclib_BEST/msdap_results/2025-03-16_17-11-57_vwmb_BEST/pdata_glob.RDS")

inp_ptn$rawdata_centered$abs <- pdata@main
inp_ptn$rawdata_centered$proteinid_all <- pdata@annotCols$protein_id
inp_ptn$rawdata_centered$interaction <- pdata@annotRows[-c(1,3)]
inp_ptn$rawdata_centered$maneuver <- pdata@annotRows[-c(1:2)]

saveRDS(inp_ptn, "inp_ptn.RDS")

#### Kruskal-Wallis test ####
inp_ptn$proc_centered$kruskal$kw_tab_interact <- cbind.data.frame(inp_ptn$rawdata_centered$interaction, t(inp_ptn$rawdata_centered$abs))
inp_ptn$proc_centered$kruskal$kw_tab_maneuver <- cbind.data.frame(inp_ptn$rawdata_centered$maneuver, t(inp_ptn$rawdata_centered$abs))

inp_ptn$proc_centered$kruskal$kw_result_interact <- lapply(inp_ptn$proc_centered$kruskal$kw_tab_interact[-1], 
                                                  function(x) kruskal.test(x ~ inp_ptn$proc_centered$kruskal$kw_tab_interact[,1]))
inp_ptn$proc_centered$kruskal$kw_result_maneuver <- lapply(inp_ptn$proc_centered$kruskal$kw_tab_maneuver[-1], 
                                                  function(x) kruskal.test(x ~ inp_ptn$proc_centered$kruskal$kw_tab_maneuver[,1]))


inp_ptn$proc_centered$kruskal$kw_pvalues_interact <- sapply(inp_ptn$proc_centered$kruskal$kw_result_interact, function(model) model$p.value)
inp_ptn$proc_centered$kruskal$kw_qvalues_interact <- p.adjust(inp_ptn$proc_centered$kruskal$kw_pvalues_interact, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

inp_ptn$proc_centered$kruskal$kw_pvalues_maneuver <- sapply(inp_ptn$proc_centered$kruskal$kw_result_maneuver, function(model) model$p.value)
inp_ptn$proc_centered$kruskal$kw_qvalues_maneuver <- p.adjust(inp_ptn$proc_centered$kruskal$kw_pvalues_maneuver, method = "BH") # ajuste redundante, pois teste de Dunn já corrige multiplas comparações

##################  Dunn post-hoc test  ###################################################################
inp_ptn$proc_centered$kruskal$DunnT_interact <- lapply(inp_ptn$proc_centered$kruskal$kw_tab_interact[-1], function(x) {
  dunn.test(x,inp_ptn$proc_centered$kruskal$kw_tab_interact[,1],method="bh",kw=TRUE,alpha=0.05)})

inp_ptn$proc_centered$kruskal$DunnT_maneuver <- lapply(inp_ptn$proc_centered$kruskal$kw_tab_maneuver[-1], function(x) {
  dunn.test(x,inp_ptn$proc_centered$kruskal$kw_tab_maneuver[,1],method="bh",kw=TRUE,alpha=0.05)})

##################  Extracting the signif pairs from Dunn poshoc  ###################################################################
inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_interact <- sapply(inp_ptn$proc_centered$kruskal$DunnT_interact, function(x) x$comparisons[x$P < 0.05 | x$P.adjusted < 0.05 ] )
inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_interact <- sapply(inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_interact, function(x) paste(x, collapse = ";"))

inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_maneuver <- sapply(inp_ptn$proc_centered$kruskal$DunnT_maneuver, function(x) x$comparisons[x$P < 0.05 | x$P.adjusted < 0.05 ] )
inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_maneuver <- sapply(inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_maneuver, function(x) paste(x, collapse = ";"))

################# Convert the results into a data frame for easier viewing ##################
inp_ptn$proc_centered$kw_dunn_tab_output <- cbind.data.frame(inp_ptn$rawdata_centered$abs, 
                                                    inp_ptn$rawdata_centered$proteinid_all,
                                                    inp_ptn$proc_centered$kruskal$kw_pvalues_interact,
                                                    inp_ptn$proc_centered$kruskal$kw_qvalues_interact,
                                                    inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_interact,
                                                    inp_ptn$proc_centered$kruskal$kw_pvalues_maneuver,
                                                    inp_ptn$proc_centered$kruskal$kw_qvalues_maneuver,
                                                    inp_ptn$proc_centered$kruskal$DunnT_signif_pairs_maneuver)
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[21] <- "protein_id"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[22] <- "KW_pvalue_interaction"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[23] <- "KW_qvalue_interaction"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[24] <- "DunnT_signif_pairs_interaction"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[25] <- "KW_pvalue_maneuver"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[26] <- "KW_qvalue_maneuver"
colnames(inp_ptn$proc_centered$kw_dunn_tab_output)[27] <- "DunnT_signif_pairs_maneuver"

saveRDS(inp_ptn, "inp_ptn.RDS")

write_tsv(inp_ptn$proc$kw_dunn_tab_output, "kw_dunn_tab_output_non-centered.tsv")
write_tsv(inp_ptn$proc_centered$kw_dunn_tab_output, "kw_dunn_tab_output_centered.tsv")
