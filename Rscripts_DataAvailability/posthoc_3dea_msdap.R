
####                   post-hoc combine results from multiple DEA models                               #############
#### recommend selecting proteins significant in at least 2 of these algorithms as a compromise        #############
#### between using multiple algorithms to maximize your search for proteins-of-interest and robustness #############
#### against false positives                                                                           #############
da_CcGc <- da_3contrast %>%
  filter(rowSums(across(c(`signif_msempire_contrast: CTRL.cent vs GAM.cent`, 
                          `signif_deqms_contrast: CTRL.cent vs GAM.cent`, 
                          `signif_msqrob_contrast: CTRL.cent vs GAM.cent`), ~ . == TRUE)) >= 2)

da_CcGc$log2FC_avg_CcGc <- rowMeans(da_CcGc[, c("foldchange.log2_msempire_contrast: CTRL.cent vs GAM.cent", 
                                                "foldchange.log2_deqms_contrast: CTRL.cent vs GAM.cent",
                                                "foldchange.log2_msqrob_contrast: CTRL.cent vs GAM.cent")], na.rm = TRUE)

da_CcGc$pvalue_avg_CcGc <- rowMeans(da_CcGc[, c("pvalue_msempire_contrast: CTRL.cent vs GAM.cent", 
                                                "pvalue_deqms_contrast: CTRL.cent vs GAM.cent",
                                                "pvalue_msqrob_contrast: CTRL.cent vs GAM.cent")], na.rm = TRUE)


