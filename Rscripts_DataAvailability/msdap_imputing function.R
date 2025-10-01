library(msdap)

# imputation using Random Forests, as implemented in the missForest R package. doi: 10.1093/bioinformatics/btr597
# reference for chosing this imputation approach @ doi: 10.1038/s41598-021-81279-4
# assuming peptide filters were applied upstream, the % missingness should be low and imputation would have relatively low impact
# note that this can be quite slow
impute_missforest = function(x_as_log2, mask_sample_groups = NA, ...) {
  cat("Now calling the custom normalization function, impute_missforest(), that imputes missing values (if any) using the missForest package\n")
  time_start = Sys.time()
  
  # report the number of missing values
  nmiss = sum(is.na(x_as_log2))
  nrowmiss = sum( rowSums(is.na(x_as_log2)) > 0 )
  
  if(nmiss > 0) {
    cat(sprintf("%d/%d datapoints were NA (%.1f%%), %d/%d peptides had at least 1 NA value (%.1f%%)\n",
                nmiss, length(x_as_log2), nmiss / length(x_as_log2) * 100,
                nrowmiss, nrow(x_as_log2), nrowmiss / nrow(x_as_log2) * 100 ))
    
    # impute the matrix using the missForest package at default settings
    x_as_log2 = missForest::missForest(x_as_log2, verbose = F, )$ximp
    # to optionally silence missForest messages to the console, wrap above statement in capture.output() and use <- instead of =
    cat(sum(is.na(x_as_log2)), "NA post imputation\n")
    print(Sys.time() - time_start)
  } else {
    cat("nothing to impute, no NA values\n")
  }
  return(x_as_log2)
}


# imputation using local least squares (LLS), as implemented in the pcaMethods R package. doi: 10.1093/bioinformatics/btm069
# reference for chosing this imputation approach @ doi: 10.1038/s41598-021-81279-4
# assuming peptide filters were applied upstream, the % missingness should be low and imputation would have relatively low impact
# note that this can be quite slow
impute_lls = function(x_as_log2, mask_sample_groups = NA, ...) {
  cat("Now calling the custom normalization function, impute_lls(), that imputes missing values (if any) using the local least squares (LLS) function of the pcaMethods package\n")
  time_start = Sys.time()
  
  # report the number of missing values
  nmiss = sum(is.na(x_as_log2))
  nrowmiss = sum( rowSums(is.na(x_as_log2)) > 0 )
  
  if(nmiss > 0) {
    cat(sprintf("%d/%d datapoints were NA (%.1f%%), %d/%d peptides had at least 1 NA value (%.1f%%)\n",
                nmiss, length(x_as_log2), nmiss / length(x_as_log2) * 100,
                nrowmiss, nrow(x_as_log2), nrowmiss / nrow(x_as_log2) * 100 ))
    
    # ! this choice of K is arbitrary, we didn't validate nor benchmark as this is just a proof-of-concept to demonstrate imputation in MS-DAP !
    transposed_imputation_result = pcaMethods::llsImpute(t(x_as_log2), allVariables = T, k = ifelse(nrow(x_as_log2) < 100, 3, floor(nrow(x_as_log2) * 0.01)) )
    x_as_log2 = t(transposed_imputation_result@completeObs)
    cat(sum(is.na(x_as_log2)), "NA post imputation\n")
    print(Sys.time() - time_start)
  } else {
    cat("nothing to impute, no NA values\n")
  }
  return(x_as_log2)
}



### use these imputation by specifying the respective function names as a "normalization" algorithm
### in the example below, we add these new function in the `norm_algorithm` parameter of `msdap::analysis_quickstart()`

f <- system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")
dataset = import_dataset_skyline(f, confidence_threshold = 0.01, return_decoys = F, acquisition_mode = "dia")
dataset = sample_metadata_custom(dataset, group_regex_array = c(A = "007|009|011", B = "008|010|012") )
dataset = setup_contrasts(dataset, contrast_list = list(c("A", "B")))
print(dataset$samples %>% select(sample_id, group))


dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 0,
  filter_min_quant = 2,
  # fraction_min_detect = 0.25, # optionally add more stringent filtering
  # fraction_min_quant = 0.5,
  # here we apply filtering and normalization to each sample group, across the entire dataset at once
  filter_by_contrast = FALSE,
  # imputation setting is mixed in the sequence application of normalizations:
  # first, normalize the peptide-level data matrix
  # next, call the random-forest imputation function (applied to peptide-level)
  # finally, normalize at protein-level to ensure the mode of log2 foldchanges between sample groups is zero
  norm_algorithm = c("vwmb", "impute_missforest", "modebetween_protein"),
  # norm_algorithm = c("vwmb", "impute_lls", "modebetween_protein"), # alternatively, use LLS imputation
  dea_algorithm = c("deqms", "msempire"), # select whichever DEA algorithms. DEqMS is a reliable default
  dea_qvalue_threshold = 0.05, # 5% FDR
  dea_log2foldchange_threshold = 0,
  output_dir = "C:/temp/",
  output_within_timestamped_subdirectory = TRUE
)