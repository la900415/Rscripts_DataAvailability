library(msdap)
library(tidyverse)
############## Loading functions of MSDAP ##########################################
ggplot_sample_detect_counts_barplots = function(samples, samples_colors_long) {
  tib = samples
  
  # if there are 'all peptide' counts (detect + MBR/quant-only), include these in plot tibble
  if("all_peptides" %in% colnames(samples)) {
    tib$quantified_peptides = tib$all_peptides - tib$detected_peptides
    tib$quantified_proteins = tib$all_proteins - tib$detected_proteins
  }
  
  # from wide to long format
  tib = tib %>%
    select(!!c("shortname", grep("(detected|quantified)_(peptides|proteins)", colnames(tib), value = T))) %>%
    pivot_longer(cols = -shortname, names_to = "type", values_to = "count")
  
  # flip levels/sorting because we use coord_flip() downstream in ggplot
  tib$shortname = factor(tib$shortname, levels = rev(samples$shortname))
  tib$pep_or_prot = ifelse(grepl("peptide", tib$type), "peptides", "proteins")
  tib$detect_or_quant = ifelse(grepl("detected", tib$type), "detect", "quant")
  
  # rename labels for plot clarity
  tib$type["detected_peptides" == tib$type] = "peptides: identif & quant"
  tib$type["quantified_peptides" == tib$type] = "peptides: only quant"
  tib$type["detected_proteins" == tib$type] = "proteins: identif & quant"
  tib$type["quantified_proteins" == tib$type] = "proteins: only quant"
  
  # color-coding
  # clr = c("detected_peptides"="#0570b0", "quantified_peptides"="#74a9cf", "detected_proteins"="#6a51a3", "quantified_proteins"="#9e9ac8")
  clr = c("peptides: identif & quant"="#0570b0", "peptides: only quant"="#74a9cf",
          "proteins: identif & quant"="#E80909", "proteins: only quant"="#ED8E8E")
  
  # coordinates for our customized dual-barplot
  tib = tib %>% group_by(pep_or_prot, shortname) %>% mutate(total=sum(count)) %>% arrange(shortname, type)
  tib_dual = tib %>% group_by(pep_or_prot) %>% mutate(count_scaled = count / max(total) * ifelse(pep_or_prot=="proteins", -1, 1),
                                                      count_scaled_max = total / max(total) * ifelse(pep_or_prot=="proteins", -1, 1))
  tib_dual$outlier_lowside = abs(tib_dual$count_scaled_max) < 0.3
  
  # custom x-axis labels, since this dimension not ranges from -1:1
  ticks_peptides = pretty(1:max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total)), n = 5, min.n = 4)
  ticks_peptides_scaled = ticks_peptides / max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total))
  ticks_proteins = pretty(1:max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total)), n = 5, min.n = 4)
  ticks_proteins_scaled = ticks_proteins / max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total))
  # combine
  ticks_coord = c(-1 * rev(ticks_proteins_scaled), ticks_peptides_scaled[-1])
  ticks_label = c(rev(ticks_proteins), ticks_peptides[-1])
  
  
  # plot code is somewhat convoluted by careful alignment of text labels. simplest QC plot of data as-is: ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) + geom_bar(stat = "identity", position = position_stack(reverse = T)) + coord_flip()
  text_size = ifelse(nrow(tib) > 50, 2, 3.5)
  p = ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T)) +
    geom_text(aes(label = count,
                  y = ifelse(detect_or_quant=="detect", ifelse(pep_or_prot=="proteins", -0.025, 0.025), count_scaled_max),
                  hjust = ifelse((type %in% c("peptides: identif & quant", "proteins: only quant") & !c(type=="proteins: only quant" & outlier_lowside)) |
                                   (type %in% c("peptides: only quant") & outlier_lowside), -0.1, 1.1)),
              colour = ifelse(tib_dual$outlier_lowside & tib_dual$detect_or_quant=="quant", "darkgrey", "white"),
              size = text_size, # 4 is default
              check_overlap = F) +
    scale_y_continuous(breaks=ticks_coord, labels=ticks_label, ) +
    scale_fill_manual(values = clr, guide = guide_legend(nrow = 2)) +
    coord_flip() +
    labs(x = "samples", y = "protein / peptide counts") +
    ggpubr::theme_pubr(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(0.75,"line"), legend.text = element_text(size=9),
          panel.grid = element_blank(), axis.line.y.right = element_blank(),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if(nrow(tib) > 50)
    p = p + theme(axis.text.y.left = element_text(size=6))
  
  return(p)
  ###### some reference code for separate plots
  # tib_pep = tib %>% filter(pep_or_prot == "peptides") %>% droplevels()
  # tib_prot = tib %>% filter(pep_or_prot == "proteins") %>% droplevels()
  #
  # plot_text_lim = 0.1 * max(tib_pep$count)
  # p_pep = ggplot(tib_pep, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
  #
  # plot_text_lim = 0.1 * max(tib_prot$count)
  # p_prot = ggplot(tib_prot, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
}
ggplot_peptide_detect_frequency = function(peptides, samples) {
  ## peptide detect counts, mapped to samples
  tib = peptides %>% filter(detect) %>% select(peptide_id, sample_id) %>% add_count(peptide_id, name = "z")
  
  tib_plot = tib %>%
    count(sample_id, z) %>%
    arrange(desc(z))
  
  ## sample sorting
  tib_samples = tib_plot %>%
    count(sample_id, wt=n) %>%
    arrange(n) %>%
    select(-n) %>%
    left_join(samples %>% select(sample_id, shortname), by="sample_id")
  
  tib_plot = tib_plot %>% left_join(tib_samples, by="sample_id")
  tib_plot$sample_id = factor(tib_plot$sample_id, levels = tib_samples$sample_id)
  tib_plot$shortname = factor(tib_plot$shortname, levels = tib_samples$shortname)
  
  # enforce integer breaks on color scale; https://stackoverflow.com/a/44886993
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
  
  p = ggplot(tib_plot, aes(x=shortname, y=n, fill=z)) +
    geom_bar(position="stack", stat="identity", ) +
    coord_flip() +
    scale_fill_distiller(palette = "Spectral", direction=-1, breaks = int_breaks) +
    # viridis::scale_fill_viridis(option = "B") +
    labs(title = "Number of samples in which a peptide is identified vs presence in individual sample", y = "Identified peptides", x = "Samples", fill = "#samples") +
    ggpubr::theme_pubr() +
    theme(plot.title = element_text(size=10), legend.position = "right")
  
  if(nrow(tib_plot) > 50)
    p = p + theme(axis.text.y.left = element_text(size=8))
  
  return(p)
}
plot_sample_pca = function(matrix_sample_intensities, samples, samples_colors, sample_label_property = "auto", pch_as_exclude = TRUE, infer_continuous_scale = TRUE, pca_dims = list(1:2, c(1, 3), 2:3)) {
  if(!sample_label_property %in% c("auto", "shortname", "index", "index_asis")) {
    append_log(paste("invalid value for parameter 'sample_label_property':", sample_label_property), type = "error")
  }
  if(!all(c("sample_index", "sample_id", "shortname") %in% colnames(samples))) {
    append_log(paste(paste(c("sample_index", "sample_id", "shortname"), collapse = ", "), "are required column in the samples table"), type = "error")
  }
  if(!(is.list(pca_dims) && all(unlist(lapply(pca_dims, function(x) length(x) == 2 && is.numeric(x) && all.equal(x, round(x)) == TRUE))) && max(unlist(pca_dims)) <= 15) ) {
    append_log("pca_dims parameter must be a list that contains integer value pairs (max value: 15)", type = "error")
  }
  
  
  PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = as.integer(max(unlist(pca_dims))), seed = 123, maxIterations = 2000)
  mat_pca = PPCA@scores # rownames = sample_id
  pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores)))
  
  props = grep("^(sample_id$|sample_index$|shortname$|contrast:)", colnames(samples_colors), ignore.case = T, value = T, invert = T)
  
  # if more than 50 samples, don't ggrepel the sample labels and instead simply plot each label at the PCA x/y location (color-coded by metadata as per usual). No point/shape is drawn in this case
  textonly_sample_labels = nrow(samples) > 50 # default
  prop_sample_labels = "shortname" # default
  if(sample_label_property == "auto") {
    # only shortname if less than 25 samples AND 80% of those are actually short strings (<=30 characters)
    prop_sample_labels = ifelse(nrow(samples) < 25 && sum(nchar(samples$shortname) > 30) < nrow(samples) * 0.2, "shortname", "sample_index") # alternatively; nrow(samples) * stats::median(nchar(samples$shortname)) < 25*10
  }
  if(sample_label_property == "index") {
    prop_sample_labels = "sample_index"
  }
  if(sample_label_property == "index_asis") {
    prop_sample_labels = "sample_index"
    textonly_sample_labels = TRUE
  }
  
  
  pcaplots = list()
  for (prop in props) { #prop=props[1]
    # don't plot 'exclude' property if there are no excluded samples
    if(prop == "exclude" && !any(samples %>% filter(sample_id %in% rownames(mat_pca)) %>% pull(exclude))) {
      next
    }
    
    plotlist = list()
    for (dims in pca_dims) { # dims=1:2
      # extract data from PCA object  &  join with sample metadata to get color-coding and label
      tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>%
        left_join(samples %>% select(sample_id, label=!!prop_sample_labels, prop=!!prop), by="sample_id") %>%
        left_join(samples_colors %>% select(sample_id, clr=!!prop), by="sample_id")
      tib = left_join(tib, samples %>% select(sample_id, exclude), by="sample_id")
      
      plot_as_numeric = FALSE
      if(infer_continuous_scale) {
        # check if current property (column in sample metadata) only contains numbers
        prop_is_numeric = suppressWarnings(all(is.finite(as.numeric(na.omit(tib$prop)))))
        # make a numeric plot only if there are more than 3 unique values. Also covers the case of booleans that should be plotted as categorical variables
        plot_as_numeric = !prop %in% c("group", "exclude") && prop_is_numeric && n_distinct(as.numeric(tib$prop), na.rm = T) > 2
      }
      
      # mutate the actual data before making ggplot object
      if(plot_as_numeric) {
        tib$prop = as.numeric(tib$prop)
      } else {
        # for categorical variables, convert NA values to a character so they show up in the legend
        tib$prop[is.na(tib$prop)] = "<NA>"
        tib$prop = as.character(tib$prop)
      }
      
      # base plot
      p = ggplot(tib, aes(x = x, y = y, label = label, colour = prop, fill = prop)) +
        guides(alpha = "none", fill = "none") +
        labs(
          title = "",
          x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
          y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
        ) +
        # coord_cartesian(clip = "off") +
        theme_bw()
      
      
      if(plot_as_numeric) {
        tib = tib %>%
          mutate(
            # add alpha to colors
            clr_fill = paste0(clr, "66"),
            clr = paste0(clr, "BB"),
            # scale color values between 0~1
            prop = scales::rescale(prop)) %>%
          # sort values from low to high
          # doesn't matter that the order of the values submitted to ggplot2::scale_color_X is changed from when we created ggplot object above
          arrange(prop)
        
        clr_na = "#BBBBBBBB"
        if(any(is.na(tib$prop))) {
          # get NA color from data
          clr_na = tib %>% filter(is.na(prop)) %>% pull(clr) %>% head(1)
          # set color code to NA or they'll show up in ggplot's color scale / legend
          tib$clr[is.na(tib$prop)] = NA
          tib$clr_fill[is.na(tib$prop)] = NA
        }
        
        
        # for gradient colors, sorting by respective value is important (already converted the values to numeric type a few lines above)
        p = p +
          scale_color_gradientn(name = gsub("[ _]+", " ", prop),
                                values = tib$prop, # set the values parameter to use the exact value-to-color mapping we computed upstream
                                colours = tib$clr,
                                aesthetics = "colour",
                                na.value = clr_na) +
          scale_color_gradientn(values = tib$prop,
                                colours = tib$clr_fill,
                                aesthetics = "fill",
                                na.value = clr_na) +
          theme(legend.text = element_text(angle=90, hjust=0.5, vjust=0.5),
                legend.title = element_text(size=10, face = "bold"))
      } else {
        uprop = gtools::mixedsort(unique(tib$prop))
        clr_map = array(tib$clr[match(uprop, tib$prop)], dimnames=list(uprop))
        
        p = p +
          scale_colour_manual(
            name = gsub("[ _]+", " ", prop),
            values = paste0(clr_map, "BB"), # named array, value=color, name=property
            breaks = names(clr_map), # sort the legend
            aesthetics = "colour") +
          scale_colour_manual(
            values = paste0(clr_map, "66"), # named array, value=color, name=property
            breaks = names(clr_map), # sort the legend
            aesthetics = "fill") +
          guides(colour = guide_legend(title.position = "top")) +
          theme(legend.text = element_text(size = ifelse(length(clr_map) < 6,
                                                         10,
                                                         ifelse(length(clr_map) < 10, 8, 6)) ),
                legend.title = element_text(size=10, face = "bold"))
      }
      
      if(textonly_sample_labels) {
        p_labeled = p + geom_text(show.legend = FALSE, size = 2)
      } else {
        p_labeled = p +
          geom_point(aes(shape = I(ifelse(exclude & pch_as_exclude, 0, 21)))) +
          ggrepel::geom_text_repel(show.legend = FALSE, size = 2, segment.alpha = .3, min.segment.length = 0, na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, box.padding = 0.2, seed = 123)
      }
      
      plotlist[[length(plotlist) + 1]] = p + geom_point(aes(shape = I(ifelse(exclude & pch_as_exclude, 0, 21))))
      plotlist[[length(plotlist) + 1]] = p_labeled
    }
    
    # finally, collapse individual PCA plots
    # importantly: if samples are labeled by index, add a warning/note
    # eg; if user labels samples 1~6 as shortname and our dataset$samples table has a different order (causing sample_index to not align), users may mistake sample identities
    # users labeling samples with shortnames like  WT:1,2,4 KO:5,6,7  happens but all too often, so this note must be visible on same page as the plot
    plotlist_merged = ggpubr::ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
    if(prop_sample_labels == "sample_index") {
      plotlist_merged = ggpubr::annotate_figure(plotlist_merged, top = ggpubr::text_grob('samples are labeled by "sample_index" (described in samples.xlsx included with MS-DAP output)', size = 10))
    }
    pcaplots[[prop]] = plotlist_merged
  }
  
  
  ### code snippet for adding PCA on any custom continuous scale. eg, some aspect that reflects sample quality such as; #detect, outliers deviation in retention time, impact on CoV
  # counts = peptides %>% group_by(sample_id) %>% summarise(detect = sum(detect), quant=n()) %>% left_join(samples %>% select(sample_id, shortname), by="sample_id")
  # plotlist = list()
  # for (dims in list(1:2, c(1, 3), 2:3)) { # dims=1:2
  #   tib = tibble(x = mat_pca[, dims[1]], y = mat_pca[, dims[2]], sample_id = rownames(mat_pca)) %>% left_join(counts, by="sample_id")
  #
  #   p = ggplot(tib, aes(x = x, y = y, label = shortname, colour = detect)) +
  #     geom_point() +
  #     scale_color_gradient(name = "detected peptides") +
  #     labs(
  #       title = "",
  #       x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
  #       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100)
  #     ) +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5),
  #           legend.text = element_text(angle=90, hjust=0.5, vjust=0.5))
  #
  #   p_labeled = p + ggrepel::geom_text_repel(show.legend = FALSE, size = 2)
  #
  #   plotlist[[length(plotlist) + 1]] = p
  #   plotlist[[length(plotlist) + 1]] = p_labeled
  # }
  # pcaplots[[length(pcaplots) + 1]] = ggarrange(plotlist = plotlist, ncol = 2, nrow = length(plotlist) / 2, common.legend = T, legend = "bottom")
  
  return(pcaplots)
}
plot_sample_pca__sample_in_contrast = function(dataset, contr, sample_label_property = "auto") {
  if(length(contr) != 1 || !(contr %in% colnames(dataset$peptides) || contr %in% colnames(dataset$samples)) ) {
    append_log("The contrast has to be a column name in the dataset's samples table, as created by the setup_contrasts function. Typical use-case: setup_contrasts(), analysis_quickstart(), then call this function.", type="error")
  }
  
  if(contr %in% colnames(dataset$peptides)) {
    intensity_col_contr = contr
  } else {
    intensity_col_contr = paste0("intensity_", contr)
    if(! intensity_col_contr %in% colnames(dataset$peptides)) {
      append_log(sprintf("The column '%s' cannot be found in the dataset's peptides table, by_contrast filtering has to be enabled for this function to work (it relies on the subsetting of samples + filtering + normalization procedures that create a data matrix specific to some statistical contrast). Typical use-case: setup_contrasts(), analysis_quickstart(), then call this function.", intensity_col_contr), type="error")
    }
  }
  
  # sample color-coding, exactly like report_as_rmarkdown.R  (should compute colors based on entire dataset so color-coding is the same as the report that features the full dataset)
  samples_colors_long = sample_color_coding__long_format(dataset$samples)
  samples_colors = samples_colors_long %>% dplyr::select(sample_id, shortname, prop, clr) %>% tidyr::pivot_wider(id_cols = c(sample_id, shortname), names_from=prop, values_from=clr)
  # case data to wide-format, then plot PCA
  tibw = dataset$peptides %>%
    dplyr::select(key_peptide, sample_id, intensity=!!intensity_col_contr) %>% dplyr::filter(!is.na(intensity)) %>%
    tidyr::pivot_wider(id_cols = key_peptide, names_from = sample_id, values_from = intensity)
  # return the results of the respective plot function
  return( suppressWarnings(plot_sample_pca(as_matrix_except_first_column(tibw), samples = dataset$samples, samples_colors = samples_colors, sample_label_property = sample_label_property)) )
}
plot_volcano = function(stats_de, log2foldchange_threshold = NA, qvalue_threshold = NA, mtitle = "", label_mode = "topn_pvalue", label_target = 25, label_avoid_overlap = TRUE, show_plots = FALSE) {
  if(!is.data.frame(stats_de) || nrow(stats_de) == 0) {
    return(list())
  }
  
  # input validation
  stopifnot(length(label_mode) == 1 && label_mode %in% c("topn_pvalue", "signif", "protein_id"))
  stopifnot(!(label_mode == "topn_pvalue" && !is.finite(label_target))) # if labeling 'top N best pvalue', must specify an integer amount
  stopifnot(!(label_mode == "protein_id" && !all(is.character(label_target))) ) # if labeling protein_id, these must all be strings
  # note; for label_mode="signif" setting the label_target parameter is ignored so we don't have to input validate it
  
  # replace invalid input with NA (0 is invalid input, as it doesn't filter/discard anything)
  if(length(log2foldchange_threshold) != 1 || !is.finite(log2foldchange_threshold) || log2foldchange_threshold == 0) {
    log2foldchange_threshold = NA
  }
  if(length(qvalue_threshold) != 1 || !is.finite(qvalue_threshold) || qvalue_threshold == 0) {
    qvalue_threshold = NA
  }
  
  # iterating over contrasts should be done upstream
  if("contrast" %in% names(stats_de) && n_distinct(stats_de$contrast) != 1) {
    append_log("stats_de parameter contains more than 1 unique contrast. If your DEA results contain multiple contrasts, iterative over these and call this volcano plot function for clean subsets of data that only contain stats for 1 contrast", type = "error")
  }
  
  
  ########### format input data
  stats_de$pvalue[!is.finite(stats_de$pvalue)] = NA
  stats_de$qvalue[!is.finite(stats_de$qvalue)] = NA
  
  # optionally, if the user wants to re-define the qvalue cutoff for this plot, do that first. Next, optionally take all significant hits and add filtering by foldchange
  if(!is.na(qvalue_threshold)) {
    stats_de = stats_de %>% mutate(signif = is.finite(qvalue) & qvalue <= qvalue_threshold)
  }
  if(!"signif" %in% colnames(stats_de) || !all(stats_de$signif %in% c(T,F))) {
    append_log("stats_de parameter must contain a column 'signif', or provide a numeric value for parameter 'qvalue_threshold'", type = "error")
  }
  if(!is.na(log2foldchange_threshold)) {
    log2foldchange_threshold = abs(log2foldchange_threshold)
    stats_de = stats_de %>% mutate(signif = signif & abs(foldchange.log2) >= log2foldchange_threshold)
  }
  
  if(!"dea_algorithm" %in% colnames(stats_de)) {
    stats_de$dea_algorithm = "statistics"
  }
  
  ## pretty-print labels
  stats_de = add_protein_prettyprint_label(stats_de)
  
  
  ########### create volcano plots
  
  result = list()
  for (algo_name in unique(stats_de$dea_algorithm)) { #algo_name = "ebayes"
    # prepare data for current contrast
    tib = stats_de %>%
      filter(dea_algorithm == algo_name) %>%
      drop_na(foldchange.log2, pvalue, qvalue) %>%
      # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
      mutate(minlog10qval = minlog10(qvalue),
             foldchange.log2_abs = abs(foldchange.log2)) %>%
      # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
      arrange(desc(pvalue)) %>%
      # reduce tibble size
      select(protein_id, label, foldchange.log2, foldchange.log2_abs, minlog10qval, signif)
    
    if(nrow(tib) == 0) {
      next
    }
    
    # which proteins should get a text label?
    tib$flag_plot_label = FALSE
    lbl_style = ""
    if(label_mode == "signif") { # all significant proteins
      tib$flag_plot_label = tib$signif == TRUE
      lbl_style = "significant"
    }
    if(label_mode == "topn_pvalue") { # topN 'best' pvalue
      tmp = min(nrow(tib), label_target)
      tib$flag_plot_label = rep(c(F,T), c(nrow(tib) - tmp, tmp)) # set last N rows to TRUE, this works because we just sorted `tib` by pvalue descending
      lbl_style = paste(tmp, "best qvalue")
    }
    if(label_mode == "protein_id") { # user-specified protein(group) IDs
      tib$flag_plot_label = tib$protein_id %in% label_target
      lbl_style = "selected proteins"
    }
    
    # classify up/down regulated
    tib$updown = ifelse(tib$foldchange.log2 < 0, "down", "up")
    tib$updown[tib$signif != TRUE] = "unchanged"
    
    ### find outliers, values that are so far away that they may skew the plot's appearance, and classify data points accordingly
    xmax_nooutlier = c(-1,1) * max(tib$foldchange.log2_abs, na.rm = T)
    ymax_nooutlier = c(0, max(2, tib$minlog10qval, na.rm=T)) # hardcoded limit; y-axis data goes to 10^-2 at least
    xmax = c(-1,1) * max(abs(quantile(tib$foldchange.log2, probs = c(.005, .995), na.rm = T)))
    ymax = c(0, max(2, quantile(tib$minlog10qval, probs = .995, na.rm = T))) # hardcoded limit; y-axis data goes to 10^-2 at least
    
    tib$isoutlier_x_low = tib$foldchange.log2 < xmax[1]
    tib$isoutlier_x_high = tib$foldchange.log2 > xmax[2]
    tib$isoutlier_y = tib$minlog10qval > ymax[2]
    
    tib$x_outlier = tib$foldchange.log2
    tib$y_outlier = tib$minlog10qval
    tib$x_outlier[tib$isoutlier_x_low] = xmax[1]
    tib$x_outlier[tib$isoutlier_x_high] = xmax[2]
    tib$y_outlier[tib$isoutlier_y] = ymax[2]
    
    tib$updown_outlier = tib$updown
    tib$updown_outlier[tib$updown == "unchanged" & (tib$isoutlier_x_low | tib$isoutlier_x_high | tib$isoutlier_y)] = "unchanged_outlier"
    tib$updown_outlier[tib$updown == "down" & (tib$isoutlier_x_low | tib$isoutlier_y)] = "down_outlier"
    tib$updown_outlier[tib$updown == "up" & (tib$isoutlier_x_high | tib$isoutlier_y)] = "up_outlier"
    
    ### construct facets
    
    plottype_labels = c(asis = "data as-is, no labels",
                        asis_lab = paste("data as-is, label", lbl_style),
                        lim = "limited x- and y-axis, no labels",
                        lim_lab = paste("limited x- and y-axis, label", lbl_style))
    
    tib_facets = bind_rows(tib %>% select(label, x=foldchange.log2, y=minlog10qval, pch=updown, flag_plot_label) %>% mutate(flag_plot_label=FALSE) %>% add_column(plottype = "asis"),
                           tib %>% select(label, x=foldchange.log2, y=minlog10qval, pch=updown, flag_plot_label)                                   %>% add_column(plottype = "asis_lab"),
                           tib %>% select(label, x=x_outlier, y=y_outlier, pch=updown_outlier, flag_plot_label) %>% mutate(flag_plot_label=FALSE)  %>% add_column(plottype = "lim"),
                           tib %>% select(label, x=x_outlier, y=y_outlier, pch=updown_outlier, flag_plot_label)                                    %>% add_column(plottype = "lim_lab"))
    tib_facets$pch = factor(tib_facets$pch, levels = c("up", "down", "up_outlier", "down_outlier", "unchanged", "unchanged_outlier"))
    
    # some mock data to enforce symmetric x-axis and expand the y-limit a bit to make room for the labels (this is a workaround because we cannot hardcode separate y-axis limits per facet)
    blank_data = bind_rows(tibble(x=xmax_nooutlier, y=ymax_nooutlier[2] * 1.2, plottype = "asis"),
                           tibble(x=xmax_nooutlier, y=ymax_nooutlier[2] * 1.2, plottype = "asis_lab"),
                           tibble(x=xmax, y=ymax[2] * 1.2, plottype = "lim"),
                           tibble(x=xmax, y=ymax[2] * 1.2, plottype = "lim_lab") ) %>%
      add_column(pch="unchanged", label="")
    
    ### volcano plot
    p = ggplot(tib_facets, aes(x, y, colour = pch, fill = pch, shape = pch, label = label)) +
      geom_point(na.rm = TRUE, show.legend = TRUE) + # show legend is needed since ggplot 3.5.0 , i.e. `drop=FALSE` in the scale is now ignored... https://github.com/tidyverse/ggplot2/issues/5728
      geom_blank(mapping = aes(x, y), data = blank_data, inherit.aes = FALSE, show.legend = FALSE) +
      scale_discrete_manual(
        aesthetics = c("colour", "fill"),
        values = c(up = "#d55e00aa", down = "#56b4e9aa", up_outlier = "#d55e00aa", down_outlier = "#56b4e9aa", unchanged = "#22222299", unchanged_outlier = "#22222299"),
        labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
                   unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
        breaks = c("up", "down", "up_outlier", "down_outlier", "unchanged", "unchanged_outlier"),
        drop = F
      ) +
      scale_shape_manual(
        values = c(up = 24, down = 25, up_outlier = 2, down_outlier = 6, unchanged = 19, unchanged_outlier = 1),
        labels = c(up = "up regulated", down = "down regulated", up_outlier = "up regulated & outside plot limits", down_outlier = "down regulated & outside plot limits",
                   unchanged = "not significant", unchanged_outlier = "not significant & outside plot limits"),
        breaks = c("up", "down", "up_outlier", "down_outlier", "unchanged", "unchanged_outlier"),
        drop = F
      ) +
      guides(colour = guide_legend(byrow = FALSE, override.aes = list(alpha = 1))) + # don't have to repease the guide_legend for fill and shape
      facet_wrap(~plottype, nrow = 2, ncol = 2, scales = "free", labeller = labeller(plottype=plottype_labels)) +
      labs(x = "log2 fold-change", y = "-log10 FDR adjusted p-value", title = paste(algo_name, "@", mtitle)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=8),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size=8))
    
    if(any(tib_facets$flag_plot_label)) {
      if(label_avoid_overlap) {
        p = p + ggrepel::geom_text_repel(data = tib_facets %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                                         show.legend = FALSE, size = 2, segment.alpha = .3, min.segment.length = 0, na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, box.padding = 0.2, seed = 123) # min.segment.length = unit(0.2, 'lines')
      } else {
        p = p + geom_text(alpha=1, data = tib_facets %>% filter(flag_plot_label == TRUE), vjust = 1.1, show.legend = FALSE, size = 2)
      }
    }
    
    if(!is.na(log2foldchange_threshold)) {
      p = p + geom_vline(xintercept = c(-1,1) * log2foldchange_threshold, colour = "darkgrey", linetype = "dashed")
    }
    if(!is.na(qvalue_threshold)) {
      p = p + geom_hline(yintercept = -log10(qvalue_threshold), colour = "darkgrey", linetype = "dashed")
    }
    
    result[[algo_name]] = list(ggplot = p, ggplot_data = tib_facets)
    
    if(show_plots) {
      print(p)
    }
  }
  
  return(invisible(result))
}
ggplot_coefficient_of_variation = function(tib_input, samples, samples_colors) {
  
  ## for CoV computation, we need natural log while intensities are log2; log_b(x) = log_d(x) / log_d(b)  @  https://www.purplemath.com/modules/logrules5.htm
  # toy example; given y = log2(x=100), we need z = log10(x);
  # x = 100
  # y = log2(x)
  # z = y / log2(10)
  # x;y;z
  
  # !! here we use the by-group filtering, and normalization, as configured by the user
  # this is important, because minpep, topN and different normalization make it different from leave-one-out CoV (there, we must use 'mode' norm for speed as we have to normalize the dataset as often as there are samples)
  tibw_abundance_naturallog = tib_input %>%
    select(peptide_id, sample_id, intensity_by_group) %>% # technically, we don't need this, but here to explicitly state input data for now. can this comment out
    filter(!is.na(intensity_by_group)) %>% # this drops all peptides not passing the filter in any sample/group, which makes downstream wide format intensity matrix much smaller
    mutate(intensity_by_group = log2_to_ln(intensity_by_group)) %>%
    pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_by_group)
  
  
  ## compute stats per group
  groups = unique(samples$group)
  mat_cov = matrix(NA, nrow=nrow(tibw_abundance_naturallog), ncol=length(groups), dimnames = list(tibw_abundance_naturallog$peptide_id, groups))
  dropcols = NULL
  for(grp in groups) {
    # samples for current group
    sid = intersect(samples %>% filter(group == grp) %>% pull(sample_id),
                    colnames(tibw_abundance_naturallog))
    # skip if <n samples in group
    if(length(sid) < 3) {
      dropcols = c(dropcols, grp)
      append_log(sprintf("no CoV computation for sample group '%s', require at least 3 replicates", grp), type = "info")
      next
    }
    
    # fast CoV computation
    m = as.matrix(tibw_abundance_naturallog %>% select(!!sid))
    rows_fail = matrixStats::rowSums2(!is.na(m)) < 3
    # less than 50 peptides have a value, not a meaningful set of datapoints for CoV analysis
    if(sum(!rows_fail) < 50) {
      dropcols = c(dropcols, sid_exclude)
      next
    }
    # m[rows_fail, ] = NA # remove rows with less than 3 values (can technically calculate on 2 values, but we chose to require at least 3)
    mat_cov[!rows_fail,grp] = coefficient_of_variation_vectorized(m[!rows_fail,])
  }
  
  # remove sample groups that have too few replicates
  mat_cov = mat_cov[ , !(colnames(mat_cov) %in% dropcols), drop=F]
  if(ncol(mat_cov) == 0) {
    append_log("No data available for CoV computation, skipping plots", type = "info")
    return(list())
  }
  
  # cov as percentages in all downstream analyses
  mat_cov = mat_cov * 100
  # debug/QC: bp=boxplot(mat_cov, outline=F, las=2); boxplot_add_text(bp, cex=.5)
  
  
  ## summary stats: boxplot. we plot these as text labels onto the ggplot downstream
  bp = boxplot(mat_cov, plot = F)
  cov_data_summ = data.frame(bp$stats)
  colnames(cov_data_summ) = bp$names
  # We split the text labels for the CoV boxplot (all minus last  vs  last) so we can apply different vertical justification
  cov_data_summ_long_14 = cov_data_summ[c(1,4), , drop = F] %>% gather(group, summ) %>% rename(cov = summ)
  cov_data_summ_long_23 = cov_data_summ[c(2,3), , drop = F] %>% gather(group, summ) %>% rename(cov = summ)
  cov_data_summ_long_5 = cov_data_summ[5,,drop = F] %>% gather(group, summ) %>% rename(cov = summ)
  
  ## summary stats: number of peptides used in CoV computation per group
  tib_cov_group_n = tibble::enframe(colSums(!is.na(mat_cov)), name="group", value = "n")
  # tib_cov_group_n$cov_median = as.numeric(cov_data_summ[3, match(colnames(cov_data_summ), tib_cov_group_n$group)]) # coordinates for plotting at boxplot median line
  
  
  # adjust text size to the number of groups
  n_groups = ncol(mat_cov)
  txt_size = 4 - 2*min(10,n_groups) / 10
  
  # prepare plot data in long format
  tib_plot = matrix_to_long(mat_cov, value_name = "cov", column_name = "group", row_name = "peptide_id") %>%
    mutate(group = factor(group, levels = colnames(mat_cov)))
  
  ## boxplot
  p_boxplot = ggplot(tib_plot, aes(group, cov, fill = group)) +
    geom_boxplot(outlier.shape = NA, na.rm = T) +
    geom_text(data = cov_data_summ_long_14, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = -.5, nudge_x = .05, size = txt_size) +
    geom_text(data = cov_data_summ_long_23, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0.5, vjust = -.5, size = txt_size) +
    geom_text(data = cov_data_summ_long_5, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = 1, nudge_x = .05, size = txt_size) +
    geom_text(data = tib_cov_group_n, aes(x = group, y = -4, label = n), hjust = 0.5, vjust = 0.5, size = txt_size * 0.8) +
    # geom_text(data = cov_data_summ_long_2, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = .5, nudge_x = .05, size = txt_size) +
    # geom_text(data = tib_cov_group_n, aes(x = group, y = cov_median, label = paste0("n=",n)), hjust = 1, vjust = -.5, nudge_x = -.05, size = txt_size) +
    scale_color_manual(values = array(unique(samples_colors$group), dimnames = list(unique(samples$group))), aesthetics = c("fill")) +
    coord_cartesian(ylim = c(-5, max(cov_data_summ)), clip = 'off') +
    scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(cov_data_summ, na.rm=T) / 10) * 10, by = 10)) +
    labs(x = "", y = "Coefficient of Variation (%)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=ifelse(n_groups<8, 11, 8)))
  
  
  ## violin plot
  p_violin = ggplot(tib_plot, aes(group, cov, fill = group)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size = .7, na.rm = T) +
    scale_color_manual(values = array(unique(samples_colors$group), dimnames = list(unique(samples$group))), aesthetics = c("fill")) +
    ylim(c(0, max(cov_data_summ))) +
    labs(x = "", y = "Coefficient of Variation (%)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=ifelse(n_groups<8, 11, 8)))
  
  
  return(list(violin = p_violin, boxplot = p_boxplot))
}
plot_differential_detect = function(dataset, zscore_threshold = 6) {
  result = list()
  if(!"dd_proteins" %in% names(dataset) || nrow(dataset$dd_proteins) == 0) {
    return(result)
  }
  
  for(contr in dataset$contrasts) {
    tib_contr = dataset$dd_proteins %>% filter(contrast == contr$label & is.finite(zscore))
    if(nrow(tib_contr) == 0) {
      next
    }
    
    tib_contr = tib_contr %>% filter(type %in% c("detect", "quant"))
    tib_contr_detect = tib_contr %>% filter(type == "detect")
    tib_contr_quant = tib_contr %>% filter(type == "quant")
    lbl_detect = sprintf("only detected peptides; #proteins tested: %d  #abs(zscore) >= %s: %d", nrow(tib_contr_detect), as.character(zscore_threshold), sum(abs(tib_contr_detect$zscore) >= zscore_threshold))
    lbl_quant = sprintf("all quantified peptides; #proteins tested: %d  #abs(zscore) >= %s: %d", nrow(tib_contr_quant), as.character(zscore_threshold), sum(abs(tib_contr_quant$zscore) >= zscore_threshold))
    
    tib_contr = tib_contr %>%
      arrange(type) %>%
      mutate(type_label = ifelse(type == "detect", lbl_detect, lbl_quant),
             type_label = factor(type_label, levels = unique(type_label)))
    
    p_hist = ggplot(tib_contr, aes(zscore)) +
      geom_histogram(bins=25, boundary = 0, colour = "black", fill="lightgrey", na.rm=T) +
      geom_vline(xintercept = c(-zscore_threshold, zscore_threshold), colour = "red") +
      facet_wrap(.~type_label, ncol = 1, scales = "free") +
      labs(x="Differential z-score for observed peptides", y="Protein count", colour = "", title = paste("contrast:", contr$label_contrast)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=9),
            plot.subtitle = element_text(hjust = 0.5, size=9),
            legend.position = "none")
    
    result[[contr$label]] = p_hist
  }
  
  return(result)
}

############# Plot identification and quantification counts ####################
dataset$samples <- dataset$samples %>%  
  mutate(shortname = recode(shortname, 
                            "Efl_Env_1" = "O-Efl-1", 
                            "Efl_Env_2" = "O-Efl-2",
                            "Efl_Env_3"="O-Efl-3",
                            "Efl_Jov_4"="Y-Efl-1",
                            "Efl_Jov_5"="Y-Efl-2") )

ggplot_sample_detect_counts_barplots(
  dataset$samples, 
  samples_colors_long=c("peptides: identif & quant"="#0570b0", 
                        "peptides: only quant"="#74a9cf",
                        "proteins: identif & quant"="#E80909", 
                        "proteins: only quant"="#ED8E8E")) +
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12) )

############# plot the data completeness #######################################

ggplot_peptide_detect_frequency(dataset$peptides, dataset$samples) + 
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12) )


############# Plot differential detection results as a histogram ################
plot_differential_detect(dataset, zscore_threshold = 6)


############ Plot the volcano plot manually ######################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))
tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Efl_Env vs Efl_Jov # condition_variable: group") %>%
  filter(dea_algorithm == "msempire") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = minlog10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(protein_id, gene_symbols,foldchange.log2, foldchange.log2_abs, minlog10qval, signif)
tib <- separate(tib, gene_symbols, into = c("gene_symbols", "G2"), sep = ";")

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.8, size=2) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = 1.301, colour = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0.86, colour = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = -0.86, colour = "darkgrey", linetype = "dashed") +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "#d55e00aa", "down-regulated" = "#56b4e9aa", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change") + ylab("-log10 FDR adj. p-value") + ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=16, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y = 12, label = "Efl_Env", hjust = 0.8, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 3, y = 12, label = "Efl_Jov",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  annotate("text", x = 0, y = 12, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="black") 



############ plot pca of peptides ################################
# peptide level data matrix, using filtered+normalized peptides across all groups
tibw_noexclude = dataset$peptides %>% 
  select(peptide_id, sample_id, intensity_all_group) %>%
  filter(!is.na(intensity_all_group)) %>%
  pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_all_group)
matrix_sample_intensities = msdap:::as_matrix_except_first_column(tibw_noexclude)

# compute PCA
PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
mat_pca = PPCA@scores # PCA coordinate matrix; rownames = sample_id and colnames = PCA dimension
pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores))) # variation per dimension

# if you only want a table with the PCA coordinates, just print the mat_pca variable and you're done at this point

# example plot code
dims = c(1,2) # plot dimensions. in this example, can be; 1:2, 2:3, c(1,3)
p = ggplot(data.frame(x=mat_pca[,dims[1]], y=mat_pca[,dims[2]], label=rownames(mat_pca)), 
           aes(x = x, y = y, label = label)) +
  geom_point() +
  geom_text() +
  labs(x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100) ) +
  theme_bw()
print(p)
