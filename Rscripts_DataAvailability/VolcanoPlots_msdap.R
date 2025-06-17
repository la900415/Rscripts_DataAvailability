library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
######################## Fig 3. Volcano plots #########################################
load("dataset.RData") # Load the dataset archive in the msdap file just clicking on it
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ctrl vs Gc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = minlog10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(protein_id, gene_symbols,foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=2) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = 2, colour = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0.556, colour = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = -0.556, colour = "darkgrey", linetype = "dashed") +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "#d55e00aa", "down-regulated" = "#56b4e9aa", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change") + ylab("-log10 FDR adj. p-value") + ggtitle("Gc_vs_Ctrl") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"),
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
  annotate("text", x = min(tib$foldchange.log2), y = 33, label = "Ctrl", hjust = 0.8, vjust = 0, size = 4, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = max(tib$foldchange.log2), y = 33, label = "Gc",   hjust = 0.5, vjust = 0, size = 4, fontface="bold", colour="#d55e00aa") +
  annotate("text", x = max(tib$foldchange.log2), y = 35, label = "(601)",   hjust = 0.5, vjust = 0, size = 4, fontface="bold", colour="#d55e00aa") +
  annotate("text", x = min(tib$foldchange.log2), y = 35, label = "(617)",   hjust = 0.5, vjust = 0, size = 4, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 35, label = "(2057)",   hjust = 0.5, vjust = 0, size = 4, fontface="bold", colour="#22222299") 

