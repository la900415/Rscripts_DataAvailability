################################################################################
#########   Protein rank script                                 ################
#########                                                       ################
################################################################################
library(ggplot2)
library(ggrepel)
library(readxl)

Prot_rank <- read_excel("Protein rank plot_astrocytes_3248quantified.xlsx", sheet = "Protein rank plot_astrocytes_32")

Prot_rank[ ] |> 
  ggplot( aes(x=Rank, 
              y=log10_intensity, label=Astrocyte.marker) ) +
  geom_point(aes(color = ifelse(Astrocyte.marker == "", "grey", "red")) ) +
  theme_classic( ) +
  theme(panel.grid.major.y = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank() ) +
  ylim(5, 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_blank( ),
        legend.title = element_blank( ),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text = element_text(size = 12) ) +
  scale_size_continuous(range = c(3, 7) )  +
  labs(title = NULL, 
       subtitle = NULL, 
       caption = NULL, 
       tag = NULL, 
       x= "Proteins rank", 
       y= "log10 (LFQ)" ) +
  geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.75, max.overlaps = 10 ) 
