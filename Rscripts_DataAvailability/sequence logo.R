library(ggseqlogo)
library(ggplot2)
library(cowplot)

############ Create a list for sequence logo ########################################
Leu <- list(S23_K43=c("SNYTVGKVGVENLVNAVPQLK", "SNYTAGKVGVENLVNAVPQLK", "SNYTAGKVGVENLVNAVPQLK"),
            G50_K71=c("GEQVVNIGSQDMNDNVWLTLAK", "GEQVVNIGSQDMNDNVWLTLAK", "GEQVVNIGSQDMNDDVWLTLAK"),
            V244_R269=c("VGNGNLYKSVFDTLATAAKTGTAVVR", "VGNGNLYKSVFDTLATAAKNGTAVVR", "VGNGNLYKTVFDTLATAAKNGTAVVR")
            )

ggseqlogo(Leu, method='prob', ncol = 1) +  
  theme_classic() + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.background=element_blank(),
        legend.box.background=element_blank())

cs1 = make_col_scheme(chars=c('V', 'A', 'N', 'D', 'S', 'T'), groups=c('gr1', 'gr2', 'gr1', 'gr2', 'gr1', 'gr2'), 
                      cols=c('purple', 'blue', 'purple', 'blue', 'purple', 'blue'))

p27 = ggseqlogo(Leu$S23_K43, method='prob', col_scheme=NULL, ncol=1) + theme_cowplot() + ggtitle("S23-K43") +
      annotate('rect', xmin = 4.5, xmax = 5.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow') +
      theme(axis.text.x = element_blank())

p64 = ggseqlogo(Leu$G50_K71, method='prob', col_scheme=NULL, ncol=1) + theme_logo() + ggtitle("G50-K71") +
      annotate('rect', xmin=14.5, xmax=15.5, ymin=-0.05, ymax=1.05, alpha=.1, col='black', fill='yellow') 

p252_263 = ggseqlogo(Leu$N244_R269, method='prob', col_scheme=NULL, ncol=1) + theme_logo() + ggtitle("N244-R269") +
           annotate('rect', xmin=8.5, xmax=9.5, ymin=-0.05, ymax=1.05, alpha=.1, col='black', fill='yellow') +
           annotate('rect', xmin=19.5, xmax=20.5, ymin=-0.05, ymax=1.05, alpha=.1, col='black', fill='yellow')

plot_grid(p27, p64, p252_263,  ncol = 1, align = 'h')


aln = data.frame(
  letter=strsplit("SNYTVGKVGVENLVNAVPQLKSNYTAGKVGVENLVNAVPQLKSNYTAGKVGVENLVNAVPQLK", "")[[1]], 
  strain = rep(c("K12", "AS1.357", "BL21"), each=21),
  x       = rep(1:21, 3) )
aln$mut = 'no'
aln$mut[ c(5,26,47) ] = 'yes'

algn = ggplot(aln, aes(x, strain)) +
  geom_text(aes(label=letter, color=mut, size=mut)) + 
  scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') + 
  scale_color_manual(values=c('black', 'red')) + 
  scale_size_manual(values=c(5, 6)) + 
  theme_logo() + 
  theme(legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14)) 

plot_grid(p27, algn,  ncol=1, align='v')
