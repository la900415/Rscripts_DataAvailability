###### Install and load packages ############################################################
library(eulerr)
library(readxl) 
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)
library(ggplot2) 

###### Load data ##############################################################################
super <- read_excel("Venn diagram HCP.xlsx", 
                    sheet = "super")
pellet <- read_excel("Venn diagram HCP.xlsx", 
                    sheet = "pellet")
Leu_super <- read_excel("Venn diagram HCP.xlsx", 
                        sheet = "Leu-super")
Leu_pellet <- read_excel("Venn diagram HCP.xlsx", 
                         sheet = "Leu-pellet")
all <- read_excel("Venn diagram HCP.xlsx", 
                               sheet = "all")

set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
  )
################################################################################################


##################### ggVennDiagram (venn + upset plots) #######################################
ggVennDiagram(Leu_pellet,
              set_size = 5,
              label = "count",
              label_alpha=0 ) + scale_fill_gradient(low="grey90",high = "red")

ggVennDiagram(Leu_super,
              set_size = 5,
              label = "count",
              label_alpha=0 ) + scale_fill_distiller(palette = "Reds", direction = 1)

ggVennDiagram(super,
              set_size = 6,
              label = "count",
              label_alpha=0 ) + scale_fill_distiller(palette = "Reds", direction = 1)

ggVennDiagram(pellet,
              set_size = 6,
              label = "count",
              label_alpha=0 ) + scale_fill_distiller(palette = "Reds", direction = 1)

###### upset plot #####################################################################
ggVennDiagram(all[],
              force_upset = TRUE,
              relative_height = 2, 
              relative_width = 0.3,
              order.set.by = "name", 
              order.intersect.by = "none") 
###################################################################################################


