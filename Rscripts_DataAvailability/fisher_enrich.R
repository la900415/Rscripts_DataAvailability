library(tidyverse)
library(readxl)

setwd("C:/Users/Luis Ariel/Downloads/Laura_regions")

TR_FISHER <- read.delim2("C:/Users/Luis Ariel/Downloads/Laura_regions/TR_FISHER.txt")
CB_FISHER <- read.delim2("C:/Users/Luis Ariel/Downloads/Laura_regions/CB_FISHER.txt")
HC_FISHER <- read.delim2("C:/Users/Luis Ariel/Downloads/Laura_regions/HC_FISHER.txt")
CX_FISHER <- read.delim2("C:/Users/Luis Ariel/Downloads/Laura_regions/CX_FISHER.txt")

write_xlsx(TR_FISHER, "TR_FISHER.xlsx")
write_xlsx(CB_FISHER, "CB_FISHER.xlsx")
write_xlsx(HC_FISHER, "HC_FISHER.xlsx")
write_xlsx(CX_FISHER, "CX_FISHER.xlsx")

Mixed <- read_excel("Mixed.xlsx")

Mixed %>%  
  ggplot(aes(Sample, Category.value) ) +
  geom_point(aes(size=Intersection.size, color=FDR)) +
  scale_color_continuous(low="blue", high="red") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold") ) 

##### v02 Figure all regions ###########

BS <- read_delim("enrichment.Process (2)TR.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
HT <- read_delim("enrichment.Process (2)HT.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
HC <- read_delim("enrichment.Process (2)HC.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
DI <- read_delim("enrichment.Process (2)DI.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
CX <- read_delim("enrichment.Process (2)CX.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
CB <- read_delim("enrichment.Process (2)CB.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
BL <- read_delim("enrichment.Process (2)BL.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)


BS2 <- BS %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="BS")
BL2 <- BL %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="BL")
CB2 <- CB %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="CB")
CX2 <- CX %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="CX")
DI2 <- DI %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="DI")
HC2 <- HC %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="HC")
HT2 <- HT %>% arrange(`false discovery rate`) %>% slice(1:5) %>% mutate(sample="HT")

BS2 <- as.data.frame(BS2)
CB2 <- as.data.frame(CB2)
BL2 <- as.data.frame(BL2)
CX2 <- as.data.frame(CX2)
DI2 <- as.data.frame(DI2)
HC2 <- as.data.frame(HC2)
HT2 <- as.data.frame(HT2)

Mixed <- rbind(BS2, BL2, CB2, 
               CX2, 
               DI2, 
               HC2 
               #HT2
               )
#Export and manually combine
writexl::write_xlsx(Mixed, "Mixed.xlsx")
writexl::write_xlsx(HT2, "HT2.xlsx")

#7 regions imported
Mixed <- read_excel("Mixed.xlsx")
Mixed <- as.data.frame(Mixed)

#Filtering with final 4 regions
Mixed <- Mixed %>% filter(sample %in% c("BS", "CB", "CX", "HC"))

# 
Mixed$term.description <- str_wrap(Mixed$term.description, width = 40)

Mixed %>%  
  ggplot(aes(sample, fct_reorder(term.description, count))) +
  geom_point(aes(size=count, color=FDR)) +
  scale_color_continuous(low="blue", high="red") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = .5),
        axis.text.y = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold") ) 









