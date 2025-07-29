# Validity of ANOVA (normality & homogeneity of variance)
library(tidyverse)
library(PerseusR)
library(readxl)
library(readr)

# Pre-processing of data
setwd("E:/RESEARCH/ygor/DN/2_tecido_FPspeclib_BEST/msdap_results/2025-03-16_17-11-57_vwmb_BEST")
dir.create("normality_test")

Tecid <- read_delim("protein_abundance__global data filter.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

samples <- read_excel("samples.xlsx")
samples <- samples[samples$exclude == "FALSE",]
#change the symbol "+" in samples$shortname by "_"
samples$shortname <- gsub("\\+", "_", samples$shortname)

T_tidy <- Tecid %>% 
  pivot_longer(cols = 4:23, names_to = "sample", values_to = "abundance")

T_tidy <- T_tidy %>% 
  left_join(samples[,2:7], by = c("sample" = "sample_id"))
saveRDS(T_tidy, "normality_test/T_tidy.RDS")

######## Normality (nortest package) ##############################################################
library(nortest)
library(dplyr)
library(broom)

ks <- 
  T_tidy %>%
  group_by(shortname, group, age, maneuver) %>% 
  do(tidy(lillie.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))
writexl::write_xlsx(ks, "normality_test/nortest_lillie_ks.xlsx")
print(ks)

ad <- 
  T_tidy %>%
  group_by(shortname, group, age, maneuver) %>% 
  do(tidy(ad.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))
writexl::write_xlsx(ad, "normality_test/nortest_ad.xlsx")
print(ad)

cvm <- 
  T_tidy %>%
  group_by(shortname, group, age, maneuver) %>% 
  do(tidy(cvm.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))
writexl::write_xlsx(cvm, "normality_test/nortest_cvm.xlsx")
print(cvm)

##### Kolmogorov-Smirnov normality test (p>0.05 normal distribution) #########################################################

"E_Efl-J_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_2"]), "pnorm", 
                       mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_2"])), 
                       sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_2"] ) ) ) #p.value=0.176
"E_Efl-J_3" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_3"]), "pnorm",
                       mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_3"])), 
                       sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_3"] ) ) ) #p.value=0.282
"E_Efl-J_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_4"]), "pnorm",
                       mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_4"])), 
                       sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_Efl-J_4"] ) ) ) #p.value=0.204
"E_IR_1" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_1"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_1"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_1"] ) ) ) #p.value=0.217
"E_IR_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_2"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_2"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_2"] ) ) ) #p.value=0.3093
"E_IR_3" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_3"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_3"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_3"] ) ) ) #p.value=0.1784
"E_IR_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_4"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_4"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_IR_4"] ) ) ) #p.value=0.2144
"E_PCI_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_2"]), "pnorm",
                     mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_2"])), 
                     sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_2"] ) ) ) #p.value=0.1813
"E_PCI_3" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_3"]), "pnorm",
                     mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_3"])), 
                     sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_3"] ) ) ) #p.value=0.4262
"E_PCI_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_4"]), "pnorm",
                     mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_4"])), 
                     sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="E_PCI_4"] ) ) ) #p.value=0.2054
"J_Efl-E_1" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_1"]), "pnorm",
                      mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_1"])), 
                      sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_1"] ) ) ) #p.value=0.6186
"J_Efl-E_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_2"]), "pnorm",
                      mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_2"])), 
                      sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_2"] ) ) ) #p.value=0.2455
"J_Efl-E_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_4"]), "pnorm",
                      mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_4"])), 
                      sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_Efl-E_4"] ) ) ) #p.value=0.7603
"J_IR_1" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_1"]), "pnorm",
                   mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_1"])), 
                   sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_1"] ) ) ) #p.value=0.06144
"J_IR_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_2"]), "pnorm",
                   mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_2"])), 
                   sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_2"] ) ) ) #p.value=0.7452
"J_IR_3" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_3"]), "pnorm",
                   mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_3"])), 
                   sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_3"] ) ) ) #p.value=0.04332
"J_IR_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_4"]), "pnorm",
                   mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_4"])), 
                   sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_IR_4"] ) ) ) #p.value=0.05066
"J_PCI_1" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_1"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_1"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_1"] ) ) ) #p.value=0.1212
"J_PCI_2" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_2"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_2"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_2"] ) ) ) #p.value=0.2507
"J_PCI_4" <- ks.test(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_4"]), "pnorm",
                    mean = mean(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_4"])), 
                    sd = sd(unique(T_tidy_na$abundance[T_tidy_na$shortname=="J_PCI_4"] ) ) ) #p.value=0.4562

#### Homogeneity of variance (p>0.05 homogeneity, p<0.05 not homogeneity) ########
library(car)
library(vartest)
library(ggpubr)

ggboxplot(T_tidy, x = "group", y = "abundance", color = "group",palette = "jco")
ggboxplot(T_tidy, x = "age", y = "abundance", color = "age",palette = "jco")
ggboxplot(T_tidy, x = "maneuver", y = "abundance", color = "maneuver",palette = "jco")

# levene's test
T_tidy$group <- as.factor(T_tidy$group)
with(T_tidy, leveneTest(abundance, group)) #Pr(>F) < 2.2e-16

T_tidy$age <- as.factor(T_tidy$age)
with(T_tidy, leveneTest(abundance, age)) #Pr(>F) = 0.5525

T_tidy$maneuver <- as.factor(T_tidy$maneuver)
with(T_tidy, leveneTest(abundance, maneuver)) #Pr(>F) < 2.2e-16

# Variance Ratio Test (F-test)
var.test(abundance ~ group, data=T_tidy) #
var.test(abundance ~ age, data=T_tidy) #p-value = 0.7672
var.test(abundance ~ maneuver, data=T_tidy) #

f.test(abundance ~ group, data=T_tidy) #p.value     : 2.650172e-26 
f.test(abundance ~ age, data=T_tidy) #p.value     : 0.3835801 
f.test(abundance ~ maneuver, data=T_tidy) #p.value     : 2.900948e-24

# Bartlett's test
bartlett.test(abundance ~ group, data=T_tidy) #p-value < 2.2e-16
bartlett.test(abundance ~ age, data=T_tidy) #p-value = 0.7672
bartlett.test(abundance ~ maneuver, data=T_tidy) #p-value < 2.2e-16

# Fligner-Killeen test
fligner.test(abundance ~ group, data=T_tidy) #p-value < 2.2e-16
fligner.test(abundance ~ age, data=T_tidy) #p-value = 0.529
fligner.test(abundance ~ maneuver, data=T_tidy) #p-value < 2.2e-16
