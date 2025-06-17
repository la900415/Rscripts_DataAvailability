########## Several test for normality ##########
library(nortest)
library(dplyr)
library(broom)

data <- rnorm(100, mean = 5, sd = 3)

# Perform the Anderson-Darling test
result <- ad.test(data)
print(result)

ad <- 
  df %>%
  group_by(age, virus, interaction, shortname) %>% 
  do(tidy(ad.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))


ks <- 
  df %>%
  group_by(age, virus, interaction, shortname) %>% 
  do(tidy(lillie.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))

cvm <- 
  df %>%
  group_by(age, virus, interaction, shortname) %>% 
  do(tidy(cvm.test(.$abundance))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))

library(PerseusR)

pdata <- read.perseus("E:/RESEARCH/astrocytes/session3perseus_20_27samples_with_exp-speclib_MBR_again_25.03.2024_matrix262.txt", check = TRUE, additionalMatrices = FALSE)
pdata2 <- pivot_longer(pdata@main,
                   cols = c(1:18), 
                   names_to = c("sample"),
                   values_to = "abundance" )
pdata2 <- merge(pdata2, pdata@annotRows, by.x = "sample")

ad_pdata <- pdata@main %>%  do(tidy(ad.test(pdata@main))) %>% 
  mutate(Normality = ifelse(p.value > 0.05, "YES", "NO"))


