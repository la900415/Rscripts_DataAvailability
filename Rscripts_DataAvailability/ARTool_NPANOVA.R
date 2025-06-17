# ARTool: R Package for the Aligned Rank Transform for Nonparametric Factorial ANOVAs 
# https://github.com/mjskay/ARTool/

library(ARTool)
data(Higgins1990Table5, package = "ARTool")


df <- pivot_longer(pdata@main %>% cbind(pdata@annotCols$protein_id),
                   cols = c(1:ncol(pdata@main)), 
                   names_to = c("sample"),
                   values_to = "abundance"
                   )
annotRows <- pdata@annotRows[, -c(2:3, )]
annotRows <- annotRows %>% rownames_to_column("sample")
df <- as.data.frame(df)
df <- merge(df, annotRows, by.x = "sample")
colnames(df)[colnames(df) == "pdata@annotCols$protein_id"] <- "protein_id"
df$sample <- as.factor(df$sample)
df$protein_id <- as.factor(df$protein_id)
saveRDS(df, "df.RDS")
#remove column protein_id

m <- art(DryMatter ~ Moisture*Fertilizer + (1|Tray), data=Higgins1990Table5)
m <- art(DryMatter ~ Moisture*Fertilizer, data=Higgins1990Table5)
summary(m)
m_anova <- anova(m)

summary.art(m)
