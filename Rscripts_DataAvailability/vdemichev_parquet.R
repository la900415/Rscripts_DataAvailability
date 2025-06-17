library(arrow)
library(reshape2)

# pivot data frame to obtain features x runs table
pivot <- function(df, file.header, id.header, quantity.header) {
  df <- unique(df[,c(file.header, id.header, quantity.header)])
  melted <- as.data.frame(melt(df, id.vars = c(file.header, id.header), measure.vars = c(quantity.header)))
  melted$value[which(melted$value == 0)] <- NA
  piv <- dcast(melted, as.formula(paste0(id.header,'~',file.header)), value.var = "value") 
  rownames(piv) <- piv[[1]]
  piv[,1] <- NULL
  as.matrix(piv)
}

data <- read_parquet('report.parquet') # put the report file name & location here
df <- data[data$Q.Value <= 0.01 & data$PG.Q.Value <= 0.05 & data$Lib.Q.Value <= 0.01 & data$Lib.PG.Q.Value <= 0.01,] # with MBR
# df <- data[data$Q.Value <= 0.01 & data$PG.Q.Value <= 0.05 & data$Global.Q.Value <= 0.01 & data$Global.PG.Q.Value <= 0.01,] # without MBR

pg <- pivot(df, 'Run', 'Protein.Group', 'PG.MaxLFQ')