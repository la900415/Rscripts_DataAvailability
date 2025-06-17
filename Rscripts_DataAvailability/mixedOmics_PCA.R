data(nutrimouse)
X <- nutrimouse$gene
MyResult.pca <- pca(X)  # 1 Run the method
plotIndiv(MyResult.pca) # 2 Plot the samples


data(multidrug)
X <- matrix2#multidrug$ABC.trans
dim(X) # Check dimensions of data
tune.pca.multi <- tune.pca(X, ncomp = 10, scale = TRUE)
plot(tune.pca.multi)
final.pca.multi <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)
final.pca.multi$prop_expl_var$X
final.pca.multi$cum.var
# final.pca.multi  # Lists possible outputs
head(selectVar(final.pca.multi, comp = 1)$value)
plotIndiv(final.pca.multi, style = '3d',
          comp = c(1, 2,3),   # Specify components to plot
          ind.names = F, # Show row names of samples
          group = metadata$condition,
          ellipse = T, ellipse.level = 0.95,
          title = 'PCA comp 1 - 2',
          legend = TRUE, legend.title = 'Cell line')

NRP1.scale <- scale(X[, 'NRP1'], center = TRUE, scale = TRUE)
boxplot(NRP1.scale ~
          metadata$condition, col = color.mixo(1:9),
        xlab = 'Cell lines', ylab = 'log2 abundance, scaled',
        par(cex.axis = 1), # Font size
        main = 'NRP1 transporter')

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

matrix2 <- t(matrix)
saveRDS(matrix2, "E:/astrocytes/1_msdap/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01/Fig2_PCA/matrix_transposed.RDS")
