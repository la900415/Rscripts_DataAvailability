#### Multivariate Normality Test #################################
library(MVN)

data(iris)
setosa <- iris[1:50, 1:4]

result = mvn(data = iris[-4], subset = "Species", mvnTest = "hz",
             univariateTest = "AD", univariatePlot = "qq",
             multivariatePlot = "qq", multivariateOutlierMethod = "adj",
             showOutliers = TRUE, showNewData = TRUE)
#### Multivariate Normality Result
result$multivariateNormality

### Univariate Normality Result
result$univariateNormality

### Descriptives
result$Descriptives

### Multivariate Outliers
result$multivariateOutliers

### New data without multivariate outliers
result$newData

df_2 <- df[,-1] %>% group_by(age, virus, interaction, shortname) 

df_2_interact = mvn(data = df_2, subset = "shortname", 
             mvnTest = "hz",
             univariateTest = "AD", univariatePlot = "qq",
             multivariatePlot = "qq", multivariateOutlierMethod = "adj",
             showOutliers = TRUE, showNewData = TRUE)
