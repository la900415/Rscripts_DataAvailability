library(usethis)
library(devtools)
install_github("guokai8/rcellmarker")
library(rcellmarker)

library(easybio)
data(pbmc.markers)
view(pbmc.markers)

library(rcellmarker)   
gene=sample(unique(humancells$SYMBOL),20)
data(human) # human cell markers for cell type identification
