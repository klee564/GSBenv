library(dplyr)
SNP <- read.csv("Data/SNPunq.csv",header = FALSE)



## Choose the number of PC scores.

dim(SNP)
snpdata <- SNP[,-1]

NC <- dim(snpdata)[2]
spca <- ClassDiscovery::SamplePCA(t(scale(snpdata)))
ag.obj <- PCDimension::AuerGervini(spca)
num_PC <- PCDimension::agDimension(ag.obj) ## chooose p^* using Auer and Gervini's method


## Derive PC scores and construct the covariate matrix.
preParam <- caret::preProcess(snpdata,method=c("pca","center","scale"),
                              pcaComp=num_PC)  ## X : principal component scores


SNP_PC <- data.frame(RID = SNP[,1],
                     data.frame(scale(predict(preParam,snpdata)))%>%
                       purrr::set_names(paste0("SNP_PC",1:num_PC)))
