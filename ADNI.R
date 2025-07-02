library(dplyr)
ROIname <- readxl::read_xlsx("Data/ROIname.xlsx",
                             col_names=c(paste0("V",1:5),"name","ind")) %>%
  dplyr::mutate(LR = substr(name,nchar(name),nchar(name)),
                region = substr(name,1,nchar(name)-2))
groupind <- split(ROIname$ind,ROIname$region)
# corp.c        pec       prec

ROIunq <- read.csv("Data/ROIunq.csv",header = FALSE) %>%
  purrr::set_names(c("RID",ROIname$name))
demo <- read.csv("Data/demounq1.csv",header = FALSE) %>%
  purrr::set_names(c("RID","gender","age","handedness","education","apoe","ICV"))
SNP <- read.csv("Data/SNPunq.csv",header = FALSE)
pc5 <- read.csv("Data/5pc.csv") %>%
  dplyr::select(-one_of(c("X","ID")))



##filtering
megaDB <- ROIunq %>% dplyr::left_join(demo,by="RID") %>%
  dplyr::left_join(SNP_PC,by="RID") %>%
  dplyr::left_join(pc5,by="RID") %>%
  dplyr::filter(handedness >0 & education >0)



#input data with scaling
Y <- log((megaDB %>% dplyr::select(ROIname$name))/megaDB$ICV) %>% scale
X <- megaDB %>% dplyr::select(c(names(demo)[2:6],names(SNP_PC)[-1])) %>% scale
#X <- megaDB %>% dplyr::select(c(names(demo)[2:6])) %>% scale

u <- 4
spest <- Bsenvlp::wrap_spenv(X,Y,u)
betahat <- spest$beta
sort(colnames(Y)[spest$fullxi==1])

envest <- Renvlp::env(X,Y,u)


permres <- Bsenvlp::perm_init(X,Y,u)

wrapmcmc <- function(geneps,groupbeta,mcmc.num=100){
  b0 <- 1
  zeta <- 1e-6
  geneps <- 1e-6
  prior_beta <- purrr::rerun(length(groupind),c(ak=groupbeta,
                                                bk=groupbeta))
  MC_Bsenv2(permres,Y,X,u,groupind,prior_beta,geneps,
                         zeta=zeta,b0= b0,mcmc.num)
}

genepsvec <- c(1e-4,1e-8,1e-12,1e-20,1e-30)
groupbetavec <- c(1e-4,1e-8,1e-12,1e-20,1e-30)

griddf <- expand.grid(geneps = genepsvec,groupbeta=groupbetavec)

mcmc.outlist <- list()
for(i in 2:nrow(griddf)){
  print(i)
  mcmc.outlist[[i]] <- wrapmcmc(griddf$geneps[i],griddf$groupbeta[i])
}

tout <- wrapmcmc(1e-08,1e-400,mcmc.num=20)

tout$fullxi.list %>% purrr::map(~length(sort(colnames(Y)[.x==1]))) %>% unlist
sort(colnames(Y)[tout$fullxi.list[[20]]==1])


zztbl <- mcmc.outlist %>% purrr::map(function(mcmc.out){
  mcmc.out$fullxi.list %>% purrr::map(~length(sort(colnames(Y)[.x==1]))) %>% unlist
}) %>% do.call("rbind",.)

image(zztbl)
summary(zztbl)

#mcmc.out1$fullxi.list[1:10] %>% purrr::map(~sort(colnames(Y)[.x==1]))



saveRDS(mcmc.out1,"mcmc.out1.rds")
geneps <- 1e-12
groupbeta <- 1e-20
prior_beta <- purrr::rerun(length(groupind),c(ak=groupbeta,
                                              bk=groupbeta))

mcmc.out2 <- MC_Bsenv(Y,X,u,groupind,prior_beta,geneps,
                     zeta=zeta,b0= b0,mcmc.num)
saveRDS(mcmc.out2,"mcmc.out2.rds")




mcmcxi1 <- do.call("rbind",mcmc.out1$fullxi.list[1001:2000]) %>% colMeans %>% round
mcmcxi2 <- do.call("rbind",mcmc.out2$fullxi.list[1001:2000]) %>% colMeans %>% round

sort(colnames(Y)[spest$fullxi==1])
sort(colnames(Y)[mcmcxi1==1])
sort(colnames(Y)[mcmcxi2==1])

#########################################
library(ggplot2)
betahat <- spest$beta
sort(colnames(Y)[spest$fullxi==1])

#betahat <- mcmc.out2$beta.list[[2000]]


ylevel_key <- colnames(Y)
names(ylevel_key) <- 1:length(ylevel_key)
xlevel_key <- colnames(X)
names(xlevel_key) <- 1:length(xlevel_key)

gg <- Matrix::summary(Matrix::Matrix(betahat,sparse=TRUE)) %>%
  dplyr::mutate(i=factor(dplyr::recode_factor(i,!!!ylevel_key),levels = sort(ylevel_key)),
                j=factor(dplyr::recode_factor(j,!!!xlevel_key),levels = xlevel_key))

ggplot(gg, aes(j, i,fill = x)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_discrete(drop=F) +
  scale_x_discrete(drop=F)


#########################################
