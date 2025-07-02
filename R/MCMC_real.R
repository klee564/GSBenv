
#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' params <- generate_GSpar(r,p,u,G)
#' ex_groupind=c(rep(G+1,u),params$groupind)
#' inputdata <- generate_data(params$param,100)
#'
#' K <- max(ex_groupind)
#' ex_group_prior <- list(ex_groupind=ex_groupind,
#' ex_Betashapes = matrix(1/1000,K,2))
#' X=inputdata$X;Y=inputdata$Y
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' init <- initial_value(X=X,Y=Y,u=u)
#' res <- preproc_MCMC(X,Y,u,ex_group_prior=ex_group_prior,ex_hyperlist=ex_hyperlist,init=init)
#'
#'
preproc_MCMC <- function(X,Y,u,ex_hyperlist,ex_group_prior,init=NULL,permres,GEidx){

  r <- ncol(Y)
  SS <- SS_hyper(X,Y,ex_hyperlist$kappa,ex_hyperlist$M)
  #if(is.null(permres)){
  #  permres <- perm_init(X,Y,u)
  #}
  permres <- perm_init(X,Y,u,init=init)
  GEidx <- permres$GEidx
  permY <- Y[,GEidx]
  #init <- Renvlp::env(X=X,Y=permY,u = u)
  init <- env.high(X=X,Y=permY,u = u)

  ex_groupind <- ex_group_prior$ex_groupind
  groupind <- factor(ex_groupind[GEidx][-(1:u)],
                     level=  unique(ex_groupind[GEidx][-(1:u)])) %>% as.integer
  ex_Betashapes <- ex_group_prior$ex_Betashapes

  subind <- unique(ex_groupind[GEidx][-(1:u)])
  Betashapes <- ex_Betashapes[subind,,drop=FALSE]
  for(i in 1:length(subind)){
    if(is.element(subind[i], ex_groupind[GEidx][1:u])){
      Betashapes[i,1] <- Betashapes[i,1]+1
    }
  }

  #Betashapes <- ex_Betashapes[unique(ex_groupind[GEidx][-(1:u)]),,drop=FALSE]

  group_prior <- list(groupind=groupind,Betashapes=Betashapes)

  hyperlist <- ex_hyperlist
  hyperlist$Psi <- hyperlist$Psi[1:u,1:u]
  hyperlist$Psi0 <- hyperlist$Psi0[1:(r-u),1:(r-u)]
  hyperlist$nu = u +1
  hyperlist$nu0 = r-u +1


  list(GEidx = GEidx, SS=SS, init=init,group_prior=group_prior,hyperlist=hyperlist,
       permres=permres)
  #group_prior
  #hyperlist
}


MCMC_A <- function(SS,init,tau0,tau1,hyperlist,group_prior,mcmc.num=200,permres){

  postlist <- list()

  r <- nrow(init$beta)
  p <- ncol(init$beta)
  n <- init$n
  u <- nrow(init$eta)
  groupindex <- group_prior$groupind
  K <- max(groupindex)
  phi <- rep(1/2,K)

  #set hyperparameter
  M <- hyperlist$M
  B0 <- hyperlist$B0
  Psi <- hyperlist$Psi
  Psi0 <- hyperlist$Psi0
  nu <- hyperlist$nu
  nu0 <- hyperlist$nu0
  kappa <- hyperlist$kappa

  Betashapes <- group_prior$Betashapes

  #set init
  #initA <- as.matrix(find_A_from_gamma(init$Gamma),col=u)
  initA <- permres$A
  gamma_gamma0 <- find_gammas_from_A(initA)

  A <- initA
  xi <- rep(0,r-u)
  rw_var <- rep(1,dim(A)[1])

  postlist <- list()
  for(iter in 1:mcmc.num){
    tausqvec <- xi *tau1^2 + (1-xi) * tau0^2

    #gen A
    lpd_val <- lpd_A_pred(n,A,SS$Omegafactor,SS$Omega0factor,Psi,Psi0,nu,nu0,tausqvec)
    MHsample <- rwmh_rowwise(lpd_val, rw_var,n,A,SS$Omegafactor,SS$Omega0factor,Psi,Psi0,nu,nu0,tausqvec)

    #tune_int <- iter-as.integer(mcmc.num/4)
    tune_int <- iter
    if(tune_int>0) rw_var <- rw_var * exp(tune_int^(-0.7)*(MHsample$alphas-0.44))
    A <- MHsample$A


    for(k in 1:K){

      xig_star <- rbinom(sum(groupindex==k),1,1/2)
      xig_pre <- xi[groupindex==k]


      tausqvec_star <- xig_star *tau1^2 + (1-xig_star) * tau0^2
      tausqvec_pre <- xig_pre *tau1^2 + (1-xig_pre) * tau0^2

      lpd_xi_star <- 1:length(tausqvec_star) %>%
        purrr::map(~mvnfast::dmvn(A[which(groupindex==k)[.x],],mu=rep(0,u),
                                  sigma = diag(tausqvec_star[.x],u),log = T)) %>% unlist %>%sum

      lpd_xi_pre <- 1:length(tausqvec_pre) %>%
        purrr::map(~mvnfast::dmvn(A[which(groupindex==k)[.x],],mu=rep(0,u),
                                  sigma = diag(tausqvec_pre[.x],u),log = T)) %>% unlist %>% sum

      lpd_xi_star <- lpd_xi_star +
        extraDistr::dbbinom(sum(xig_star),length(xig_star),
                            Betashapes[k,1],Betashapes[k,2],log = TRUE)
      lpd_xi_pre <- lpd_xi_pre +
        extraDistr::dbbinom(sum(xig_pre),length(xig_pre),
                            Betashapes[k,1],Betashapes[k,2],log = TRUE)

      aprob <- min(exp(lpd_xi_star-lpd_xi_pre),1)
      if (runif(1) < aprob) {
        #  # accept
        xi[groupindex==k] <- xig_star
      }
    }


    #for(j in 1:(r-u)){
    #  logpropratio <- u*log(tau0/tau1) + LaplacesDemon::logit(max(min(0.999,phi[groupindex[j]]),0.0001))+
        #logpropratio <- u*log(tau0/tau1) + max(min(LaplacesDemon::logit(phi[groupindex[j]]),10^(30)),-10^(30))+
    #    sum(A[j,]^2)/(2*tau0^2)*(1- (tau0/tau1)^2)
    #  propxi <- LaplacesDemon::invlogit(logpropratio)
    #  xi[j] <- rbinom(1,1,propxi)
    #}
    for(k in 1:K){
      phi[k] <- rbeta(1,Betashapes[k,1]+sum(xi[groupindex==k]),
                      Betashapes[k,2]+sum(1-xi[groupindex==k]))
    }
    postlist[[iter]] <- list(A=A,xi=xi,phi=phi,alphas=MHsample$alphas)
  }
  return(postlist)
}




load_ADNI <- function(){
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


  MARRY <- read.csv("Data/RID_MARRY.csv")[,-1] %>%
    purrr::set_names(c("RID","marry")) %>%
    fastDummies::dummy_cols(select_columns = "marry") %>%
    dplyr::select(-c("marry"))

  ##filtering
  megaDB <- ROIunq %>% dplyr::left_join(demo,by="RID") %>%
    dplyr::left_join(SNP_PC,by="RID") %>%
    #dplyr::left_join(pc5,by="RID") %>%
    dplyr::left_join(MARRY,by="RID") %>%
    dplyr::filter(handedness >0 & education >0) %>%
    dplyr::filter(marry_5==0) %>%
    na.omit



  #input data with scaling
  Y <- log((megaDB %>% dplyr::select(ROIname$name))/megaDB$ICV) %>% scale
  #X <- megaDB %>% dplyr::select(c(names(demo)[2:6],paste0("PC",1:5),names(SNP_PC)[2:6])) %>% scale
  #X <- megaDB %>% dplyr::select(c(names(demo)[2:6],names(SNP_PC)[2:6])) %>% scale
  X <- megaDB %>% dplyr::select(c(names(demo)[2:6],names(MARRY)[3:5],names(SNP_PC)[-1])) %>% scale
  #X <- megaDB %>% dplyr::select(c(names(demo)[2:6])) %>% scale

  list(megaDB=megaDB,X=X,Y=Y,ROIname=ROIname)
}


ADNI_hypergroup <- function(ADNIdata,alpha=0.1,beta=0.1){
  X<- ADNIdata$X
  Y<- ADNIdata$Y
  ROIname <- ADNIdata$ROIname
  p <- ncol(X)
  r <- ncol(Y)

  group_list <- split(ROIname$ind,ROIname$region)
  ex_groupind <- rep(0,length(ROIname$ind))

  for(g in 1:length(group_list)){
    ex_groupind[group_list[[g]]] <- g
  }


  ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
                       Psi0 = diag(1,r),kappa = 10^(-6))
  G <- max(ex_groupind)
  ex_Betashapes <- matrix(c(alpha,beta),G,2,byrow=TRUE)
  ex_group_prior <- list(ex_groupind=ex_groupind,ex_Betashapes=ex_Betashapes)

  list(ex_group_prior=ex_group_prior,ex_hyperlist=ex_hyperlist)

}
