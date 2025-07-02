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
#' res <- wrap_envGS_MCMC(X,Y,u+1,ex_group_prior=ex_group_prior,
#' tau0=0.05,tau1=1,ex_hyperlist=ex_hyperlist,mcmc.num=20)
#'
wrap_envGS_MCMC <- function(X,Y,u,ex_group_prior=ex_group_prior,
                          tau0,tau1,ex_hyperlist,mcmc.num=20){

  r <- ncol(Y)
  permres <- perm_init(X,Y,u)
  GEidx <- permres$GEidx
  permY <- Y[,GEidx]

  ex_groupind <- ex_group_prior$ex_groupind
  groupind <- factor(ex_groupind[GEidx][-(1:u)],
                     level=  unique(ex_groupind[GEidx][-(1:u)])) %>% as.integer
  ex_Betashapes <- ex_group_prior$ex_Betashapes
  Betashapes <- ex_Betashapes[unique(ex_groupind[GEidx][-(1:u)]),,drop=FALSE]

  hyperlist <- ex_hyperlist
  hyperlist$Psi <- hyperlist$Psi[1:u,1:u]
  hyperlist$Psi0 <- hyperlist$Psi0[1:(r-u),1:(r-u)]
  hyperlist$nu = u +1
  hyperlist$nu0 = r-u +1

  group_prior <- list(groupind=groupind,Betashapes=Betashapes)
  envres <- envGS_MCMC(X,permY,u,group_prior=group_prior,tau0=tau0,tau1=tau1,
                       hyperlist=hyperlist,mcmc.num=mcmc.num)

  invGEidx <- Matrix::invPerm(GEidx)
  envres$xi <- c(rep(1,u),envres$xi)[invGEidx]
  envres$beta <- envres$beta[invGEidx,]
  envres$betaxi <- envres$betaxi[invGEidx,]

  list(GEidx = GEidx,envres=envres )
}
