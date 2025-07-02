
#input : train data, uvec, tau0vec, tau1vec

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
#' trainX=inputdata$X;trainY=inputdata$Y
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' hypergroup <- list(ex_group_prior=ex_group_prior,ex_hyperlist=ex_hyperlist)
#' init <- initial_value(X=trainX,Y=trainY,u=u)
#' uvec <- u; tau0vec <- 0.1; tau1vec <- 10
#'
fit_winit <- function(trainX,trainY,uvec,tau0vec,tau1vec,hypergroup,ncore=1,mcmc.num=100,
                     init,permres,GEidx){

  if(ncore>1){
    library(future)
    plan(multisession, workers = ncore)

    preproc_list <- uvec %>%
      furrr::future_map(~preproc_MCMC(X=trainX,Y=trainY,.x,
                               ex_group_prior=hypergroup$ex_group_prior,
                               ex_hyperlist=hypergroup$ex_hyperlist,
                               init=init,GEidx=GEidx,permres=permres))

  } else{

    preproc_list <- uvec %>%
      purrr::map(~preproc_MCMC(X=trainX,Y=trainY,.x,
                               ex_group_prior=hypergroup$ex_group_prior,
                               ex_hyperlist=hypergroup$ex_hyperlist,
                               init=init,GEidx=GEidx,permres=permres))

  }


  taudf <- expand.grid(uind = 1:length(uvec),tau0=tau0vec,tau1=tau1vec)

  if(ncore>1){
    library(future)
    plan(multisession, workers = ncore)
    MCMCres <- 1:nrow(taudf) %>%
      furrr::future_map(
        function(ind){
          preproc_u <- preproc_list[[taudf$uind[ind]]]
          MCMC_A(preproc_u$SS,preproc_u$init,tau0=taudf$tau0[ind],
                 tau1=taudf$tau1[ind],hyperlist=preproc_u$hyperlist,
                 group_prior=preproc_u$group_prior,mcmc.num=mcmc.num,permres=preproc_u$permres) %>%
            generate_others(X=trainX,Y=trainY,hyperlist=preproc_u$hyperlist)%>%
            c(list(GEidx =preproc_u$GEidx ,u=uvec[ind]))
        }
        )
  } else{
    MCMCres <-  1:nrow(taudf) %>%
      purrr::map(function(ind){
        preproc_u <- preproc_list[[taudf$uind[ind]]]
        MCMC_A(preproc_u$SS,preproc_u$init,tau0=taudf$tau0[ind],
               tau1=taudf$tau1[ind],hyperlist=preproc_u$hyperlist,
               group_prior=preproc_u$group_prior,mcmc.num=mcmc.num,permres=preproc_u$permres) %>%
          generate_others(X=trainX,Y=trainY,hyperlist=preproc_u$hyperlist) %>%
          c(list(GEidx =preproc_u$GEidx ,u=uvec[ind]))
      })


  }

  BICres <- data.frame(taudf,BIC= MCMCres %>% purrr::map(~.x$BIC) %>% unlist)
  fitlist <- MCMCres %>% purrr::map(~summary_MCMC(.x))
  #predperf <- fitlist %>% purrr::map(~norm(testX %*% t(.x$beta) - testY,type="F"))

  list(BICres=BICres,fitlist=fitlist)
  #sel perf

}

summary_MCMC <- function(mcmc.list,u){
  xi.mean <- mcmc.list$postlist %>% purrr::map(~.x$xi) %>% do.call("rbind",.) %>% colMeans
  beta.mean <- mcmc.list$postlist %>%
    purrr::map(~.x$beta) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)

  invGEidx <- Matrix::invPerm(mcmc.list$GEidx)
  xi_extend <- c(rep(1,mcmc.list$u),xi.mean)[invGEidx]
  beta_extend <- beta.mean[invGEidx,]
  list(xi=xi_extend,beta=beta_extend,betaxi =beta_extend*round(xi_extend))
}


