

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
#' hypergroup <- list(ex_group_prior=ex_group_prior,ex_hyperlist=ex_hyperlist)
#' uvec <- u; tau0vec <- c(0.01); tau1vec <- 10
#' tt <- run_boot(X,Y,u,hypergroup,tau0vec,tau1vec,N=5)
#' tt2 <- run_boot(X,Y,3,hypergroup,tau0vec,tau1vec,N=5)
#' tt1 <- run_boot(X,Y,1,hypergroup,tau0vec,tau1vec,N=5)
#' eval_boot(tt);
#' eval_boot(tt2)
#' eval_boot(tt1)
#'
run_boot <- function(X,Y,u,hypergroup,tau0vec,tau1vec,N=10){
  preinit <- initial_value(X,Y,u)
  #permres <- perm_init(X,Y,u,init=preinit)
  #GEidx <- permres$GEidx
  #permY <- Y[,GEidx]
  #init <- Renvlp::env(X=X,Y=permY,u = u)
  init <- preinit
  permres <- NULL
  GEidx <- NULL



  restot <- fit_winit(trainX=X,trainY=Y,uvec=u,
                      tau0vec=tau0vec,tau1vec=tau1vec,hypergroup,ncore=1,mcmc.num=2000,
                      init=init,GEidx=GEidx,permres=permres)
  srestot <- senv_selu(Y=Y,X=X,uvec=u,init=preinit)


  result <- list()
  for(i in 1:N){
    bootind_list <- purrr::map(1:10,~ sample(1:nrow(Y),size = as.integer(nrow(Y)*1),replace = T))

    library(future)
    plan(multisession, workers = 10)
    Bsenv_list <- bootind_list %>%
      furrr::future_map(~fit_winit(trainX=X[.x,],trainY=Y[.x,],uvec=u,
                                  tau0vec=tau0vec,tau1vec=tau1vec,hypergroup,ncore=1,mcmc.num=2000,
                                  init=init,GEidx=GEidx,permres=permres))
    plan(multisession, workers = 10)

    senv_list <- bootind_list %>%  furrr::future_map(~senv_selu(Y=Y[.x,],X=X[.x,],uvec=u,init=preinit,
                                                                lambda1=srestot$result$lambda[1],
                                                                lambda2=srestot$result$lambda[2]))

    result[[i]] <- list(bootind_list=bootind_list,Bsenv_list=Bsenv_list,senv_list=senv_list)
  }

  list(resboot = result,restot = restot, srestot=srestot,init=init)
}


eval_boot <- function(reslist){

  bootlist <- reslist$resboot %>% purrr::map(~.x$Bsenv_list) %>% do.call("c",.)
  sbootlist <- reslist$resboot %>% purrr::map(~.x$senv_list) %>% do.call("c",.)
  tres <- bootlist %>%
    purrr::map(~round(as.vector(.x$fitlist[[1]]$xi))) %>% do.call("rbind",.)

  calxi <- function(x){
    vec <- rep(0,length(x$result$alpha))
    vec[x$result$where1] <- 1
    vec
  }
  stres <- sbootlist %>% purrr::map(~calxi(.x)) %>% do.call("rbind",.)


  tot_mean <- reslist$restot$fitlist[[1]]$xi
  bias <- apply(tres, 2,mean) - tot_mean
  vars <-  apply(tres, 2,var)

  sbias <-  apply(stres, 2,mean) - calxi(reslist$srestot)
  svars <-  apply(stres, 2,var)


  list(list(bias=bias,vars=vars,rmse= sqrt(mean(bias^2 + vars))),
  list(sbias=sbias,svars=svars,srmse= sqrt(mean(sbias^2 + svars))))
}


eval_boot_beta <- function(reslist){

  bootlist <- reslist$resboot %>% purrr::map(~.x$Bsenv_list) %>% do.call("c",.)
  sbootlist <- reslist$resboot %>% purrr::map(~.x$senv_list) %>% do.call("c",.)
  tres <- bootlist %>%
    purrr::map(~as.vector(.x$fitlist[[1]]$beta)) %>% do.call("rbind",.)

  stres <- sbootlist %>% purrr::map(~as.vector(.x$result$beta)) %>% do.call("rbind",.)


  tot_mean <- reslist$restot$fitlist[[1]]$beta
  bias <- apply(tres, 2,mean) - as.vector(tot_mean)
  vars <-  apply(tres, 2,var)

  sbias <-  apply(stres, 2,mean) - as.vector(reslist$srestot$result$beta)
  svars <-  apply(stres, 2,var)


  c(sqrt(mean(bias^2 + vars)), sqrt(mean(sbias^2 + svars)))
}
