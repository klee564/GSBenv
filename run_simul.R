

run_simul <- function(params,nvec,uvec=1:4,tau0vec=c(0.01,0.05,0.1,0.2),ncore=4,
                      alpha=10^(-6),beta=10^(-6)){

  library(future)
  plan(multisession, workers = ncore)


  r <- nrow(params$param$beta.tru)
  p <- ncol(params$param$beta.tru)

  ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
                       Psi0 = diag(1,r),kappa = 10^(-6))
  G <- max(params$groupind)
  ex_Betashapes <- matrix(c(alpha,beta),G+1,2,byrow=TRUE)
  ex_group_prior <- list(ex_groupind=c(rep(G+1,u),params$groupind),ex_Betashapes=ex_Betashapes)


  data_list <- nvec %>%
    purrr::map(~generate_data(params$param,.x))

  plan(multisession, workers = ncore)
  Bsenvres <- data_list %>%
    furrr::future_map(~envGS_selu(.x$Y,.x$X,ex_group_prior=ex_group_prior,uvec=uvec,
                                  tau0vec = tau0vec,tau1=10,ex_hyperlist=ex_hyperlist))
  senvres <- data_list %>%
    furrr::future_map(~senv_selu(Y=.x$Y,X=.x$X,uvec=uvec))

  uhat <- c(purrr::map2(senvres,nvec,~data.frame(u=.x$u,n=.y,method="senv")),
        purrr::map2(Bsenvres,nvec,~data.frame(u=.x$u,n=.y,method="Bsenv"))) %>%
    do.call("rbind",.)

  tau0hat <- purrr::map2(Bsenvres,nvec,~data.frame(tau0=.x$tau0,n=.y,method="Bsenv")) %>%
    do.call("rbind",.)

  list(uhat=uhat,Bsenvres=Bsenvres,senvres=senvres,
       tau0hat=tau0hat,params=params,data_list=data_list)
  #senv_summ <- purrr::map2(senvres,nvec,~data.frame(perf_senv(.x$result %>% summary_senv,params),n=.y,method="senv"))
  #Bsenv_summ <- purrr::map2(Bsenvres,nvec,~data.frame(perf_Bsenv(.x$result$envres,params),n=.y,method="Bsenv"))

  #list(uhat=uhat,perf=do.call("rbind",c(senv_summ,Bsenv_summ)),tau0hat=tau0hat,params=params,data_list=data_list)

}





wrap_simulrun <- function(r,p,q){
  nvec <- rep(c(50,100,200),50)
  tau0vec=c(0.01,0.05,0.1)
  uvec <- 1:3
  u <- 2
  params <- generate_qspar(r,p,u,q)
  simulres <- run_simul(params,nvec,uvec=uvec,tau0vec=tau0vec,ncore=10,alpha = 10^(-200),beta = 10^(-200))
  saveRDS(simulres,paste0("simulres",r,"_",p,"_",q,".rds"))
}


wrap_simulrun(10,20,4)
