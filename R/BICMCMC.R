lpdf <- function(A, Ymat, Xmat,eta,mu,Omega0,Omega) {

  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0

  Sigma <- gamma %*% Omega %*% t(gamma) +
    gamma0 %*% Omega0 %*% t(gamma0)
  beta <- gamma %*% eta

  sum(mvnfast::dmvn(Ymat-t(beta %*% t(Xmat)), mu= mu,sigma = Sigma,log=TRUE))
}


BIC_MCMC <- function(mcmc.list,Xmat,Ymat){

  r <- ncol(Ymat)
  p <- ncol(Xmat)
  u <- ncol(mcmc.list[[1]]$A)
  s <- mcmc.list %>% purrr::map(~.x$xi) %>% do.call("rbind",.) %>%
    colMeans %>% round %>% sum
  #logLuhat <- purrr::map(mcmc.list,~lpdf(.x$A * (.x$xi %*% t(rep(1,u))),Ymat,Xmat,.x$eta,.x$mu,
  logLuhat <- purrr::map(mcmc.list,~lpdf(.x$A ,Ymat,Xmat,.x$eta,.x$mu,
                             .x$Omega0,.x$Omega)) %>% unlist %>% max
  #s <- r-u

  -2*logLuhat + log(nrow(Ymat)) * numberofparam(r,p,u,s)
}
