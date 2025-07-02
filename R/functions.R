
#'
#' @examples
#' r <- 20
#' p <- 10
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' eta <- all_pars$eta.tru
#' Omega <- all_pars$Omega.tru
#' Omega0 <- all_pars$Omega0.tru
#' A <- all_pars$A.tru
#'
param_trans <- function(eta,Omega,Omega0,A){
  gammas <- find_gammas_from_A(A)
  sqrtCAtCA <- sqrtmat(gammas$CAtCA)
  sqrtDAtDA <- sqrtmat(gammas$DAtDA)
  etatilde <- sqrtCAtCA %*% eta
  Omegatilde <- sqrtCAtCA %*% Omega %*% sqrtCAtCA
  Omega0tilde <- sqrtDAtDA %*% Omega0 %*% sqrtDAtDA

  list(etatilde=etatilde,Omegatilde=Omegatilde,Omega0tilde=Omega0tilde)
}

#'
#' @examples
#' r <- 20
#' p <- 10
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' eta <- all_pars$eta.tru
#' Omega <- all_pars$Omega.tru
#' Omega0 <- all_pars$Omega0.tru
#' A <- all_pars$A.tru
#'
param_trans_inv <- function(post){
  gammas <- find_gammas_from_A(post$A)
  sqrtCAtCA_inv <- sqrtmatinv(gammas$CAtCA)
  sqrtDAtDA_inv <- sqrtmatinv(gammas$DAtDA)
  eta <- sqrtCAtCA_inv %*% post$etatilde
  Omega <- sqrtCAtCA_inv %*% post$Omegatilde %*% sqrtCAtCA_inv
  Omega0 <- sqrtDAtDA_inv %*% post$Omega0tilde %*% sqrtDAtDA_inv

  list(eta=eta,Omega=Omega,Omega0=Omega0, beta=gammas$gamma %*%eta )
}

#'
#' @examples
#' r <- 20
#' p <- 10
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' eta <- all_pars$eta.tru
#' Omega <- all_pars$Omega.tru
#' Omega0 <- all_pars$Omega0.tru
#' A <- all_pars$A.tru
#'
param_invtrans <- function(mu,etatilde,Omegatilde,Omega0tilde,A){
  gammas <- find_gammas_from_A(A)
  sqrtCAtCAinv <- sqrtmatinv(gammas$CAtCA)
  sqrtDAtDAinv <- sqrtmatinv(gammas$DAtDA)
  eta <- sqrtCAtCAinv %*% etatilde

  Omega <- sqrtCAtCAinv %*% Omegatilde %*% sqrtCAtCAinv
  Omega0 <- sqrtDAtDAinv %*% Omega0tilde %*% sqrtDAtDAinv

  list(mu=mu,eta=eta,Omega=Omega,Omega0=Omega0,A=A,beta= gammas$gamma %*% eta,
       Sigma = gammas$gamma %*% Omega %*% t(gammas$gamma)+gammas$gamma0 %*% Omega0 %*% t(gammas$gamma0))
}


#functions.R


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- rWishart(1,11,diag(1,10))[,,1]
#' condparam_W(Sigma)
#'
condparam_W <- function(Sigma){
  library(magrittr)
  param <- 2:dim(Sigma)[1] %>%
    purrr::map(function(i){
      coef <- solve(Sigma[1:(i-1),1:(i-1)],Sigma[1:(i-1),i])
      var <-  Sigma[i,i]-sum(coef * Sigma[1:(i-1),i])
      list(var=var,coef=coef)
    }
    )
  c(list(list(var=Sigma[1,1],coef=0)),param)
}


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- rWishart(1,11,diag(1,10))[,,1]
#' yval <- sample(0:1,10,replace=T)
#' lbvec  <- -1/yval +1
#' ubvec <- (1/!yval)-1
#' condparam <- condparam_W(Sigma)
#' gen_W(mu,condparam,lbvec,ubvec)
#'
#'
gen_W <- function(mu,condparam,lbvec,ubvec){
  W <- numeric(length(lbvec))
  W[1] <- extraDistr::rtnorm(1,mean=mu[1],sd=sqrt(condparam[[1]]$var),
                             a=lbvec[1],b=ubvec[1])
  for(i in 2:length(lbvec)){
    W[i] <- extraDistr::rtnorm(1,mean=mu[i] +sum(condparam[[i]]$coef *( W[1:(i-1)]-mu[1:(i-1)])),
                               sd=sqrt(condparam[[i]]$var),a=lbvec[i],b=ubvec[i])

  }
  W
  #  mytruncNormal(mu,Sigma,lbvec,ubvec)
}

#'
#' @export
#'
gen_Wmat <- function(mumat,condparam,lbmat,ubmat){
  1:dim(lbmat)[1] %>%
    purrr::map(~gen_W(mu=mumat[.x,],condparam=condparam,
                      lbvec=lbmat[.x,],ubvec=ubmat[.x,]))%>%
    do.call("rbind",.)
}


#'
#' @export
#'
#' @examples
#' mu <- rnorm(10)
#' Sigma <- diag(runif(10)+1,10)
#' yval <- sample(0:1,10,replace=T)
#' lbvec  <- -1/yval +1
#' ubvec <- (1/!yval)-1
#'
#'
mytruncNormal <- function(mu,Sigma,lbvec,ubvec,Winit=NULL){
  1:length(mu) %>%
    purrr::map(~extraDistr::rtnorm(1,mean = mu[.x],
                                   sd=sqrt(Sigma[.x,.x]),
                                   a = lbvec[.x], b = ubvec[.x])) %>%
    unlist
}


sqrtmatinv <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*% diag(1/sqrt(eig$values), nrow = nrow(mat)) %*% t(eig$vec)

}

sqrtmat <- function(mat) {
  eig <- eigen(mat)
  eig$vec %*%
    diag(sqrt(eig$values), nrow = nrow(mat)) %*%
    t(eig$vec)

}

find_gammas_from_A <- function(A) {
  dims <- dim(A)
  u <- dims[2]
  r <- sum(dims)
  CA <- matrix(0, nrow = r, ncol = u)
  DA <- matrix(0, nrow = r, ncol = r-u)
  CA[(u+1):r, ] <- A
  CA[1:u, 1:u] <- diag(1, u)
  DA[1:u, ] <- -t(A)
  DA[-(1:u), ] <- diag(1, r-u)
  CAtCA <- crossprod(CA)
  DAtDA <- crossprod(DA)
  gamma <- CA %*% sqrtmatinv(CAtCA)
  gamma0 <- DA %*% sqrtmatinv(DAtDA)

  list(gamma = gamma, gamma0 = gamma0,
       CA = CA, CAtCA = CAtCA,
       DA = DA, DAtDA = DAtDA)
}

find_A_from_gamma <- function(gamma) {
  m <- ncol(gamma)
  G1 <- as.matrix(gamma[1:m, ])
  # check if G1 is invertible - else reorganize the predictors
  if (abs(det(G1)) < 1e-7) {
    gamma.t <- t(gamma)
    X.order <- qr(gamma.t, tol = 1e-7)$pivot
    X <- X[, X.order]
    gamma <- gamma[X.order, ]
  }

  G2 <- gamma[-(1:m), ]
  G2 %*% solve(G1)
}




#'
#' @examples
#' coefvec <- c(1,0.9,0.9)
#' m <- 500
#' zz <- mvp(rep("x",3), power=0:2, abs(coefvec)/sum(coefvec))^m
#' zz$coef %>% plot
#' expandpol_fft(abs(coefvec)/sum(coefvec),m) %>% lines
#'
#'
expandpol_fft <- function(coefvec,m){
  Re(fft(fft(c(coefvec,numeric(2*m+1-3)))^m,inverse = T))/(2*m+1)
}



numberofparam <- function(r,p,u,s){
  u*(u+1)/2 + (r-u)*(r-u+1)/2 + r + p*u + u*s
}
