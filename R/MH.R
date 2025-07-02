

#'
#' @export
#' @examples
#'
#' r <- 5;u <- 2;p <- 10;
#' xi.tru <- c(rep(1,2),rep(0,1))
#' groupind <- c(rep(1,2),rep(2,1))
#' n <- 200
#' paramlist <- generate_sparsepar(r,p,u, xi.tru)
#' data <- generate_data(paramlist,n)
#' Xmat <- data$X
#' Ymat <- data$Y
#' kappa <- 0.01
#' M <- diag(0.01,p)
#'
SS_hyper <- function(Xmat,Ymat,kappa,M){
  n <- nrow(Xmat)

  barX <- colMeans(Xmat)
  Xc <- scale(Xmat,scale=FALSE)
  barY <- colMeans(Ymat)
  Yc <- scale(Ymat,scale=FALSE)

  XctXc <- crossprod(Xc)
  barXbarXt <- crossprod(t(barX))
  YctYc <- crossprod(Yc)
  barYbarYt <- crossprod(t(barY))
  XctYc <- crossprod(Xc,Yc)
  barXbarYt <- barX %*% t(barY)

  nkf <- n*kappa/(n+kappa)
  CXY <- solve(XctXc + nkf * barXbarXt + M,
               XctYc + nkf * barXbarYt )

  Omegafactor <- crossprod(Xc %*% CXY - Yc) + t(CXY) %*% M %*% CXY +
    nkf * tcrossprod(t(CXY) %*% barX - barY)
  Omega0factor <- YctYc + nkf * barYbarYt

  list(Omegafactor=Omegafactor,Omega0factor=Omega0factor,CXY=CXY)
}


#'
#' @export
#' @examples
#'
#' r <- 5;u <- 2;p <- 3;
#' xi.tru <- c(rep(1,2),rep(0,1))
#' groupind <- c(rep(1,2),rep(2,1))
#' n <- 200
#' paramlist <- generate_sparsepar(r,p,u, xi.tru)
#' data <- generate_data(paramlist,n)
#' Xmat <- data$X
#' Ymat <- data$Y
#' kappa <- 0.0001
#' M <- diag(0.0001,p)
#' Psi <- diag(1,u)
#' Psi0 <- diag(1,r-u)
#' nu <- u+1
#' nu0 <- r-u+1
#'
#' tau1 <- 1; tau0 <- 0.1
#' xi <- paramlist$xi.tru
#' tausqvec <- xi *tau1^2 + (1-xi) * tau0^2
#'
#' SS <- SS_hyper(Xmat,Ymat,kappa,M)
#' A <- paramlist$A.tru +matrix(runif((r-u)*u)*0.1,r-u,u)
#'
#' lpd_A_pred(n,A,SS$Omegafactor,SS$Omega0factor,Psi,Psi0,nu,nu0,tausqvec)
#'
lpd_A_pred <- function(n,A, Omegafactor, Omega0factor,Psi,Psi0,nu,nu0,tausqvec){
  u <- ncol(A)

  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0


  lpdf <- -(n+nu)*determinant(t(gamma) %*% Omegafactor %*% gamma + Psi,logarithm = TRUE)$modulus -
    (n+nu0)*determinant(t(gamma0) %*% Omega0factor %*% gamma0 + Psi0,logarithm = TRUE)$modulus

  Aprior <- 1:length(tausqvec) %>%
    purrr::map(~mvnfast::dmvn(A[.x,],mu=rep(0,u),
                              sigma = diag(tausqvec[.x],u),log = T)) %>%
    unlist() %>% sum

  return(lpdf +Aprior)
}


#'
#' @examples
#'
#' r <- 3
#' p <- 4
#' u <- u.tru <- 2
#' all_pars <- generate_par(r, p, u)
#' simuldata <- generate_data(all_pars,50)
#' A0 <- matrix(0,r-u,u)
#' K <- diag(10^4, r-u)
#' L <- diag(10^4, u)
#' M <- diag(0.001,p)
#' K.half.inv <- sqrtmatinv(K)
#' L.half.inv <- sqrtmatinv(L)
#' M.half <- sqrtmat(M)
#'
#' Wmat <- simuldata$Z
#' Xmat <- simuldata$X
#' A <- all_pars$A.tru
#' teta <- all_pars$eta.tru
#' tmu <- all_pars$mu.tru
#' Omega.inv <- solve(all_pars$Omega.tru)
#' Omega0.inv <- solve(all_pars$Omega0.tru)
#' lpd_val <- lpd_A_pred(A,Wmat,Xmat,teta,tmu,Omega0.inv,Omega.inv,K.half.inv, L.half.inv, A0,M.half)
#'
#' rw_var <- rep(0.4,dim(A)[2])
#' rw_var <- c(0.6,0.3)
#' rwmh_colwise(lpd_val, rw_var,A,Wmat, Xmat,teta,tmu,Omega0.inv,
#' Omega.inv,K.half.inv, L.half.inv, A0,M.half)
#'
#'
rwmh_rowwise <- function(lpd_val, rw_var,n,A,...){

  #if(is.null(lpd_val)){
  #  lpd_val <- lpd_fun(A=A,...)
  #}
  u <- dim(A)[2]
  r_u <- dim(A)[1]
  alphas <- numeric(u)

  for(j in sample(1:(r_u))) {
    # if (autotune_size == "single") tau_curr <- tau
    A_j_star <- A[j, ] + rnorm(u, 0, sqrt(rw_var[j]))
    A_star <- A
    A_star[j, ] <- A_j_star
    lpd_A_star <- lpd_A_pred(n=n,A = A_star, ...)
    alphas[j] <- min(exp(lpd_A_star-lpd_val),1)
    if (runif(1) < alphas[j]) {
      # accept
      A <- A_star
      lpd_val <- lpd_A_star
    }
  }

  list(A=A,lpd_val=lpd_val,alphas=alphas)
}
