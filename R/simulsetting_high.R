
#'
#' @examples
#' r <- 100
#' u <- 2
#' p <- 6
#' G <- 10
#' params <- generate_GSpar_high(r,p,u,G)
#'
generate_GSpar_high <- function(r,p,u,G){

  ind_partition <- caret::createFolds(1:(r-u),k = G)
  xi.tru <- rep(0,r-u)
  xi.tru[ind_partition[[1]]] <- 1

  groupind <- rep(0,r-u)
  for(g in 1:G){
    groupind[ind_partition[[g]]] <- g
  }

  param <- generate_sparsepar_high(r,p,u, xi.tru)

  list(param=param,groupind=groupind)
}


#'
#' @export
#'
#' @examples
#' r <- 5
#' u <- 2
#' p <- 6
#' xi.tru <- sample(c(0,1),r-u,replace=TRUE)
#'
#' paramlist <- generate_sparsepar(r,p,u, xi.tru)
#' paramlist$A.tru
#'
#'
generate_sparsepar_high <- function(r, p, u, xi.tru) {

  mu.tru <- runif(r, 0, 0)
  eta.tru <- matrix(runif(u * p, min = 0.1, max = 1), nrow = u, ncol = p)
  A <- A.tru <- matrix(runif(u * (r - u), min =0.2, max = 1.5), nrow = r-u, ncol = u)*
    matrix(xi.tru,r-u,u)
  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0
  CA <- gamma_gamma0$CA

  # CA <- rbind(diag(1, u), A)
  # DA <- rbind(-t(A), diag(1, r-u))
  #
  # calculate beta
  beta.tru <- gamma %*% eta.tru
  #
  # CAtCA <- t(CA) %*% CA
  # gamma <- CA %*% sqrtmatinv(CAtCA)
  #
  # DAtDA <- t(DA) %*% DA
  # gamma0 <- DA %*% sqrtmatinv(DAtDA)

  # omega.tru <- runif(0.5, 0.8)
  # Omega.tru <- rinvwish(dim = u, Phi = diag(0.8, u), nu = 1)
  # Omega0.tru <- rinvwish(dim = r-u, Phi = diag(250, r-u), nu = 1)
  # Omega.tru <- diag(sort(runif(u, 0, 1), decreasing = TRUE),
  #                   ncol = u, nrow = u)
  # Omega0.tru <- diag(sort(runif(r-u, 2, 10), decreasing = TRUE),
  #                    ncol = r-u, nrow = r-u)

  V.tru <- diag(0.04,ncol = u, nrow = u)
  V0.tru <- diag(25,ncol = r-u, nrow = r-u)

  R.tru <- diag(1,u)
  R0.tru <- diag(1,r-u)

  Omega.tru <- sqrt(V.tru) %*% R.tru %*% sqrt(V.tru)
  Omega0.tru <- sqrt(V0.tru) %*% R0.tru %*% sqrt(V0.tru)


  #Omega.tru <- diag(sort(runif(u, 0,1), decreasing = TRUE),
  #                  ncol = u, nrow = u)
  #Omega0.tru <- diag(sort(runif(r-u, 5, 10), decreasing = TRUE),
  #                   ncol = r-u, nrow = r-u)


  Sigma1 <- gamma %*% Omega.tru %*% t(gamma)
  Sigma2 <- gamma0 %*% Omega0.tru %*% t(gamma0)

  Sigma.tru <- Sigma1 + Sigma2

  list(mu.tru = mu.tru,
       eta.tru=eta.tru,
       beta.tru = beta.tru,
       Omega.tru = Omega.tru,
       Omega0.tru = Omega0.tru,
       A.tru = A.tru,
       gamma.tru = gamma,
       gamma0.tru = gamma0,
       Sigma.tru = Sigma.tru,
       xi.tru = xi.tru,
       mux = 0,
       sigmax = 1)
}

