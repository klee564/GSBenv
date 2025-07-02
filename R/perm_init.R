


#'
#' @export
#'
#' @example
#' X <- data$X; Y <- data$Y
#' permres <- perm_init(X,Y,u)
#'
#'
perm_init <- function(X, Y, u, lambda1=NULL,lambda2=NULL, ftol = 1e-4, maxiter = 1e2,verbose=0,init) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  r = ncol(Y)
  if(missing(lambda1)) lambda1 <- exp(seq(log(1),log(1e-5),len = 15))
  if(missing(lambda2)) lambda2 <- exp(seq(log(1),log(1e-5),len = 15))
  if(missing(init))  init =  initial_value(X,Y,u)

  GEidx = GE(init)
  newY = Y[, GEidx]
  newinit = init[GEidx,,drop = FALSE]

  m1 <- LassoLambda.spenv(X, newY, u, lambda=lambda1, ftol = ftol,
                                   maxiter = maxiter, weight = rep(1,r-u),
                                   init = newinit,verbose=verbose)
  Gammahat <- m1$Gamma
  w <- Gammahat %*% solve(Gammahat[1:u, ])
  w_norm <- 1/(rowSums(w^2)^2)[(u+1):r]

  m2 <- LassoLambda.spenv(X, newY, u, lambda=lambda2, ftol = ftol,
                                   maxiter = maxiter, weight = w_norm,
                                   init = newinit,verbose=verbose)

  xi <- rep(0,r-u)
  xi[setdiff(m2$where1,1:u)-u] <- 1

  list(GEidx = GEidx,X=X,Y=newY,A=find_A_from_gamma(m2$Gamma),xi=xi)
}
