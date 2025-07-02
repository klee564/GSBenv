
#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' uvec <- 1:3
#' params <- generate_GSpar(r,p,u,G)
#' inputdata <- generate_data(params$param,50)
#' X <- inputdata$X
#' Y <- inputdata$Y
#' sCVfun(X,Y,u)
#'
sCVfun <- function(X,Y,u){

  n <- nrow(X)
  validx <- caret::createMultiFolds(1:n,k=10,times=2)

  library(furrr)
  plan(multisession, workers = 10)

  trainfit <- validx %>%
    furrr::future_map(~spenv(X=X[.x,],Y=Y[.x,],u=u))


  purrr::map(1:length(trainfit),
             ~stesterr(trainfit[[.x]],newX=X[-validx[[.x]],],newY=Y[-validx[[.x]],])) %>%
    unlist %>% mean

}

stesterr <- function(resenv,newX,newY){

  mean((newY - rep(1,nrow(newX)) %*% t(resenv$alpha) - newX %*% t(resenv$beta))^2)
}



#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' uvec <- 2:3
#' params <- generate_GSpar(r,p,u,G)
#' inputdata <- generate_data(params$param,50)
#' X <- inputdata$X
#' Y <- inputdata$Y
#' zz <- wrap_sCV(X=X,Y=Y,uvec = uvec)
#'
wrap_sCV <- function(X,Y,uvec){
  library(magrittr)
  tdf <- expand.grid(u=uvec)

  value <- 1:nrow(tdf) %>%
    purrr::map(~sCVfun(X=X,Y=Y,u=tdf$u[.x])) %>% unlist

  data.frame(tdf,value=value)

}


