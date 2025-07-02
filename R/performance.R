

#'
#' @examples
#' groupind <- params$groupind
#' xi <- params$param$xi.tru
#' xi_gxi(xi,groupind)
#'
xi_to_gxi <- function(xi,groupind){
  G <- max(groupind)
  gxi <- rep(0,G)

  for(g in 1:G){
    gxi[g] <- max(xi[groupind==g])
    #gxi[g] <- mean(xi[groupind==g])
    #gxi[g] <- as.integer(mean(xi[groupind==g])==1)
  }
  gxi
}


#'
#'
#' @examples
#' r <- 10
#' u <- 2
#' p <- 6
#' G <- 4
#' params <- generate_GSpar(r,p,u,G)
#' perf_sel(params$param$xi.tru,pmax(1,params$param$xi.tru),params$groupind)
#'
perf_sel <- function(xi.tru, xi.est,groupind){
  library(magrittr)
  measure <- ROCit::measureit(score = xi.est, class = xi.tru,
                       measure = c("ACC", "TPR","TNR"))

  Gmeasure <- ROCit::measureit(score = round(xi_to_gxi(xi.est,groupind)),
                               class = round(xi_to_gxi(xi.tru,groupind)),
                              measure = c("ACC", "TPR","TNR"))
  data.frame(value= c(measure$TPR[2],measure$TNR[2],measure$ACC[2],
                      Gmeasure$TPR[2],Gmeasure$TNR[2],Gmeasure$ACC[2]),
             measure= c("TPR","TNR","ACC","GTPR","GTNR","GACC"))
}


#'
#' @examples
#'
#'
perf_Bsenv <- function(postmean,params){
  G <- max(params$groupind)
  u <- ncol(params$param$A.tru)

  xi.tru <- c(rep(1,u),params$param$xi.tru)
  xi.est <- round(postmean$xi)

  ex_groupind=c(rep(G+1,u),params$groupind)

  rbind(data.frame(value=c(norm(postmean$beta - params$param$beta.tru,type="F")/norm(params$param$beta.tru,type="F"),
                           norm(postmean$betaxi - params$param$beta.tru,type="F")/norm(params$param$beta.tru,type="F")),
                   measure=c("betaFnorm","betaFnorm2")),perf_sel(xi.tru, xi.est,ex_groupind))
}


#'
#' @examples
#'
#'
perf_senv <- function(senv.out,params){
  G <- max(params$groupind)
  u <- ncol(params$param$A.tru)

  xi.tru <- c(rep(1,u),params$param$xi.tru)
  xi.est <- senv.out$xi

  ex_groupind=c(rep(G+1,u),params$groupind)

  rbind(data.frame(value=norm(senv.out$beta - params$param$beta.tru,type="F")/norm(params$param$beta.tru,type="F"),
                   measure="betaFnorm"),perf_sel(xi.tru, xi.est,ex_groupind))
}

