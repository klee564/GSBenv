# sechol.R
#
# Schnabel-Eskow generalized cholesky.
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2004-6  Jeff Gill, Micah Altman
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

######################################################
#       
# sechol
#       
# Schnabel-Eskow generalized cholesky.
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

"sechol" <- function(A, tol = .Machine$double.eps, silent= TRUE )  {
  if (is.complex(A))  {
    warning("complex matrices not permitted at present")
    return(NULL)
  } else if (!is.numeric(A))  {
    warning("non-numeric argument to 'sechol'")
    return(NULL)
  }
  if (is.matrix(A)) {
    if (nrow(A) != ncol(A)) {
      warning("non-square matrix in 'sechol'")
      return(NULL)
    }
    if (nrow(A)==1&&A>0) return(as.matrix(sqrt(A)))
  } else {
    if (length(A) != 1) {
      warning("non-matrix argument to 'sechol'")
      return(NULL)
    }
    if (A>0) {
      return(as.matrix(sqrt(A)))
    } 
    warning("the leading minor of order 1 is not positive definite")
    return(NULL)
  }
  n <- nrow(A)
  L <- matrix(rep(0,n*n),ncol=ncol(A))
  tau <- tol ^(1/3)  # made to match gauss
  gamm <- max(A)
  deltaprev <- 0
  Pprod <- diag(n)
  if (n > 2)  {
    for (k in 1:(n-2))  {
      if( (min(diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm) 
          && (min(svd(A[(k+1):n,(k+1):n])$d)) < 0) {
        dmax <- order(diag(A[k:n,k:n]))[(n-(k-1))]
        if (A[(k+dmax-1),(k+dmax-1)] > A[k,k])  {
          if (!silent) {
            print(paste("iteration:",k,"pivot on:",dmax,"with absolute:",(k+dmax-1)))
          }
          P <- diag(n)
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        g <- rep(0,length=(n-(k-1)))
        for (i in k:n)  {
          if (i == 1) sum1 <- 0
          else sum1 <- sum(abs(A[i,k:(i-1)]))
          if (i == n) sum2 <- 0
          else sum2 <- sum(abs(A[(i+1):n,i]))
          g[i-(k-1)] <- A[i,i] - sum1 - sum2
        }
        gmax <- order(g)[length(g)]
        if (gmax != k)  {
          if (!silent) {
            print(paste("iteration:",k,
                        "gerschgorin pivot on:",gmax,"with absolute:",(k+gmax-1)))
          }
          P <- diag(ncol(A))
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        normj <- sum(abs(A[(k+1):n,k]))
        delta <- max(0,deltaprev,-A[k,k]+max(normj,tau*gamm))
        if (delta > 0)  {
          A[k,k] <- A[k,k] + delta
          deltaprev <- delta
        }
      }
      L[k,k] <- A[k,k] <- sqrt(A[k,k])
      for (i in (k+1):n)  {
        L[i,k] <- A[i,k] <- A[i,k]/L[k,k]
        A[i,(k+1):i] <- A[i,(k+1):i] - L[i,k]*L[(k+1):i,k]
        if(A[i,i] < 0) A[i,i] <- 0
      }
    }
  }
  A[(n-1),n] <- A[n,(n-1)]
  eigvals <- eigen(A[(n-1):n,(n-1):n])$values
  delta <- max(0,deltaprev,
               -min(eigvals)+tau*max((1/(1-tau))*(max(eigvals)-min(eigvals)),gamm))
  if (delta > 0)  {
    if (!silent) {
      print(paste("delta:",delta))
    }
    A[(n-1),(n-1)] <- A[(n-1),(n-1)] + delta
    A[n,n] <- A[n,n] + delta
    deltaprev <- delta
  }
  L[(n-1),(n-1)] <- A[(n-1),(n-1)] <- sqrt(A[(n-1),(n-1)])
  L[n,(n-1)] <- A[n,(n-1)] <- A[n,(n-1)]/L[(n-1),(n-1)]
  L[n,n] <- A[n,n] <- sqrt(A[n,n] - L[n,(n-1)]^2)
  
  r = t(Pprod)%*%t(L)%*%t(Pprod)
  attr(r,"delta")=delta
  return(r)
}

######################################################
#       
# secholTest
#       
# [Internal Function]
#
# self test of sechol, sanity checks
#       
# Parameters:
#
# silent - print debugging output
# 
######################################################


"secholTest"<-function(silent=TRUE) {
  rv = TRUE
  # non singular
  S <- matrix(c(2,0,2.4,0,2,0,2.4,0,3),ncol=3)
  rv = (sum( signif(chol(S),digits=14) == signif(sechol(S),digits=14)) ==9)
  if (!rv && !silent) {
    warning("sechol alters PD matrix")
  }
  S <- matrix(c(2,0,10,0,2,0,10,0,3),ncol=3)
  
  if (!is.R()){
    isTRUE<-function (x) identical(TRUE, x)
  }
  
  t = isTRUE(
    all.equal(matrix(signif(sechol(S), digits=10),ncol=3), 
              matrix(c(1.414213562,0,0,0,1.414234971,0,7.071067812,0,0.007781680058), ncol=3)))
  
  if (!t && !silent) {
    warning("sechol results don't match benchmark")
  }
  rv= rv && t
  return(rv)
}








