# grams  Gram-Schmidt orthogonalization of the columns of A.
# The columns of A are assumed to be linearly independent.
#
# Q = grams(A) returns an m by n matrix Q whose columns are 
# an orthonormal basis for the column space of A.
#
# [Q, R] = grams(A) returns a matrix Q with orthonormal columns
# and an invertible upper triangular matrix R so that A = Q * R.
#
# Warning: For a more stable algorithm, use [Q, R] = qr(A, 0) .

grams =function(A){
  A=as.matrix(A)
  m=nrow(A)
  n=ncol(A)
  if(n==1) {R=sqrt(sum(A^2)); return(Q=A/R)} #return(list(Q=A/R,R=R))
  Asave = A
  Q=A
  eps=2.2204e-16
  for (j in 2:n){
    for (k in 1:(j-1)){    
      mult = (t(A[, j]) %*% A[,k]) / (sum(A[,k]^2))
      A[,j] = A[,j]-mult*A[,k];
    }
  }
  for (j in 1:n){
    if (sum(A[,j]^2)<sqrt(eps))  stop('Columns of A are linearly dependent.')
    Q[,j] = A[,j]/ sqrt(sum(A[,j]^2))
  }
  R = t(Q) %*% Asave
  return(Q)
  #return(list(Q=Q,R=R))
}