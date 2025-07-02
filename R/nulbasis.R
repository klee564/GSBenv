# nulbasis  Basis for nullspace.
#
# N = nulbasis(A) returns a basis for the nullspace of A
# in the columns of N. The basis contains the n-r special 
# solutions to Ax=0.  freecol is the list of free columns.
#
# Example:
#
#
#A=matrix(c(1,2,0,3,0,0,1,4),2,4,T)
# nulbasis(A)
## See also fourbase.
nulbasis<-function (A){
  temp= rref(A);
  R=temp$A
  pivcol=temp$jb
  m = nrow(A)
  n = ncol(A)
  r = length(pivcol);
  freecol = 1:n;
  freecol=freecol[-pivcol]
  N = matrix(0,n,n-r)
  N[freecol,  ] = diag(n - r)
  N[pivcol, ] = - R[1 : r, freecol]
  return(N)
}
  
  
  
# {
#   A=t(A)
#   tmp <- qr(A)
#   set <- if(tmp$rank == 0) 1:ncol(M) else  - (1:tmp$rank)
#   qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
# }