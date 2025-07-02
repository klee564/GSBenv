#RREF   Reduced row echelon form.
#   R = RREF(A) produces the reduced row echelon form of A.
#
#   [R,jb] = RREF(A) also returns a vector, jb, so that:
  #       r = length(jb) is this algorithm's idea of the rank of A,
#       x(jb) are the bound variables in a linear system, Ax = b,
#       A(:,jb) is a basis for the range of A,
#       R(1:r,jb) is the r-by-r identity matrix.
#
#   [R,jb] = RREF(A,TOL) uses the given tolerance in the rank tests.
#
#   Roundoff errors may cause this algorithm to compute a different
#   value for the rank than RANK, ORTH and NULL.
#
#   Class support for input A:
#      float: double, single
#
#   See also RANK, ORTH, NULL, QR, SVD.

#   Copyright 1984-2005 The MathWorks, Inc. 
# A=matrix(c(16,2,3,13,5,11,10,8,9,7,6,12,4,14,15,1),4,4,T)
rref<-function(A){
  m = nrow(A)
  n = ncol(A)
  # Loop over the entire matrix.
  i = 1
  j = 1
  jb = c()
  tol = max(m,n)*2.2204e-016*norm(A,'I')
  while((i<= m)&(j<= n)){
    # Find value and index of largest element in the remainder of column j.
    p = max(abs(A[i:m,j]))
    k = which.max(abs(A[i:m,j]))
    k=k+i-1
    if (p<=tol){
       #The column is negligible, zero it out.
       A[i:m,j] = 0
       j = j + 1
    }
    else{
      # Remember column index
      jb = c(jb,j)
      # Swap i-th and k-th rows.
      A[c(i,k),j:n] = A[c(k,i),j:n]
      # Divide the pivot row by the pivot element.
      A[i,j:n] = A[i,j:n]/A[i,j];
      # Subtract multiples of the pivot row from all the other rows.
      for (k in (1:m)[-i]){
        A[k,j:n] = A[k,j:n] - A[k,j]*A[i,j:n];
      }  
      i = i + 1;
      j = j + 1;
    } #end else
  }#end while  
  result<-list()
  result$A=A
  result$jb=jb
  return(result)
}