logdet =function(A){
  eigtem=eigen(A,only.values=TRUE,symmetric=TRUE)$values
  y=sum(log(eigtem[eigtem > 1e-10])) #[eigtem > 1e-10]
  return(y)
}

