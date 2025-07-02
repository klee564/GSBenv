GE <- function(A) {

# Gaussian elimination, p must be less than or equal to n
	a <- dim(A)
	n <- a[1]
	p <- a[2]
	idx <- rep(0, p)
	res.idx <- 1:n
	
	i <- 1
	while (i <= p) {
		tmp <- max(abs(A[res.idx, i]))
		Stmp <- setdiff(which(abs(A[, i]) == tmp), idx)
		idx[i] <- Stmp[1]
		res.idx <- setdiff(res.idx, idx[i])
		for (j in 1:(n-i)) {
			A[res.idx[j], ] <- A[res.idx[j], ] - A[res.idx[j], i] / A[idx[i], i] * A[idx[i], ]
		}
		i <- i + 1			
	}
	c(idx, res.idx)
}