regfac.expand.1par <- function(beta, X, y, fbase) {
  # obtain base distribution derivatives
  ret <- fbase(X%*%beta,y)
  # expand base derivatives
  f <- sum(ret$f)
  g <- t(X)%*%ret$g
	xtw <- 0*X
	for (k in 1:ncol(X)) xtw[,k] <- X[,k]*ret$h # TODO: convert for loop to sapply
	h <- t(xtw)%*%X
  
  return (list(f=f, g=g, h=h))
}

regfac.expand.2par <- function(coeff, X, Z=matrix(1.0, nrow=nrow(X), ncol=1), y, fbase, block.diag=FALSE) {
  # extracting coefficients of X and Z
  beta <- coeff[1:ncol(X)]
  gamma <- coeff[(ncol(X)+1):(ncol(X)+ncol(Z))]
  
  # obtain base distribution derivatives
  ret <- fbase(X%*%beta,Z%*%gamma,y)

  # expand base derivatives; TODO: vectorize for loops
  # function
  f <- sum(ret$f)
  # first derivative
  g <- c(t(X)%*%ret$g[,1],t(Z)%*%ret$g[,2])
  # second derivative
	h <- array(0, dim=c(ncol(X)+ncol(Z),ncol(X)+ncol(Z)))
	# XX block
	xtw <- 0*X
	for (k in 1:ncol(X)) xtw[,k] <- X[,k]*ret$h[,1]
	h[1:ncol(X),1:ncol(X)] <- t(xtw)%*%X
	# ZZ block
	ztw <- 0*Z
	for (k in 1:ncol(Z)) ztw[,k] <- Z[,k]*ret$h[,2]
	h[(ncol(X)+1):(ncol(X)+ncol(Z)),(ncol(X)+1):(ncol(X)+ncol(Z))] <- t(ztw)%*%Z
	# XZ and ZX blocks
	if (!block.diag) {
	  ztw2 <- 0*Z
	  for (k in 1:ncol(Z)) ztw2[,k] <- Z[,k]*ret$h[,3]
	  h[(ncol(X)+1):(ncol(X)+ncol(Z)),1:ncol(X)] <- t(ztw2)%*%X
	  h[1:ncol(X),(ncol(X)+1):(ncol(X)+ncol(Z))] <- t(h[(ncol(X)+1):(ncol(X)+ncol(Z)),1:ncol(X)])
	}
	
  return (list(f=f,g=g,h=h))
}




