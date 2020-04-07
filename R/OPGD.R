OPGD <- function(X, y, V0 = NULL, ndim = NULL, lambda = 0){
  n <- nrow(X)
  d <- ncol(X)
  yy <- as.numeric(as.factor(y))
  nc <- max(yy)
  if(is.null(ndim)){
    if(is.null(V0)) ndim <- max(yy)-1
    else ndim <- ncol(V0)
  }
  else if(ndim > d) stop('cannot compute projection with more dimensions than the data')
  PI <- sapply(1:nc, function(k) sum(yy==k))/n
  MU <- t(sapply(1:nc, function(k) colMeans(X[which(yy==k),])))
  M <- sapply(1:nc, function(k) yy==k)%*%MU
  R <- X-M
  S <- lapply(1:nc, function(k) (1-lambda)*cov(X[which(yy==k),]) + lambda*cov(R))
  if(is.null(V0)){
    Sw <- cov(X-M)
    Sb <- cov(X)-Sw
    V0 <- eigen(solve(Sw + diag(d)*1e-5)%*%Sb + cov(X)*1e-5)$vectors[,1:ndim]
  }

  tX <- t(X)
  V <- V0
  V <- optim(V0, fnb_G_2, dfnb_G_2, tX, yy, MU, S, PI, n, method = 'L-BFGS-B', control = list(fnscale = -1))$par
  Vold <- V
  V <- matrix(0, nrow(Vold), 0)
  for(ii in 1:(ndim-1)){
    ix <- which.max(sapply(1:ncol(Vold), function(rm) fnb_G_2(cbind(V, Vold[,rm]), tX, yy, MU, S, PI, n)))
    V <- cbind(V, Vold[,ix])
    Vold <- Vold[,-ix]
  }
  V <- cbind(V, Vold)
  sol <- list(V=V, X = X, y = yy, MU=MU, S=S, PI=PI, nc = max(yy))
  class(sol) <- "OPGD"
  sol
}




predict.OPGD <- function(sol, Xte){
  posterior <- matrix(0,nrow(Xte),sol$nc)
  colnames(posterior) <- sol$y_labels
  mu <- sol$MU%*%sol$V
  s <- lapply(sol$S, function(s) diag(t(sol$V)%*%s%*%sol$V))
  for(k in 1:sol$nc){
    posterior[,k] <- exp(-.Internal(colSums((t(Xte%*%sol$V) - mu[k,])^2/s[[k]]/2, ncol(sol$V), nrow(Xte), TRUE)))/prod(s[[k]]^.5)*sol$PI[k]
  }
  list(class = apply(posterior, 1, which.max), posterior = posterior/rowSums(posterior))
}
