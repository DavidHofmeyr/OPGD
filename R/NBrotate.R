OPGD <- function(X, y, V0 = NULL, ndim = NULL, Joint = FALSE, commonS = 0, cv = FALSE, cv_folds = NULL){
  n <- nrow(X)
  d <- ncol(X)
  yy <- numeric(n)
  u <- unique(y)
  for(k in 1:length(u)) yy[which(y==u[k])] = k
  if(is.null(ndim)){
    if(is.null(V0)) ndim <- max(yy)-1
    else ndim <- ncol(V0)
  }
  else if(ndim > d) stop('cannot compute projection with more dimensions than the data')
  PI <- sapply(1:length(u), function(k) sum(yy==k))/n
  MU <- t(sapply(1:max(yy), function(k) colMeans(X[which(yy==k),])))
  if(is.null(V0)){
    M <- sapply(1:max(yy), function(k) yy==k)%*%MU
    Sw <- cov(X-M)
    Sb <- cov(X)-Sw
    V0 <- eigen(solve(Sw + diag(d)*1e-10)%*%Sb + cov(X)*1e-10)$vectors[,1:ndim]
  }
  M <- sapply(1:length(u), function(k) yy==k)%*%MU
  R <- X
  S <- lapply(1:length(u), function(k) (1-commonS)*cov(X[which(yy==k),]) + commonS*cov(R))

  #else S <- lapply(1:length(u), function(k) cov(X[which(yy==k),]))
  tX <- t(X)
  V <- V0
  #V <- optim(V0, fnb_G_2, dfnb_G_2, tX, yy, MU, S, PI, n, Joint, method = 'L-BFGS-B', control = list(fnscale = -1))$par
  if(1){
    for(it in 1:1000){
      fval <- fnb_G_2(V, tX, yy, MU, S, PI, n)
      dir <- dfnb_G_2(V, tX, yy, MU, S, PI, n)
      ndir2 <- sum(dir^2)
      stp <- 1
      repeat{
        fnew <- fnb_G_2(V + stp*dir, tX, yy, MU, S, PI, n)
        if(fnew > (fval + .95*stp*ndir2)) break
        else stp <- stp/2
        if(stp < 1e-9) break
      }
      V <- V + stp*dir
      if(stp < 1e-9) break
      plot(X%*%V, col = yy)
    }
  }
  #V <- sweep(V, 2, apply(V, 2, norm_vec), '/')
  print(fnb_G_2(V, tX, yy, MU, S, PI, n))
  sol <- list(V=V, X = X, y = yy, MU=MU, S=S, PI=PI, nc = max(yy), y_labels = u, method = 'Gaussian')
  class(sol) <- "OPGD"
  sol
}




predict.OPGD <- function(sol, Xte){
  posterior <- matrix(0,nrow(Xte),sol$nc)
  colnames(posterior) <- sol$y_labels
  mu <- sol$MU%*%sol$V
  s <- lapply(sol$S, function(s) diag(t(sol$V)%*%s%*%sol$V))
  for(k in 1:sol$nc) posterior[,k] <- exp(-distmat(Xte%*%t(t(sol$V)/sqrt(s[[k]])),mu[k,]/sqrt(s[[k]]))^2/2)/prod(s[[k]]^.5)*sol$PI[k]
  posterior <- posterior[,order(sol$y_labels)]
  list(class = sort(sol$y_labels)[apply(posterior, 1, which.max)], posterior = posterior/rowSums(posterior))
}
