
OPGDcluster <- function(X, k, V0 = NULL, ndim = NULL, omega = NULL, standardise = FALSE){
  if(is.null(ndim)){
    if(k==2) ndim <- 2
    else ndim <- min(k-1, ncol(X))
  }

  if(standardise){
    X <- sweep(X, 2, colMeans(X), '-')
    for(i in 1:ncol(X)) if(sd(X[,i])>1e-50) X[,i] <- X[,i]/sd(X[,i])
  }

  n <- nrow(X)


  if(is.null(omega)) omega <- n

  MU <- matrix(0, k, ncol(X))
  S <- list()
  PI <- numeric(k)

  gmm <- Mclust(X, G = c(k), modelNames = 'VVV')
  if(is.null(gmm$z)){
    km <- kmeans(X, k)
    p <- sapply(1:k, function(kk) km$cluster==kk)
    y00 <- apply(p, 1, which.max)
    for(upd in 1:20){
      pold <- p
      for(kk in 1:k){
        PI[kk] <- sum(p[,kk])/sum(p)
        MU[kk,] <- colSums(X*p[,kk])/sum(p[,kk])
        trX <- sweep(X, 2, MU[kk,], '-')
        S[[kk]] <- t(trX*p[,kk])%*%trX/sum(p[,kk])
      }
      for(kk in 1:k){
        trX <- sweep(X, 2, MU[kk,], '-')
        E <- eigen(S[[kk]])
        E$values <- Re(E$values)
        E$vectors <- Re(E$vectors)
        E$values[E$values<.Machine$doule.xmin] <- Inf
        d2 <- rowSums((trX%*%E$vectors%*%diag(1/E$values^.5))^2)
        p[,kk] <- exp(-d2/2)*PI[kk]/(prod(E$values[which(is.finite(E$values))])+.Machine$double.xmin)^.5 + .Machine$double.xmin
      }
      p <- p/rowSums(p)
      if(max(abs(p-pold))<1e-6) break
    }
  }
  else{
    PI <- gmm$parameters$pro
    MU <- t(gmm$parameters$mean)
    p <- gmm$z
    S <- list()
    for(kk in 1:k) S[[kk]] <- gmm$parameters$variance$sigma[,,kk]
  }

  yy <- apply(p, 1, which.max)

  tX <- t(X)

  if(is.null(V0)){
    Sb <- t(MU*PI)%*%MU
    Sw <- cov(X)-Sb
    V0 <- Re(eigen(solve(Sw + diag(ncol(X))*1e-5)%*%Sb + cov(X)*1e-5)$vectors[,1:ndim])
    V0[,1] = V0[,1]/sqrt(sum(V0[,1]^2))
    for(i in 2:ncol(V0)){
      for(j in 1:(i-1)) V0[,i] = V0[,i]-(V0[,i]%*%V0[,j])[1]*V0[,j]
      V0[,i] = V0[,i]/sqrt(sum(V0[,i]^2))
    }
  }

  y0 <- yy

  V <- optim(V0, fnb_G_3, dfnb_G_3, tX, MU, S, PI, n, omega, method = 'L-BFGS-B', control = list(fnscale = -1))$par

  PIo <- .8*PI + .2/k
  MUo <- .8*MU
  smn <- S[[1]]*PI[1]
  for(kk in 2:k) smn <- smn + S[[kk]]*PI[kk]
  So <- list()
  for(kk in 1:k) So[[kk]] <- .8*S[[kk]]+.2*smn

  ll <- 0
  for(upd in 1:30){
    pold <- p
    p <- matrix(0, nrow(X), k)
    for(kk in 1:k){
      trX <- sweep(X, 2, MU[kk,], '-')%*%V
      s <- sqrt(diag(t(V)%*%S[[kk]]%*%V))
      d2 <- rowSums((trX/matrix(s, nrow(X), ncol(V), byrow = TRUE))^2)
      p[,kk] <- exp(-d2/2)*PI[kk]/(prod(s)+.Machine$double.xmin)
    }
    ll <- sum(log(rowSums(p)))
    p <- p/rowSums(p)
    for(kk in 1:k){
      PI[kk] <- sum(p[,kk])/sum(p)
      MU[kk,] <- colSums(X*p[,kk])/sum(p[,kk])
      trX <- sweep(X, 2, MU[kk,], '-')
      S[[kk]] <- t(trX*p[,kk])%*%trX/sum(p[,kk])
    }
    if(max(abs(p-pold))<1e-6) break
  }
  popt <- p


  llo <- 0
  for(upd in 1:30){
    pold <- p
    p <- matrix(0, nrow(X), k)
    for(kk in 1:k){
      trX <- sweep(X, 2, MUo[kk,], '-')%*%V
      s <- sqrt(diag(t(V)%*%So[[kk]]%*%V))
      d2 <- rowSums((trX/matrix(s, nrow(X), ncol(V), byrow = TRUE))^2)
      p[,kk] <- exp(-d2/2)*PIo[kk]/(prod(s)+.Machine$double.xmin)
    }
    llo <- sum(log(rowSums(p)))
    p <- p/rowSums(p)
    for(kk in 1:k){
      PIo[kk] <- sum(p[,kk])/sum(p)
      MUo[kk,] <- colSums(X*p[,kk])/sum(p[,kk])
      trX <- sweep(X, 2, MUo[kk,], '-')
      So[[kk]] <- t(trX*p[,kk])%*%trX/sum(p[,kk])
    }
    if(max(abs(p-pold))<1e-7) break
  }

  if(llo > ll){
    popt <- p
    MU <- MUo
    PI <- PIo
    S <- So
  }
  Vold <- V
  V <- matrix(0, nrow(Vold), 0)
  for(ii in 1:(ndim-1)){
    ix <- which.max(sapply(1:ncol(Vold), function(rm) fnb_G_3(cbind(V, Vold[,rm]), tX, MU, S, PI, n, 0)))
    V <- cbind(V, Vold[,ix])
    Vold <- Vold[,-ix]
  }
  V <- cbind(V, Vold)

  y <- apply(popt, 1, which.max)

  list(V=V, X = X, y = y, MU=MU, S=S, PI=PI, k = k, y0 = y0, posterior = p, V0 = V0)
}

