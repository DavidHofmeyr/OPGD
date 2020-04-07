### function evaluation for projection matrix V
## V = projection matrix
## X = data matrix
## Y = vector of class labels
## MU = matrix of class means
## S = list of class covariance matrices
## PI = vector of class priors
## N = total number of data (i.e., nrow(X))

fnb_G_2 <- function(V, tX, Y, MU, S, PI, N){

  nc <- length(PI)

  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)


  mu <- matrix(MU%*%V, ncol = d)
  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)

  x <- matrix(t(V)%*%tX, nrow = d)

  T1 <- sum(N*PI*log(PI)-sapply(1:nc, function(k) N*PI[k]*sum(log(s[,k]))/2+sum(.Internal(rowSums((x[,Y==k] - mu[k,])^2, d, N*PI[k], TRUE))/s[,k])/2))

  p <- numeric(N)
  for(k in 1:nc) p <- p + exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/prod(s[,k]^.5)*PI[k]

  T1 - sum(log(p+.Machine$double.xmin*(p==0)))
}

dfnb_G_2 <- function(V, tX, Y, MU, S, PI, N, Joint = FALSE){

  nc <- length(PI)

  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)

  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)
  sv <- lapply(S, function(s) s%*%V)


  ssv <- lapply(1:nc, function(k) t(sv[[k]])/s[,k])
  dv <- t(V*0)
  for(k in 1:nc) dv <- dv - ssv[[k]]*N*PI[k]

  dets <- apply(s, 2, prod)^.5

  mu <- matrix(MU%*%V, ncol = d)

  x <- matrix(t(V)%*%tX, nrow = d)

  p <- matrix(0,N,nc)
  for(k in 1:nc) p[,k] <- exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/dets[k]*PI[k]

  p[p==0] <- .Machine$double.xmin

  ps <- rowSums(p)

  for(k in 1:nc){
    Xk <- tX-MU[k,]
    Xkv <- t(V)%*%Xk
    dv <- dv - (Xkv[,Y==k]/s[,k])%*%t(Xk[,Y==k])
    dv <- dv + ssv[[k]]*.Internal(rowSums(Xkv[,Y==k]^2/s[,k], d, N*PI[k], TRUE))
    dv <- dv - ssv[[k]]*.Internal(colSums((p[,k]/ps)*t(Xkv^2/s[,k]-1), N, d, TRUE))
    dv <- dv + (Xkv/s[,k])%*%((p[,k]/ps)*t(Xk))
  }
  t(dv)
}



fnb_G_3 <- function(V, tX, MU, S, PI, N, om = 0){

  nc <- length(PI)

  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)

  mu <- matrix(MU%*%V, ncol = d)
  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)

  x <- matrix(t(V)%*%tX, nrow = d)

  p <- matrix(0, N, nc)
  for(k in 1:nc) p[,k] <- exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/prod(s[,k]^.5)*PI[k]

  p[p==0] <- .Machine$double.xmin

  sum(log(apply(p, 1, function(x) max(x)/sum(x)))) - om*sum((t(V)%*%V-diag(d))^2)
}


dfnb_G_3 <- function(V, tX, MU, S, PI, N, om){

  nc <- length(PI)

  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)

  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)
  sv <- lapply(S, function(s) s%*%V)

  dets <- apply(s, 2, prod)^.5

  mu <- matrix(MU%*%V, ncol = d)

  x <- matrix(t(V)%*%tX, nrow = d)

  p <- matrix(0,N,nc)
  for(k in 1:nc) p[,k] <- exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/dets[k]*PI[k]

  p[p==0] <- .Machine$double.xmin

  Y <- apply(p, 1, which.max)

  mypi <- sapply(1:nc, function(k) sum(Y==k))/N

  ssv <- lapply(1:nc, function(k) t(sv[[k]])/s[,k])
  dv <- t(V*0)
  for(k in 1:nc) dv <- dv - ssv[[k]]*N*mypi[k]

  ps <- rowSums(p)

  for(k in 1:nc){
    Xk <- tX-MU[k,]
    Xkv <- t(V)%*%Xk
    dv <- dv - (Xkv[,Y==k]/s[,k])%*%t(Xk[,Y==k])
    dv <- dv + ssv[[k]]*.Internal(rowSums(Xkv[,Y==k]^2/s[,k], d, N*mypi[k], TRUE))
    dv <- dv - ssv[[k]]*.Internal(colSums((p[,k]/ps)*t(Xkv^2/s[,k]-1), N, d, TRUE))
    dv <- dv + (Xkv/s[,k])%*%((p[,k]/ps)*t(Xk))
  }
  t(dv) - 4*om*V%*%(t(V)%*%V - diag(d))
}




dfapp3 = function(V, tX, MU, S, PI, N, om){
  out = V*0
  for(i in 1:nrow(V)){
    for(j in 1:ncol(V)){
      e = V*0
      e[i,j] = 1e-5
      out[i,j] = (fnb_G_3(V+e, tX, MU, S, PI, N, om)-fnb_G_3(V-e, tX, MU, S, PI, N, om))/2e-5
    }
  }
  out
}




