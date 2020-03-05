### function evaluation for projection matrix V
## V = projection matrix
## X = data matrix
## Y = vector of class labels
## MU = matrix of class means
## S = list of class covariance matrices
## PI = vector of class priors
## N = total number of data (i.e., nrow(X))

fnb_G_2 <- function(V, tX, Y, MU, S, PI, N, Joint = FALSE){

  nc <- length(PI)

  ### normalise V and turn into a matrix in case it is given as a vector
  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)

  ### compute statistics on projection
  #mu <- MU%*%V
  #s <- sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V))

  mu <- matrix(MU%*%V, ncol = d)
  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)

  #x <- t(V)%*%tX
  x <- matrix(t(V)%*%tX, nrow = d)

  ### The first term in the log-likelihood essentially only depends on the variance terms:
  #T1 <- sum(N*PI*(log(PI) - sapply(1:nc, function(k) sum(log(s[,k])))/2 - d/2))

  T1 <- sum(N*PI*log(PI)-sapply(1:nc, function(k) N*PI[k]*sum(log(s[,k]))/2+sum(.Internal(rowSums((x[,Y==k] - mu[k,])^2, d, N*PI[k], TRUE))/s[,k])/2))

  if(Joint) return(T1)

  ### compute projected data (in each class. This could be removed if we always use the actual covariances in S,
  ### or scaled versions of it. It might be useful to try other alternatives, though, such as adding a ridge
  ### term to decrease the flexibility)




  ### compute the (scaled) density of all points in all classes
  p <- numeric(N)
  for(k in 1:nc) p <- p + exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/prod(s[,k]^.5)*PI[k]

  T1 - sum(log(p))
}

dfnb_G_2 <- function(V, tX, Y, MU, S, PI, N, Joint = FALSE){
  ### The first steps in the gradient are the same as in the function evaluation

  nc <- length(PI)

  ### normalise V and turn into a matrix in case it is given as a vector
  V <- matrix(V, nrow = ncol(MU))

  d <- ncol(V)

  #s <- sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V))
  s <- matrix(sapply(1:nc, function(k) diag(t(V)%*%S[[k]]%*%V)), ncol = nc)
  sv <- lapply(S, function(s) s%*%V)


  ### The derivative of the first term in the log-lkelihood also only depends on the derivatives of the variance terms
  ssv <- lapply(1:nc, function(k) t(sv[[k]])/s[,k])
  dv <- t(V*0)
  for(k in 1:nc) dv <- dv - ssv[[k]]*N*PI[k]

  #T1 <- sum(N*PI*log(PI)-sapply(1:nc, function(k) N*PI[k]*sum(log(s[,k]))/2+sum(.Internal(rowSums((x[,Y==k] - mu[k,])^2, d, N*PI[k], TRUE))/s[,k])/2))

  if(Joint) return(t(dv))

  dets <- apply(s, 2, prod)^.5

  ### compute statistics on projection
  #mu <- MU%*%V
  mu <- matrix(MU%*%V, ncol = d)

  #x <- t(V)%*%tX
  x <- matrix(t(V)%*%tX, nrow = d)

  ### compute the (scaled) density of all points in all classes
  p <- matrix(0,N,nc)
  for(k in 1:nc) p[,k] <- exp(-.Internal(colSums((x - mu[k,])^2/s[,k]/2, d, N, TRUE)))/dets[k]*PI[k]

  ps <- rowSums(p)




  ### Next compute the derivative of the second term in the log-likelihood
  for(k in 1:nc){
    Xk <- tX-MU[k,]
    Xkv <- t(V)%*%Xk
    #dv <- dv - ssv[[k]]*colSums((p[,k]/ps)*t(Xkv^2/s[,k]-1))
    #dv <- dv - ssv[[k]]*rowSums(sweep(Xkv^2/s[,k]-1, 2, p[,k]/ps, '*'))
    #dv <- dv - ssv[[k]]*sapply(1:d, function(j) sum((Xkv[j,]^2/s[j,k]-1)*p[,k]/ps))
    #  rowSums(sweep(Xkv^2/s[,k]-1, 2, p[,k]/ps, '*'))
    dv <- dv - (Xkv[,Y==k]/s[,k])%*%t(Xk[,Y==k])
    dv <- dv + ssv[[k]]*.Internal(rowSums(Xkv[,Y==k]^2/s[,k], d, N*PI[k], TRUE))
    dv <- dv - ssv[[k]]*.Internal(colSums((p[,k]/ps)*t(Xkv^2/s[,k]-1), N, d, TRUE))
    dv <- dv + (Xkv/s[,k])%*%((p[,k]/ps)*t(Xk))
  }
  t(dv)
}



dfapp = function(V, tX, Y, MU, S, PI, N, Joint = FALSE){
  out = V*0
  for(i in 1:nrow(V)){
    for(j in 1:ncol(V)){
      e = V*0
      e[i,j] = 1e-5
      out[i,j] = (fnb_G_2(V+e, tX, Y, MU, S, PI, N, Joint)-fnb_G_2(V-e, tX, Y, MU, S, PI, N, Joint))/2e-5
    }
  }
  out
}
