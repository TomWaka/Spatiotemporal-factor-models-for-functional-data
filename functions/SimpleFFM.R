
FFM <- function(Y, Ws, M, mc, bn) {
  
  # Y  ; observed data              ; N \times T \times K
  # yt ; observed data              ; N \times K
  # zt ; state parameters           ; N \times K
  # B  ; factor loading matrix      ; N \times M
  # X  ; factor                     ; M \times T \times K
  # G  ; AR coefficicient matrix    ; MK \times MK
  # R  ; kernel of GP
  # W  ; adjacent matrix
  
  {
    
    
    R <- function(a, b , phi) {
      return(exp( -abs(a-b)/2/phi ) )
    }
    
    
    Q <- function(psi, W) {
      if (is.matrix(W)) {
        MM <- ncol(W)
        I  <- diag(1, MM)  
        return( Matrix((I - psi*W) %*% t(I - psi*W)) )
      } else {
        return(0)
      }
    } 
    
    
    post_phi <- function(phi, phi_new, beta, Z, X, B, TR, eta, s) {
      re   <- phi_new^(-3)*exp(-beta/phi_new) / phi^(-3) / exp(-beta/phi) 
      TRN  <- tildeR(eta, phi_new)
      STRN <- solve(TRN)
      STR  <- solve(TR)
      DTRN <- det(TRN)
      DTR  <- det(TR)
      ZZZ <- (Z[s, 1:T, 1:K]- t(t(kronecker(B[s,],diag(1,K)) )  %*%array( aperm(X[1:M,1:T,1:K], c(3,1,2)) , dim=c(M*K,T) ) ))
      re <- re * as.matrix( DTRN^{-1/2} / DTR^{-1/2} * exp(- sum(diag(  ZZZ %*% (STRN*eta^2)  %*%  t(ZZZ)    ))  /2   +  sum(diag(  ZZZ %*% (STR*eta^2)  %*%  t(ZZZ)    ))   /2))
      return(re)
    }
    
    post_phid <- function(phid, phid_new, beta, mu, etad, s) {
      re   <- phid_new^(-3)*exp(-beta/phid_new) / phid^(-3) / exp(-beta/phid) 
      TR  <-  tildeR(etad, phid)
      TRN  <- tildeR(etad, phid_new)
      STRN <- solve(TRN)
      STR  <- solve(TR)
      DTRN <- det(TRN)
      DTR  <- det(TR)
      re <- re * matrix( DTRN^{-1/2} / DTR^{-1/2} * exp(- t(mu) %*% (STRN*etad^2) %*%  mu/2 +  t(mu) %*% (STR*etad^2)  %*% mu  /2) )
      return(re)
    }
    
    
    post_psi <- function(psi, psi_new, th, B, Ws) {
      re <- psi_new^(17)*(1-psi_new)/psi^(17)/(1-psi)
      for (s in 1:M) {
        #if (s==1){
        #  re <- re * det(Q(psi,Ws[[s]]))^{1/2} * exp( Q(psi,Ws[[s]]) * (B[s,]) %*% B[s,]  /2 /th[s]^2   )
        #} else if(s-1>M) {
        #  re <- re * det(Q(psi,Ws[[M]]))^{1/2} * exp( (B[s,1:M]) %*% Q(psi,Ws[[M]]) %*% (B[s,1:M]) /2 /th[s]^2   )
        #} else {
        function_psi <- function(psi) {
          ret <- det(Q(psi,Ws[[s]]))^{1/2} * exp(- (B[(s+1):N,s]) %*% Q(psi,Ws[[s]]) %*% (B[(s+1):N,s]) /2 /tau^2/th[s]^2 )
          if (as.numeric(ret)==0) {
            ret<- 10^(-50)
          }
          return(ret)
        }
        if (N-s==1){
          re <- re * exp(- (B[(s+1):N,s])  * (B[(s+1):N,s]) /2/tau^2 /th[s]^2 )
        } else {
          re <- re * function_psi(psi_new)/ function_psi(psi)
        }
        #}
      }
      return(as.numeric(re))
    }
    
    N <- dim(Y)[1]
    T <- dim(Y)[2]
    K <- dim(Y)[3]
    ups <- c(1:K)
    
    
    
    # initial values
    beta  <- (K-1)/(-2*log(0.05))
    lam   <- rep(1, M)  
    e     <- rep(1, N)  
    eta   <- rep(1, N)  
    etad  <- rep(1, N)  
    th    <- c(rep(1,M))
    tau   <- 0.5
    ze    <- rep(1/2, M)  
    nu    <- 1
    phi   <- rep(1/2, N)  
    phid  <- rep(1/2, N)  
    psi   <- 1/2
    gam   <- rep(0.9, M)  
    Z     <- Y 
    X     <- Y[1:M,,] 
    mu    <- apply(X[,DD==1,], c(1,3), mean) - apply(X[,DD==0,], c(1,3), mean)
    for (s in 1:M) {
      mu[s,]    <- t(D[(D!=0)]) %*% (X[s,c(1:T)[(D!=0)],]- X[s,c(1:T)[(D!=0)]-1,])/ length(D[(D!=0)])
    }
    for (i in 1:M) {
      if (i==1) {
        B <- diag(0, N, M)
      }
      B[i,i] <- 1
    }
    
    TR <-  array(NA, dim = c(N, K, K))   
    STR <- array(0, dim = c(N*K, N*K))   
    
    ## preparation
    
    #mg <- rep(0.95, M)
    #sg <- rep(1, M)
    
    #nlam <- 1
    #slam <- 1
    
    #neta <- 1
    #seta <- 1
    
    #ne <- 1
    #se <- 1/10
    
    #nth <- 1
    #sth <- 1
    
    
    # new version
    in_B <- function(X, B, Z, t, s, TR) { 
      TTZ <- as.vector(t(Z[(s+1):N,t,] - B[(s+1):N,-s] %*% X[-s,t,]))
      #CCC <- t(X[s,t,])%*%sTR
      #TTZ <- Z[(s+1):N,t,] - B[(s+1):N,-s] %*% X[-s,t,]
      CCC <- matrix(0, nr=N-s, nc=K*(N-s))
      for (j in 1:(N-s)) {
        CCC[j,(24*(j-1)+1):(24*(j))] <- X[s,t,] %*% solve(TR[s+j,,])
      }
      #AAb <- Matrix(kronecker.prod(diag(N-s) , CCC))
      return( CCC %*% TTZ) 
    }
    in_ccov <- function(X, TR, s, j) {
      ret <- 0
      INVTR <- solve(TR[s+j,,])
      for (t in 1:T) {
        ret <- ret + X[s,t,] %*% INVTR %*% X[s,t,]
      }
      return(ret)
    }
    inin_B <- function(s,Z,X,B,TR,psi,Ws,th,tau) {
      #sTR <- solve(TR[s,,])
      #sTR <- TR[((s+1):N),,]
      #ccov<- kronecker.prod(diag(N-s), sum( diag(X[s,,] %*% sTR %*% t(X[s,,]) ) ))
      #ccov <- diag(0, N-s)
      #for(t in 1:T) {
      #  for (j in 1:(N-s)) {
      #    ccov[j,j] <- ccov[j,j] + X[s,t,] %*% solve(TR[s+j,,]) %*% X[s,t,]
      #  }
      #}
      #for (j in 1:(N-s)) {
      #  print(in_ccov(X, TR, s, j))
      #}
      ccov <- bdiag(mclapply(1:(N-s), function(j) { in_ccov(X, TR, s, j) }, mc.cores = 20))
      Red <- Reduce("+", mclapply(1:T, function(t) { in_B(X, B, Z, t, s, TR) }, mc.cores = 20))
      #Sb  <- solve(as(ccov + Q(psi, Ws[[s]])/tau^2/th[s]^2 , "sparseMatrix"))  
      Sb  <- solve(ccov + Q(psi, Ws[[s]])/tau^2/th[s]^2 )  
      #return( rtmvnorm(1, mu=Sb %*% Red, sigma=Sb))#,lb=rep(-1,length(Sb %*% Red)), ub=rep(1.2,length(Sb %*% Red)) ))
      return(as.array(rmvn.sparse(1, mu=Sb %*% Red, CH=Cholesky(Sb), FALSE)))
    } 
    
    colB <- function(Z, X, B, TR, psi, Ws, th,tau) {
      B   <- as.array(B)
      Red <-  mclapply(1:M, function(s) { inin_B(s,Z,X,B,TR,psi,Ws,th,tau)  }, mc.cores = 2)
      for (s in 1:M) {
        B[(s+1):N,s] <- Red[[s]]
      }
      return(B)
    }
    
    para_post <- function (s,X,lam,B,Z,ze,th,mu,eta,e,phi,tau,nu) {
      if (s < M+1) {
        # update G
        mg <- rep(0.95, M)
        sg <- rep(1, M)
        sg[s] <-  ( sum(diag(t(X[s,(1:T-1),1:K]) %*% X[s,(1:T-1),1:K])) /lam[s]^2 + 1 )^(-1/2)
        mg[s] <-  sg[s]^2 * (sum(diag(t(X[s,(2:T),1:K] - D[2:T]%o%mu[s,] ) %*% X[s,(1:T-1),1:K])) /lam[s]^2 + 1 )
        gam[s] <- rtruncnorm(1, a=-0.999, b=0.999, mg[s], sd=sg[s]) #0.98
        
        # update lambda
        sum_lam　<-　sum( (X[s, 2:T, 1:K] - gam[s]*X[s, 1:(T-1),1:K] -  D[2:T] %o% mu[s,] ) ^2 )
        alam   <- (1 + (T-1)*K) /2
        blam   <- (1 + sum_lam )/2
        lam[s] <- sqrt(rinvgamma(1, alam, blam))
        
        # update theta
        a_th <- (N-s) / 2
        b_th <- as.numeric((B[(s+1):N,s]) %*% Q(psi,Ws[[s]]) %*% (B[(s+1):N,s]) /tau^2/2 + 1/ze[s])
        th[s] <- sqrt(rinvgamma(1, a_th, b_th))
        
        # update zeta
        a_ze  <- 1
        b_ze  <- 1 + 1/th[s]^2
        ze[s] <- rinvgamma(1, a_ze, b_ze)
        
        # update mu
        Smu <- solve(sum(D^2)/lam[s]^2 *diag(K) )
        Mmu <- c(t( t(D[(D!=0)]) %*% (X[s,c(1:T)[(D!=0)],]- gam[s]*X[s,c(1:T)[(D!=0)]-1,]) /lam[s]^2))
        mu[s,] <- mvrnorm(1, Smu %*% Mmu, Smu)    
        
        # tau
        a_tau <- (N-s)/2
        b_tau <- as.numeric((B[(s+1):N,s]) %*% Q(psi,Ws[[s]]) %*% (B[(s+1):N,s]))/th[s]^2 /2 
      }
      
      #sum_eta <- sum(diag( ( Z[s, 1:T, 1:K]- t(t(kronecker(B[s,],diag(1,K)) )  %*%array( aperm(X[1:M,1:T,1:K], c(3,1,2)) , dim=c(M*K,T) ) )) %*% (solve(TR[s,,])*eta[s]^2)  %*%  t(Z[s, 1:T, 1:K]- t(t(kronecker(B[s,],diag(1,K)) )  %*%array( aperm(X[1:M,1:T,1:K], c(3,1,2)) , dim=c(M*K,T) ) ))    ))
      sum_eta <- 0
      for (t in 1:T) {
        sum_eta <- sum_eta + t(Z[s, t, 1:K] - kronecker(t(B[s,]),diag(1,K)) %*% as.vector(t(X[1:M,t,1:K])) ) %*% (solve(TR[s,,])*eta[s]^2) %*% (Z[s, t, 1:K] - kronecker(t(B[s,]),diag(1,K)) %*% as.vector(t(X[1:M,t,1:K])) )    
      }
      
      # update eta
      a_eta  <- (1 + T*K) /2
      b_eta  <- as.matrix((1 + sum_eta )/2)
      eta[s] <- sqrt(rinvgamma(1, a_eta, b_eta))
      
      # update e
      ae    <- (1+T*K) /2
      be    <- (1 + sum((Y[s,,]-Z[s,,])^2) )/2
      e[s]  <- sqrt(rinvgamma(1, ae, be))
      
      # update phi
      phi_new <- phi[s] + rnorm(1, mean = 0, sd = 1/5)     
      if (phi_new <0) {
        phi_new <- - phi_new
      }
      if (phi_new/phi[s] >10) {
        phi_new <- 10*phi[s]
      } else if (phi_new/phi[s] < 0.1) {
        phi_new <- 0.1*phi[s]
      }
      r_phi <- post_phi(phi[s], phi_new, beta, Z, X, B, TR[s,,], eta[s], s)     # acceptance rate
      if (is.na(r_phi)) {
        r_phi <- 0
      }
      u_phi <- runif(1)
      if (r_phi > u_phi) {
        phi[s] <- phi_new
      }
      if (s<M+1){
        GAM<-gam[s]
        LAM<-lam[s]
        ZE<-ze[s]
        TH<-th[s]
        MU<-mu[s,]
      } else {
        GAM<-LAM<-ZE<-TH<-a_tau<-b_tau<-1
        MU<-rep(1,K)
      }
      return(cbind(GAM,LAM,ZE,TH,eta[s],e[s],phi[s],MU,a_tau,b_tau))
    }
    
    tildeR <-  function(eta, phi) {
      for (i in 1:K) {
        if(i==1){
          TR <- diag(K)      # tilde R                      
        }
        for (j in 1:K) {
          TR[i,j] <- TR[j,i] <- eta^2 * R(ups[i], ups[j] , phi)
        }
      }
      return(TR)
    }
    
    in_X <- function(BSTR, SLG, SL,  GL, CSx2t, Sx2t, t, X, D, mu) { 
      mx <- Sx2t %*% ( BSTR %*% as.vector(t(Z[1:N,t,1:K]))  +  SLG %*% as.vector(t(X[1:M,t-1,1:K]))  + GL %*% as.vector(t(X[1:M,t+1,1:K])) + D[t] * SL %*% as.vector(t(mu)) - D[t+1] * GL %*% as.vector(t(mu))  )
      return( t(array( rmvn.sparse(1, mx, CSx2t, FALSE) , dim=c(K,M))) )
    }
    
    in_Z <- function(Y, X, SzA, SzB, SSz, t) { 
      mz <- SzA %*% as.vector(t(Y[1:N,t,1:K])) + SzB %*% as.vector(t(X[1:M,t,1:K]))
      return( t(array( rmvn.sparse(1, mz, SSz, FALSE),dim=c(K,N))) )
    }
    
  } 　
  
  
  # MCMC box
  B.pos <- array(NA, dim = c(N, M, (mc-bn)/3 ))     
  #Z.pos  <- array(NA, dim = c(N, T, K, (mc-bn)/3))
  X.pos  <- array(NA, dim = c(M, T, K, (mc-bn)/3))
  gam.pos <- array(NA, dim = c(M, (mc-bn)/3))
  mu.pos <- array(NA, dim = c(M, K, (mc-bn)/3))
  th.pos  <- array(NA, dim = c(M, (mc-bn)/3))
  
  ## MCMC 
  j <- 1
  pb <- progress_bar$new(total = mc, format = "[:bar] :percent remain: :eta", clear = TRUE)
  
  for (p in 1:mc) {
    pb$tick()
    
    
    # update  tilde R 
    TTTT <- mclapply(1:N, function(s) { tildeR(eta[s], phi[s]) }, mc.cores = 40)
    TR   <- aperm(array(unlist(TTTT) ,dim =c(K,K,N)) , c(3,1,2))    
    
    
    # update Z
    SE  <- Matrix(kronecker.prod( diag(e^(-2)) , diag(1,K) ))
    Sz  <- solve( SE + STR )
    SzA <- Sz %*% SE
    SzB <- Sz %*% STR %*% kronecker(B,diag(1,K))
    SSz <- Cholesky(Sz)
    Z <- aperm( array(unlist( mclapply(1:T, function(t) { in_Z(Y, X, SzA, SzB, SSz, t)}, mc.cores = 40)), dim=c(N,K,T) ), c(1,3,2))
    
    
    # update X    
    SL   <- as(kronecker.prod(diag(lam^(-2)), diag(K)), "sparseMatrix")
    GL   <- t(G) %*% SL
    GLG  <- GL %*% G
    #STR  <- as(solve(kronecker(diag(1,N),TR)), "sparseMatrix")
    STR  <- solve(bdiag(TTTT))
    Sx2t <- solve( t(kronecker(B,diag(1,K))) %*% STR %*% kronecker(B,diag(1,K))  +  SL + GLG  ) 
    #kSx2t  <- kronecker(diag(1,T-2), Sx2t)
    SLG <- SL %*% G
    BSTR <- t(kronecker(B,diag(1,K))) %*% STR
    CSx2t <- Cholesky(Sx2t)
    
    X[1:M,2:(T-1),1:K] <- aperm( array(unlist( mclapply(2:(T-1), function(t) { in_X(BSTR, SLG, SL, GL, CSx2t, Sx2t, t, X, D, mu) }, mc.cores = 40) ), dim=c(M,K,T-2) ), c(1,3,2))
    
    Sx <- solve(  BSTR %*% kronecker(B,diag(1,K))  + GLG  ) 
    mx <- Sx %*%( BSTR %*% as.vector(t(Z[1:N,1,1:K]))   +  GL %*% as.vector(t(X[1:M,2,1:K])) + D[2]* SL %*% as.vector(t(mu)) )
    X[1:M,1,1:K] <- t(array( mvrnorm(1, mx, Sx) , dim=c(K,M)))
    
    Sx <- solve(  BSTR %*% kronecker(B,diag(1,K))  +  SL  ) 
    mx <- Sx %*%( BSTR %*% as.vector(t(Z[1:N,T,1:K]))  +  SLG %*% as.vector(t(X[1:M,T-1,1:K])) - D[T]* GL %*% as.vector(t(mu)) )
    X[1:M,T,1:K] <- t(array( mvrnorm(1, mx, Sx) , dim=c(K,M)))
    print(X[1,T,])
    
    
    # update  G
    G  <- Matrix(kronecker(diag(gam), diag(K)))    
    
    Paras<- mclapply(1:N, function(s) { para_post(s,X,lam,B,Z,ze,th,mu,eta,e,phi,tau)}, mc.cores = 40)
    
    a_tau <- 1/2
    b_tau <- 1/nu
    for (s in 1:N) {
      if (s<M+1){
        gam[s] <- Paras[[s]][1,1]
        lam[s] <- Paras[[s]][1,2]
        ze[s]  <- Paras[[s]][1,3]
        th[s]  <- Paras[[s]][1,4]
        mu[s,] <- Paras[[s]][,8]
        a_tau  <- a_tau + Paras[[s]][1,9]
        b_tau  <- b_tau + Paras[[s]][1,10]
      }
      eta[s] <- Paras[[s]][1,5]
      e[s]   <- Paras[[s]][1,6]
      phi[s] <- Paras[[s]][1,7]
    }
    
    tau <- sqrt(rinvgamma(1, a_tau, b_tau))
    tau <- 1
    
    # update nu
    nu  <- sqrt(rinvgamma(1, 1, tau^(-2)+1))
    
    
    # update psi
    psi_new <- psi + rnorm(1, mean = 0, sd = 1/5)   
    if (psi_new >1) {
      psi_new <- 2 - psi_new
    } else if (psi_new <0) {
      psi_new <- - psi_new
    }
    r_psi <- post_psi(psi, psi_new, th, B, Ws)  # acceptance rate
    u_psi <- runif(1)
    if (r_psi > u_psi) {
      psi <- psi_new
    }
    
    
    # update  B(col type)
    B <- Matrix(colB(Z, X, B, TR, psi, Ws, th, tau))  
    
    
    # store samples
    if (p > bn && p%%3==0) {
      B.pos[,,j]  <- as.matrix(B)
      #Z.pos[,,,j] <- as.matrix(Z)
      X.pos[,,,j] <- as.matrix(X)
      mu.pos[,,j] <- mu
      gam.pos[,j] <- gam
      th.pos[,j]  <- th
      j <- j+1
    }
  }
  
  
  
  Zmed <- apply(Z.pos[,,,] , c(1,2,3), median)
  uci <- apply(Z.pos, c(1,2,3), quantile, probs=0.975)
  lci <- apply(Z.pos, c(1,2,3), quantile, probs=0.025)
  result <- list()
  result$B <- B.pos
  result$Z <- Z.pos
  result$G <- gam.pos
  result$RMSE <- sqrt(mean(abs(Zmed-Z_true)^2))       
  result$CP <- sum(lci<Z_true & Z_true<uci)/T/N/K   
  return(result)
  
}


