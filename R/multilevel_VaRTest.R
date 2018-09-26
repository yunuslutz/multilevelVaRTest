multilevel_VaRTest <- function(alphas, actual, VaR, confidence, m = 5, B = 2000){
  
  actual <- as.matrix(actual)
  K <- length(alphas)
  T <- nrow(VaR)
  
  if(nrow(actual)!=T) {stop("VaR Forecast and actual Returns not of same length")}
  if(ncol(VaR)!=K) {stop("Number of Forecasts does not equal number of alpha values")}
  
  if(alphas[1]<alphas[2]) {warning("VaR forecasts and alphas must be ordered from largest to smallest quantile and were reordered.")
    VaR <- VaR[,seq(from = K, to = 1)]
    alphas <- alphas[seq(from = K, to = 1)]}
  
  J <- matrix(0, nrow = T, ncol = length(alphas))
  I <- matrix(0, nrow = T, ncol = length(alphas)+1)
  
  for(i in 1:K) {
    for(j in 1:K) {
      I[,i] <- if_else(actual<VaR[,i], 1, 0)              #Violation Indicator
      J[,j] <- if_else(actual<VaR[,j], I[,j] - I[,j+1], 0)#Violation Indicator, only highestqr
    }}
  
  J0 <- if_else(apply(J, 1, sum)==0, 1, 0) 
  J <- cbind(J0, J) #add column J0 to matrix, indicating no violations in t
  I <- I[,1:(ncol(I)-1)] #delete last column of 0s
  N <- apply(I, 1, sum) #total number of violations for t=1,...,T
  T_i <- apply(J, 2, sum) #number of observations for which N_t=i
  T_table <- as.matrix(table(head(N, -1), tail(N, -1))) #matrix of violations in t(columns) and t-1(rows)
  
  pi_hat_ij <- matrix(0, nrow = K+1, ncol = K+1) #initializing matrix
  
  for(i in 1:(K+1)){ #transition matrix of conditional probabilities
    for(j in 1:(K+1)){
      pi_hat_ij[i,j] <- T_table[i,j] / T_i[i]
    }
  }
  
  pi_hat_i <- c()
  thetas <- c()
  
  for(i in 1:(K+1)) {
    pi_hat_i[i] <- T_i[i]/T #observed probabilities of falling in between VaR quantiles
  }
  
  for(i in 1:(K-1)) {
    thetas[i] <- alphas[i]-alphas[i+1] 
  }
  thetas <- c((1-alphas[1]), thetas, alphas[K]) #Bernoulli prob. of J under null
  
  thetas_matrix <- matrix(thetas, ncol = 1) %*% thetas
  
  pearson <- function(N, m, K, thetas_matrix){
    T_table <- list()
    X_j <- c()
    
    for (j in 1:m) {
      T_table[[j]] <- as.matrix(table(head(N, -j), tail(N, -j)))
      
      #in case some J are not observed, add empty column/row to matrix to ensure right dimensions
      while (ncol(T_table[[j]])<K+1) {
        T_table[[j]] <- cbind(T_table[[j]], 0)
      }
      while (nrow(T_table[[j]])<K+1) {
        T_table[[j]] <- rbind(T_table[[j]], 0)
      }
      X_j[j] <- sum((T_table[[j]]-(T-j)*thetas_matrix)^2/((T-j)*thetas_matrix))
    }
    X_m <- sum(X_j)
    U <- runif(1)
    return(list(X_m = X_m, U = U))
  }
  
  #pearson test statistics
  X_m_0 <- pearson(N = N, m = m, K = K, thetas_matrix = thetas_matrix)$X_m
  U0 <- pearson(N = N, m = m, K = K, thetas_matrix = thetas_matrix)$U
  
  ##MC simulation afer Dufour (2006)
  sample_space <- list()
  for (i in seq_along(thetas)) {
    sample_space[[i]] <- rep(K-4+i, thetas[i]*10000) #multiplication by a sufficiently large number to ensure integer 
  }
  sample_space <- unlist(sample_space)
  
  X_mj <- c()
  N_mc <- list()
  U_mc <- c()
  
  #generate B statistics from randomly sampled data
  for (i in 1:B) {
    N_mc[[i]] <- sample(sample_space, size = T)
    X_mj[i] <- pearson(N = N_mc[[i]], m = m, K = K, thetas_matrix = thetas_matrix)$X_m
    U_mc[i] <- pearson(N = N_mc[[i]], m = m, K = K, thetas_matrix = thetas_matrix)$U
  }
  #compare simulated test statistic to "real" one. 
  G_hat <- 1 - (1/B)*sum(X_m_0 >= X_mj) + (1/B)*sum((X_m_0 == X_mj)*(U0 <= U_mc))
  pvalue.pearson <- (B*G_hat+1)/(B+1)
  decision.pearson <- ifelse(pvalue.pearson < (1-confidence), "Reject H0", "Can not reject H0")
  
  LR.ind <- 2*(sum(T_table*log(pi_hat_ij), na.rm = TRUE)-sum(T_i*log(pi_hat_i), na.rm = TRUE))
  LR.cc <- 2*(sum(T_table*log(pi_hat_ij), na.rm = TRUE) - sum(T_i*log(thetas), na.rm = TRUE))
  LR.uc <- LR.cc - LR.ind
  
  critvalue_cc <- qchisq(confidence, df = K^2)
  critvalue_uc <- qchisq(confidence, df = K)
  
  pvalue.cc <- 1 - pchisq(LR.cc, df = K^2)
  pvalue.uc <- 1 - pchisq(LR.uc, df = K)
  
  decision.cc <- if_else(pvalue.cc < (1-confidence), "Reject H0", "Can not reject H0")
  decision.uc <- if_else(pvalue.uc < (1-confidence), "Reject H0", "Can not reject H0")
  
  result <- list(LR.cc, critvalue_cc, pvalue.cc, decision.cc, LR.uc, critvalue_uc, pvalue.uc, decision.uc, X_m_0, pvalue.pearson, decision.pearson)
  names(result) <- c("likelihood ratio CC", "critical value CC", "p-value CC", "decision CC", 
                     "likelihood ratio UC", "critical value UC", "p-value UC", "decision UC",
                     "X_m", "p-value pearson", "decision.person")
  return(result)
}






