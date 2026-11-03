start_time <- Sys.time()
M <- 1000 # Number of simulation runs
n <- 100 #Sample size
BB <-50 #bootstrap Sample
c0 <-5
G<-3
############
#cure rate
rho1=0.2
rho2=0.2
eta=0
phi01 <- rho1*(1-rho2) - eta
phi10 <- (1-rho1)*rho2 - eta
phi11 <- (1-rho1)*(1-rho2) + eta
phi00 <- rho1 * rho2 + eta
#True value
true_params  = c(beta1 =0.6, beta2 =0.6, beta3 =0.1,rho1=rho1,rho2=rho2)  
##############Matrix definition############# 
estimates_matrix= matrix(0, M, 5) 
B_estimates_matrix= matrix(0, BB, 5) 
ese.repe <- matrix(0,M, 5)
cp.delta.repe <- matrix(0, M, 5)
###########Vector initialization#######################
convergence <- numeric(M)  
beta1_estimates <- numeric(M)  
beta2_estimates <- numeric(M)  
beta3_estimates <- numeric(M)  
phi00_estimates <- numeric(M)  
phi01_estimates <- numeric(M)  
phi10_estimates <- numeric(M)  
phi11_estimates <- numeric(M) 
B_beta1_estimates <- numeric(BB)  
B_beta2_estimates <- numeric(BB)  
B_beta3_estimates <- numeric(BB)  
B_phi00_estimates <- numeric(BB)  
B_phi01_estimates <- numeric(BB)  
B_phi10_estimates <- numeric(BB)  
B_phi11_estimates <- numeric(BB) 

convergence_count <- 0  
max_iter <- M+1000 
###########Main loop################################
for (iter in 1:max_iter) {  
  if (convergence_count >= M) break  
  
  Z1 <- rnorm(n, mean = 0, sd = 1)  #Covariate 
  Z2 <- rnorm(n, mean = 0, sd = 1)
  u <- runif(n, min = 0, max = 1)  
  v <- runif(n, min = 0, max = 1)  
  w <- runif(n, min = 0, max = 1)
  
  S1 <- log(u) / (-G * exp(Z1 * true_params[1]))  
  S2 <- log(v) / (-G * exp(Z2 * true_params[2]))  
  S3 <- log(w) / (-true_params[3])  
  T1 <- sapply(1:n, function(i) min(S1[i], S3[i]))  
  T2 <- sapply(1:n, function(i) min(S2[i], S3[i])) 
  # mixture cure model #
  #Cure indicator 
  V1 = rbinom(n, 1, 1-rho1)  
  V2 = rbinom(n, 1, 1-rho2)  
  #Censoring time
  C1 <- runif(n, 0, c0)   
  C2 <- runif(n, 0, c0)
  #Failure time 
  T1_1=ifelse(V1==1,T1,999)
  T2_2=ifelse(V2==1,T2,999)
  delta1=as.integer(T1_1<=C1) 
  delta2=as.integer(T2_2<=C2)
  #obversed time
  t1=pmin(T1_1,C1)
  t2=pmin(T2_2,C2)
  ##################
  cure1 <- sum(V1 == 0)  
  cure2 <- sum(V2 == 0)  
  #Calculate the cure proportion
  cure_rate1 <- cure1 / n  
  cure_rate2 <- cure2 / n  
  print(paste("Cure proportion (Group 1):", cure_rate1))  
  print(paste("Cure proportion (Group 2):", cure_rate2))
  #Censoring proportion
  rate1 <- (n - sum(delta1)) / n
  rate2 <- (n - sum(delta2)) / n
  rate1
  rate2
  
  
  ###############################Sieve##########
  #Basis functions of the Bernstein polynomial 
  b.km = function(t, k, m, cj, uj){  
    if ((t - cj) / (uj - cj) < 0 || (t - cj) / (uj - cj) > 1) {  
      return(0)  
    }  
    y = choose(m, k) * ((t- cj) / (uj - cj))^k * (1 - (t - cj) / (uj - cj))^(m - k)  
    return(y)}  
  
  n <- length(t1)  
  cj <- min(t1)-0.01  
  uj <- max(t1)+0.01 
  v=0.3
  m =floor(n^v) 
  A = matrix(0, nrow = n, ncol = m + 1)  
  
  for(i in 1:n){  
    for(k in 0:m){  
      A[i, k + 1] = b.km(t1[i], k, m, cj, uj)  
    }  
  }  
  
  db.km = function(t, k, m, cj, uj){  
    if ((t - cj) / (uj - cj) < 0 || (t - cj) / (uj - cj) > 1) {  
      return(0)  
    }  
    
    dt = 1 / (uj - cj)  
    term1 = choose(m, k) * k * ((t - cj) / (uj - cj))^(k - 1) * dt * (1 - (t - cj) / (uj - cj))^(m - k)  
    term2 = choose(m, k) * ((t - cj) / (uj - cj))^k * (m - k) * (-dt) * (1 - (t - cj) / (uj - cj))^(m - k - 1)
    y = term1 + term2  
    return(y)  
  }  
  
  D_A = matrix(0, nrow = n, ncol = m + 1)  
  
  for(i in 1:n){  
    for(k in 0:m){  
      D_A[i, k + 1] = db.km(t1[i], k, m, cj, uj)  
    }  
  }  
  
  
  #####lamda2#######
  ####Basis functions of the Bernstein polynomial #####
  b.km1 = function(t, k, m, cj1, uj1){  
    if ((t - cj1) / (uj1 - cj1) < 0 || (t - cj1) / (uj1 - cj1) > 1) {  
      return(0)  
    }  
    y1 = choose(m, k) * ((t - cj1) / (uj1 - cj1))^k * (1 - (t - cj1) / (uj1 - cj1))^(m - k)  
    return(y1)  
  }  
  
  n <- length(t2)  
  cj1 <- min(t2)- 0.01  
  uj1 <- max(t2)+0.01  
  
  v=0.3
  m =floor(n^v) 
  B = matrix(0, nrow = n, ncol = m + 1)  
  
  for(i in 1:n){  
    for(k in 0:m){  
      B[i, k + 1] = b.km1(t2[i], k, m, cj1, uj1)  
    }  
  }  
  
  db.km1 = function(t, k, m, cj1, uj1){  
    if ((t - cj1) / (uj1 - cj1) < 0 || (t - cj1) / (uj1 - cj1) > 1) {  
      return(0)  
    }  
    
    scale_factor <- 1 / (uj1 - cj1)  
    term1 <- choose(m, k) * k * ((t - cj1) / (uj1 - cj1))^(k - 1) * scale_factor *  
      (1 - (t - cj1) / (uj1 - cj1))^(m - k)  
    term2 <- choose(m, k) * ((t - cj1) / (uj1 - cj1))^k * (m - k) *   
      (-scale_factor) * (1 - (t - cj1) / (uj1 - cj1))^(m - k - 1)  
    
    y1_prime <- term1 + term2  
    return(y1_prime)  
  }  
  
  D_B = matrix(0, nrow = n, ncol = m + 1)  
  for(i in 1:n){    
    for(k in 0:m){    
      D_B[i, k + 1] = db.km1(t2[i], k, m, cj1, uj1)    
    }    
  }  
  #####################Likelihood#########################
  logL <- function(par){  
    beta1<-par[1]
    beta2<-par[2]
    beta3 <- exp(par[3])
    phi00<-exp(par[4])/(1+exp(par[4]))
    phi01<-1-phi00-exp(par[5])*(1-phi00)/(1+exp(par[5]))
    phi10 <- 1 - phi00 - phi01 - exp(par[6]) * (1 - phi00 - phi01) / (1 + exp(par[6]))  
    phii1<-par[7:(7+m)]
    phii2<-par[(8+m):(8+2*m)]
    phi11=1 - phi00 - phi10 - phi01 
    phi1=cumsum(exp(phii1))
    phi2=cumsum(exp(phii2))
    Dlamda1=rowSums(A%*%phi1)
    Dlamda2=rowSums(B%*%phi2)
    lamda2=rowSums(D_B%*%phi2)
    lamda1=rowSums(D_A%*%phi1) 
    
    pian1  <- phi11*exp(-Dlamda1*exp(Z1*beta1)-Dlamda2*exp(Z2*beta2)-beta3*pmax(t1,t2))*(-lamda1 *exp(Z1 * beta1 )- beta3*(t1 >= t2)) + phi10 * exp(-Dlamda1 *exp(Z1 * beta1)  - beta3 * t1)*(-lamda1 *exp(Z1 * beta1 )- beta3) 
    
    
    pian2 <- phi11*exp(-Dlamda1*exp(Z1*beta1)-Dlamda2*exp(Z2*beta2)-beta3*pmax(t1,t2))*(-lamda2 *exp(Z2 * beta2 )- beta3*(t2 >= t1))+phi01 * exp(-Dlamda2*exp(Z2*beta2)-beta3*t2)*(-lamda2 *exp(Z2* beta2 )- beta3)
    
    
    St1t2<- phi11*exp(-Dlamda1*exp(Z1*beta1)-Dlamda2*exp(Z2*beta2)-beta3*pmax(t1,t2))+phi10*exp(-Dlamda1*exp(Z1*beta1)-beta3*t1)+phi01*exp(-Dlamda2*exp(Z2*beta2)-beta3*t2)+phi00
    ft1t2 <- phi11*exp(-Dlamda1*exp(Z1*beta1)-Dlamda2*exp(Z2*beta2)-beta3*pmax(t1,t2))*(-lamda1 *exp(Z1 * beta1)- beta3*(t1 >= t2))*(-lamda2 *exp(Z2 * beta2 )- beta3*(t2 >= t1))  
    
    piandao1_value<- -pian1
    piandao2_value <- -pian2
    piandao1_value<-ifelse(piandao1_value==0,1e-4,-pian1) 
    piandao2_value <-ifelse(piandao2_value==0,1e-4,-pian2)
    
    ft1t2_shifted <-  ft1t2
    likelognorm=sum(delta1*delta2*log(ft1t2_shifted )+delta1*(1-delta2)*log(piandao1_value)+delta2*(1-delta1)*log(piandao2_value)
                    +(1-delta1)*(1-delta2)*log(St1t2))
    return(-likelognorm)  
  }
  
  #Initial value
  initial_values<- c(true_params[1]+0.2,true_params[2]+0.2,log(true_params[3]+0.2),0,0,0.07, rep(0,2*(m+1))) 
  #Optimize the log-likelihood function
  result <- try(optim(par = initial_values, fn = logL, control = list(maxit = 12000)), silent = TRUE)
  if ('try-error' %in% class(result) || result$convergence != 0) {
    next  
  }
  ###Store the current convergence result
  convergence_count <- convergence_count + 1
  convergence[convergence_count]= result$convergence
  beta1_estimates[convergence_count] <- result$par[1]
  beta2_estimates[convergence_count] <- result$par[2]
  beta3_estimates[convergence_count] <- exp(result$par[3])
  phi00_estimates[convergence_count] <- exp(result$par[4]) / (1 + exp(result$par[4]))
  phi01_estimates[convergence_count] <- 1 - phi00_estimates[convergence_count] - exp(result$par[5]) * (1 - phi00_estimates[convergence_count]) / (1 + exp(result$par[5]))
  phi10_estimates[convergence_count] <- 1 - phi00_estimates[convergence_count] - phi01_estimates[convergence_count] - exp(result$par[6]) * (1 - phi00_estimates[convergence_count] - phi01_estimates[convergence_count]) / (1 + exp(result$par[6]))
  phi11_estimates[convergence_count] <- 1 - phi00_estimates[convergence_count] - phi10_estimates[convergence_count] - phi01_estimates[convergence_count]
  
  ##############rho1.rho2#######
  rho1_estimates <- phi00_estimates + phi01_estimates
  rho2_estimates <- phi00_estimates + phi10_estimates
  ############Parameter vector##############
  est=c(
    beta1_estimates[convergence_count],
    beta2_estimates[convergence_count],
    beta3_estimates[convergence_count],
    rho1_estimates[convergence_count],
    rho2_estimates[convergence_count] )
  ##############Parameter matrix#######
  estimates_matrix[convergence_count, ] <- c(
    beta1_estimates[convergence_count],
    beta2_estimates[convergence_count],
    beta3_estimates[convergence_count],
    rho1_estimates[convergence_count],
    rho2_estimates[convergence_count] )
  
  print(paste("Number of successful convergences:", convergence_count, "beta1 =", beta1_estimates[convergence_count], "beta2 =", beta2_estimates[convergence_count], "beta3 =", beta3_estimates[convergence_count]))
  
  ################### bootstrap  ##############################
  #########bootstrap loop############
  for (ii in 1:BB) {
    
    ################# Sampling with replacement( t1, t2, Z1, Z2, delta1, delta2 )  
    Sam <- sample(n, n, replace = TRUE)
    St1 <- t1[Sam]
    St2 <- t2[Sam]
    SZ1 <- Z1[Sam]
    SZ2 <- Z2[Sam]
    Sdelta1 <- delta1[Sam]
    Sdelta2 <- delta2[Sam]
    
    #############################
    ###############################Sieve##########
    b.km = function(t, k, m, cj, uj){  
      if ((t - cj) / (uj - cj) < 0 || (t - cj) / (uj - cj) > 1) {  
        return(0)  } 
      
      y = choose(m, k) * ((t- cj) / (uj - cj))^k * (1 - (t - cj) / (uj - cj))^(m - k)  
      return(y)  
    }  
    
    n <- length(St1)  
    cj <- min(St1)-0.01  
    uj <- max(St1)+0.01 
    v=0.3
    m =floor(n^v) 
    A = matrix(0, nrow = n, ncol = m + 1)  
    for(i in 1:n){  
      for(k in 0:m){  
        A[i, k + 1] = b.km(St1[i], k, m, cj, uj)  
      }  
    }  
    
    db.km = function(t, k, m, cj, uj){  
      if ((t - cj) / (uj - cj) < 0 || (t - cj) / (uj - cj) > 1) {  
        return(0)  
      }  
      dt = 1 / (uj - cj)  
      term1 = choose(m, k) * k * ((t - cj) / (uj - cj))^(k - 1) * dt * (1 - (t - cj) / (uj - cj))^(m - k)  
      term2 = choose(m, k) * ((t - cj) / (uj - cj))^k * (m - k) * (-dt) * (1 - (t - cj) / (uj - cj))^(m - k - 1)
      y = term1 + term2  
      return(y)  
    }  
    
    D_A = matrix(0, nrow = n, ncol = m + 1)  
    
    for(i in 1:n){  
      for(k in 0:m){  
        D_A[i, k + 1] = db.km(St1[i], k, m, cj, uj)  
      }  
    }  
    
    
    #####lamda2 #####
    b.km1 = function(t, k, m, cj1, uj1){  
      if ((t - cj1) / (uj1 - cj1) < 0 || (t - cj1) / (uj1 - cj1) > 1) {  
        return(0)  
      }  
      y1 = choose(m, k) * ((t - cj1) / (uj1 - cj1))^k * (1 - (t - cj1) / (uj1 - cj1))^(m - k)  
      return(y1)  
    }  
    
    n <- length(St2)  
    cj1 <- min(St2)- 0.01  
    uj1 <- max(St2)+0.01  
    v=0.3
    m =floor(n^v) 
    B = matrix(0, nrow = n, ncol = m + 1)  
    for(i in 1:n){  
      for(k in 0:m){  
        B[i, k + 1] = b.km1(St2[i], k, m, cj1, uj1)  
      }  
    }  
    
    db.km1 = function(t, k, m, cj1, uj1){  
      if ((t - cj1) / (uj1 - cj1) < 0 || (t - cj1) / (uj1 - cj1) > 1) {  
        return(0)  
      }  
      
      scale_factor <- 1 / (uj1 - cj1)  
      term1 <- choose(m, k) * k * ((t - cj1) / (uj1 - cj1))^(k - 1) * scale_factor *  
        (1 - (t - cj1) / (uj1 - cj1))^(m - k)  
      term2 <- choose(m, k) * ((t - cj1) / (uj1 - cj1))^k * (m - k) *   
        (-scale_factor) * (1 - (t - cj1) / (uj1 - cj1))^(m - k - 1)  
      
      y1_prime <- term1 + term2  
      return(y1_prime)  
    }  
    
    D_B = matrix(0, nrow = n, ncol = m + 1)  
    for(i in 1:n){    
      for(k in 0:m){    
        D_B[i, k + 1] = db.km1(St2[i], k, m, cj1, uj1)    
      }    
    } 
    
    ###############Likelihood###############################
    logL <- function(par){  
      beta1<-par[1]
      beta2<-par[2]
      beta3 <- exp(par[3])
      phi00<-exp(par[4])/(1+exp(par[4]))
      phi01<-1-phi00-exp(par[5])*(1-phi00)/(1+exp(par[5]))
      phi10 <- 1 - phi00 - phi01 - exp(par[6]) * (1 - phi00 - phi01) / (1 + exp(par[6]))  
      phii1<-par[7:(7+m)]
      phii2<-par[(8+m):(8+2*m)]
      phi11=1 - phi00 - phi10 - phi01 
      phi1=cumsum(exp(phii1))
      phi2=cumsum(exp(phii2))
      Dlamda1=rowSums(A%*%phi1)
      Dlamda2=rowSums(B%*%phi2)
      lamda2=rowSums(D_B%*%phi2)
      lamda1=rowSums(D_A%*%phi1) 
      
      pian1  <- phi11*exp(-Dlamda1*exp(SZ1*beta1)-Dlamda2*exp(SZ2*beta2)-beta3*pmax(St1,St2))*(-lamda1 *exp(SZ1 * beta1 )- beta3*(St1 >= St2)) + phi10 * exp(-Dlamda1 *exp(SZ1 * beta1)  - beta3 * St1)*(-lamda1 *exp(SZ1 * beta1 )- beta3) 
      
      pian2 <- phi11*exp(-Dlamda1*exp(SZ1*beta1)-Dlamda2*exp(SZ2*beta2)-beta3*pmax(St1,St2))*(-lamda2 *exp(SZ2 * beta2 )- beta3*(St2 >= St1))+phi01 * exp(-Dlamda2*exp(SZ2*beta2)-beta3*St2)*(-lamda2 *exp(SZ2* beta2 )- beta3)
      
      St1t2<- phi11*exp(-Dlamda1*exp(SZ1*beta1)-Dlamda2*exp(SZ2*beta2)-beta3*pmax(St1,St2))+phi10*exp(-Dlamda1*exp(SZ1*beta1)-beta3*St1)+phi01*exp(-Dlamda2*exp(SZ2*beta2)-beta3*St2)+phi00
      ft1t2 <- phi11*exp(-Dlamda1*exp(SZ1*beta1)-Dlamda2*exp(SZ2*beta2)-beta3*pmax(St1,St2))*(-lamda1 *exp(SZ1 * beta1)- beta3*(St1 >= St2))*(-lamda2 *exp(SZ2 * beta2 )- beta3*(St2 >= St1))  
      
      piandao1_value<- -pian1
      piandao2_value <- -pian2
      piandao1_value<-ifelse(piandao1_value==0,1e-4,-pian1) 
      piandao2_value <-ifelse(piandao2_value==0,1e-4,-pian2)
      
      ft1t2_shifted <-  ft1t2
      likelognorm=sum(Sdelta1*Sdelta2*log(ft1t2_shifted )+Sdelta1*(1-Sdelta2)*log(piandao1_value)+Sdelta2*(1-Sdelta1)*log(piandao2_value)
                      +(1-Sdelta1)*(1-Sdelta2)*log(St1t2))
      return(-likelognorm)  
    }
    
    B_out <- try(optim(par = initial_values, fn = logL, control = list(maxit = 12000)), silent = TRUE)
    if ('try-error' %in% class(B_out) || B_out$convergence != 0) next
    
    B_beta1_estimates[ii] <- B_out$par[1]  
    B_beta2_estimates[ii] <- B_out$par[2]  
    B_beta3_estimates[ii] <- exp(B_out$par[3])
    B_phi00_estimates[ii] <- exp(B_out$par[4]) / (1 + exp(B_out$par[4]))  
    B_phi01_estimates[ii] <- 1 - B_phi00_estimates[ii] - exp(B_out$par[5]) * (1 - B_phi00_estimates[ii]) / (1 + exp(B_out$par[5]))  
    B_phi10_estimates[ii] <- 1 - B_phi00_estimates[ii] - B_phi01_estimates[ii] - exp(B_out$par[6]) * (1 - B_phi00_estimates[ii] - B_phi01_estimates[ii]) / (1 + exp(B_out$par[6])) 
    B_phi11_estimates[ii] <- 1 - B_phi00_estimates[ii] - B_phi10_estimates[ii] - B_phi01_estimates[ii] 
    B_rho1_estimates <- B_phi00_estimates + B_phi01_estimates
    B_rho2_estimates <- B_phi00_estimates + B_phi10_estimates
    ####
    B_estimates_matrix[ii, ] <- c(
      B_beta1_estimates[ii],  
      B_beta2_estimates[ii],  
      B_beta3_estimates[ii],  
      B_rho1_estimates[ii],  
      B_rho2_estimates[ii]  )
  }######################
  ##############ese
  ese <- apply(B_estimates_matrix, 2, sd)
  ese.repe[convergence_count, ] <- ese
  upper <- est + qnorm(0.975) * ese
  lower <- est - qnorm(0.975) * ese
  interval <- data.frame(value=true_params, upper=upper, lower=lower) 
  cp.delta <- ((lower < true_params) & (true_params < upper))
  cp.delta.repe[convergence_count, ] <- cp.delta
}
###################

# bias
bias <- c(  
  mean(beta1_estimates) - true_params["beta1"],  
  mean(beta2_estimates) - true_params["beta2"],  
  mean(beta3_estimates) - true_params["beta3"],
  mean(rho1_estimates) - true_params["rho1"],
  mean(rho2_estimates) - true_params["rho2"])

#Sample standard deviation
SSE <- apply(estimates_matrix[1:convergence_count, ], 2, sd)
###################
ESE<- colMeans(ese.repe[1:convergence_count, ])
CP <- colMeans(cp.delta.repe[1:convergence_count, ])
end_time <- Sys.time()
time=end_time-start_time
##########Print results###########
#(bias,Sample standard deviation,Theoretical standard deviation,CP)
cbind(bias, SSE, ESE, CP)
cat("The running time of the code is：","\n")
time
print(paste("Number of iterations until successful convergence:", convergence_count))
#######
jieguo<- cbind(bias, SSE, ESE, CP)
write.csv(jieguo, file = "table1.csv", row.names = FALSE)
