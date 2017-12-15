require("VGAM")
require("bbmle")

LANGS<-c("no-growth", "pref-attach", "rand-attach")
NROWS<-length(LANGS)
# Truncated zeta (very slow)
rzeta <- function(gamma, kmax){
  ret = 0.0
  for(i in seq(kmax)){
    ret = ret + i**(-gamma)
  }
  return(ret)
}

# Compute the AIC corrected for small samples
get_AIC <- function(m2logL, K, N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

f1 <- function(x, indicator=1) {
  N = length(x)
  C = 0
  for (j in seq(N)) {
    C = C + sum(log(seq(2, x[j])))
  }
  mll_dpois <- function(lambda) {
    -(sum(x)*log(lambda)-N*(lambda+log(1-exp(-lambda)))-C)
  }
  ll <- mle2(
    mll_dpois,
    start = list(lambda = 2),
    method = "L-BFGS-B",
    lower = c(1e-5)
  )
  attr = attributes(summary(ll))
  m2ll = attr$m2logL
  best_lambda = attr$coef[1]
  print(sprintf('f1: lambda = %f', best_lambda))
  if(indicator==0){
    return(best_lambda)
  }else{
    return(get_AIC(m2ll, 1, N))
  }
  
  
}

f2 <- function(x, indicator=1) {
  N = length(x)
  M = sum(x)
  mll_dgeo <- function(q) {
    -((M-N)*log(1-q)+N*log(q))
    
  }
  ll <- mle(
    mll_dgeo,
    start = list(q = 2),
    method = "L-BFGS-B",
    lower = c(1e-5),
    upper = c(1-1e-5)
  )
  attr = attributes(summary(ll))
  m2ll = attr$m2logL
  best_q = attr$coef[1]
  #t3$q[i] = best_q
  #t3 <<- t3
  print(sprintf('f2: q = %f', best_q))
  if(indicator==0){
    return(best_q)
  }else{
    return(get_AIC(m2ll, 1, N))
  }
}



f3 <- function(x) {
  N = length(x)
  m2ll = -2*(-2*sum(log(x)) - N*log(pi**2/6))
  
  return(get_AIC(m2ll, 0, N))
}

get_AIC <- function(m2logL, K, N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

f4 <- function(x, indicator=1) {
  N = length(x)
  mll_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
  }
  ll <- mle2(
    mll_zeta,
    start = list(gamma = 3),
    method = "L-BFGS-B",
    lower = c(1.0000001),
    upper = c(10)
  )
  attr = attributes(summary(ll))
  m2ll = attr$m2logL
  best_gamma = attr$coef[1]
  #t3$gamma1[i] = best_gamma
  #t3 <<- t3
  print(sprintf('f4: gamma = %f', best_gamma))
  if(indicator==0){
    return(best_gamma)
  }else{
    return(get_AIC(m2ll, 1, N))
  }
  
}

f5 <- function(x, indicator=1) {
  N = length(x)
  mll_rzeta <- function(gamma, kmax) {
    length(x) * log(rzeta(gamma, kmax)) + gamma * sum(log(x))
  }
  ll <- mle2(
    mll_rzeta,
    start = list(gamma = 2, kmax = 50),
    method = "L-BFGS-B",
    lower = c(1.0000001, 30),
    upper = c(10, 300)
  )
  attr = attributes(summary(ll))
  m2ll = attr$m2logL
  best_gamma = attr$coef[1]
  best_kmax = attr$coef[2]
  print(sprintf('f5: gamma = %f kmax = %f', best_gamma, best_kmax))
  #t3$gamma2[i] = best_gamma
  #t3$kmax[i] = best_kmax
  #t3 <<- t3
  if(indicator==0){
    return(cbind(best_kmax,best_gamma))
  }else{
    return(get_AIC(m2ll, 2, N))
  }
}


#No-Growth
version = "No vertex growth"
name = "no-growth-10000-5000-50-degree-dist"
full.name = paste0("data/", name, ".txt")
data <- read.degree.dist.raw(full.name)
x<-data$degree
remove_zeros<-x!=0
#x<-x[remove_zeros]
new_data<-data[remove_zeros,]
AIC_ng<-AIC_value(new_data$degree)
best_param[1,]<-best_parameters(new_data$degree)


#Pref Attachment
#Preferential Attachment
version = "Preferential Attachment"
name = "pref-attach-10000-15-8-degree-dist"
full.name = paste0("data/", name, ".txt")
data <- read.degree.dist.raw(full.name)

x<-data$degree
remove_zeros<-x!=0
#x<-x[remove_zeros]
new_data<-data[remove_zeros,]
AIC_pref_attach<-AIC_value(new_data$degree)
best_param[2,]<-best_parameters(new_data$degree)



#Random Attachment
version = "Random Attachment"
name = "rand-attach-10000-15-8-degree-dist"
full.name = paste0("data/", name, ".txt")
data <- read.degree.dist.raw(full.name)

x<-data$degree
remove_zeros<-x!=0
#x<-x[remove_zeros]
new_data<-data[remove_zeros,]
AIC_rand_attach<-AIC_value(new_data$degree)
best_param[3,]<-best_parameters(new_data$degree)


AIC_Score<-function(){
  score<-rbind(AIC_ng,AIC_pref_attach,AIC_rand_attach)
  row.names(score)<-LANGS
  colnames(score)<-c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
  return(as.data.frame(score))
}

AIC_delta<-function(){
  best_AIC_NG<-min(AIC_Score()[1,])
  best_AIC_PA<-min(AIC_Score()[2,])
  best_AIC_RA<-min(AIC_Score()[3,]) 
  AIC_score_data<-AIC_Score()
  AIC_score_data[1,]<-AIC_score_data[1,]-best_AIC_NG
  AIC_score_data[2,]<-AIC_score_data[2,]-best_AIC_PA
  AIC_score_data[3,]<-AIC_score_data[3,]-best_AIC_RA
  return(AIC_score_data)
}
AIC_value<-function(x){
  AIC_value<-numeric(5)
  AIC_value[1]<-f1(x)
  AIC_value[2]<-f2(x)
  AIC_value[3]<-f3(x)
  AIC_value[4]<-f4(x)
  AIC_value[5]<-f5(x)
  return(AIC_value)
}

best_param<- as.data.frame(rbind(best_NG<-numeric(5), best_PA<-numeric(5), best_RA<-numeric(5)))
best_parameters<-function(x){
  best_param<-numeric(5)
  best_param[1]<-f1(x,0)
  best_param[2]<-f2(x,0)
  best_param[3]<-f4(x,0)
  best_param[4]<-f5(x,0)[1]
  best_param[5]<-f5(x,0)[2]
  return(best_param)
}

get_best_parameters_table<-function(){
  row.names(best_param)<-c("No Growth","Preferential Attachment","Random Attachment")
  colnames(best_param)<-c("lambda", "q", "gamma_1","gamma_2", "K_max")
  return(as.data.frame(best_param))
}

