require("VGAM")
require("bbmle")

LANGS=data
NROWS=length(LANGS)
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

f1 <- function(x, i) {
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
  t3$lambda[i] = best_lambda
  print(sprintf('f1: lambda = %f', best_lambda))
  t3 <<- t3
  return(get_AIC(m2ll, 1, N))
}
f2 <- function(x, i) {
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
  t3$q[i] = best_q
  t3 <<- t3
  #print(sprintf('f2: q = %f', best_q))
  return(get_AIC(m2ll, 1, N))
}



f3 <- function(x) {
  N = length(x)
  m2ll = -2*(-2*sum(log(x)) - N*log(pi**2/6))
  return(get_AIC(m2ll, 0, N))
}

get_AIC <- function(m2logL, K, N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

f4 <- function(x, i) {
  N = length(x)
  mll_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
  }
  ll <- mle2(
    mll_zeta,
    start = list(gamma = 2),
    method = "L-BFGS-B",
    lower = c(1.0000001),
    upper = c(10)
  )
  attr = attributes(summary(ll))
  m2ll = attr$m2logL
  best_gamma = attr$coef[1]
  t3$gamma1[i] = best_gamma
  t3 <<- t3
  #print(sprintf('f4: gamma = %f', best_gamma))
  return(get_AIC(m2ll, 1, N))
}
f5 <- function(x, i) {
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
  #print(sprintf('f5: gamma = %f kmax = %f', best_gamma, best_kmax))
  t3$gamma2[i] = best_gamma
  t3$kmax[i] = best_kmax
  t3 <<- t3
  return(get_AIC(m2ll, 2, N))
}

t3 = data.frame(
  #lang=LANGS,
  lambda=numeric(NROWS),
  q=numeric(NROWS),
  gamma1=numeric(NROWS),
  gamma2=numeric(NROWS),
  kmax=numeric(NROWS)
)
t4 = data.frame(
  #lang=LANGS,
  a1=numeric(NROWS),
  a2=numeric(NROWS),
  a3=numeric(NROWS),
  a4=numeric(NROWS),
  a5=numeric(NROWS)
)


# Fill a row
table_row <- function(lang, ds, i) {
  #d$lang[i] = lang # Doesn't work (?)
  fill_row(ds$V1, i)
}

# Build the table
# for(i in seq(NROWS)){
#   lang = LANGS[i]
#   DB = sprintf('data/%s_in-degree_sequence.txt', lang)
#   ds = read.table(DB, header = FALSE)
#   table_row(lang, ds, i)
# }

#No-Growth
version = "No vertex growth"
name = "no-growth-10000-5000-50-degree-dist"
full.name = paste0("data/", name, ".txt")
data <- read.degree.dist.raw(full.name)
AIC_no_growth<-function(x){
  AIC_ng<-numeric(5)
  AIC_ng[1]<-f1(x)
  AIC_ng[2]<-f2(x)
  AIC_ng[3]<-f3(x)
  AIC_ng[4]<-f4(x)
  AIC_ng[5]<-f5(x)
  return(AIC_ng)
}

AIC_ng<-AIC_no_growth(data$degree)
