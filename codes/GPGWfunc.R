library(AdequacyModel)
library(GenSA)
#library(gsl)
# library(reliaR)

gammainc <- function(x, a, b){
  r1 <- gamma(x)
  if (missing(b)){
    if (isTRUE(all.equal(a, 0)))
      return(r1)
    return(r1*pgamma(a, x, lower=FALSE))
  }
  gammainc(x, a)-gammainc(x, b)
}

# generalized power weibull pdf
f_GPW <- function(par,x){
  alpha = par[1]
  #beta = par[2]
  lambda = par[2]
  gamma = par[3]
  beta = 1
  const = gamma*alpha*beta*lambda
  numerator = x^(gamma-1) * (1+lambda*x^gamma)^(alpha-1) * exp(1-(1+lambda*x^gamma)^alpha)
  denominator = (1-exp(1-(1+lambda*x^gamma)^alpha))^(1-beta)
  const*numerator/denominator
}

## generalized power weibull cdf
F_GPW <- function(par,x){
  alpha = par[1]
  #beta = par[2]
  lambda = par[2]
  gamma = par[3]
  beta = 1
  (1 - exp(1-(1+lambda*x^gamma)^alpha))^beta
}
# GPGW pdf
fdp_gam <-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = par[1] # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[2] # c eh o segundo parametro da distribuicao expweibull
  gamma = par[3] # a eh o terceiro parametro da distribuicao expweibull
  a = par[4] # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}

## GPGW cdf
Fda_gam <-function(par,x){
    alpha = par[1]
    lambda = par[2]
    gamma = par[3]
    a = par[4]
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}

# GNH pdf
fdp_gam1<-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = par[1] # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[2] # c eh o segundo parametro da distribuicao expweibull
  gamma = 1 # a eh o terceiro parametro da distribuicao expweibull
  a = par[3] # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}
# GNH cdf
Fda_gam1 <-function(par,x){
  alpha = par[1]
  lambda = par[2]
  gamma = 1
  a = par[3]
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}

# Exponencial pdf
fdp_EXP <-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = 1 # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[1] # c eh o segundo parametro da distribuicao expweibull
  gamma = 1 # a eh o terceiro parametro da distribuicao expweibull
  a = 1 # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}
# Exponencial cdf
Fda_EXP <-function(par,x){
  alpha = 1
  lambda = par[1]
  gamma = 1
  a = 1
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}

# Gama-Exponencial pdf
fdp_GEXP <-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = 1 # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[1] # c eh o segundo parametro da distribuicao expweibull
  gamma = 1 # a eh o terceiro parametro da distribuicao expweibull
  a = par[2] # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}
# Gama-Exponencial cdf
Fda_GEXP <-function(par,x){
  alpha = 1
  lambda = par[1]
  gamma = 1
  a = par[2]
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}


# Gama-rayleigh pdf
fdp_GR <-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = 1 # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[1] # c eh o segundo parametro da distribuicao expweibull
  gamma = 2 # a eh o terceiro parametro da distribuicao expweibull
  a = par[2] # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}
# Gama-rayleigh cdf
Fda_GR <-function(par,x){
  alpha = 1
  lambda = par[1]
  gamma = 2
  a = par[2]
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}


# Gama-Weibull pdf
fdp_GW <-function(par, x){
  #ATENCAO!! o nome desses vetores precisa ser par e x. caso contrario o pacote nao roda
  alpha = 1 # beta eh o primeiro parametro da distribuicao expweibull
  lambda = par[1] # c eh o segundo parametro da distribuicao expweibull
  gamma = par[2] # a eh o terceiro parametro da distribuicao expweibull
  a = par[3] # a eh o quarto parametro da distribuicao expweibull
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  #a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1) # aqui a fda da expweibull
  ((alpha*lambda*gamma*x^(gamma-1))*((1+lambda*x^gamma)^(alpha-1))*exp(1-(1+lambda*x^gamma)^alpha)/(gamma(a))*(z^(a-1)))
}
# Gama-weibull cdf
Fda_GW <-function(par,x){
  alpha = 1
  lambda = par[1]
  gamma = par[2]
  a = par[3]
  z<- -log(1-exp(1-(1+lambda*x^gamma)^alpha))
  1-gammainc(a,0,z)/gamma(a)  
}

## Burr-X pdf
f_BurrX <- function(par,x){
  alpha = 1#par[1]
  beta = par[1]
  lambda = par[2]
  gamma = 2#par[3]
  const = gamma*alpha*beta*lambda
  numerator = x^(gamma-1) * (1+lambda*x^gamma)^(alpha-1) * exp(1-(1+lambda*x^gamma)^alpha)
  denominator = (1-exp(1-(1+lambda*x^gamma)^alpha))^(1-beta)
  const*numerator/denominator
}
## Burr-X cdf
F_BurrX <- function(par,x){
  alpha = 1#par[1]
  beta = par[1]
  lambda = par[2]
  gamma = 2#par[3]
  (1 - exp(1-(1+lambda*x^gamma)^alpha))^beta
}
## Rayleigh pdf
f_Ray <- function(par,x){
  alpha = 1#par[1]
  beta = 1#par[1]
  lambda = par[1]
  gamma = 2#par[3]
  const = gamma*alpha*beta*lambda
  numerator = x^(gamma-1) * (1+lambda*x^gamma)^(alpha-1) * exp(1-(1+lambda*x^gamma)^alpha)
  denominator = (1-exp(1-(1+lambda*x^gamma)^alpha))^(1-beta)
  const*numerator/denominator
}
## Rayleigh cdf
F_Ray <- function(par,x){
  alpha = 1#par[1]
  beta = 1#par[1]
  lambda = par[1]
  gamma = 2#par[3]
  (1 - exp(1-(1+lambda*x^gamma)^alpha))^beta
}
## exp generalized power weibull pdf
f_EGPW <- function(par,x){
  alpha = par[1]
  beta = par[2]
  lambda = par[3]
  gamma = par[4]
  const = gamma*alpha*beta*lambda
  numerator = x^(gamma-1) * (1+lambda*x^gamma)^(alpha-1) * exp(1-(1+lambda*x^gamma)^alpha)
  denominator = (1-exp(1-(1+lambda*x^gamma)^alpha))^(1-beta)
  const*numerator/denominator
}
## exp generalized power weibull cdf
F_EGPW <- function(par,x){
  alpha = par[1]
  beta = par[2]
  lambda = par[3]
  gamma = par[4]
  (1 - exp(1-(1+lambda*x^gamma)^alpha))^beta
}

## Nadarajah-Haghighi pdf
fdp_nh <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  alpha*lambda*(1+lambda*x)^(alpha-1) * exp(1-(1+lambda*x)^alpha)
}
##Nadarajah-Haghighi cdf
fda_nh <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  1-exp(1-(1+lambda*x)^alpha)
}
##Exp-Nadarajah-Haghighi pdf Lemonte(2013).
fdp_enh <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  beta = par[3]
  const = alpha*beta*lambda
  numerator = (1+lambda*x)^(alpha-1) * exp(1-(1+lambda*x)^alpha)
  denominator = (1-exp(1-(1+lambda*x)^alpha))^(1-beta)
  const*numerator/denominator
}
##Exp-Nadarajah-Haghighi cdf Lemonte(2013).
fda_enh <- function(par,x){
  alpha = par[1]
  lambda = par[2]
  beta = par[3]
  (1-exp(1-(1+lambda*x)^alpha))^beta
}
## Weibull pdf
fdp_we <- function(par,x){
  beta = par[1]
  c = par[2]
  beta*c*exp(-beta*x^c)*x^(c-1)
}
## Weibull cdf
fda_we <- function(par,x){
  beta = par[1]
  c = par[2]    
  1 - exp(-beta*x^c)
}


## Burr XII pdf
f_BXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  c*k*s^(-c)*(1+(x/s)^c)^(-k-1)*x^(c-1)
}
## Burr XII cdf
F_BXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  1-(1+(x/s)^c)^(-k)
}  

## Log-Logistic pdf
f_LL<- function(par, x){
  c = par[1]
  k = 1
  s = par[2]
  c*k*s^(-c)*(1+(x/s)^c)^(-k-1)*x^(c-1)
}
## Log-Logistic cdf
F_LL<- function(par, x){
  c = par[1]
  k = 1
  s = par[2]
  1-(1+(x/s)^c)^(-k)
}  

## Kw Weibull pdf
f_kww <- function(par,x){
  a = par[1]
  b = par[2]
  c = par[3]
  beta = par[4]
  a*b*beta*c*exp(-(beta*x)^c)*(beta*x)^(c-1) * (1-exp(-(beta*x)^c))^(a-1) * (1-(1-exp(-(beta*x)^c))^a)^(b-1)
}

## Kw Weibull cdf
F_kww <- function(par,x){
  a = par[1]
  b = par[2]
  c = par[3]
  beta = par[4]
  1 - (1 - (1 - exp(-(beta*x)^c))^a)^b
}

## Exp-Weibull pdf
f_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  (a * c * (beta * x)^c * (1 - exp(-(beta*x)^c))^a) /(x * (exp((beta*x)^c)-1))
}


## Exp-Weibull cdf
F_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  (1 - exp(-(beta*x)^c))^a
}

  
fdp_we <- function(par,x){
  lambda = par[1]
  k = par[2]
  k/lambda*(x/lambda)^(k-1)*exp(-(x/lambda)^k)
}

fda_we <- function(par,x){
  lambda = par[1]
  k = par[2]    
  1 - exp(-(x/lambda)^k)
}

## zb Burr XII pdf
f_ZB<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  f_BXII(par,x)/gsl_sf_gamma(a)*(-log(1-F_BXII(par,x)))^(a-1)
}
## zb Burr XII cdf
F_ZB<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  z = -log(1-F_BXII(par,x))
  gammainc(a,0,z)/gamma(a)  
}  

## RB Burr XII pdf
f_RB<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  f_BXII(par,x)/gamma(a)*(-log(F_BXII(par,x)))^(a-1)
}
## RB Burr XII cdf
F_RB<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  z = -log(F_BXII(par,x))
  1-gammainc(a,0,z)/gamma(a)  
}  

## Beta Burr XII pdf
f_BetaBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  f_BXII(par,x)/beta(a,b)*F_BXII(par,x)^(a-1)*(1-F_BXII(par,x))^(b-1)
}
## Beta Burr XII cdf
F_BetaBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  pbeta(F_BXII(par,x),a,b)
}  

## Kw Burr XII pdf
f_KwBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  a*b*f_BXII(par,x)*F_BXII(par,x)^(a-1)*(1-F_BXII(par,x)^a)^(b-1)
}
## Kw Burr XII cdf
F_KwBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  1-(1-F_BXII(par,x)^a)^b
}  

## Weibull Burr XII pdf
f_WeiBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  a*b*c*k*s^(-c)*x^(c-1)*exp(-a*((1+(x/s)^c)^k-1)^b)*((1+(x/s)^c)^k-1)^(b-1)*(1+(x/s)^c)^(k-1)#/(1+(x/s)^c)
}  

## Weibull Burr XII cdf
F_WeiBXII<- function(par, x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = par[4]
  b = par[5]
  1-exp(-a*((1+(x/s)^c)^k-1)^b)
}  

## Weibull LL pdf
f_WeiLL<- function(par, x){
  c = par[1]
  k = 1
  m = par[2]
  a = par[3]
  b = par[4]
  s=m^(-1)
  a*b*c*k*s^(-c)*x^(c-1)*exp(-a*((1+(x/s)^c)^k-1)^b)*((1+(x/s)^c)^k-1)^(b-1)*(1+(x/s)^c)^(k-1)#/(1+(x/s)^c)
}  

## Weibull LL cdf
F_WeiLL<- function(par, x){
  c = par[1]
  k = 1
  m = par[2]
  a = par[3]
  b = par[4]
  s= m^(-1)
  1-exp(-a*((1+(x/s)^c)^k-1)^b)
}  

## Nikulin-Haghighi pdf
fdp_n_h <- function(par,x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = 1
  b = 1
  a*b*c*k*s^(-c)*x^(c-1)*exp(-a*((1+(x/s)^c)^k-1)^b)*((1+(x/s)^c)^k-1)^(b-1)*(1+(x/s)^c)^(k-1)#/(1+(x/s)^c)
}

## Nikulin-Haghighi cdf

fda_n_h <- function(par,x){
  c = par[1]
  k = par[2]
  s = par[3]
  a = 1
  b = 1
  1-exp(-a*((1+(x/s)^c)^k-1)^b)
}

## Weibull Lomax pdf
f_WeiL<- function(par, x){
  c = 1
  k = par[1]
  s = par[2]
  a = par[3]
  b = par[4]
  a*b*c*k*s^(-c)*x^(c-1)*exp(-a*((1+(x/s)^c)^k-1)^b)*((1+(x/s)^c)^k-1)^(b-1)*(1+(x/s)^c)^(k-1)#/(1+(x/s)^c)
}  

## Weibull Lomax cdf
F_WeiL<- function(par, x){
  c = 1
  k = par[1]
  s = par[2]
  a = par[3]
  b = par[4]
  1-exp(-a*((1+(x/s)^c)^k-1)^b)
}  

## Weibull Nadarajah-Haghighi cdf
F_WNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  1-exp(-a*(exp((1+l*x)^alpha-1)-1)^b)
}

## Weibull Nadarajah-Haghighi pdf
f_WNH <- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  a*b*l*alpha*(1+l*x)^(alpha-1)*(1-exp(1-(1+l*x)^alpha))^(b-1)*
    exp(-b*(1-(1+l*x)^alpha)-a*(exp((1+l*x)^alpha-1)-1)^b)
}

## Weibull exponential cdf
F_WExp<- function(par, x){
  a = par[1]
  b = par[2]
  l = par[3]
  1-exp(-a*(exp(l*x)-1)^b)
}

## Weibull exponential pdf
f_WExp <- function(par, x){
  a = par[1]
  b = par[2]
  alpha = 1
  l = par[3]
  a*b*l*(1-exp(-l*x))^(b-1)*
    exp(b*l*x-a*(exp(l*x)-1)^b)
}

## Gompertz cdf
F_Go<- function(par, x){
  #if(prod(as.numeric(par>0))){
    theta = par[1]
    b = 1
    alpha = 1
    l = par[2]
    1-exp(-(theta/l)*(exp((1+l*x)^alpha-1)-1)^b)
  } #else 0}


## Gompertz pdf
f_Go <- function(par, x){
  #if(prod(as.numeric(par>0))){
  theta = par[1]
  b = 1
  alpha = 1
  l = par[2]
  (theta/l)*b*l*alpha*(1+l*x)^(alpha-1)*(1-exp(1-(1+l*x)^alpha))^(b-1)*
    exp(-b*(1-(1+l*x)^alpha)-(theta/l)*(exp((1+l*x)^alpha-1)-1)^b)
} #else 0}

## EG Nadarajah-Haghighi cdf OK
F_EGNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  (1-exp(1-(1+l*x)^alpha)^a)^b
}

## EG Nadarajah-Haghighi pdf  OK
f_EGNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  a*b*l*alpha*(1+l*x)^(alpha-1)*exp(1-(1+l*x)^alpha)^a/
    (1-exp(1-(1+l*x)^alpha))^(1-b)
}

## Kw Nadarajah-Haghighi cdf  OK
F_KwNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  1-(1-(1-exp(1-(1+l*x)^alpha))^a)^b
}

## Kw Nadarajah-Haghighi pdf OK
f_KwNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  a*b*l*alpha*(1+l*x)^(alpha-1)*exp(1-(1+l*x)^alpha)*
    (1-exp(1-(1+l*x)^alpha))^(a-1)/(1-(1-exp(1-(1+l*x)^alpha))^a)^(1-b)
}

## Beta Nadarajah-Haghighi cdf  
F_BNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  pbeta(fda_nh(c(alpha,l),x),a,b)
}

## Beta Nadarajah-Haghighi pdf  
f_BNH<- function(par, x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  l*alpha*(1+l*x)^(alpha-1)*(1-exp(1-(1+l*x)^alpha))^(a-1)*
    exp(1-(1+l*x)^alpha)^b/beta(a,b)
}


## ZB Nadarajah-Haghighi cdf  
F_ZBNH<- function(par, x){
  a = par[1]
  alpha = par[2]
  l = par[3]
  gammainc(a,0,(1+l*x)^(alpha)-1)/gamma(a)  
}

## ZB Nadarajah-Haghighi pdf  
f_ZBNH<- function(par,x){
 a = par[1]
 alpha = par[2]
 l = par[3]
  l*alpha*(1+l*x)^(alpha-1)*((1+l*x)^alpha-1)^(a-1)*
    exp(1-(1+l*x)^alpha)/gamma(a)
}

## Logistic Nadarajah-Haghighi cdf  
F_LNH<- function(par, x){
  g = par[1]
  a = par[2]
  l = par[3]
  ((1+l*x)^a-1)^g/(1+((1+l*x)^a-1)^g)
  }

## Logistic Nadarajah-Haghighi pdf  
f_LNH<- function(par,x){
  g = par[1]
  a = par[2]
  l = par[3]
  g*a*l*(1+l*x)^(a-1)*((1+l*x)^a-1)^(g-1)/(1+((1+l*x)^a-1)^g)^2
}

## Logistic exponential cdf  
F_Lexp<- function(par, x){
  g = par[1]
  a = 1
  l = par[2]
  ((1+l*x)^a-1)^g/(1+((1+l*x)^a-1)^g)
}

## Logistic exponential pdf  
f_Lexp<- function(par,x){
  g = par[1]
  a = 1
  l = par[2]
  g*a*l*(1+l*x)^(a-1)*((1+l*x)^a-1)^(g-1)/(1+((1+l*x)^a-1)^g)^2
}
## modified Weibull pdf
f_mw <- function(par,x){
    alpha = par[1]
    beta = par[2]
    lambda = par[3]
    lambda*beta*(x/alpha)^(beta-1)*exp((x/alpha)^beta+lambda*alpha*(1-exp((x/alpha)^beta)))
}
## modified Weibull cdf
F_mw <- function(par,x){
    alpha = par[1]
    beta = par[2]
    lambda = par[3]
    1 - exp(lambda*alpha*(1-exp(x^beta/alpha^beta)))
}
## Chen pdf
f_chen <- function(par,x){
    lambda = par[1]
    beta = par[2]
    exp(x^beta+(1-exp(x^beta))*lambda)*x^(-1+beta)*beta*lambda
}
## Chen cdf
F_chen <- function(par,x){
    lambda = par[1]
    beta = par[2]
    1-exp(lambda*(1-exp(x^beta)))
}
## flexible Weibull pdf
f_flexw <- function(par,x){
    alpha = par[1]
    beta = par[2]
    (alpha+beta/x^2)*exp(alpha*x-beta/x)*exp(-exp(alpha*x-beta/x))
}
## flexible Weibull cdf
F_flexw <- function(par,x){
    alpha = par[1]
    beta = par[2]
    1-exp(-exp(alpha*x-beta/x))
}

## Marsahll Olkin Nadarajah-Haghighi
fda_monh <- function(par,x){
  l=par[1]
  a=par[2]
  b=par[3]
  (1-exp(1-(1+l*x)^a))/(1-(1-b)*exp(1-(1+l*x)^a))
}

fdp_monh <- function(par,t){
  l=par[1]
  a=par[2]
  b=par[3]
  a*b*l*((1+l*x)^(a-1)*exp(1-(1+l*x)^a))/(1-(1-b)*exp(1-(1+l*x)^a))^2
}

## Beta Weibull
F_BW<- function(par,x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  pbeta(1-exp(-(x/l)^alpha),a,b)
}


f_BW<- function(par,x){
  a = par[1]
  b = par[2]
  alpha = par[3]
  l = par[4]
  c=alpha
  g=l
  1/beta(a,b)*(c/g)*(x/g)^(c-1)*(1-exp(-(x/g)^c))^(a-1)*exp(-b*(x/g)^c)
}

fit.sa2 <- function(data,density) {
  minusllike <- function(x) -sum(log(density(c(x[1],x[2]),data)))
  lower <- c(0.001,0.001) #may need some changes here
  upper <- c(10,10)
  out <- GenSA(lower = lower, upper = upper,
               fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}

fit.sa3 <- function(data,density) {
  minusllike <- function(x) -sum(log(density(c(x[1],x[2],x[3]),data)))
  lower <- rep(.1,3) #may need some changes here
  upper <- rep(10,3)
  out <- GenSA(lower = lower, upper = upper,
               fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}

fit.sa4 <- function(data,density) {
  minusllike <- function(x) -sum(log(density(c(x[1],x[2],x[3],x[4]),data)))
  lower <- c(0.001,0.001,0.001,0.001) #may need some changes here
  upper <- c(10,10,10,10)
  out <- GenSA(lower = lower, upper = upper,
               fn = minusllike, control=list(verbose=TRUE,max.time=2))
  return(out[c("value","par","counts")])
}