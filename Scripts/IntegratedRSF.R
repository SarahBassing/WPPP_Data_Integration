#  Date Integration RSF model 
#  Combining SIMULATED camera trap & telemetry data under one habitat use model
#  Beth Gardner
#  November 2020
#  ---------------------------------------------------------

library(rjags)

###simulate some 'data'
ngrid=200   #number of grid cells
ncam=100     #number of camera traps
n=25        #number of telemetered animals
beta=.3     #covariate parameter

#create a standardized covariate (forest)
forest=rnorm(ngrid, 0, 1)

#Simulate telemetry data one of two ways
if(1 == 1){ #set 1==2 to turn off this 
  #simulate gps data by grid (this is a clunky way to do it)

M<-matrix(NA, n, ngrid) # rows = individual animal, columns = individual grid cells
int<-rnorm(n, 0, .5)
for(i in 1:n){
lambda=exp(int[i] + beta*forest)
M[i,]<-rpois(ngrid,lambda)
}

R<- rowSums(M) # total number of telemetry locations by individual
}


if(1==2){  #set 1==1 to turn this one on if you want
  #Basically simulates data the same as above
M<-matrix(NA, n, ngrid)
pi=matrix(NA, n, ngrid)
int<-rnorm(n, 0, .5)
R<-round(rnorm(n, 250, 50))  # total number of telemetry locations by individual

int<-rnorm(n, 0, .5)
for(i in 1:n){
  lambda=exp(int[i] + beta*forest)
  pi[i,]=lambda/sum(lambda)
  M[i,]<-rmultinom(1, R[i], pi[i,])
}
}

#Simulte camera data
betac<-.35   #you can set this to whatever you want, but should be close
            #to the other beta, you could test how it works when it's different?
intc<-rnorm(ncam, 1, 1)
y<-NULL
for(i in 1:ncam){
  y[i]<-rpois(1, exp(intc[i] + betac*forest[i]))
} 


#Combined model for both gps and camera data
cat("
model{

for(i in 1:n){

alpha[i]~dnorm(int, tau)
M[i,1:ngrid]~dmulti(pi[i,1:ngrid], R[i])

for(j in 1:ngrid){
log(lam[i,j]) = alpha[i] + b1*X[j]
pi[i,j]=lam[i,j]/sum(lam[i,1:ngrid])

}
}


for(k in 1:ncam){
a0[k]~dnorm(intc, tauc)
y[k]~dpois(lamc[k])
log(lamc[k])<-a0[k] + b1*X[k]
}

int~dnorm(0,.1)
tau~dunif(0,10)
b1~dnorm(0, .1)

intc~dnorm(0,.1)
tauc~dunif(0,4)


}

", fill=TRUE, file="combo.txt")



# arguments for jags()
data = list(M=M, R=R, X=forest, ngrid=ngrid, n=n, ncam=ncam, y=y)
parameters = c('alpha','b1', 'int','tau', 'intc', 'tauc')

inits = function() {list(b1=rnorm(1))}


# call to jags
mod <- jags.model("combo.txt", data, inits, n.chains=3, n.adapt=500)
fit <- coda.samples(mod,parameters,n.iter=5000)
summary(fit)
plot(fit)


#################################

#RSF only
cat("
model{

for(i in 1:n){

alpha[i]~dnorm(int, tau)
M[i,1:ngrid]~dmulti(pi[i,1:ngrid], R[i])

for(j in 1:ngrid){

lam[i,j] = alpha[i] + b1*X[j]
hold[i,j]<-exp(lam[i,j])
pi[i,j]=hold[i,j]/sum(hold[i,1:ngrid])

}
}
int~dnorm(0,.1)
tau~dunif(0,10)
b1~dnorm(0, .1)

}

", fill=TRUE, file="rsf.txt")



# arguments for jags()
data = list(M=M, R=R, X=forest, ngrid=ngrid, n=n)
parameters = c('alpha','b1', 'int','tau')

inits = function() {list(b1=rnorm(1))}


# call to jags
mod <- jags.model("rsf.txt", data, inits, n.chains=3, n.adapt=500)
fit <- coda.samples(mod,parameters,n.iter=5000)
summary(fit)
plot(fit)

############################################

#Camera only
cat("
model{

for(k in 1:ncam){
a0[k]~dnorm(intc, tauc)
y[k]~dpois(lamc[k])
log(lamc[k])<-a0[k] + b1*X[k]
}


b1~dnorm(0, .1)

intc~dnorm(0,.1)
tauc~dunif(0,4)


}

", fill=TRUE, file="cam.txt")



# arguments for jags()
data = list( X=forest, ncam=ncam, y=y)
parameters = c('b1', 'intc', 'tauc')

inits = function() {list(b1=rnorm(1))}


# call to jags
mod <- jags.model("cam.txt", data, inits, n.chains=3, n.adapt=500)
fit <- coda.samples(mod,parameters,n.iter=15000)
summary(fit)
plot(fit)