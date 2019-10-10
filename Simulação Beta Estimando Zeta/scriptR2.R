rm(list=ls(all=TRUE))
require(rstan)

set.seed(2018)

# Gerando dados de uma distribuição de uma regressão:
n = 1000 # tamanho amostral
x1 = rbinom(n,1,0.5) # covariável 1
x2 = runif(n,-1,1) # covariável 2
beta0_real = 1 # intercepto real
beta1_real = 1.5 # efeito da covariável x1
beta2_real = -2 # efeito da covariável x2

lpred_real = beta0_real + beta1_real*x1 + beta2_real*x2 
theta_real = exp(lpred_real)/(1+exp(lpred_real))

w1 = rbinom(n,1,0.5)
w2 = runif(n,-0.5,0.5)
gama0_real = 2
gama1_real = 1
gama2_real = -0.5
lpredg_real = gama0_real + gama1_real*w1 + gama2_real*w2 
zeta_real = exp(lpredg_real)

# Parâmetros originais da Beta no software R.
A = theta_real*zeta_real
B = zeta_real-A

# pesquisador tem acesso somente a esta informação.
y = rep(0,n)
for(i in 1:n){
 y[i] = rbeta(1,A[i],B[i]) 
}

dat = cbind(y,x1,x2)

# Hitograma para avaliar a distribuição dos dados gerados
hist(y,20,prob=TRUE,main="Histograma dos dados artificiais",ylab="densidade",cex.lab=2,cex.axis=2)
lines(density(y),lwd=2,col="red")

##############################
# Preparando dados para o Stan

# especificações a priori
# beta0 ~ N(m0,v0)
# beta1 ~ N(m1,v1)
# beta2 ~ N(m2,v2)
# sigma2 ~ Gama(a,b)
# gama0 ~ N(mg0,vg0)
# gama1 ~ N(mg1,vg1)
# gama2 ~ N(mg2,vg2)
m0 = 0; v0 = 10
m1 = 0; v1 = 10
m2 = 0; v2 = 10
mg0 = 10; vg0 = 0.1
mg1 = 10; vg1 = 0.1
mg2 = 10; vg2 = 0.1
# Densidades de probabilidade a priori.
# Expressa a opinião do pesquisador sobre mu e sigma2
# antes de observar os dados.
beta0 = seq(-10,10,0.1)
beta1 = seq(-10,10,0.1)
beta2 = seq(-10,10,0.1)
gama0 = seq(-10,10,0.1)
gama1 = seq(-10,10,0.1)
gama2 = seq(-10,10,0.1)
dm = dnorm(beta0,m0,sqrt(v0))
ds = dnorm(gama0,mg0,sqrt(vg0))
par(mfrow=c(1,2))
plot(beta0,dm,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2)
plot(gama0,type="l",ds,lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2)

# Organizando em lista todos os elementos que o Stan vai usar para trabalhar.
data = list(n=n, y=y, x1=x1, x2=x2, m0=m0, v0=v0, m1=m1, v1=v1, m2=m2, v2=v2,
            w1=w1, w2=w2, mg0=mg0, vg0=vg0, mg1=mg1, vg1=vg1, mg2=mg2, vg2=vg2)
pars = c("beta0","beta1","beta2","gama0","gama1","gama2","zeta","theta") # nome dos parâmetros alvo.

# Chutes iniciais para mu e sigma2
init = list()
init[[1]] = list(beta0=0,beta1=0,beta2=0,gama0=0,gama1=0,gama2=0) 


# Aspectos relacionados ao algorítmo
iter = 10000 # total de iterações (incluindo warm-up).
warmup = 5000
chains = 1

output = stan(file = "scriptSTAN2.stan", data=data, 
              iter=iter, warmup=warmup, chains=chains, 
              pars=pars, init=init, verbose=FALSE)

# Explorando resultados
print(output, pars=c("beta0","beta1","beta2","gama0","gama1","gama2","zeta[1]","theta[1]"))

samp = extract(output)
sbeta0 = samp$beta0; 
szeta = samp$zeta;
sbeta1 = samp$beta1;
sbeta2 = samp$beta2;
sgama0 = samp$gama0; 
sgama1 = samp$gama1; 
sgama2 = samp$gama2; 

# Traceplots das cadeias (sequência de amostras a posteriori)
par(mfrow=c(1,2))
plot(sbeta0,type="l",ylab="beta0",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=beta0_real,lwd=2,col="red")
plot(sbeta1,type="l",ylab="beta1",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=beta1_real,lwd=2,col="red")
plot(sbeta2,type="l",ylab="beta2",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=beta2_real,lwd=2,col="red")
plot(sgama0,type="l",ylab="gama0",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=gama0_real,lwd=2,col="red")
plot(sgama1,type="l",ylab="gama1",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=gama1_real,lwd=2,col="red")
plot(sgama2,type="l",ylab="gama2",xlab="iterações",cex.lab=1.5,cex.axis=1.5)
abline(h=gama2_real,lwd=2,col="red")

# Densidades a posteriori.
par(mfrow=c(1,2))
plot(density(sbeta0),main="Desnidade a posteriori beta0",col="red")
plot(density(sbeta1),main="Desnidade a posteriori beta1",col="red")
plot(density(sbeta2),main="Desnidade a posteriori beta2",col="red")
plot(density(sgama0),main="Desnidade a posteriori gama0",col="red")
plot(density(sgama1),main="Desnidade a posteriori gama1",col="red")
plot(density(sgama2),main="Desnidade a posteriori gama2",col="red")
plot(density(szeta),main="Densidade a posteriori zeta",col="red")

# Posteriori versus Priori.
par(mfrow=c(1,2))
plot(beta0,dm,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sbeta0),lwd=2,main="Desnidade a posteriori mu",col="red")
plot(beta1,dm,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sbeta1),lwd=2,main="Desnidade a posteriori mu",col="red")
plot(beta2,dm,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sbeta2),lwd=2,main="Desnidade a posteriori mu",col="red")
plot(gama0,ds,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sgama0),lwd=2,main="Desnidade a posteriori mu",col="red")
plot(gama1,ds,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sgama1),lwd=2,main="Desnidade a posteriori mu",col="red")
plot(gama2,ds,type="l",lwd=2,col="blue",ylab="densidade",cex.lab=2,cex.axis=2,ylim=c(0,3.5))
lines(density(sgama2),lwd=2,main="Desnidade a posteriori mu",col="red")

require(coda)
# Estimativa intervalar (intervalo de credibilidade HPD - High Posterior Density)
beta0_mcmc = as.mcmc( c(sbeta0) )
beta1_mcmc = as.mcmc( c(sbeta1) )
beta2_mcmc = as.mcmc( c(sbeta2) )
gama0_mcmc = as.mcmc( c(sgama0) )
gama1_mcmc = as.mcmc( c(sgama1) )
gama2_mcmc = as.mcmc( c(sgama2) )
zeta_mcmc = as.mcmc( c(szeta) )
hpd1 = HPDinterval(beta0_mcmc, prob=0.95)
hpd2 = HPDinterval(beta1_mcmc, prob=0.95)
hpd3 = HPDinterval(beta2_mcmc, prob=0.95)
hpd4 = HPDinterval(gama0_mcmc, prob=0.95)
hpd5 = HPDinterval(gama1_mcmc, prob=0.95)
hpd6 = HPDinterval(gama2_mcmc, prob=0.95)
hpd7 = HPDinterval(zeta_mcmc, prob=0.95)

# Estimativas pontuais
est1 = c(mean(sbeta0),beta0_real,median(sbeta0),var(sbeta0),sd(sbeta0),hpd1[1:2])
est2 = c(mean(sbeta1),beta1_real,median(sbeta1),var(sbeta1),sd(sbeta1),hpd2[1:2])
est3 = c(mean(sbeta2),beta2_real,median(sbeta2),var(sbeta2),sd(sbeta2),hpd3[1:2])

est4 = c(mean(sgama0),gama0_real,median(sgama0),var(sgama0),sd(sgama0),hpd4[1:2])
est5 = c(mean(sgama1),gama1_real,median(sgama1),var(sgama1),sd(sgama1),hpd5[1:2])
est6 = c(mean(sgama2),gama2_real,median(sgama2),var(sgama2),sd(sgama2),hpd6[1:2])

est7 = c(mean(szeta),median(szeta),var(szeta),sd(szeta),hpd7[1:2])
tab = rbind(est1,est2,est3,est4,est5,est6)
colnames(tab) = c("M�dia","VR","Mediana","Vari�ncia","Desvio_Padr�o","HPD_inf","HPD_sup")
rownames(tab) = c("beta0","beta1","beta2","gama0","gama1","gama2")
tab


