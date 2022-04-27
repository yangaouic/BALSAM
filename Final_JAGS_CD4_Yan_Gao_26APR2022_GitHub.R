######################################################################
## Project: Bayesian Arc Length Survival Analysis Model           ####
## Script purpose: JAGS Code for Model I in the HIV/AIDS Study    ####
## Date: Apr 26th, 2022                                           ####
## Author: Yan Gao, Assistant Professor in Biostatistics          ####
## Division: Division of Biostatistics, MCW                       ####
## Note: this codes were used to conduct for Model I              ####
######################################################################

######################################################
####      Clean the temporary files          #########
######################################################

#Clean the console
cat("\014") 
dir()
#clean all objects
rm(list = ls())
#System information
sessionInfo()
getwd()
setwd("C:/Users/yagao/Documents")

#######################################################
## Install R packages                                ##
#######################################################
library(dplyr)
#install.packages("tibble")
library(tibble)
#install.packages("reshape")
library(reshape)
#library(lattice)
library("jagsUI")
library(ggplot2)
library(JMbayes)
#install.packages("knitr")
library(knitr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
#install.packages("moments")
library(moments)
library(ggplot2)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library("ggpubr")
#install.packages("ggpmisc")
library(ggpmisc)
library(lattice)
library(survival)
library(nlme)
require(lattice)

######################################################
## Data Importation                                  #
######################################################
data(aids)
data(aids.id)
pbc2=aids
pbc2.id=aids.id
#View(aids)
names(pbc2)

### Calculate the frequency table per time points.
n_obs_per_pat<-data.frame(table(pbc2$patient))
#View(n_obs_per_pat)
names(n_obs_per_pat)
table(n_obs_per_pat$Freq)

par(mfrow=c(2,3),mai=c(1, 1, 1,0.8),
    cex.main=1.1, cex.lab=1.1, cex.axis=1.1, cex=1)
xyplot(CD4 ~ obstime |drug, group = patient, data = pbc2,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, type = "l", col = 1, ...)
       }, xlab = "Time (months)", ylab = "sqrt(CD4)")

#######################################################
## Data manipulation                                 ##
#######################################################
ADS_surv=pbc2.id
id<-as.integer(ADS_surv$patient)
ADS_surv<-pbc2.id %>% arrange(patient)
#View(ADS_surv)
#Calculate failure time per patient
surv_dat<-ADS_surv 
obs.t<-as.numeric(ADS_surv[,"Time"])
#ddC:0; ddI:1 for X
X<-as.numeric(surv_dat[,"drug"]) -1
#Male:1 vs. female:0;
gender<-as.numeric(surv_dat[,"gender"]) -1
#AZT: failure 1 vs. Intolerance: 0;
AZT<-as.numeric(surv_dat[,"AZT"]) -1
#prevOI: AIDS: 1 vs. 0 no AIDS
prevOI<-as.numeric(surv_dat[,"prevOI"]) -1
event<-surv_dat$death
n_sub<-length(unique(id))
n_obs<-nrow(pbc2)
event_Dat<- surv_dat %>%  filter(death==1)
n.event<-nrow(event_Dat)
ADS_long_cov<-pbc2 %>% arrange(patient, obstime)
X_long<-as.numeric(ADS_long_cov[,"drug"]) -1
prevOI_long<-as.numeric(ADS_long_cov[,"prevOI"]) -1
#View(ADS_long_cov)
sub<-as.vector(as.numeric(ADS_long_cov$patient))
#View(sub)
obs_t_mat<-as.vector(obs.t)
oij<-as.vector(ADS_long_cov$obstime)
#View(oij)
Vij<-as.vector(ADS_long_cov$CD4)
#View(Vij)

######################################################
##  Initial Values                                  ##
######################################################
init_values <-function(){
  list( prec.sigma2=1,
        mu = c(1,1), prec.mat =matrix(c(1,0,0,1), byrow=T, nrow=2),
        lambda=0.1,beta=0,alpha=0.1,gamma=0,
        b_gender=0,b_AZT=0, b_prevOI=0, ga_prevOI=0
  )}

######################################################
##  parameters                                      ##
######################################################
parms <- c("lambda", "beta","alpha","gamma",
           "b_gender", "b_AZT","b_prevOI",
           "sigma2","sigma","mu", 
           "HT", "Survival","hT","c_bs","G.T")

###########################################w###########
##  Input constants: covariates, fixed values       ##
######################################################
jags_Dat<-list(Vij= Vij, 
               oij= oij,
               n_sub=n_sub,
               n_obs=n_obs, 
               obs_t_mat=obs_t_mat,
               X=X,
               gender=gender,
               AZT=AZT,
               prevOI=prevOI,
               X_long=X_long,
               zeros=rep(0,n_sub),
               sub=sub,
               event=event)

###########################################w##########
##  Model I of BALSAM                               ##
######################################################
writeLines("
           model{
           #######################################################
           ##   Priors for the longitudinal part                ##
           #######################################################
           
           prec.sigma2~dunif(0.1,10)
           sigma2<-1/prec.sigma2
           lambda~dunif(0,10)
           beta~dnorm(0,10)
           alpha~dunif(-10,10)
           gamma~dnorm(0,10)
       
           #Add three covarites
           b_gender~dunif(-10,10)
           b_AZT~dunif(-10,10)
           b_prevOI~dunif(-10,10)
           
           mu[1] ~ dunif(-10,10)
           mu[2] ~ dunif(-10,10)
           Omega[1,1]<-1
           Omega[2,2]<-0.5
           Omega[1,2]<-0
           Omega[2,1]<-0
           sigma<-inverse(prec.mat)
           prec.mat~dwish(Omega,2)
           
           for (i in 1:n_sub) {

               c_bs[i,1:2]~dmnorm(mu, prec.mat)

               #######################################################
               ##   Survival submodel                              ##
               #######################################################
               G.T[i]<-sqrt(1+c_bs[i,2]^2)*obs_t_mat[i]
               
               ###########################################################
               #Constant baseline hazard and time-independent covariate ##
               ###########################################################
               hT[i]<-lambda*exp(X[i]*beta+b_gender*gender[i]+b_AZT*AZT[i]+b_prevOI*prevOI[i]+alpha*G.T[i])
               HT[i]<-lambda*exp(X[i]*beta+b_gender*gender[i]+b_AZT*AZT[i]+b_prevOI*prevOI[i])*(1/(alpha*sqrt(1+c_bs[i,2]^2) ) )*(exp( sqrt(1+c_bs[i,2]^2)*obs_t_mat[i]*alpha)-1)
               Survival[i]<-exp(-HT[i])
               Surv.Log.Lik[i] <- -100000+event[i]*log(hT[i])-HT[i]
               ###########################################################
               #Zeros trick to work with the likelihood of the JM       ##
               ###########################################################
               zeros[i] ~ dpois(-Surv.Log.Lik[i]) 
           }

           #######################################################
           #Long part: Vij=c1 + c2*oij +gamma*x+ epsilon        ##
           #######################################################
           #Normal distribution
           for (j in 1:n_obs) {
               U[j] <- c_bs[sub[j],1] +c_bs[sub[j],2]*oij[j]+gamma*X_long[j]
               Vij[j]~dnorm(U[j], prec.sigma2) 
           }

           }", con="AIDS_JAGS.txt")  

######################################################
##  Run the modeling BALSAM                         ##
######################################################
set.seed(78910)
surv_jag <- jags(data=jags_Dat,
                 model.file = "AIDS_JAGS.txt",
                 inits = init_values,
                 #seed = 12345,
                 parallel=TRUE,
                 parameters.to.save = parms, 
                 n.adapt = 90000,
                 n.chains = 3, 
                 n.iter = 100000, 
                 n.burnin = 10000,
                 n.thin = 10,
                 DIC = T)

#print(Sys.time()-start.time)
#surv_jag
save(surv_jag, file="CD4_surv_jag.RData")
load("CD4_surv_jag.RData")
summary(surv_jag)
plot(surv_jag)

###########################################w##########
##  Extract the posterior data                      ##
######################################################
post_dat<-surv_jag$sims.list
#View(post_dat)

###########################################w##########
##  Cumulative Hazard                               ##
######################################################
post_HT_dat<-as.data.frame(post_dat$HT)
#summary(post_HT_dat)
min(post_HT_dat)
#mean(post_HT_dat)
HT_mean<-apply(post_HT_dat,2, mean)

###########################################w##########
##  Hazard Rate                                     ##
######################################################
post_hT_dat<-as.data.frame(post_dat$hT)
ht_mean<-apply(post_hT_dat,2, mean)
min(post_hT_dat)

###########################################w##########
##  Survival Probability                            ##
######################################################
post_sur_dat<-as.data.frame(post_dat$Survival)
surv_mean<-apply(post_sur_dat,2, mean)
min(post_sur_dat)

###########################################w##########
##  Cumulative Variation                            ##
######################################################
post_Gt_dat<-as.data.frame(post_dat$G.T)
#View(post_Gt_dat)
Gt_mean<-apply(post_Gt_dat,2, mean)
summary(Gt_mean)
q_95_Gt<-quantile(Gt_mean, 0.95)
q_obs.t<-quantile(obs.t, 0.95)
Gt_pat_list<-Gt_mean[Gt_mean>=q_95_Gt&obs.t<=q_obs.t]

###########################################w##########
##  Mixed effects Model for Longitudinal data       ##
######################################################
c_bs_1<-as.data.frame(post_dat$c_bs[,,1])
c_bs_1_mean<-apply(c_bs_1,2, mean)
c_bs_2<-as.data.frame(post_dat$c_bs[,,2])
c_bs_2_mean<-apply(c_bs_2,2, mean)

dat_post_surv<-cbind(ADS_surv, obs.t, HT_mean, ht_mean, surv_mean, Gt_mean, c_bs_1_mean, c_bs_2_mean) 
#View(dat_post_surv)
dat_post_all<-as.data.frame(dat_post_surv) %>% arrange(patient) 
#View(dat_post_all)
dim(dat_post_all)
#print(dat_post_all)
#tail(dat_post_all)

obs_t_seq<-sort(unique(obs_t_mat))
c_bs_1_pre<-mean(c_bs_1_mean)
c_bs_2_pre<-mean(c_bs_2_mean)
long_mean_outcome<-c_bs_1_pre+c_bs_2_pre*obs_t_seq

###########################################w##########
## Figure 2: Spaghetti plot for each treatment group #
######################################################
theme_set(
  theme_classic() +
    theme(legend.position = c(.5, .5),
          text = element_text(size = 13))
)

gamma_est<-0.182
ddI_intercept<-c_bs_1_pre+gamma_est
pat_ddc<-  pbc2 %>% filter(drug=="ddC")
pat_ddI<-  pbc2 %>% filter(drug=="ddI")

p_ddC <- ggplot(data = pat_ddc, aes(x =obstime , y = CD4, group = patient))+
  geom_line(linetype=4,size=0.1)+
  geom_abline(intercept = c_bs_1_pre, slope = c_bs_2_pre, color="red", 
                 linetype="dashed", size=2)+
  ggtitle(expression(paste("   ", italic(z)," = 7.101 - 0.152", italic(s))))+xlim(0,20)+
  ylab(expression(sqrt(CD4))) + xlab("Measurement Time (months)")+ylim(0,25)

p_ddI <- ggplot(data = pat_ddI, aes(x =obstime , y = CD4, group = patient))+
  geom_line(linetype=4, size=0.1)+
  geom_abline(intercept = ddI_intercept, slope = c_bs_2_pre, color="red", 
              linetype="dashed", size=2)+
  ggtitle(expression(paste(italic(z)," = 7.101 - 0.152", italic(s), " + 0.182")))+xlim(0,20)+
  ylab(expression(sqrt(CD4))) + xlab("Measurement Time (months)") +ylim(0,25)

figure2 <- ggarrange(p_ddI, p_ddC, 
                    labels = c("ddI", "ddC"),
                    ncol = 2, nrow = 1)

figure2


#######################################################
## Figure 4: Fitting curve for four patients         ##
#######################################################

id_list=c("46", "87", "451","390")
names(pbc2)
#View(pbc2)
post_gamma_dat<-as.data.frame(post_dat$gamma)
#View(post_gamma_dat)
gamma_mean<-apply(post_gamma_dat,2, mean)

pbc_pats_extrem<-pbc2 %>% filter(patient %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$patient))

par(mfrow=c(2,2),mai=c(1, 1, 1,0.8),
    cex.main=1.1, cex.lab=1.1, cex.axis=1.1, cex=1)

#Patient 46: ddI, death
pat_46<-pbc_pats_extrem %>% filter(patient %in% "46")
t_46<-round(pat_46$Time[1],digits = 2)
t_seq<-seq(0,t_46, by=0.2)
pred_46<-c_bs_1_mean[46]+t_seq*c_bs_2_mean[46]+gamma_mean*X[46]
drug_46=pat_46$drug[1]
plot(pat_46$obstime, pat_46$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 46"," (",drug_46,")"))
lines(t_seq,pred_46, col="blue",lty = 4,lwd=3)
segments(t_46,-1,t_46,tail(pred_46,1),lwd=2, lty=1, col="red")
Gt_mean_46<-round(Gt_mean[46],digits=2)
Gt_mean_46 #11.53 
t_46 #11.07
text(6, 13, expression(paste(italic(G(t)), " = 11.53", "; ", italic(t), " = 11.07",  " (death)") ))

#Patient 451: ddI, censor
pat_451<-pbc_pats_extrem %>% filter(patient %in% "451")
t_451<-round(pat_451$Time[1],digits = 2)
t_seq<-seq(0,t_451, by=0.2)
pred_451<-c_bs_1_mean[451]+t_seq*c_bs_2_mean[451]+gamma_mean*X[451]
drug_451=pat_451$drug[1]
plot(pat_451$obstime, pat_451$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 451"," (",drug_451,")"))
lines(t_seq,pred_451, col="blue",lty = 4,lwd=3)
segments(t_451,-1,t_451,tail(pred_451,1),lwd=2, lty=5, col="black")
Gt_mean_451<-round(Gt_mean[451],digits=2)
Gt_mean_451 #8.28 
t_451 #7.93
text(10, 13, expression(paste(italic(G(t)), " = 8.28", "; ", italic(t), " = 7.93",  " (censor)") ))

#Paitent 87: death
pat_87<-pbc_pats_extrem %>% filter(patient %in% "87")
t_87<-round(pat_87$Time[1],digits = 2)
t_seq<-seq(0,t_87, by=0.2)
pred_87<-c_bs_1_mean[87]+t_seq*c_bs_2_mean[87]+gamma_mean*X[87]
drug_87=pat_87$drug[1]
plot(pat_87$obstime, pat_87$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 87"," (",drug_87,")"))
lines(t_seq,pred_87, col="blue",lty = 4,lwd=3)
segments(t_87,-1,t_87,tail(pred_87,1),lwd=2, lty=1, col="red")
Gt_mean_87<-round(Gt_mean[87],digits=2)
Gt_mean_87 #10.65
t_87 #10.4
text(6, 13, expression(paste(italic(G(t)), " = 10.65", "; ", italic(t), " = 10.40",  " (death)") ))


#Paitent 390: ddC, censor
pat_390<-pbc_pats_extrem %>% filter(patient %in% "390")
t_390<-round(pat_390$Time[1],digits = 2)
t_seq<-seq(0,t_390, by=0.2)
pred_390<-c_bs_1_mean[390]+t_seq*c_bs_2_mean[390]+gamma_mean*X[390]
drug_390=pat_390$drug[1]
plot(pat_390$obstime, pat_390$CD4, ylim=c(0,15), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 390"," (",drug_390,")"))
lines(t_seq,pred_390, col="blue",lty = 4,lwd=3)
segments(t_390,-1,t_390,tail(pred_390,1),lwd=2, lty=5, col="black")
Gt_mean_390<-round(Gt_mean[390],digits=2)
Gt_mean_390 #12.74
t_390 #12.2
text(6, 8, expression(paste(italic(G(t)), " = 12.74", "; ", italic(t), " = 12.20",  " (censor)") ))



######################################################
## Figure 5; high risk patients illustration         #
######################################################  
par(mfrow=c(2,3),mai=c(1, 1, 0.5,0.5),
    cex.main=1.1, cex.lab=1.1, cex.axis=1.1, cex=1)
id_list=c("130","258","318")

pbc_pats_extrem<-pbc2 %>% filter(patient %in% id_list)
pbc_pats_extrem$ID_cha<-paste("Subject", as.character(pbc_pats_extrem$patient))

#Patient 130
pat_130<-pbc_pats_extrem %>% filter(patient %in% "130")
t_130<-round(pat_130$Time[1],digits = 2)
t_seq<-seq(0,t_130, by=0.2)
pred_130<-c_bs_1_mean[130]+t_seq*c_bs_2_mean[130]+gamma_mean*X[130]
drug_130=pat_130$drug[1]
plot(pat_130$obstime, pat_130$CD4, ylim=c(0,20), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 130"," (",drug_130,")"))
lines(t_seq,pred_130, col="blue",lty = 4,lwd=3)
t_130<-round(pat_130$Time[1],digits = 2)
segments(t_130,-1,t_130,tail(pred_130,1),lwd=2, lty=5, col="black")
Gt_mean_130<-round(Gt_mean[130],digits=2)
Gt_mean_130 #20.30
t_130 #19.20
text(8, 3, expression(paste(italic(G(t)), " = 20.30", "; ", italic(t), " = 19.20") ))
#Patient 258
pat_258<-pbc_pats_extrem %>% filter(patient %in% "258")
t_258<-round(pat_258$Time[1],digits = 2)
t_seq<-seq(0,t_258, by=0.2)
pred_258<-c_bs_1_mean[258]+t_seq*c_bs_2_mean[258]+gamma_mean*X[258]
drug_258=pat_258$drug[1]
plot(pat_258$obstime, pat_258$CD4, ylim=c(0,20), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 258"," (",drug_258,")"))
lines(t_seq,pred_258, col="blue",lty = 4,lwd=3)
segments(t_258,-1,t_258,tail(pred_258,1),lwd=2, lty=5, col="black")
Gt_mean_258<-round(Gt_mean[258],digits=2)
Gt_mean_258 #20.03 
t_258 # 19.27
text(8, 3, expression(paste(italic(G(t)), " = 20.03", "; ", italic(t), " = 19.27") ))
#Patient 318
pat_318<-pbc_pats_extrem %>% filter(patient %in% "318")
t_318<-round(pat_318$Time[1],digits = 2)
t_seq<-seq(0,t_318, by=0.2)
pred_318<-c_bs_1_mean[318]+t_seq*c_bs_2_mean[318]+gamma_mean*X[318]
drug_318=pat_318$drug[1]
plot(pat_318$obstime, pat_318$CD4, ylim=c(0,20), xlim=c(0,20),
     xlab="Time (months)", ylab=expression(sqrt(CD4)), main=paste0("Subject 318"," (",drug_318,")"))
lines(t_seq,pred_318, col="blue",lty = 4,lwd=3)
segments(t_318,-1,t_318,tail(pred_318,1),lwd=2, lty=5, col="black")
Gt_mean_318<-round(Gt_mean[318],digits=2)
Gt_mean_318 #20.35 
t_318 # 18.37
text(8, 3, expression(paste(italic(G(t)), " = 20.35", "; ", italic(t), " = 18.37") ))

###########################################################
# Figure 3: Survival fitting curves                      ##
###########################################################

#Posterior mean from the model 
lambda<-0.009
alpha<-0.061
beta<-0.166
G.T<-sqrt(1+c_bs_2_pre^2)*obs_t_seq
#Make all the other parameters at the reference group;

#Treatment ddI
X1=1
X=X1
hT<-lambda*exp(X*beta+alpha*G.T)
HT<-lambda*exp(X*beta)*(1/(alpha*sqrt(1+c_bs_2_pre^2) ) )*(exp( sqrt(1+c_bs_2_pre^2)*obs_t_seq*alpha)-1)
Survival_X1<-exp(-HT)
hT_1=hT
#Treatment ddC
X0=0
X=X0
hT<-lambda*exp(X*beta+alpha*G.T)
HT<-lambda*exp(X*beta)*(1/(alpha*sqrt(1+c_bs_2_pre^2) ) )*(exp( sqrt(1+c_bs_2_pre^2)*obs_t_seq*alpha)-1)
Survival_X0<-exp(-HT)
hT_0<-hT

par(mfrow=c(1,2),mai=c(1, 1, 1,0.8),
    cex.main=1.1, cex.lab=1.1, cex.axis=1.1, cex=1.1)

plot(obs_t_seq,Survival_X1, type="l",col="black", lwd=2, lty=1, 
     xlab="Time (months)", ylab="Survival Probability", 
     xlim=c(0,22), ylim=c(0.6,1))
lines(obs_t_seq,Survival_X0, type="l",col="red", lwd=2, lty=2, 
     xlim=c(0,22), ylim=c(0.6,1))
legend(0, 0.7, legend=c("Didanosine (ddI)", "Zalcitabine (ddC)"),
       col=c("black", "red"), lty=1:2, lwd=2, bty="n")

plot(obs_t_seq,hT_1, type="l",col="black", lwd=2, lty=1, 
     xlab="Time (months)", ylab="Hazard Rate", 
     xlim=c(0,22), ylim=c(0,0.15))
lines(obs_t_seq,hT_0, type="l",col="red", lwd=2, lty=2, 
      xlim=c(0,22), ylim=c(0,0.15))
legend(5, 0.1, legend=c("Didanosine (ddI)", "Zalcitabine (ddC)"),
       col=c("black", "red"), lty=1:2, lwd=2, bty="n")

######################################################
## End of Model I                                    #
######################################################  
















