#packages
library(survival)
library(ranger)
library(ggplot2)
library(ggfortify)
library(tidyverse)
#=======================================================================================================================================================
#creating the data
Patient_name<-c("John","James","Linda","Arnold","keter","peter","kelvin","bryan","Esset","Devis","William")
cell_type<-c("squamous","adeno","small cell","small cell","adeno","adeno","small cell","squamous","adeno","small cell",
             "adeno")
survival_time<-c(78,46,100,120,125,67,91,65,70,52,95)
censoring_status<-c(1,1,0,1,0,0,1,1,1,1,1)
Age<-c(70,38,45,23,28,67,49,53,23,19,23)
prior_therapy<-c(0,1,0,1,1,0,1,1,1,1,1)
treatment<-c(1,1,2,2,1,2,1,2,2,1,1)
#data frame
data_survival<-data.frame(Patient_name,cell_type,survival_time,censoring_status,Age,prior_therapy,treatment)
View(data_survival)
#=======================================================================================================================================================
#Kaplan meir Analysis
#censoring_status=1 indicates whether death is observed
#a “+” after the survival time in the print out of km indicates censoring.
km<-with(data_survival,Surv(survival_time,censoring_status))
km
#To begin our analysis we use the formula surv(survival_time,censoring_status)~1
#survfit()function is used to produce the km estimates of probability of survival over time
#summary()function gives some control over which times to print
#fitting the km
km_fit<-survfit(Surv(survival_time,censoring_status)~1,data = data_survival)
km_fit
summary(km_fit)
#plotting the kmfit
plot(km_fit,xlab = "Survival Time in days",ylab = "survival",main ="KM_fit graph",col = "red")
autoplot(km_fit,main = "KM_fit graph")
#============================================================================================================================================================
#survival curve by treatment
km_treatment_fit<-survfit(Surv(survival_time,censoring_status)~treatment,data = data_survival)
km_treatment_fit
summary(km_treatment_fit)
#plot
autoplot(km_treatment_fit)
#============================================================================================================================================================
#Now lets check survival by age,create a new variable from a categorical variables that has LT70 and  GT 70
#make treatment and prior therapy into factor variables
data_survival$treatment<-as.factor(data_survival$treatment)
data_survival$prior_therapy<-as.factor(data_survival$prior_therapy)
data_survival$cell_type<-as.factor(data_survival$cell_type)
data_surv<-mutate(data_survival,AG =ifelse((Age < 70),"LT70","GT70"),
                  AG =factor(AG),treatment,prior_therapy,cell_type)
#fit
km_AG_fit<-survfit(Surv(survival_time,censoring_status)~AG,data = data_surv)
summary(km_AG_fit)
#plot
plot(km_AG_fit,col = "red",main ="Survival by Age",xlab = "time",ylab = "surv")
#=============================================================================================================================================================
#cox proportional model
cox<-coxph(Surv(survival_time,censoring_status) ~ cell_type + Age + prior_therapy + treatment,data =data_surv )
summary(cox)
#cox fit
cox_fit<-survfit(cox)
cox_fit
#plot
autoplot(cox_fit)
plot(cox_fit, main = "cph model", xlab="Days")
#=============================================================================================================================================================
aa_fit<-aareg(Surv(survival_time,censoring_status) ~ Age + prior_therapy + treatment,data =data_surv)
aa_fit
summary(aa_fit)  # provides a more complete summary of results
autoplot(aa_fit)
#==============================================================================================================================================================
# ranger model
r_fit <- ranger(Surv(survival_time, censoring_status) ~ treatment + cell_type + 
                  Age + prior_therapy,
                data = data_surv,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)
summary(r_fit)
#=============================================================================================================================================================
# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)
# Plot the survival models for each patient
plot(r_fit$unique.death.times,r_fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Days",
     ylab = "survival",
     main = "Patient Survival Curves")
#=================================================================================================================================================================
cols <- colors()
for (n in sample(c(2:dim(data_surv)[1]), 10)){
  lines(r_fit$unique.death.times, r_fit$survival[n,], type = "l", col = cols[n])
}
lines(death_times, avg_prob, lwd = 2)
legend(500, 0.7, legend = c('Average = black'))















