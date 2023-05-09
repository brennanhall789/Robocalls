
############ Required Packages #############
############################################
library("sas7bdat")
library("survival")
library("ggplot2")
library("survminer")
library("ggfortify")
library("timereg")
library("gridExtra")
library("coxphw")
library("dplyr")
library("survRM2")
############################################

DAT<-read.sas7bdat("N:/press/moy/sas/sas_analysis.sas7bdat")
names(DAT)


DAT[1,]

table(DAT$RSPCDE00)
table(DAT$X_LI00)
length(DAT$NAICS00)




#############################################
### Initial variables
#############################################
# trt_A: Panel A
# trt_B: Panel B
# trt_C: control (ie. everyone else) - data available? 

initial_mail <- as.Date("02/10/22","%m/%d/%y") 
due_date <- as.Date("03/22/22","%m/%d/%y") 
second_mail <- as.Date("05/09/22","%m/%d/%y")
tfu_start <- as.Date("05/17/22","%m/%d/%y")
tfu_end <- as.Date("06/03/22","%m/%d/%y")
third_email <-as.Date("06/16/22","%m/%d/%y")
end_date <- as.Date("06/15/22","%m/%d/%y")

## check actual dates
trt_A<-as.Date(c("05/18/22"),"%m/%d/%y")
trt_B<-as.Date(c("05/31/22"),"%m/%d/%y")


## Set start date (run program from beg)
#start <- initial_mail 
# start <- due_date
start <- second_mail


order0<-c(0,0,0,0)
order1<-c(1,2,1,2)
order2<-c(2,1,2,1)


#############################################
### Clean up data
#############################################

# option to use response_code vs check_date 
use_response <- T
# option to restrict timeframe to stop at third_email
restrict_time <- F
# option to add tfu to trt column for 'control' group
use_tfu <- F

surv<-data.frame(as.Date(DAT$CKNDTE00,origin="1960-01-01"), 
                 DAT$X_LI00, DAT$WGT00, substr(DAT$NAICS00,1,3), DAT$EMOS_E00,DAT$X_TFU00,DAT$RSPCDE00)
names(surv)<-c("check_date","trt","wgt","naics3","mos","tfu","response")
if(use_response){
  surv$resp <- ifelse(surv$response=='Y',1,0)
  surv$resp[is.na(surv$check_date)]<-0
  surv$check_date[is.na(surv$check_date)]<-max(surv$check_date,na.rm=T)
  surv$check_date[which(surv$response != 'Y' & surv$check_date > 0)] <- max(surv$check_date,na.rm=T)
} else{   
  surv$resp<-1
  surv$resp[is.na(surv$check_date)]<-0
  surv$check_date[is.na(surv$check_date)]<-max(surv$check_date,na.rm=T)
}


surv$time<-surv$check_date-start       
surv$cert<-1
surv$cert[surv$wgt>1]<-0
surv$cert<-as.factor(surv$cert)

if(use_tfu){
  levels(surv$trt) <- c(levels(surv$trt), "T")
  surv$trt[which(surv$tfu == 'Y')] <- "T"
}

end_time <- ifelse(restrict_time, end_date - start, max(surv$check_date,na.rm=T)-start)

#remove any without panel assignment
dim(surv)
surv<-surv[surv$trt!="",]
dim(surv)

#remove any with negative response time
dim(surv)
surv<-surv[surv$time>0,]
dim(surv)

# restrict analysis time frame to before 06/16/22 (day 38)
if (restrict_time){
  surv$time <- ifelse(surv$time >= 38, 38, surv$time)
  surv$resp <- ifelse(surv$time >= 38, 0, surv$resp)
  Surv(surv$time, surv$resp)
}

surv$id<-1:nrow(surv)


#############################################
### Survial models
#############################################

bottom <- -.05+(length(surv$resp)-sum(surv$resp))/length(surv$resp)
top <- .05+(sum(surv$resp))/length(surv$resp)


Surv(surv$time, surv$resp)

total_fit <- survfit(Surv(time, resp) ~ 1, data = surv, conf.int=.9)
total_fit_wgt <- survfit(Surv(time, resp) ~ 1, data = surv, weights = wgt, conf.int=.9)

trt_fit <- survfit(Surv(time, resp) ~ trt, data = surv, conf.int=.9)
wgt_fit <- survfit(Surv(time, resp) ~ trt, data = surv, weights=wgt, conf.int=.9)

(total_chkin <- sum(total_fit$n.event)/total_fit$n)
(trtA_chkin <- sum(trt_fit[1]$n.event)/trt_fit[1]$n)
(trtB_chkin <- sum(trt_fit[2]$n.event)/trt_fit[2]$n)
if(use_tfu){
  (trtT_chkin <- sum(trt_fit[3]$n.event)/trt_fit[3]$n)
}

# important dates for graphs
trtAstart = data.frame(x_value=as.numeric(trt_A-start))
trtBstart = data.frame(x_value=as.numeric(trt_B-start))
mail2 = data.frame(x_value=as.numeric(second_mail-start))
tfu1 <- data.frame(x_value=as.numeric(tfu_start-start))
tfu2 <- data.frame(x_value=as.numeric(tfu_end-start))
email3 = data.frame(x_value=as.numeric(third_email-start))


p1 <- ggsurvplot(
  fit=total_fit,
  fun="event",
  data = surv, 
  conf.int = TRUE,         
  # point estimaes of survival curves.
  xlim = c(0,end_time),        # present narrower X axis, but not affect
  ylim = c(0,top),
  ylab = "Check-in Probability",
  # survival estimates.
  break.time.by = 7,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
) 
p1$plot + 
  geom_vline(xintercept=as.numeric(trt_A-start),col="red",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(trt_B-start),col="turquoise",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(second_mail-start),col="green",lwd=1,lty=c(2)) +
  geom_vline(xintercept=as.numeric(third_email-start),col="black",lwd=1,lty=c(2))


p2<-ggsurvplot(
  fit=trt_fit,
  fun="event",
  data = surv, 
  conf.int = TRUE,         
  # point estimaes of survival curves.
  xlim = c(0,end_time),        # present narrower X axis, but not affect
  ylim = c(0,top),
  ylab = "Check-in Probability",
  # survival estimates.
  break.time.by = 7,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
) 
p2$plot + 
  geom_vline(xintercept=as.numeric(trt_A-start),col="red",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(trt_B-start),col="turquoise",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(second_mail-start),col="green",lwd=1,lty=c(2)) +
  geom_vline(xintercept=as.numeric(third_email-start),col="black",lwd=1,lty=c(2)) +
  geom_vline(xintercept=as.numeric(tfu_start-start),col="blue",lwd=1,lty=c(6)) +
  geom_vline(xintercept=as.numeric(tfu_end-start),col="blue",lwd=2,lty=c(3))


# Difference in response rate? #
tot_rates <- summary(total_fit)
dat_total_rates <- data.frame(tot_rates$time, 'All',
                      round(100*cbind(1-tot_rates$surv,1-tot_rates$upper,1-tot_rates$lower),1))
names(dat_total_rates)<-c("time","trt","rr_rate","90lower","90upper")
rownames(dat_total_rates)<-NULL
print(dat_total_rates[dat_total_rates$time>=end_time,], row.names = F)


rates <- summary(trt_fit)
dat_rates <- data.frame(rates$time, rates$strata,
                      round(100*cbind(1-rates$surv,1-rates$upper,1-rates$lower),1))
names(dat_rates)<-c("time","trt","rr_rate","90lower","90upper")
rownames(dat_rates)<-NULL
print(dat_rates[dat_rates$time>=end_time-1,], row.names = F)


wgt_rates <- summary(wgt_fit)
dat_wgt_rates <- data.frame(wgt_rates$time, wgt_rates$strata,
                      round(100*cbind(1-wgt_rates$surv,1-wgt_rates$upper,1-wgt_rates$lower),1))
names(dat_wgt_rates)<-c("time","trt","rr_rate","90lower","90upper")
rownames(dat_wgt_rates)<-NULL
print(dat_wgt_rates[dat_wgt_rates$time>=end_time-1,], row.names = F)


##### mean/median time to respond #####
# total response
print(total_fit, print.rmean=T)

# treatment responses 
print(trt_fit, print.rmean=T)
print(wgt_fit, print.rmean=T)


# response time CI
means = c(65.3,66.4)
se = c(0.0270,0.0255)
(ciA = cbind(means[1] - qt(0.9,74-1) * se[1], means[1] + qt(0.9,74-1) * se[1]))
(ciB = cbind(means[2] - qt(0.9,74-1) * se[2], means[2] + qt(0.9,74-1) * se[2]))
(ciT = cbind(means[3] - qt(0.9,74-1) * se[3], means[3] + qt(0.9,74-1) * se[3]))

wmeans = c(11.9,13.9)
wse = c(0.035, 0.041)
(wciA = cbind(wmeans[1] - qt(0.9,38-1) * wse[1], wmeans[1] + qt(0.9,38-1) * wse[1]))
(wciB = cbind(wmeans[2] - qt(0.9,38-1) * wse[2], wmeans[2] + qt(0.9,38-1) * wse[2]))

# response time difference using RMST2 package
trt01 = if_else(surv$trt=='A',1,0)
rmst2(surv$time, surv$resp, trt01, alpha=0.9)

#############################################
### Test for expected differnces 
#############################################

sd1 <- survdiff(Surv(time, resp) ~ trt, data = surv)#,subset = which(surv$trt %in% c('B','T')))
sd1


obs<-sd1$obs
exp<-sd1$exp
n<-sd1$n

p12<-prop.test(x = obs[c(1,2)], n = n[c(1,2)],conf.level=.90, alternative = 'two.sided')
p13<-prop.test(x = obs[c(1,3)], n = n[c(1,3)],conf.level=.90, alternative = 'two.sided')
p23<-prop.test(x = obs[c(2,3)], n = n[c(2,3)],conf.level=.90, alternative = 'two.sided')

stuff<-rbind(p12$conf.int,
             p13$conf.int,
             p23$conf.int
             )

rownames(stuff)<-c("Diff 1-2","Diff 1-3","Diff 2-3")
colnames(stuff)<-c("90lower","90upper")

round(stuff*100,2)


## weighted logrank test



#############################################
### Cox proportional hazard models
#############################################

#A HR < 1 indicates reduced hazard of death whereas a HR > 1 indicates an increased hazard of death
#So our HR = 0.90 implies that around 0.90 times as many C are responding as A and C

cox1<-coxph(Surv(time, resp) ~ trt, data = surv)
summary(cox1)
ggforest(cox1)

cox1<-coxph(Surv(time, resp) ~ trt, data = surv, weights=wgt)
summary(cox1)
ggforest(cox1)

cox1<-coxph(Surv(time, resp) ~ trt + strata(cert) + strata(naics3), data = surv)
summary(cox1)


cox1<-coxph(Surv(time, resp) ~ trt + strata(cert) + strata(naics3), data = surv, weights=wgt)
summary(cox1)


cox1<-coxph(Surv(time, resp) ~ trt + cert + naics3, data = surv)
summary(cox1)
survminer::ggforest(cox1)

cox1<-coxph(Surv(time, resp) ~ trt + cert + naics3, data = surv, weights=wgt)
summary(cox1)
ggforest(cox1)

cox1<-coxph(Surv(time, resp) ~ trt + tt(trt) + strata(cert) + strata(naics3), data = surv)
summary(cox1)


#proportional assumption test
(test.ph <- cox.zph(cox1,terms = F))
plot(test.ph)
ggcoxzph(test.ph)


