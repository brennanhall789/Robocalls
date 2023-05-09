
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

table(DAT$PCFLG00)
table(DAT$STATUS00)
table(DAT$RSPCDE00)
length(DAT$NAICS00)
length(DAT$SICRCD00)
table(substr(DAT$SICRCD00,1,3),substr(DAT$NAICS00,1,3))



#############################################
### Initial variables
#############################################
# trt_A: Panel A
# trt_B: Panel B
# trt_C: control (ie. everyone else) - data available? 

initial_mail <- as.Date("02/10/22","%m/%d/%y") 
due_date <- as.Date("03/22/22","%m/%d/%y") 
second_mail <- as.Date("05/09/22","%m/%d/%y")
third_email <-as.Date("06/16/22","%m/%d/%y")
end_date <- as.Date("06/16/22","%m/%d/%y")

## check actual dates
trt_A<-as.Date(c("05/18/22"),"%m/%d/%y")
trt_B<-as.Date(c("05/31/22"),"%m/%d/%y")

# option to restrict timeframe to stop at third_email
restrict_time <- F

order0<-c(0,0)
order1<-c(1,2)
order2<-c(2,1)

#############################################
### Set start date (run program from beg)
#############################################
# star t<- initial_mail 
# start <- due_date
start <- second_mail

#############################################
### Clean up data
#############################################

# option to use response_code vs check_in 
use_response <- FALSE
# option to restrict timeframe to stop at third_email
restrict_time <- F

surv2<-data.frame(as.Date(DAT$CKNDTE00,origin="1960-01-01"), 
                  DAT$X_LI00, DAT$WGT00, substr(DAT$NAICS00,1,3), DAT$EMOS_E00,DAT$X_TFU00,DAT$RSPCDE00)
names(surv2)<-c("check_date","trt","wgt","naics3","mos","tfu","response")
if(use_response){
  surv$resp <- ifelse(surv$response=='Y',1,0)
  surv$resp[is.na(surv$check_date)]<-0
  surv$check_date[is.na(surv$check_date)]<-max(surv$check_date,na.rm=T)
  surv$check_date[which(surv$response != 'Y' & surv$check_date > 0)] <- max(surv$check_date,na.rm=T)
} else{   
  surv2$check_in<-1
  surv2$check_in[is.na(surv2$check_date)]<-0
  surv2$check_date[is.na(surv2$check_date)]<-max(surv2$check_date,na.rm=T)
}

surv$time<-surv$check_date-start       
surv$cert<-1
surv$cert[surv$wgt>1]<-0
surv$cert<-as.factor(surv$cert)

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
Surv(surv$time, surv$resp)

bottom <- -.05+(length(surv$resp)-sum(surv$resp))/length(surv$resp)
top <- .05+(sum(surv$resp))/length(surv$resp)

total_fit <- survfit(Surv(time, resp) ~ 1, data = surv)
total_fit_wgt <- survfit(Surv(time, resp) ~ 1, data = surv, weights = wgt)

trt_fit <- survfit(Surv(time, resp) ~ trt, data = surv)
wgt_fit <- survfit(Surv(time, resp) ~ trt, data = surv, weights=wgt)

(total_chkin <- sum(total_fit$n.event)/total_fit$n)
(trtA_chkin <- sum(trt_fit[1]$n.event)/trt_fit[1]$n)
(trtB_chkin <- sum(trt_fit[2]$n.event)/trt_fit[2]$n)

# important dates for graphs
trtAstart = data.frame(x_value=as.numeric(trt_A-start))
trtBstart = data.frame(x_value=as.numeric(trt_B-start))
mail2 = data.frame(x_value=as.numeric(second_mail-start))
email3 = data.frame(x_value=as.numeric(third_email-start))

###### unused plots ######
# colors <- c("Robocall trtA" = "red", "Robocall trtB" = "turquoise",
#            "Due Date" = "black", "2nd Mail Follow-up" = "green")

# p1 <- ggsurvplot(
#   fit=total_fit,
#   fun="event",
#   data = surv, 
#   conf.int = F,         
#   # point estimaes of survival curves.
#   xlim = c(0,end_time), 
#   ylim = c(0,top),
#   ylab = "Response Probability",
#   title = 'Overall Response Rate',
#   break.time.by = 7,     
#   ggtheme = theme_minimal()
# ) 
# 
# 
# p1$plot + 
#   geom_vline(data = trtAstart, mapping = aes(xintercept=x_value, color='Robocall trtA'),lwd=1,lty=c(1)) +
#   geom_vline(data = trtBstart, mapping = aes(xintercept=x_value, color='Robocal trtB'),lwd=1,lty=c(1)) +
#   geom_vline(data = due, mapping = aes(xintercept=x_value, color='Due Date'),lwd=1,lty=c(2)) +
#   geom_vline(data = mail2, mapping = aes(xintercept=x_value, color='2nd Mail Follow-up'),lwd=1,lty=c(2)) +
#   labs(color="") +
#   theme(legend.position = 'right') +
#   scale_size_manual(values = colors)
# # +
#   # scale_x_date(date_labels = "%b %d")
# 
# 
# 
# p4 <- ggsurvplot(
#   fit=total_fit_wgt,
#   fun="event",
#   data = surv, 
#   conf.int = F,         
#   # point estimaes of survival curves.
#   xlim = c(0,end_time), 
#   ylim = c(0,top),
#   ylab = "Response Probability",
#   title = 'Overall Response Rate',
#   break.time.by = 7,    
#   ggtheme = theme_minimal()
# ) 
# p4$plot + geom_vline(xintercept=as.numeric(trt_A-start),col="red",lwd=1,lty=c(1,2)) +
#   geom_vline(xintercept=as.numeric(trt_B-start),col="turquoise",lwd=1,lty=c(2,1)) +
#   geom_vline(xintercept=as.numeric(due_date-start),col="black",lwd=1,lty=c(1))# +
#  # scale_x_date(date_labels = "%b %d")
#######################

p2<-ggsurvplot(
  fit=trt_fit,
  fun="event",
  data = surv, 
  conf.int = T,         
  # point estimaes of survival curves.
  xlim = c(0,end_time),        
  ylim = c(0,top),
  ylab = "Response Probability",
  title = 'Response Rate by Treatment',
  break.time.by = 7,    
  ggtheme = theme_minimal()
) 

p2$plot + 
  geom_vline(xintercept=as.numeric(trt_A-start),col="red",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(trt_B-start),col="turquoise",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(second_mail-start),col="green",lwd=1,lty=c(2)) +
  geom_vline(xintercept=as.numeric(third_email-start),col="black",lwd=1,lty=c(2))



p3<-ggsurvplot(
  wgt_fit,                     
  fun="event",
  data = surv, 
  conf.int = TRUE,         
  # point estimaes of survival curves.
  xlim = c(0,end_time),       
  ylim = c(0,top),
  ylab = "Weighted Response Probability",
  title = 'Weighted Response Rate by Treatment',
  break.time.by = 7,     
  ggtheme = theme_minimal() 
) 

p3$plot + 
  geom_vline(xintercept=as.numeric(trt_A-start),col="red",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(trt_B-start),col="turquoise",lwd=1,lty=c(1)) +
  geom_vline(xintercept=as.numeric(second_mail-start),col="green",lwd=1,lty=c(2)) +
  geom_vline(xintercept=as.numeric(third_email-start),col="black",lwd=1,lty=c(2))


##### mean/median time to respond #####
# total response
print(total_fit, print.rmean=T)

# treatment responses 
print(trt_fit, print.rmean=T)
print(wgt_fit, print.rmean=T)


# response time CI
means = c(13.7,14.5)
se = c(0.256, 0.274)
(ciA = cbind(means[1] - qt(0.9,38-1) * se[1], means[1] + qt(0.9,38-1) * se[1]))
(ciB = cbind(means[2] - qt(0.9,38-1) * se[2], means[2] + qt(0.9,38-1) * se[2]))

wmeans = c(11.9,13.9)
wse = c(0.035, 0.041)
(wciA = cbind(wmeans[1] - qt(0.9,38-1) * wse[1], wmeans[1] + qt(0.9,38-1) * wse[1]))
(wciB = cbind(wmeans[2] - qt(0.9,38-1) * wse[2], wmeans[2] + qt(0.9,38-1) * wse[2]))

#response time difference using RMST2 package
trt01 = if_else(surv$trt=='A',1,0)
rmst2(surv$time, surv$resp, trt01, alpha=0.9)


# Difference? #
tot_rates <- summary(total_fit)
dat_total_rates <- data.frame(tot_rates$time, 'All',
                      round(100*cbind(1-tot_rates$surv,1-tot_rates$upper,1-tot_rates$lower),1))
names(dat_total_rates)<-c("time","trt","rr_rate","95lower","95upper")
rownames(dat_total_rates)<-NULL
print(dat_total_rates[dat_total_rates$time>=end_time,], row.names = F)


rates <- summary(trt_fit)
dat_rates <- data.frame(rates$time, rates$strata,
                      round(100*cbind(1-rates$surv,1-rates$upper,1-rates$lower),1))
names(dat_rates)<-c("time","trt","rr_rate","95lower","95upper")
rownames(dat_rates)<-NULL
print(dat_rates[dat_rates$time>=end_time,], row.names = F)


rates_wgt <- summary(wgt_fit)
dat_rates_wgt <- data.frame(rates_wgt$time, rates_wgt$strata,
                      round(100*cbind(1-rates_wgt$surv, 1-rates_wgt$upper, 1-rates_wgt$lower),1))
names(dat_rates_wgt)<-c("time","trt","rr_rate","95lower","95upper")
rownames(dat_rates_wgt)<-NULL
print(dat_rates_wgt[dat_rates_wgt$time>=end_time,], row.names = F)

#############################################
### Test for expected differnces 
#############################################

##### response rate against trt #####
sd1 <- survdiff(Surv(time, resp) ~ trt, data = surv)
sd1


obs <- sd1$obs
exp <- sd1$exp
n <- sd1$n

(p12 <- prop.test(x = obs[c(1,2)], n = n[c(1,2)], conf.level = .9, alternative = 'two.sided'))
# p13<-prop.test(x = obs[c(1,3)], n = n[c(1,3)])
# p23<-prop.test(x = obs[c(2,3)], n = n[c(2,3)])

stuff <- rbind(p12$conf.int)
             # p13$conf.int,
             # p23$conf.int)

rownames(stuff) <- c("Diff A-B")
colnames(stuff) <- c("90lower","90upper")

round(stuff*100,2)


# weighted #
# see weighted coxph model summary


# ##### response rate vs trt + cert + naics #####
# sd1 <- survdiff(Surv(time, resp) ~ trt + strata(cert) + strata(naics3), data = surv)
# sd1
# 
# 
# obs <- rowSums(sd1$obs)
# exp <- rowSums(sd1$exp)
# n <- sd1$n
# 
# p12 < -prop.test(x = obs[c(1,2)], n = n[c(1,2)])
# # p13<-prop.test(x = obs[c(1,3)], n = n[c(1,3)])
# # p23<-prop.test(x = obs[c(2,3)], n = n[c(2,3)])
# 
# stuff <- rbind(p12$conf.int)
#              # p13$conf.int,
#              # p23$conf.int)
# 
# rownames(stuff) <- c("Diff 1-2")
# colnames(stuff) <- c("95lower","95upper")
# 
# round(stuff*100,2)
# 
# prop.test(x = obs[c(1,2)], n = n[c(1,2)],conf.level=.95)



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
(test.ph <- cox.zph(cox1))
plot(cox.zph(cox1))
ggcoxzph(test.ph)



#############################################
### Cost?
#############################################


cost_A<-as.Date(c("05/06/21","06/15/21"),"%m/%d/%y")
cost_B<-as.Date(c("05/13/21","06/22/21"),"%m/%d/%y")
cost_C<-as.Date(c("05/06/21 ","06/15/21"),"%m/%d/%y")

cost_D1<-as.Date("08/20/21","%m/%d/%y")



start

f1 <- survfit(Surv(time, resp) ~  trt, data = surv)

rates<-summary(f1)
dat_risk<-data.frame(rates$time,rates$strata,rates$n.risk)
names(dat_risk)<-c("time","trt","n_risk")

cost_out<-c(sum(dat_risk$n_risk[dat_risk$time %in% c(cost_A-start) & dat_risk$trt=="trt=A"]),
sum(dat_risk$n_risk[dat_risk$time %in% c(cost_B-start) & dat_risk$trt=="trt=B"]),
sum(dat_risk$n_risk[dat_risk$time %in% c(cost_C-start) & dat_risk$trt=="trt=C"]))



names(cost_out)<-c("trt_A","trt_B","trt_C")
cost_out

cost_start<-dat_risk$n_risk[dat_risk$time %in% c(due_date-start+1)]
cost3<-dat_risk$n_risk[dat_risk$time %in% c(cost_D1-start)]


cost_start
#mailing per nonrespondent
stuff<-round(rbind(cost_start,
      cost_out,
cost_out/cost_start),2)

rownames(stuff)<-c("Nonrespondents","2nd FU Mailings","Mailings/nonrespondents")
stuff

stuff<-round(rbind(cost_start,
             cost_out,
             cost3,
             (cost_out+cost3)/cost_start),2)

rownames(stuff)<-c("Nonrespondents","2nd FU Mailings","3rd FU Mailings","Mailings/nonrespondents")
stuff
# 
#############################################
### Balance and Distance
#############################################

mean_r<-rep(NA,end_time)
mean_nr<-rep(NA,end_time)
mean_samp<-rep(NA,end_time)
mean_rr<-rep(NA,end_time)

surv$wmos<-surv$mos*surv$wgt


#pick which test group

p<-'B'

for(i in 1:(end_time)){
  mean_r[i]<-mean(surv$wmos[surv$time<i & surv$trt==p],na.rm=T)
  mean_nr[i]<-mean(surv$wmos[surv$time>=i & surv$trt==p],na.rm=T)
  mean_samp[i]<-mean(surv$wmos[surv$trt==p],na.rm=T)
  mean_rr[i]<-sum(!is.na(surv$wmos[surv$time<i & surv$trt==p]))/sum(!is.na(surv$wmos[surv$trt==p]))
}


#p<-"C"
#for(i in 1:(end_time)){
#  mean_r[i]<-mean(surv$wmos[surv$time<i],na.rm=T)
#  mean_nr[i]<-mean(surv$wmos[surv$time>=i],na.rm=T)
#  mean_samp[i]<-mean(surv$wmos,na.rm=T)
#  mean_rr[i]<-sum(!is.na(surv$wmos[surv$time<i ]))/sum(!is.na(surv$wmos))
#}

par(mfrow=c(2,1))
par(mgp=c(2,1,0))
par(mar=c(3.1,3.1,1.1,1.1))

plot(1:(end_time),mean_rr,type="l",
     ylab="Response Rate",
     xlab="Time",
     lwd=2,
     main=paste("Treatment",p)
)
abline(v=as.numeric(trt_B-start),lty=order1)

plot(1:(end_time),mean_r,type="l",
     ylim=c(min(mean_r,mean_samp,na.rm=T),max(mean_r,mean_samp,na.rm=T)),
     ylab="Weighted Mean MOS of Respondents",
     xlab="Time",
     lwd=2,
     main=paste("Treatment",p)
)
lines(1:(end_time),mean_samp,lwd=2,lty=2,col="blue")
abline(v=as.numeric(trt_B-start),lty=order1)


#############################################
### Sandbox (extra stuff)
#############################################


aa_fit<-aalen(Surv(time, resp) ~ trt + const(cert)*const(naics3)  , data = surv)
summary(aa_fit)

aa_fit <-aareg(Surv(time, resp) ~ trt + cert + naics3, data = surv)
aa_fit

coef(aa_fit)

p1<-autoplot(aa_fit[2]) 
p2<-autoplot(aa_fit[3])

grid.arrange(p1,p2, nrow = 1, ncol=2)

# test diff of total vs tot_wgt by prop test
A_testx = c(2351,104615)
A_testn= c(11107,624027)
prop.test(A_testx,A_testn, alternative = 'greater')

B_testx = c(2230,92002)
B_testn= c(11084,603490)
prop.test(B_testx,B_testn, alternative = 'greater')


(total_wgt_chkin <- sum(total_fit_wgt$n.event)/total_fit_wgt$n)
sum(surv$wgt)

