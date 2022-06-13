
setwd("")

library(rms)
library(foreign)
library(survival)

nomogram<-read.table("lassoRisk2.txt",header=T,sep="\t")
View(nomogram)

cox <- coxph(Surv(OS.time,OS) ~ age + Grade + risk,data=nomogram)
cox$score
cox=summary(cox)
cox
cox$concordance  (c-index为0.7013671)

cox1 <- cph(Surv(OS.time,OS) ~ age + Grade + risk,surv=T,x=T, y=T,time.inc = 1*365*1,data=nomogram) 
cox3 <- cph(Surv(OS.time,OS) ~ age + Grade + risk,surv=T,x=T, y=T,time.inc = 1*365*3,data=nomogram)
cox5 <- cph(Surv(OS.time,OS) ~ age + Grade + risk,surv=T,x=T, y=T,time.inc = 1*365*5,data=nomogram)
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*1, m= 131, B=500)
cal3 <- calibrate(cox3, cmethod="KM", method="boot", u=1*365*3, m= 131, B=500)
cal5 <- calibrate(cox5, cmethod="KM", method="boot", u=1*365*5, m= 131, B=500)

pdf("校准曲线.pdf")
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.7,1),ylim=c(0.5,1),
     xlab="Nomogram-Predicted Probability of 1-Year OS",
     ylab="Actual 1-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))
plot(cal3,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.2,1),ylim=c(0.2,1),
     xlab="Nomogram-Predicted Probability of 3-Year OS",
     ylab="Actual 3-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))
plot(cal5,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 5-Year OS",
     ylab="Actual 5-Year OS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()
