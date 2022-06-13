
library(rms)
library(foreign)
library(survival)
#setwd("")
nomogram<-read.table("12-m6A-LPS-ROC.txt",header=T,sep="\t")

ddist <- datadist(nomogram)
options(datadist='ddist')
View(head(nomogram))

cox <- cph(Surv(OS.time,OS) ~age + Grade + riskScore ,surv=T,x=T, y=T,data=nomogram) 

surv <- Survival(cox)
sur_1_year<-function(x)surv(365*1,lp=x)
sur_3_year<-function(x)surv(365*3,lp=x)
sur_5_year<-function(x)surv(365*5,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-year survival','3-year survival','5-year survival'),maxscale=100,fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))
nom_sur$`3-year survival`
pdf("nomogram plot.pdf")
plot(nom_sur,xfrac=0.25)
dev.off()

