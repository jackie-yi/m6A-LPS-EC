library(survivalROC)
library(pROC)
library(ROCR)
#BiocManager::install("AUC")
library(AUC)
library(timeROC)
library(plotROC)
library(survival)
library(RColorBrewer)
rocCol=brewer.pal(8,"Set1")
aucText=c()

mayo=read.table("12-m6A-LPS-ROC.txt",sep = "\t",header=T,check.names=F)
roc1=survivalROC(Stime=mayo$OS.time, ##生存时间
                status=mayo$OS, ## 终止事件
                marker = mayo$riskScore, ## risk score
                predict.time =356, ## 预测时间截点
                method="KM")
roc3=survivalROC(Stime=mayo$OS.time, ##生存时间
                 status=mayo$OS, ## 终止事件
                 marker = mayo$riskScore, ## risk score
                 predict.time =1068, ## 预测时间截点
                 method="KM")
roc5=survivalROC(Stime=mayo$OS.time, ##生存时间
                 status=mayo$OS, ## 终止事件
                 marker = mayo$riskScore, ## risk score
                 predict.time =1780, ## 预测时间截点
                 method="KM")
#span = 0.25*nobs^(-0.20))##span,NNE法的namda)
roc1$AUC
roc3$AUC
roc5$AUC
#0.6710215, 0.7274743, 0.7754893

plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue", 
     xlab="1-Specificity", ylab="Sensitivity",
     lwd = 2, cex.main=1, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("1-year survival: ",sprintf("%.3f",roc1$AUC)))
aucText=c(aucText,paste0("3-year survival: ",sprintf("%.3f",roc3$AUC)))
abline(0,1)
aucText=c(aucText,paste0("5-year survival: ",sprintf("%.3f",roc5$AUC)))
lines(roc3$FP, roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="yellow",lwd=3)
lines(roc5$FP, roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red",lwd=3)
legend("bottomright", aucText,lwd=3,bty="n",col=c("blue","yellow","red"),cex =1.2)




library(survival)
library(timeROC)
setwd("")
data<-read.table("12-m6A-LPS-ROC-nomogram.txt",header=T,sep="\t")
predict_1_year<- 1*365
predict_3_year<- 3*365
predict_5_year<- 5*365

pdf("ROC1_plot.pdf")
#1年OS
ROC1<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$riskScore,cause=1,
              weighting="marginal",
              times=c(predict_1_year),ROC=TRUE)
ROC2<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$age,cause=1,
              weighting="marginal",
              times=c(predict_1_year),ROC=TRUE)
ROC3<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$Grade,cause=1,
              weighting="marginal",
              times=c(predict_1_year),ROC=TRUE)
ROC4<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$nomogram,cause=1,
              weighting="marginal",
              times=c(predict_1_year),ROC=TRUE)
plot(ROC1,time=predict_1_year,col="blue",title=FALSE,lwd=3)
plot(ROC2,time=predict_1_year,col="yellow",add=TRUE,title=FALSE,lwd=3)
plot(ROC3,time=predict_1_year,col="red",add=TRUE,title=FALSE,lwd=3)
plot(ROC4,time=predict_1_year,col="green",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("Risk score: ",round(ROC1$AUC[2],3)),
         paste("Age: ",round(ROC2$AUC[2],3)),
         paste("Grade: ",round(ROC3$AUC[2],3)),
         paste("Nomogram: ",round(ROC4$AUC[2],3))),col=c("blue","yellow","red","green"),lwd=3)
dev.off()

#3年OS
pdf("ROC3_plot.pdf")
ROC1<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$riskScore,cause=1,
              weighting="marginal",
              times=c(predict_3_year),ROC=TRUE)
ROC2<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$age,cause=1,
              weighting="marginal",
              times=c(predict_3_year),ROC=TRUE)
ROC3<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$Grade,cause=1,
              weighting="marginal",
              times=c(predict_3_year),ROC=TRUE)
ROC4<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$nomogram,cause=1,
              weighting="marginal",
              times=c(predict_3_year),ROC=TRUE)
plot(ROC1,time=predict_3_year,col="blue",title=FALSE,lwd=3)
plot(ROC2,time=predict_3_year,col="yellow",add=TRUE,title=FALSE,lwd=3)
plot(ROC3,time=predict_3_year,col="red",add=TRUE,title=FALSE,lwd=3)
plot(ROC4,time=predict_3_year,col="green",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("Risk score: ",round(ROC1$AUC[2],3)),
         paste("Age: ",round(ROC2$AUC[2],3)),
         paste("Grade: ",round(ROC3$AUC[2],3)),
         paste("Nomogram: ",round(ROC4$AUC[2],3))),col=c("blue","yellow","red","green"),lwd=3)
dev.off()

#5年OS
pdf("ROC5_plot.pdf")
ROC1<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$riskScore,cause=1,
              weighting="marginal",
              times=c(predict_5_year),ROC=TRUE)
ROC2<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$age,cause=1,
              weighting="marginal",
              times=c(predict_5_year),ROC=TRUE)
ROC3<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$Grade,cause=1,
              weighting="marginal",
              times=c(predict_5_year),ROC=TRUE)
ROC4<-timeROC(T=data$OS.time,delta=data$OS,
              marker=data$nomogram,cause=1,
              weighting="marginal",
              times=c(predict_5_year),ROC=TRUE)

plot(ROC1,time=predict_5_year,col="blue",title=FALSE,lwd=3)+
plot(ROC2,time=predict_5_year,col="yellow",add=TRUE,title=FALSE,lwd=3)
plot(ROC3,time=predict_5_year,col="red",add=TRUE,title=FALSE,lwd=3)
plot(ROC4,time=predict_5_year,col="green",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("Risk score: ",round(ROC1$AUC[2],3)),
         paste("Age: ",round(ROC2$AUC[2],3)),
         paste("Grade: ",round(ROC3$AUC[2],3)),
       paste("Nomogram: ",round(ROC4$AUC[2],3))),col=c("blue","yellow","red","green"),lwd=3)
dev.off()


