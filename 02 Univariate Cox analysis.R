library(survival)
rt=read.table("inputfile.txt",header=T,sep="\t",check.names=F,row.names=1)
View(head(rt[,3:ncol(rt)]))
View(colnames(rt[,3:ncol(rt)]))
colnames(rt)[2]<- "status"

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox_results.xls",sep="\t",row.names=F,quote=F)
#p<0.005，最终有187个m6A-related lncRNA与UCEC患者的预后相关。

#绘制森林图
#加载这两个R包
#BiocManager::install("forestplot")
library(forestplot)
library(haven)

#加载肺癌这套数据
setwd("../../../")
ForestPlot=read.table("12-m6A-LPS-forestplot.txt",header=T,sep="\t",check.names=F)
attach(ForestPlot)
ForestPlot[,1:3]
      forestplot(as.matrix(ForestPlot[,1:3]), HR, HR.95L, HR.95H, graph.pos=4, zero=1,xlab = "Hazard ratio", graphwidth=unit(50,"mm"), lineheight="auto", boxsize=0.2, xticks=(c(0.7,1.0,1.5,2.0)), col= fpColors(box=c("red","red","red","red","green","red","red","red","red","red","red","red"),summary="#8B008B",lines = 'black',zero = 'black'))
