
library(survival)

setwd("")
rt=read.table("Multivariate Cox.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(OS, Censor) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="MultiCox.xls",sep="\t",row.names=F,quote=F)
