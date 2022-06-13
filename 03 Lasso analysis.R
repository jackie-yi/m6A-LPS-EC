library("glmnet")
library("survival")

getwd()               
rt=read.table("Lassoinput.txt",header=T,sep="\t",row.names=1)      

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
fit
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
#两条虚线分别指示了两个特殊的λ值,一个是lambda.min,一个是lambda.1se,这两个值之间的lambda都认为是合适的。lambda.1se构建的模型最简单，即使用的基因数量少，而lambda.min则准确率更高一点，使用的基因数量更多一点。

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)
#最终获得了12个m6A-LPS

riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
median(riskScore)
outCol=c("OS.time","OS",lassoGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
View(risk)
table(risk) #高风险273，低风险273
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore),risk)
write.table(cbind(id=rownames(outTab),outTab),
            file="lassoRisk.txt",
            sep="\t",
            quote=F,
            row.names=F)

