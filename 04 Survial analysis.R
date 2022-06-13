library(survival)
library(survminer)

rt=read.table("lassoRisk.txt",header=T,sep="\t")
View(head(rt))

diff=survdiff(Surv(OS.time, OS) ~ risk,data = rt)

#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #2.766e-10
fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)
summary(fit)   
pdf(file="survivalTrain.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

#风险因子联动图和热图
#BiocManager::install("ggrisk")
library(ggrisk)
library(rms)
exp=read.table("lassoRisk1.txt",header=T,sep="\t",check.names=F)
exp[,1:5]
View(exp)

#2、构建COX模型
exp1=exp[,1:14]
fit <- cph(Surv(OS.time,OS)~rt$riskScore,exp1)

fit <- cph(Surv(OS.time,OS)~RP11.296O14.3+AC009501.4+RP11.657O9.1+LINC00501+RP11.6E9.4+NADK2.AS1+RP11.420L9.5+AC144652.1+RP4.778K6.3+POT1.AS1+CTD.2007L18.5+AC010761.9,exp1)
fit <- cph(Surv(time,status)~ANLN+CENPA+GPR182+BCO2,LIRI)
#fit
ggrisk(fit,
       cutoff.value='median', #可选‘median’, ’roc’ or ’cutoff’
       cutoff.x = 150,  #“cutoff”文本的水平位置
       cutoff.y = -1,  #“cutoff”文本的垂直位置
       )
#cutoff.value 为划分风险值cutoff的方式：
#cutoff.value = "median"：默认方式，使用风险得分的中位数作为切点值；
#cutoff.value = "roc"：将风险得分和结局时间进行roc分析，将约登点作为最佳切点值；
#cutoff.value = "cutoff"：将会使用cutoff包，通过最小p值法计算出最佳切点。
#cutoff.value = 赋值数值：根据切点值将风险得分分为高危组和低危组。

#调整细节以及颜色
?ggrisk
ggrisk(fit,
       cutoff.value='median', #可选‘median’, ’roc’ or ’cutoff’
       cutoff.x = 150,  #“cutoff”文本的水平位置
       cutoff.y = -1,  #“cutoff”文本的垂直位置
       code.highrisk = 'high risk',#高风险标签，默认为 ’High’
       code.lowrisk = 'low risk', #低风险标签，默认为 ’Low’
       title.A.ylab='Risk score', #A图 y轴名称
       title.B.ylab='Survival time(years)', #B图 y轴名称，注意区分year month day
       title.A.legend='Risk group', #A图图例名称
       title.B.legend='Status',     #B图图例名称
       title.C.legend='Expression', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度    
       color.A=c(low='green',high='red'),#A图中点的颜色
       color.B=c(code.0='green',code.1='red'), #B图中点的颜色
       color.C=c(low='green',median='white',high='red'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2  #B图中y轴标签到y坐标轴的距离,默认是2
)            
#图A为风险得分按照从小到大的顺序排列 （此示例为根据中值分组）；
#图B为风险得分与生存时间的散点图，并按照结局将散点图分成红色和蓝色；
#图C为基因表达量热图；

library(pheatmap)
library(ggplot2)
library(ggplotify)
library(cowplot)
exp1$OS=factor(exp1$OS,labels = c("Alive","Dead"))
p2=ggplot(data=exp1)+geom_point(aes(x=seq(0:545),y=OS.time,color=OS))+scale_x_continuous(breaks = seq(0,545,50))+scale_y_continuous(breaks = seq(0,12,3))
p3=p2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+labs(x="",y="Survival time(years)")+theme(axis.line = element_line(size = 1,colour = "black"))+geom_vline(aes(xintercept=273),colour="#BB0000",linetype="dashed")
p3

biomarker_risk=exp1[order(as.numeric(rt$riskScore)),]
biomarker_risk$group=ifelse(rt$riskScore>1.5,"high risk","low rsk")
table(rt$risk)
p4=ggplot(data=exp1)+geom_point(aes(x=seq(0:545),y=rt$riskScore,color=rt$risk))+scale_x_continuous(breaks = seq(0,545,50))+scale_y_continuous(breaks = seq(-2,5,1))+geom_vline(aes(xintercept=273),colour="#BB0000",linetype="dashed")
p5=p4+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+labs(x="",y="Risk score")+theme(axis.line = element_line(size = 1,colour = "black"))
p5
exp2=exp1[,3:14]
pheatmap(exp2,show_colnames = F,scale = "row",fontsize_row = 11,color = colorRampPalette(c("blue","white","red"))(100),legend = F,show_rownames = T,cluster_rows = F)

#确定不同临床特征指标亚组是否可以用于预测OS，年龄分>=60和<60，Grade分G1、G2、G3和High_Grade
rt=read.table("lassoRisk2.txt",header=T,sep="\t")
View(head(rt))

pdf("clin-risk-OS.pdf")
#Grade亚组G1
G1=rt[rt[,6]=="G1",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = G1)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #3.778e-02
fit <- survfit(Surv(OS.time, OS) ~ risk, data = G1)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("Grade G1"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))

#Grade亚组G2
G2=rt[rt[,6]=="G2",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = G2)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #3.215e-01
fit <- survfit(Surv(OS.time, OS) ~ risk, data = G2)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("Grade G2"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))

#Grade亚组G3
G3=rt[rt[,6]=="G3",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = G3)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #8.99e-04
fit <- survfit(Surv(OS.time, OS) ~ risk, data = G3)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("Grade G3"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))

#Grade亚组G4
G4=rt[rt[,6]=="High_Grade",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = G4)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #2.549e-01
fit <- survfit(Surv(OS.time, OS) ~ risk, data = G4)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("High Grade"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))

#age亚组>=60
age1=rt[rt[,7]==">=60",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = age1)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #1.638e-06
fit <- survfit(Surv(OS.time, OS) ~ risk, data = age1)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("age >= 60"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))

#age亚组<60
age2=rt[rt[,7]=="<60",]
diff=survdiff(Surv(OS.time, OS) ~ risk,data = age2)
#计算p value
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
pValue  #9.049e-04
fit <- survfit(Surv(OS.time, OS) ~ risk, data = age2)
summary(fit)   
plot(fit, lwd=2,col=c("red","blue"),xlab="Time (year)",ylab="Survival rate",  main=paste("age < 60"),mark.time=T)
legend("topright", c("high risk", "low risk"),lwd=2,col=c("red","blue"))
dev.off()
