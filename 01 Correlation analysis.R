
library(limma)

corFilter=0.5          
pvalueFilter=0.001       

rt = read.table("lncRNA_expression.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,2]
view(head(rt))
exp=rt[,3:ncol(rt)]
view(head(exp))
dimnames=list(rownames(exp),colnames(exp))
lncRNA=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
view(head(lncRNA))
lncRNA=avereps(lncRNA)
lncRNA=lncRNA[rowMeans(lncRNA)>0.5,]
dim(lncRNA)
sapply(strsplit(colnames(lncRNA),"\\-"),"[",4)
group=sapply(strsplit(colnames(lncRNA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
group
lncRNA=lncRNA[,group==0]
dim(lncRNA)

rt = read.table("m6A-related gene expression.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
view(head(rt))
rownames(rt)=rt[,2]
exp=rt[,3:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
m6AGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
m6AGene=avereps(m6AGene)
m6AGene=m6AGene[rowMeans(m6AGene)>0.5,]
group=sapply(strsplit(colnames(m6AGene),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
m6AGene=m6AGene[,group==0]
view(head(m6AGene,100))
dim(m6AGene)

outTab=data.frame()
for(i in row.names(lncRNA)){
	  if(sd(lncRNA[i,])>0.5){
			  for(j in row.names(m6AGene)){
				     x=as.numeric(lncRNA[i,])
				     y=as.numeric(m6AGene[j,])
				     corT=cor.test(x,y)
				     cor=corT$estimate
				     pvalue=corT$p.value
				     if((cor>corFilter) & (pvalue<pvalueFilter)){
				         outTab=rbind(outTab,cbind(m6AGene=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				     }
				     if((cor< -corFilter) & (pvalue<pvalueFilter)){
				         outTab=rbind(outTab,cbind(m6AGene=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				     }
			  }
		}
}
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F) 
view(head(outTab))
dim(outTab)#一共有13,181对m6A-lncRNA关系对,|相关系数|>0.5，p<0.001
m6ALncRNA=unique(outTab[,"lncRNA"])
m6ALncRNAexp=lncRNA[m6ALncRNA,]
m6ALncRNAexp=rbind(ID=colnames(m6ALncRNAexp),m6ALncRNAexp)
dim(m6ALncRNAexp)#一共有1,560个与m6A相关的lncRNA
write.table(m6ALncRNAexp,file="m6ALncRNAexp.txt",sep="\t",quote=F,col.names=F)


