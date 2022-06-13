rm(list = ls())
library(tidyverse)

#读取数据
ALLdata <- data.table::fread("../../TCGA-UCEC.htseq_counts.tsv",data.table = F)
ALLdata[1:5,1:5]
dim(ALLdata)
view(head(ALLdata))

#读取临床信息
clin <- data.table::fread("../../survival_UCEC_survival.txt",data.table = F)
dim(clin)
view(head(clin))

#读取基因注释文件（详见 TCGA数据下载与ID转换）
#基因注释文件不同，基因名可能就不同！！
gtf <- rtracklayer::import('../../Homo_sapiens.GRCh38.105.chr.gtf.gz')
#转换为数据框
gtf <- as.data.frame(gtf)
view(head(gtf,100))
#保存
save(gtf,file = "gtf.Rdata")

#去掉基因名小数点及后面的数字方便下一步转换
#colnames(ALLdata)[1]
colnames(ALLdata)[1]<-'gene_id'
ALLdata1<-separate(ALLdata,gene_id,into = c("gene_id"),sep="\\.") 
view(head(ALLdata1))


#提取lncRNA，与表达谱进行合并以转换基因名。
lncRNA<-dplyr::filter(gtf,type=="gene",gene_biotype=="lncRNA")%>%#选择编码蛋白
  select(gene_name,gene_id,gene_biotype)%>%#选择有用的三列
  inner_join(ALLdata1,by ="gene_id")%>%#与表达谱合并
  select(-gene_name,-gene_biotype)%>%
  distinct(gene_id,.keep_all = T)
dim(lncRNA)
view(head(lncRNA,100))
write.table(lncRNA$gene_id,file = "lncRNA.txt",sep = "\t")

gtf$gene_biotype
#提取lncRNA，与表达谱进行合并以转换基因名。
mRNA<-dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding")%>%#选择编码蛋白
  select(gene_name,gene_id,gene_biotype)%>%#选择有用的三列
  inner_join(ALLdata1,by ="gene_id")%>%#与表达谱合并
  select(-gene_name,-gene_biotype)%>%
  distinct(gene_id,.keep_all = T)
mRNA=t(mRNA)
write.table(mRNA,file = "mRNA.txt",sep = "\t")

#提取m6A调控因子的表达
m6A1=read.table("21_m6A.txt",sep = "\t",header = T)
m6A1=as.data.frame(m6A1)
view(head(m6A1,100))
m6A<-m6A1%>%#选择编码蛋白
  select(gene_id,gene_name)%>%#选择有用的三列
  inner_join(ALLdata1,by ="gene_id")%>%#与表达谱合并
  distinct(gene_id,.keep_all = T)
view(head(m6A,100))
write.table(m6A,file = "m6A-related gene expression.txt",sep = "\t",quote=F,row.names=F)

dim(lncRNA)#14,032个lncRNA,583个样本，35个正常样本，548个UCEC样本

mRNA<-dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding")%>%#选择编码蛋白

mRNA[1:5,1:5]

#把基因名一列作为行名
row.names(mRNA)<-mRNA[,1]
view(head(mRNA))
mRNA<-mRNA[,-1]
view(head(mRNA))
mRNA[1:5,1:5]

#电脑内存不够，先保存数据
save(lncRNA,design,group_list,file = 'UCEC_lncRNA.Rdata')
write.table(lncRNA,file = "lncRNA_expression.txt",sep = "\t",quote=F,row.names=F)
save(mRNA,file = 'UCEC_mRNA.Rdata')
save(clin,file = 'UCEC_clin.Rdata')

rm(list = ls())
load('UCEC_lncRNA.Rdata')

#将表达谱倒置，方便后续与临床资料的合并
lncRNA<-as.data.frame(t(lncRNA))
view(head(lncRNA))

#接下来是设计Normal和Cancer的分组，根据TCGA数据ID的特性可区分Normal和Cancer TCGA差异分析及ggplot作图验证
substr(rownames(lncRNA),14,15)
group_list=ifelse(as.numeric(substr(rownames(lncRNA),14,15)) < 10,'tumor','normal')
table(group_list)#35个正常样本，548个UCEC样本

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(lncRNA)
rownames(design)=rownames(lncRNA)
table(design)
save(group_list,file = 'group.Rdata')

#给design数据,数据变成数据框，添加一列ID
class(design)
design<-as.data.frame(design)
row.names(design)
design$ID<-row.names(design)
design[1:3,1:3]

#将design分组中tumor这一列的1换成Tumor，0换成Normal
design$tumor[design$tumor==1]  <-  "Tumor"#1换成Tumor
design$tumor[design$tumor==0]  <-  "Normal"#0换成Normal
design<-design[,-1]#去掉normal这一列
colnames(design)[1]<-'Type'#cancer这一列列名改为Type
head(design)

#合并数据
lncRNA$ID<-row.names(lncRNA)
lncRNA=inner_join(lncRNA,design,by ="ID",copy=T)
lncRNA[1:5,1:5]
view(head(lncRNA))

#调整列的顺序便于观察
lncRNA<-select(lncRNA,ID,Type,everything())
lncRNA[1:5,1:5]
dim(lncRNA)

#加载临床资料
load('UCEC_clin.Rdata')
clin[1:5,1:5]

#合并
clin1<-rename(clin,ID=sample)
clin1[1:5,1:5]
rownames(clin1)
write.table(clin1$ID,file = "clin.txt",sep = "\t")
write.table(lncRNA$ID,file = "lncRNA.txt",sep = "\t")

drawdata<-dplyr::inner_join(lncRNA,clin1,by ="ID",copy=T)
drawdata[1:5,1:5]
view(head(drawdata))
