
setwd("")      
  rt=read.table("12-m6A-LPS-exp.txt",sep="\t",header=T, row.names = 1, check.names=F) 
  View(rt)
  rt=t(rt)
outpdf="heatmap_exp_clin.pdf"

library(pheatmap)
library(RColorBrewer)
library(colorspace)

pdf(outpdf)
#pheatmap(rt[6:17,], annotation_col = rt[1:3,], cellwidth = 0.5, cellheight = 20, cluster_rows = T,
#         color = colorRampPalette(c("blue", "white", "orange"))(100),
#         border_color = "white",
#         annotation_height=13,
#         annotation_colors = list(
#           riskScore= c("white", "black")), 
##           Grade = c(G1 = "lightcyan2", G2 = "yellow2", G3 ="orange",High_Grade="red"),
# #          risk=c(low="blue",high="red")),
#         cluster_cols =F,
#         fontsize=10,
#         fontsize_row=10,
#         scale="row",
#         show_colnames=F,
#         show_rownames = T,
#         fontsize_col=3)
#dev.off()
#pheatmap(rt,show_colnames = F)

Type=read.table("clin.txt",sep="\t",header=T,check.names=F)
# 构建列注释信息
annotation_col = data.frame(Type$Grade,Type$age,Type$risk,Type$riskScore)
rownames(annotation_col) = Type$id
head(annotation_col)
colnames(annotation_col)[1] = 'Grade'
colnames(annotation_col)[2] = 'age'
colnames(annotation_col)[4] = 'riskScore'
colnames(annotation_col)[3] = 'risk'

pheatmap(rt,show_colnames = F,annotation_col = annotation_col,cluster_cols = F,
        annotation_colors = list(
        riskScore= c("white", "black"), 
        risk=c(low="green",high="red"),
        age=c("lightcyan2","yellow2","orange"),
        Grade = c(G1 = "lightcyan2", G2 = "yellow2", G3 ="orange",High_Grade="red")
        ),
         )
dev.off()

