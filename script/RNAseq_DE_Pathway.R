#!/usr/bin/Rscript
setwd("/home/chenlei_8th/CYZ_RNAseq_analysis/ZouQiang_compare_two_sequencer") #you need to modify this line according to where you put the folder
library('DESeq2')
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

######read in files Notice human/mouse ensembl/gencode
#y=read.table("gene_id_type_name_human",header=T,row.names=1) #human
#y=read.table("gene_id_type_name_ensembl",header=T,row.names=1) #mouse
y=read.table("gene_id_type_name",header=T,row.names=1) #mouse

meta=read.table("meta",header=T,row.names=1)
#meta$Sample=rownames(meta)
#meta$Condition_Time=paste(meta$Condition,meta$Time,sep="_")
filelist=list.files(pattern="count")
#Treg=as.character(read.table("Treg")$V1)
#Th=as.character(read.table("Th")$V1)

######
x=read.table(filelist[[1]])
ids=x$V1
x=x$V2
for(i in filelist[2:length(filelist)]){
  print(i)
  tmp=read.table(i)
  x=cbind(x,tmp$V2)
}
###mod this to suit your file suffix !and meta!
colnames(x)=gsub(".count.union.reverse","",filelist) 
meta$File<-gsub('.count.union.reverse','',filelist)
x=x[,meta$File]
colnames(x)=rownames(meta)
rownames(x)=ids
x=x[grep("ENS",rownames(x)),] #valid gene only
y=y[rownames(x),]
y=y[rowSums(x)>0,]
x=x[rowSums(x)>0,]#at least one sample has it

rownames(x)=gsub('\\.\\d+$',"",rownames(x))
rownames(y)=gsub('\\.\\d+$',"",rownames(y))

####basic stats
bigmeta=meta
bigmeta$Sum=colSums(x)
bigmeta$Gene1=colSums(x>0)
bigmeta$Gene10=colSums(x>9)

#number of genes in samples
barplot(bigmeta$Gene10,names=rownames(bigmeta))
# ggplot(bigmeta,aes(x=ILC,y=Gene1))+
#   geom_boxplot()+
#   facet_wrap(~Run)
# 
# sg1=rowSums(x[,1:3])
# sg2=rowSums(x[,4:6])
# g1=sum(sg2>0)
# g10=sum(sg2>9)
# 
# hist(log2(sg2+1),br=50,xlim=c(0,20))
# 
# x11()
# boxplot(bigmeta$Gene1,ylim=c(0,25000))
# plot(x[,1],x[,7])
# cor(x[,1],x[,7],method="spearman")
# 
x11()
par(mfcol=c(3,5))
for(i in c(1:14)){
  hist(log2(x[,i]+1),br=50,xlim=c(0,20))
#  plot(x[,i],x[,i+3],main=cor(x[,i],x[,i+3],method="spearman"))
}
# dev.set(4)
# a=9
# b=10
# plot(x[,a],x[,b],main=cor(x[,a],x[,b],method="spearman"))
# 
# 
# hist(log2(x[,3]+1),br=50,xlim=c(0,20))
#####################


#DESeq2 works on raw counts, tell DESeq2 which sample belong to which group and which groups to compare
#g_list=as.factor(c("WT","WT","WT","KO","KO","KO"))
#g_list=as.factor(c("WT","KO"))
#mycol=data.frame(row.names=colnames(myx),group=g_list)
#meta$myCondition=c("KO","WT","U","U","WT","KO")

####change the fields and which groups to compare here
tag="Condition2"
group1="KO"
group2="WT"
##########################
#x=x[,c(-5,-11,-12)]
#meta=meta[c(-5,-11,-12),]

dds=DESeqDataSetFromMatrix(countData=x,colData = meta,design=as.formula(paste("~",tag,sep="")))
colData(dds)[,tag]=as.factor(colData(dds)[,tag]) ###
dds2=DESeq(dds,fitType="parametric") # fitType parametric is default

rld=rlog(dds2)
x11()
#####simple ordination plot, need improvement
#tmpx=plotPCA(rld,intgroup=c(tag),ntop=10000,returnData=T)
xxx=plotPCA(rld,intgroup=c("Condition1"),ntop=2000)
xxx
###label dots with texts
x11()
plot(xxx$data$PC1, xxx$data$PC2,col="red")
text(xxx$data$PC1, xxx$data$PC2,labels=xxx$data$name,col="blue")

#mycolor=meta$Condition
#levels(mycolor)=c("red","green")

#################




#######diversity stuff from vegan
div=diversity(t(x))


########

res_o=results(dds2,contrast=c(tag,group1,group2),alpha=0.05)
res=as.data.frame(res_o)   #res is the result from DESeq2

#add normalized count to res so you can see the actual expression levels alongside p values etc.
norm_x=counts(dds2,normalized=T)
########grab gene length and calculate fpkm
fpkm=as.data.frame(norm_x/colSums(norm_x)*1000000)
#rownames(fpkm)=gsub('\\.\\d+$',"",rownames(fpkm))
glens=getlength(rownames(fpkm),"mm9","ensGene")
glens[is.na(glens)]=-1
glens_name=data.frame(len=glens,name=y$gene_name)
fpkm=fpkm/glens*1000


# LHB=as.data.frame(read.csv("LHB_list.csv",header=T,row.names=1))
# myfpkm=fpkm
# myfpkm$ens=rownames(fpkm)
# myfpkm$gene_name=y$gene_name
# myfpkm=myfpkm[myfpkm$gene_name %in% rownames(LHB),]
# rownames(myfpkm)=myfpkm$gene_name
# myfpkm=myfpkm[rownames(LHB),]
# myfpkm$gene_name=NULL
# rownames(myfpkm)=rownames(LHB)
# WT=rowMeans(myfpkm[,meta$Condition_Time=="WT_T3"])
# KO=rowMeans(myfpkm[,meta$Condition_Time=="KO_T3"])
# myout=data.frame(WT=WT,KO=KO)
# write.csv(myout, file="myfpkm.csv")
# 
# LHB=as.data.frame(read.csv("LHB.csv",header=T,row.names=1))
# MY=as.data.frame(read.csv("myfpkm.csv",header=T,row.names=1))
# LHB=log2(LHB+1)
# MY=log2(MY+1)
# LHB=LHB-rowMeans(LHB)
# MY=MY-rowMeans(MY)
# COMBINE=cbind(LHB,MY)
# COL=rep("black",nrow(COMBINE))
# COL[grep("Socs",rownames(COMBINE))]="red"
# COL[grep("Cish",rownames(COMBINE))]="red"
# COL[grep("Asb",rownames(COMBINE))]="red"
# COMBINE=COMBINE[COL=="red",]
# x11()
# pdf("socs.pdf",height=6,width=5)
# heatmap.2(as.matrix(COMBINE),
#           col=bluered(100),
#           Rowv=F,Colv=F,dendrogram="none",
#           key=T,cexRow=1,cexCol=1,
#           trace="none",#sepwidth=c(0,0),sepcolor="black",
#           key.par=list(mgp=c(1.5, 0.5, 0),
#                        mar=c(4.5, 0.5, 4.5, 0)),
#           density.info="none",
# #          colRow=COL,
#           colsep=2,sepcolor="white",sepwidth=c(0.1,0.1)
#           )
# 
# dev.off()  

res=cbind(res,y,norm_x) #add gene names and raw counts
res=cbind(res,fpkm)  #add fpkm
res=res[order(res$padj),]  # order according to adjusted p value
write.csv(res,file=paste(group1,"vs",group2,".csv",sep=""))


####filter res to your desire!
#myres=res[res$gene_type=="lincRNA",]
res1000=res[1:1000,]
myres=res[abs(res$log2FoldChange)>1,]
myres=myres[myres$padj<0.01,]
myres=myres[order(myres$log2FoldChange, decreasing = T),]
myres=myres[!is.na(myres$padj),]

#####volcano plot
plot(res$log2FoldChange,-log10(res$padj),xlim=c(-10,10),ylim=c(0,100))
up=(res$padj<0.01)&(res$log2FoldChange>2)
down=(res$padj<0.01)&(res$log2FoldChange<(-2))
points(res$log2FoldChange[up],-log10(res$padj)[up],col="red")
points(res$log2FoldChange[down],-log10(res$padj)[down],col="blue")
###############







#some fancy plot come with DESeq2
#plotMA(res_o,main="DESeq2",ylim=c(-2,2))  ##out of frame dot plot as triangles
#plotDispEsts(dds2)
#plotCounts(dds2,gene="ENSMUSG00000027833.16",intgroup="group")
#vsd=getVarianceStabilizedData(dds2)


#rownames(res_o)=y$gene_name
#plot(res_o$log2FoldChange,-log(res_o$padj),pch=15)
#sig=res_o[which(res_o$padj  < 0.01),]
#points(sig$log2FoldChange,-log(sig$padj),col="red",pch=15)
#identify(res_o$log2FoldChange,-log(res_o$padj))




######add db column stuff
#columns(org.Hs.eg.db)
#supportedOrganisms()
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]
go_up=as.numeric((res$log2FoldChange>1)&(res$padj<0.05))
#names(go_up)=gsub('\\.\\d+$',"",rownames(res))
go_up=go_up[!is.na(go_up)]
table(go_up)

pwf_up=nullp(go_up,"hg19","ensGene",plot.fit=F)
plotPWF(pwf_up,binsize=200)
GO.wall_up=goseq(pwf_up,"hg19","ensGene")
#GO.sample_up=goseq(pwf_up,"hg19","ensGene")

#######KEGG
en2eg=as.list(org.Hs.egENSEMBL2EG)
eg2kegg=as.list(org.Hs.egPATH)
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
kegg=lapply(en2eg,grepKEGG,eg2kegg)

my_up=goseq(pwf_up,gene2cat=kegg)

KEGG_up=goseq(pwf_up,"hg19","ensGene",test.cats="KEGG")
KEGG_up$PathName=as.list(KEGGPATHID2NAME)[KEGG_up$category]


go_dn=as.numeric((res$log2FoldChange<-1)&(res$padj<0.05))
names(go_dn)=gsub('\\.\\d+$',"",rownames(res))
table(go_dn)
#KEGGPATHNAME2ID
# sig_id=rownames(res)[1:50]
# sig_id=gsub('\\.\\d+$',"",sig_id)
# tmp=mapIds(org.Hs.eg.db, keys=sig_id, keytype="ENSEMBL",multiVals="first",
#        column="ONTOLOGYALL")

#########

####plot and stuff
#x11()
#plot(rowMeans(log2(counts(dds2,normalized=T)[,1:3] +1)), rowMeans(log2(counts(dds2,normalized=T)[,4:6] +1)), pch=".",col="gray")
#textxy(rowMeans(log2(counts(dds2,normalized=T)[1:3,1:3] +1)), rowMeans(log2(counts(dds2,normalized=T)[1:3,4:6] +1)), c("p1","p2","p3"),cex=1,offset=2,m=c(10,10))
#points(rowMeans(log2(counts(dds2,normalized=T)[1:3,1:3] +1)), rowMeans(log2(counts(dds2,normalized=T)[1:3,4:6] +1)), col="red")


#norm_counts=counts(dds2,normalized=T)
#rld=rlog(dds2)
norm_counts=assay(rld)
norm_counts1000=norm_counts[rownames(res1000),]
norm_countsmy=norm_counts[rownames(myres),]
#norm_counts=norm_counts[y$gene_type,]
rownames(norm_counts)=y$gene_name
#norm_counts=counts(dds2,normalized=T)
rownames(norm_counts)=y[rownames(norm_counts),]$gene_name
rownames(norm_counts1000)=y[rownames(norm_counts1000),]$gene_name
rownames(norm_countsmy)=y[rownames(norm_countsmy),]$gene_name

gmeanK=rowMeans(norm_counts[,which(meta[,tag]==group1)])
gmeanW=rowMeans(norm_counts[,which(meta[,tag]==group2)])
#gmean1=log2(gmean1+1)
#gmean2=log2(gmean2+1)
norm_counts_o=norm_counts  ###for correlation plot use non-centered 
norm_counts=norm_counts-rowMeans(norm_counts)
norm_counts1000=norm_counts1000-rowMeans(norm_counts1000)
norm_countsmy=norm_countsmy-rowMeans(norm_countsmy)

#####simple correlation matrix plot
mycor=cor(norm_counts_o,method="spearman")
corrplot(mycor,addCoef.col = "grey")
#corrplot(mycor,add=T,type="lower",method="number",diag=F)


##norm_counts can be used to do overall clustering and ordination
#my_counts=norm_counts[Th,]
#row_ann=data.frame(cell_signature=rep("Treg",18))
#row_ann=data.frame(cell_signature=c(rep("Th1",6),rep("Th2",7),rep("Th17",5),rep("Tfh",3),rep("Th9",1)))
#rownames(row_ann)=rownames(my_counts)
pheat_meta=meta
pheat_meta$File=NULL
pheat_meta$Sample=NULL
pheat_meta$Condition_Time=NULL
pheat_meta$Sex=NULL
#pheat_meta$Run=NULL
tmpx=norm_counts
tmpx=norm_counts[c("Socs1","Socs2","Socs3","Cish"),]
rownames(tmpx)=rep("",dim(tmpx)[[1]])
rownames(tmpx)[rownames(norm_counts)=="lnc17"]="lnc17"

####pheatmap
x11()
pheatmap(tmpx,cluster_row=F,legend=T,
         show_rownames = T,labels_row = F,
#         annotation_row=row_ann,
         breaks=c(-4,-3,-2,-1,-0.7,-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.4,0.7,1,2,3,4),
         color=redblue(18),
         annotation_col=pheat_meta)
#############

####pheatmap for genes from file
pdf("Th_heatmap1.pdf",height=6,width=4)
pheatmap(norm_counts[Th,],cluster_row=F,legend=T,
         show_rownames = T,
         annotation_col=pheat_meta,
         border_color=NA,
         colorRampPalette(c("red", "black", "green"))(50))
dev.off()
#####heatmap.2 some may like it better
x11()
pdf("Treg_heatmap2_forcolorblind.pdf",height=6,width=5)
heatmap.2(norm_counts[Treg,],
         rowsep=c(0),sepwidth=c(1,1),sepcolor="white",
          col=bluered(30),
          margins=c(10,12),
          Rowv=F,Colv=F,dendrogram="none",
          key=T,cexRow=1,cexCol=1,
          trace="none",#sepwidth=c(0,0),sepcolor="black",
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(4.5, 0.5, 4.5, 0)),
          density.info="none")
dev.off()
#######heatmap.2 of top1000 etc
x11()
pdf("cutoff.pdf",height=6,width=5)
heatmap.2(norm_countsmy,
          col=bluered(30),
          # margins=c(10,12),
          Rowv=F,Colv=F,
          key=T,cexRow=1,cexCol=1,
          labRow="",
          trace="none",
          density.info = "none",
          key.par=list(mgp=c(1.5, 0.5, 0),
                  mar=c(4.5, 0.5, 4.5, 0))
          )
dev.off()  
  
#####pheatmap for top down
res$padj[is.na(res$padj)]=1
res$sig=res$padj<0.05
sig_res=res[res$sig,]
up=as.character(sig_res[which(sig_res$log2FoldChange>0),]$gene_name[1:15])
down=as.character(sig_res[which(sig_res$log2FoldChange<0),]$gene_name[1:15])
pheatmap(norm_counts[up,],cluster_row=F,legend=T,
         show_rownames = T,
         annotation_col=pheat_meta)
pheatmap(norm_counts[down,],cluster_row=F,legend=T,
         show_rownames = T,
         annotation_col=pheat_meta)
##################

####scatterplot with labels
tmp=rep("",length(norm_counts))
names(tmp)=rownames(norm_counts)
#tmp[up]=up
#tmp[down]=down
tmp[Treg]=Treg
tmp[Th]=Th
ggdata=data.frame(WT=gmeanW,KO=gmeanK,label=tmp)
tmp=rep("other",length(norm_counts))
names(tmp)=rownames(norm_counts)
#tmp[up]="up"
#tmp[down]="down"
tmp[Treg]="Treg"
tmp[Th]="Th"
ggdata$type=tmp
#ggdata_up=ggdata[which(ggdata$type=="up"),]
#ggdata_down=ggdata[which(ggdata$type=="down"),]
ggdata_Treg=ggdata[which(ggdata$type=="Treg"),]
ggdata_Th=ggdata[which(ggdata$type=="Th"),]
pdf("scatterplot3.pdf",height=6,width=6)
ggplot(ggdata,aes(WT,KO,label=label))+
  geom_point(size=1,color="gray")+
  geom_point( data=ggdata_Treg,size = 1, color = "blue")+
  geom_point(data=ggdata_Th,size=1,color="red")+
#  geom_label_repel(aes(label=label),
#                   point.padding=NA,
#                   box.padding = unit(0.5, 'lines'),
#                   force=15)+
  # geom_label_repel(data=ggdata_Treg,
  #                  fill="red",
  #                  force=10,
  #                  nudge_x=3,
  #                  nudge_y=-5)+
  # geom_label_repel(data=ggdata_Th,
  #                  fill="blue",
  #                  force=10,
  #                  nudge_x=-3,
  #                  nudge_y=5)+
  theme_classic()
dev.off()
#norm_counts=norm_counts[rownames(res[1:100,]),]
#rownames(norm_counts)=res$gene_name[1:100]
#norm_counts=log2((norm_counts+1)/(rowMeans(norm_counts)+1))
#up=norm_counts[rowSums(norm_counts[,1:3])>0,]
#down=norm_counts[rowSums(norm_counts[,1:3])<0,]

#pheatmap(up)
#pheatmap(down)


#x <- rnorm(50)
#y <- rnorm(50)
#plot(x,y,asp=1)
#textxy(x,y,1:50,m=c(mean(x),mean(y)))
rownames(norm_x)=y$gene_name
norm_my=norm_x[c("Socs1","Socs2","Socs3","Cish"),]
norm_my=t(norm_my)
norm_my=as.data.frame(norm_my)
norm_my$group=meta$Condition_Time
xxx=ddply(norm_my,~group,summarise,s1m=mean(Socs1),
          s2m=mean(Socs2),
          s3m=mean(Socs3),
          cm=mean(Cish))

T1=xxx[1,2:5]/xxx[4,2:5]
T2=xxx[2,2:5]/xxx[5,2:5]
T3=xxx[3,2:5]/xxx[6,2:5]
xx=rbind(T1,T2,T3)
xx=log2(xx)
rownames(xx)=c("T1","T2","T3")
x11()
matplot(xx,type="l",col=1:4,pch=1,xlab="Time point",ylab="log2(KO/WT)")
legend("topleft",legend=c("Socs1","Socs2","Socs3","Cish"),pch=1,col=1:4)
