#!/usr/bin/Rscript
setwd("/home/samples") #you need to modify this line according to where you put the folder
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# meta table should record sample condition (gender, treatment, control)
meta=read.table("meta",header=T,row.names=1)
# import count files from prervious RNAseq pipeline
filelist=list.files(pattern="count")

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

par(mfcol=c(3,5))
for(i in c(1:14)){
  hist(log2(x[,i]+1),br=50,xlim=c(0,20))
#  plot(x[,i],x[,i+3],main=cor(x[,i],x[,i+3],method="spearman"))
}

#####################
Differential expression
####change the fields and which groups to compare here
tag="Condition2"
group1="KO"
group2="WT"
##########################
# make DEseq object, add group info 
dds=DESeqDataSetFromMatrix(countData=x,colData = meta,design=as.formula(paste("~",tag,sep="")))
colData(dds)[,tag]=as.factor(colData(dds)[,tag]) ###
dds2=DESeq(dds,fitType="parametric") # fitType parametric is default

rld=rlog(dds2)
#####simple ordination plot, need improvement
#tmpx=plotPCA(rld,intgroup=c(tag),ntop=10000,returnData=T)
xxx=plotPCA(rld,intgroup=c("Condition1"),ntop=2000)
xxx
###label dots with texts
plot(xxx$data$PC1, xxx$data$PC2,col="red")
text(xxx$data$PC1, xxx$data$PC2,labels=xxx$data$name,col="blue")

#mycolor=meta$Condition
#levels(mycolor)=c("red","green")

########

res_o=results(dds2,contrast=c(tag,group1,group2),alpha=0.05)
res=as.data.frame(res_o)   #res is the result from DESeq2

#add normalized count to res so you can see the actual expression levels alongside p values etc.
norm_x=counts(dds2,normalized=T)
# make file results table
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




