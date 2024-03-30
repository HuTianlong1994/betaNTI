#本函数依赖picant、ape、iCAMP三个包
#threads：使用的线程数，可用detectCores()函数查看可使用的CPU数。
#使用的CPU数并不是越多越好，还要顾及自己的内存是否足够，因为每开一个线程都会占据相同的内存。
#otu_niche表示行为样点，列为OTU的表；otu_tree表示OTU的遗传发育树；
#reps表示构建的零模型的个数；
library(picante)
library(ape)
library(iCAMP)
#读取OTU表
otu_biom2 <- read.delim("./data/otu table.txt", row.names=1)
#查看样品序列条数
summary(colSums(otu_biom2))
otu_biom2  = otu_biom2[,colSums(otu_biom2) > 8500]
#将OTU表转置为行名为样品名称，列名为OTU名称的表
otu<-t(otu_biom2)
head(rowSums(otu))
#对OTU表进行抽平
otu = as.data.frame(rrarefy(otu, 8500))
head(rowSums(otu))
#去除序列条数较少的OTU
otu_niche<-otu[,log(colSums(otu)/sum(otu),10)>-4.5]
#otu_niche<-otu[,1:200]
#读取树文件
otu_tree<-read.tree("./data/tree")
#去除发育树中不存在于OTU表中的OTU
prune_tree<-prune.sample(otu_niche,otu_tree) 
#获取发育树中OTU的名称
tip<-prune_tree$tip.label  
coln<-colnames(otu_niche)
m<-NULL
#打印OTU表中不存在于发育树中的OTU
for(i in 1:length(coln)){
  if(!coln[i]%in%tip){
    #print(coln[i])
    m<-cbind(m,coln[i])
  }
}
m<-as.vector(m)
#去除OTU表中不存在于发育树中的OTU
otu_niche<-otu_niche[,!colnames(otu_niche)%in%m]
#计算OTU两两间的遗传距离
otu_phydist <- cophenetic(prune_tree)
beta_nti <- bNTIn.p(otu_niche,otu_phydist,nworker = 10)
rc <- RC.pc(otu_niche,nworker = 10)
#输出表格
write.table(beta_nti$index,"./beta-nti.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)
write.table(rc$index,"./beta-nti.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)
