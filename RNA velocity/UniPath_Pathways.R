library(UniPath)
##Load all data files 
data("human_null_model")
data("c2.cp.v6.1.symbols")
data<-read.table(file.choose(new = T), sep="\t", header=T)
data2<-make.unique(data$X)
rownames(data)<-c(as.data.frame(data2)$data)
data<-subset(data,select = -c(X))
View(data)

##Converting mouse null data into p-values
Pval = binorm(human_null_data)
##Converting gene expression data into p-values
Pval1 = binorm(data)

##Combining of p-values for null model data matrix
combp_ref = combine(c2.cp.v6.1.symbols,human_null_data,rownames(human_null_data),Pval,thr=2)
##Combining of p-values for gene expression data matrix
combp = combine(c2.cp.v6.1.symbols,data,rownames(data),Pval1,thr=2)
scores = adjust(combp,combp_ref)
View(scores)



my_data<-data = read.csv(file.choose(new = F))
my_data <- df[order(apply(df, 1, max), decreasing = TRUE)[1:100],]
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(my_data, scale="row", col = coul))









