library(velocyto.R)
ldat <- read.loom.matrices("merge.loom")
emat <- ldat$spliced
emat <- emat[,colSums(is.na(emat))<nrow(emat)]
# this dataset has already been pre-filtered, but this is where one woudl do some filtering
emat <- emat[,colSums(emat)>=1e3]
library(pagoda2)
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Xist"],main='Xist') 

#Velocity estimation
#emat <- emat[,colSums(is.na(emat))<nrow(emat)]
emat <- ldat$spliced; nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- sccore::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE


cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

vel <- rvel.cd; arrow.scale=3; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
show.velocity.on.embedding.cor(emb,vel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)




