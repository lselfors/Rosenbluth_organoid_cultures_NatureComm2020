
options(stringsAsFactors = FALSE)
library('heatmap.plus')
library(RColorBrewer)
library(gplots)

directory=("~/Downloads")
setwd(directory) 

dataFile=c('csv_ID2096_heatmap_ORGANOIDS_Run1.csv') 

dataFromXshift=read.table(dataFile,sep=",",header=TRUE) 

colnames(dataFromXshift)=gsub("^.+?_","",colnames(dataFromXshift)) # cleans up the column names.

medianAbValueByCluster=aggregate(dataFromXshift[,7:ncol(dataFromXshift)],by=list(dataFromXshift$ClusterID),median) 

rownames(medianAbValueByCluster)=medianAbValueByCluster[,1]

numberOfClusters=nrow(medianAbValueByCluster)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colorsForCluster=sample(col_vector, numberOfClusters, TRUE) 
colorMatrix=cbind.data.frame('cluster'=unique(dataFromXshift$ClusterID),'color'=colorsForCluster)

dataWithColors=merge(colorMatrix,medianAbValueByCluster,by.x='cluster',by.y=1)
rownames(dataWithColors)=dataWithColors$cluster

dissimilarity <- 1 - cor(dataWithColors[,3:ncol(dataWithColors)],method='pearson') 
distance <- as.dist(dissimilarity)
cluster=hclust(distance, method = "average") 

dissimilarityDim2 <- 1 - cor(t(dataWithColors[,3:ncol(dataWithColors)]),method='pearson')
distanceDim2 <- as.dist(dissimilarityDim2)
clusterDim2=hclust(distanceDim2, method = "average")
AbOrderForHM=cluster$label[cluster$order]

my_palette <- colorRampPalette(c("#0000FF", "white", "#FF0000"))(n = 59)
col_breaks = c(seq(-4,-1,length=20),  # for blue
              seq(-0.99,0.99,length=20),           # for white
              seq(1,4,length=20))      #red

medianAbValueByClusterCentered=scale(dataWithColors[,3:ncol(dataWithColors)],center=TRUE,scale=TRUE)
medianAbValueByClusterCentered[medianAbValueByClusterCentered > 4] <- 4
medianAbValueByClusterCentered[medianAbValueByClusterCentered < -4 ] <- -4

medianAbValueByClusterCentered=cbind.data.frame('cluster'=dataWithColors$cluster,'cluster1'=dataWithColors$color,'cluster2'=dataWithColors$color,medianAbValueByClusterCentered)

pdf("figure_S4.pdf", width = 6, height = 8.5, family = "Helvetica", paper = "special") 
heatmap.plus(as.matrix(medianAbValueByClusterCentered[,4:ncol(medianAbValueByClusterCentered)]),RowSideColors=as.matrix(medianAbValueByClusterCentered[,2:3]),scale="none",Rowv=as.dendrogram(clusterDim2),Colv=as.dendrogram(cluster),cexRow=0.5,cexCol=1,margins = c(10,10),col=my_palette, breaks=col_breaks)
legend('topright',bty = "n",inset=c(-0.2,-0.1),legend = medianAbValueByClusterCentered$cluster, text.width = max(sapply(text, strwidth)),col=as.matrix(medianAbValueByClusterCentered[,2]), lwd=5, cex=0.4,xpd=TRUE, ncol = 1) 
dev.off()