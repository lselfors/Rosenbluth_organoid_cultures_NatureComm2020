options(stringsAsFactors = FALSE)
library('heatmap.plus')
library(RColorBrewer)
library(gplots)


###################################
# Reads in csv files from Xshift for cytof runs which can be downloaded 
# from https://flowrepository.org/id/FR-FCM-Z2H2:
#	Heatmap_ID1892_HMEC41_WT29.csv
# Generates heatmaps in Figure 5a
###################################


directory=("~/Downloads")
setwd(directory) 

dataFile=c('20190619_HOP_HMEC24_ORG24_HMEC25.csv') 

dataFromXshift=read.table(dataFile,sep=",",header=TRUE) 

colnames(dataFromXshift)=gsub("^.+?_","",colnames(dataFromXshift)) 
colnames(dataFromXshift)=gsub("c.Parp","cParp",colnames(dataFromXshift)) 
colnames(dataFromXshift)=gsub("\\..+?$","",colnames(dataFromXshift)) 

dataFromXshift=dataFromXshift[,-c(45:54,56:ncol(dataFromXshift))]
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

colorsForHM2Df=as.data.frame(medianAbValueByClusterCentered[,1:3])
rownames(colorsForHM2Df)=colorsForHM2Df$cluster

orderForHMs=cbind.data.frame(colorsForHM2Df[clusterDim2$label[clusterDim2$order],],'order'=seq(from = 1, to =numberOfClusters , by = 1))
samples=unique(dataFromXshift$File)
dataFromXshiftCentered=scale(dataFromXshift[,5:ncol(dataFromXshift)],center=TRUE,scale=TRUE)

dataFromXshiftCentered[dataFromXshiftCentered > 4] <- 4
dataFromXshiftCentered[dataFromXshiftCentered < -4 ] <- -4
dataFromXshiftCentered=cbind.data.frame('File'=dataFromXshift$File,'ClusterID'=dataFromXshift$ClusterID,dataFromXshiftCentered)

for (sample in c("export_Org24_concat_1_20190619 unprocessed_195Pt_Viabilty", "export_Primary24_concat_1_20190619 unprocessed_195Pt_Viabilty","export_HMEC24_concat_1_20190619 unprocessed_195Pt_Viabilty")){

	type=gsub("(export_|24.+$)","",sample)
	testData=dataFromXshiftCentered[dataFromXshiftCentered$File %in% sample,]
	testData=merge(orderForHMs,testData,by.x='cluster',by.y='ClusterID')
	testData=testData[order(testData$order),]
	colorsForHM=as.matrix(cbind.data.frame(cluster=testData$cluster1,cluster2=testData$cluster1))
	colorsForLegend=testData[!duplicated(testData$cluster),]
	colorsForLegend=colorsForLegend[order(colorsForLegend$order,decreasing=TRUE),]

	pdf(paste("figure_5_",type,".pdf",sep=""))
	heatmap.plus(as.matrix(testData[,AbOrderForHM]),scale="none",Rowv=NA,RowSideColors=colorsForHM, Colv=NA,labRow=NA,cexRow=0.5,cexCol=1,margins = c(12,12),col=my_palette, breaks=col_breaks)
	legend('topright',bty = "n",inset=c(-0.1,-0.1),legend = colorsForLegend$cluster, text.width = max(sapply(text, strwidth)),col=colorsForLegend$cluster1, lwd=5, cex=0.4,xpd=TRUE)
	dev.off()
	
	png(paste("figure_5_",type,".png",sep=""))
	heatmap.plus(as.matrix(testData[,AbOrderForHM]),scale="none",Rowv=NA,RowSideColors=colorsForHM, Colv=NA,labRow=NA,cexRow=0.5,cexCol=1,margins = c(12,12),col=my_palette, breaks=col_breaks)
	legend('topright',bty = "n",inset=c(-0.1,-0.1),legend = colorsForLegend$cluster, text.width = max(sapply(text, strwidth)),col=colorsForLegend$cluster1, lwd=5, cex=0.4,xpd=TRUE)
	dev.off()
}
