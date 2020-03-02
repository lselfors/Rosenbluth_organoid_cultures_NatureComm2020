options(stringsAsFactors = FALSE)
library('heatmap.plus')
library(RColorBrewer)
library(gplots)

###################################
# Reads in csv files from Xshift for 3 separate cytof runs which can be downloaded 
# from https://flowrepository.org/id/FR-FCM-Z2H2:
#	Heatmap_ID1892_HMEC41_WT29.csv
#	20190619_HOP_HMEC24_ORG24_HMEC25.csv
#	Heatmap_ID1922_AllTheRest.csv
# Generates expression signatures for each major epithelial cluster which are defined as the
# median value for each cyTOF marker within each cluster from each primary tissue sample.
# Pearson's R for each cell compared to matched primary tissue signatures are calculated
# outputs max pearson's R for each cell to a file to generate boxplots in Figure 4e and 5b in JMP 13
###################################

directory="~/Downloads"
setwd(directory)
dataFile=c('Heatmap_ID1892_HMEC41_WT29.csv') 

epithelialClusters_manual=c('59261','59254','59271','59257','59264','59265','59266')
dataFromXshift=read.table(dataFile,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift)=gsub("^.+?_","",colnames(dataFromXshift)) 
colnames(dataFromXshift)=gsub("\\..+$","",colnames(dataFromXshift)) 

dataFromXshift=dataFromXshift[,-c(41:ncol(dataFromXshift))]

tissue = c("20170602_WT29.clean")
culture=c("20170602_Hmec41-wt29.clean")
dataFromXshift_tissue=dataFromXshift[dataFromXshift$File %in% tissue,]
dataFromXshift_tissue=dataFromXshift_tissue[dataFromXshift_tissue$ClusterID %in% epithelialClusters_manual,]

dataFromXshift_culture=dataFromXshift[dataFromXshift$File %in% culture,]
medianAbValueByCluster=aggregate(dataFromXshift_tissue[,5:ncol(dataFromXshift_tissue)],by=list(dataFromXshift_tissue$ClusterID),median) 

myData=rbind.data.frame(medianAbValueByCluster[,2:ncol(medianAbValueByCluster)],dataFromXshift_culture[,5:ncol(dataFromXshift_culture)])
rvalues_all=data.frame()
for (cluster in 1:nrow(medianAbValueByCluster)){	
	sample=culture
	rvalue=cbind.data.frame('r'=apply(myData[(nrow(medianAbValueByCluster)+1):nrow(myData),],1, function(x){try(cor.test(as.numeric(x), as.numeric(myData[cluster,]),
	 alternative = c("two.sided"),
	method = c("pearson"),
	exact = FALSE, conf.level = 0.95,na.action='na.omit')$estimate[1],silent = TRUE)}))
	colnames(rvalue)=c(medianAbValueByCluster[cluster,'Group.1'])
	if (nrow(rvalues_all)>0){
		rvalues_all=cbind.data.frame(rvalues_all,rvalue)
	}
	else {rvalues_all=cbind.data.frame("sample"=sample,rvalue)}
}

rvalues_all=cbind.data.frame(rvalues_all,'max'= apply(rvalues_all[,2:ncol(rvalues_all)],1,max))
maxR=rvalues_all[,c('sample','max')]

########################################
dataFile=c('20190619_HOP_HMEC24_ORG24_HMEC25.csv')  
epithelialClusters_manual=c('7024','7023','7028','7018','7031','7022')

dataFromXshift=read.table(dataFile,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift)=gsub("^.+?_","",colnames(dataFromXshift)) 
colnames(dataFromXshift)=gsub("\\..+$","",colnames(dataFromXshift)) 

dataFromXshift=dataFromXshift[,-c(45:ncol(dataFromXshift))]

tissues=c('export_Primary24_concat_1_20190619 unprocessed_195Pt_Viabilty','export_Primary25_concat_1_20190619 unprocessed_195Pt_Viabilty')
cultures=c('export_HMEC24_concat_1_20190619 unprocessed_195Pt_Viabilty','export_HMEC25_concat_1_20190619 unprocessed_195Pt_Viabilty','export_Org24_concat_1_20190619 unprocessed_195Pt_Viabilty','export_Org25_concat_1_20190619 unprocessed_195Pt_Viabilty')  
rvalues_all=data.frame()
for (tissue in tissues){
	dataFromXshift_tissue=dataFromXshift[dataFromXshift$File %in% tissue,]
	dataFromXshift_tissue=dataFromXshift_tissue[dataFromXshift_tissue$ClusterID %in% epithelialClusters_manual,]

	tissueSample=gsub("(export_Primary|_concat.+$)","",tissue)
	print(tissueSample)
	medianAbValueByCluster=aggregate(dataFromXshift_tissue[,7:ncol(dataFromXshift_tissue)],by=list(dataFromXshift_tissue$ClusterID),median) 
	for (culture in cultures){
		if (grepl(tissueSample, culture)){
			print(culture)
			dataFromXshift_culture=dataFromXshift[dataFromXshift$File %in% culture,]
			myData=rbind.data.frame(medianAbValueByCluster[,2:ncol(medianAbValueByCluster)],dataFromXshift_culture[,7:ncol(dataFromXshift_culture)])
			rvalues_all=data.frame()
			for (cluster in 1:nrow(medianAbValueByCluster)){	
				sample=culture
				rvalue=cbind.data.frame('r'=apply(myData[(nrow(medianAbValueByCluster)+1):nrow(myData),],1, function(x){try(cor.test(as.numeric(x), as.numeric(myData[cluster,]),
		 		alternative = c("two.sided"),
				method = c("pearson"),
				exact = FALSE, conf.level = 0.95,na.action='na.omit')$estimate[1],silent = TRUE)}))
				colnames(rvalue)=medianAbValueByCluster[cluster,'Group.1']
				if (nrow(rvalues_all)>0){
					rvalues_all=cbind.data.frame(rvalues_all,rvalue)
				}
				else {rvalues_all=cbind.data.frame("sample"=sample,rvalue)}
			}
			rvalues_all=cbind.data.frame(rvalues_all,'max'= apply(rvalues_all[,2:ncol(rvalues_all)],1,max))
			maxR=rbind.data.frame(maxR,rvalues_all[,c('sample','max')])

		}
	}
}




########################################
dataFile=c('Heatmap_ID1922_AllTheRest.csv')  
epithelialClusters_manual=c('60287','60296','60290','60306','60244','60242','60251','60305','60291','60243','60256','60286','60258','60288','60293','60292','60302','60270','60304','60298','60301','60272')
dataFromXshift=read.table(dataFile,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift)=gsub("^.+?_","",colnames(dataFromXshift)) 
colnames(dataFromXshift)=gsub("\\..+$","",colnames(dataFromXshift))

tissues=c('20170602_Het35-B1.clean','20170602_Het36-B2.clean','20170602_WT14.clean','20170602_WT15.clean')
cultures=c('20170602_O30-wt14.clean','20170602_O32-wt15.clean','20170602_O43-B1het35.clean','20170602_O48-B2het36.clean')  
rvalues_all=data.frame()
for (tissue in tissues){
	dataFromXshift_tissue=dataFromXshift[dataFromXshift$File %in% tissue,]
	dataFromXshift_tissue=dataFromXshift_tissue[dataFromXshift_tissue$ClusterID %in% epithelialClusters_manual,]

	tissueSample=gsub("(^.+_|-.+$|.clean)","",tissue)
	print(tissueSample)
	medianAbValueByCluster=aggregate(dataFromXshift_tissue[,5:ncol(dataFromXshift_tissue)],by=list(dataFromXshift_tissue$ClusterID),median) 
	for (culture in cultures){
		if (grepl(tissueSample, culture,ignore.case=T)){
			dataFromXshift_culture=dataFromXshift[dataFromXshift$File %in% culture,]
			print(culture)
			myData=rbind.data.frame(medianAbValueByCluster[,2:ncol(medianAbValueByCluster)],dataFromXshift_culture[,5:ncol(dataFromXshift_culture)])
			rvalues_all=data.frame()
			for (cluster in 1:nrow(medianAbValueByCluster)){	
				sample=culture
				rvalue=cbind.data.frame('r'=apply(myData[(nrow(medianAbValueByCluster)+1):nrow(myData),],1, function(x){try(cor.test(as.numeric(x), as.numeric(myData[cluster,]),
		 		alternative = c("two.sided"),
				method = c("pearson"),
				exact = FALSE, conf.level = 0.95,na.action='na.omit')$estimate[1],silent = TRUE)}))
				colnames(rvalue)=medianAbValueByCluster[cluster,'Group.1']
				if (nrow(rvalues_all)>0){
					rvalues_all=cbind.data.frame(rvalues_all,rvalue)
				}
				else {rvalues_all=cbind.data.frame("sample"=sample,rvalue)}
			}
			rvalues_all=cbind.data.frame(rvalues_all,'max'= apply(rvalues_all[,2:ncol(rvalues_all)],1,max))
			maxR=rbind.data.frame(maxR,rvalues_all[,c('sample','max')])

		}
	}
}
maxR$sample=gsub("^(20170602_|export_)","",maxR$sample)
maxR$sample=gsub("(-|_).+$","",maxR$sample)
maxR$sample=gsub("O","Org",maxR$sample)
maxR$sample=gsub("Orgrg","Org",maxR$sample)
maxR$sample=gsub("Hmec","HMEC",maxR$sample)

### 5b
pdf("figure_5b.pdf")
stripchart(as.numeric(max)~factor(sample),data=maxR[maxR$sample %in% c("Org24","HMEC24"),],
            vertical = TRUE, method = "jitter",pch=19,col='gray',group.names=NA,ylab="Max Correlation (Pearson's R)",main="",cex=.5,cex.axis=0.8)  
boxplot(as.numeric(max)~factor(sample),data=maxR[maxR$sample %in% c("Org24","HMEC24"),],las=3,ylab="",xlab="",outpch = NA,lwd = 2,whisklty=1, outcex=0.1,boxwex=0.4,add=T,cex.axis=0.8)
dev.off()
### 4e
pdf("figure_4e.pdf")
stripchart(as.numeric(max)~factor(sample),data=maxR[maxR$sample %in% c("Org30","Org32","Org43","Org48","HMEC41"),],
            vertical = TRUE, method = "jitter",pch=19,col='gray',group.names=NA,ylab="Max Correlation (Pearson's R)",main="",cex=.5,cex.axis=0.8)  
boxplot(as.numeric(max)~factor(sample),data=maxR[maxR$sample %in% c("Org30","Org32","Org43","Org48","HMEC41"),],las=3,ylab="",xlab="",outpch = NA,lwd = 2,whisklty=1, outcex=0.1,boxwex=0.4,add=T,cex.axis=0.8)
dev.off()


#pvalues for 4e
wilcox.test(as.numeric(maxR[maxR$sample %in% c("HMEC41"),'max']),as.numeric(maxR[maxR$sample %in% c('Org30'),'max']), alternative = "two.sided")
wilcox.test(as.numeric(maxR[maxR$sample %in% c("HMEC41"),'max']),as.numeric(maxR[maxR$sample %in% c('Org32'),'max']), alternative = "two.sided")
wilcox.test(as.numeric(maxR[maxR$sample %in% c("HMEC41"),'max']),as.numeric(maxR[maxR$sample %in% c('Org43'),'max']), alternative = "two.sided")
wilcox.test(as.numeric(maxR[maxR$sample %in% c("HMEC41"),'max']),as.numeric(maxR[maxR$sample %in% c('Org48'),'max']), alternative = "two.sided")

# pvalue for 5b
wilcox.test(as.numeric(maxR[maxR$sample %in% c("HMEC24"),'max']),as.numeric(maxR[maxR$sample %in% c('Org24'),'max']), alternative = "two.sided")

write.table(maxR,file='dataForFigure4e_5b.txt',sep='\t',row.names=F)
