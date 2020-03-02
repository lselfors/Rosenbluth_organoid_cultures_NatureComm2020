options(StringsAsFactors=FALSE)

######################################### 
# Reads in four CyTOF files from Xshift:
#	csv_ID2096_heatmap_ORGANOIDS_Run1.csv
# 	Org23 Org25 for Laura_ORGANOIDS_Run2.csv
# 	Heatmap_ID1922_AllTheRest.csv
# 	20190619_HOP_HMEC24_ORG24.csv 
# These files are available at https://flowrepository.org/id/FR-FCM-Z2H2 
# Generates expression signatures for each major epithelial cluster which are defined as the
# median value for each cyTOF marker within each cluster across unmatched primary 
# tissue samples.
# Correlates data from 10 unmatched organoid epithelial samples to these signatures.
# Calculates Max Pearson's R
#########################################


intersection=c("ClusterID","EventID","File","CD44","CD54","AndrogenReceptor","K14","EGFR","EPCR","Muc1","EpCAM","anpep","CD90","hsp27","CD49f","GlucocorticoidR","CD24","SMA","anxa8","laminin5","CD10","CD95","CD73","Galectin1","ERa","CD133","Her2","CD47","K8","ProgesteronR","CK17","HLA","BRCA1")          
directory=("~/Downloads")
setwd(directory) 

dataFile=c('csv_ID2096_heatmap_ORGANOIDS_Run1.csv') 
dataFile2=c('Org23 Org25 for Laura_ORGANOIDS_Run2.csv')  

omitCluster=c("62937")
dataFromXshift_org1=read.table(dataFile,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift_org1)=gsub("^.+?_","",colnames(dataFromXshift_org1)) # cleans up the column names.
colnames(dataFromXshift_org1)=gsub('cd45',"CD45",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('cd31',"CD31",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('K8_18',"K8",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('cd140b',"CD140b",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('ProgesteronR_B',"ProgesteronR",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('K17',"CK17",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('HLA_abc',"HLA",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('Brca1',"BRCA1",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('Rank',"RANK",colnames(dataFromXshift_org1))
colnames(dataFromXshift_org1)=gsub('File.Name',"File",colnames(dataFromXshift_org1))


dataFromXshift_org2=read.table(dataFile2,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift_org2)=gsub("\\..+$","",colnames(dataFromXshift_org2)) 
colnames(dataFromXshift_org2)=gsub("^.+?_","",colnames(dataFromXshift_org2)) 
colnames(dataFromXshift_org2)=gsub('K17',"CK17",colnames(dataFromXshift_org2))
colnames(dataFromXshift_org2)=gsub('Brca1',"BRCA1",colnames(dataFromXshift_org2))

dataFromXshift_org=rbind.data.frame(dataFromXshift_org1[,intersection], dataFromXshift_org2[,intersection])
dataFromXshift_org=dataFromXshift_org[!(dataFromXshift_org$ClusterID %in% omitCluster),]

##### tissues
dataFile=c('Heatmap_ID1922_AllTheRest.csv')  
dataFile2=c('20190619_HOP_HMEC24_ORG24.csv')  

epithelialClusters_manual1=c('60287','60296','60290','60306','60244','60242','60251','60305','60291','60243','60256','60286','60258','60288','60293','60292','60302','60270','60304','60298','60301','60272')
epithelialClusters_manua12=c('7024','7023','7028','7018','7031','7022')

dataFromXshift1=read.table(dataFile,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift1)=gsub("^.+?_","",colnames(dataFromXshift1)) 
colnames(dataFromXshift1)=gsub("\\..+$","",colnames(dataFromXshift1)) 
dataFromXshift2=read.table(dataFile2,sep=",",row.names=NULL,header=T) 
colnames(dataFromXshift2)=gsub("^.+?_","",colnames(dataFromXshift2)) 
colnames(dataFromXshift2)=gsub("\\..+$","",colnames(dataFromXshift2)) 

colnames(dataFromXshift1)=gsub('cd45',"CD45",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('cd31',"CD31",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('K8_18',"K8",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('cd140b',"CD140b",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('ProgesteronR_B',"ProgesteronR",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('K17',"CK17",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('HLA_abc',"HLA",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('Brca1',"BRCA1",colnames(dataFromXshift1))
colnames(dataFromXshift1)=gsub('Rank',"RANK",colnames(dataFromXshift1))

dataFromXshift=rbind.data.frame(dataFromXshift1[,intersection], dataFromXshift2[,intersection])
tissues1=c('20170602_Het35-B1.clean','20170602_Het36-B2.clean','20170602_WT14.clean','20170602_WT15.clean')
tissues2=c('export_Primary24_concat_1_20190619 unprocessed_195Pt_Viabilty','export_Primary25_concat_1_20190619 unprocessed_195Pt_Viabilty')


dataFromXshift_tissue=dataFromXshift[dataFromXshift$File %in% c(tissues1, tissues2),]
dataFromXshift_tissue$File=as.character(dataFromXshift_tissue$File)
dataFromXshift_tissue=dataFromXshift_tissue[dataFromXshift_tissue$ClusterID %in% c(epithelialClusters_manual1,epithelialClusters_manua12),]
medianAbValueByCluster=aggregate(dataFromXshift_tissue[,4:ncol(dataFromXshift_tissue)],by=list(dataFromXshift_tissue$ClusterID),median)

rvalues_all=data.frame()

dataFromXshift_org$File=as.character(dataFromXshift_org$File)

cultures=c("20170602_O42-p53.clean","20170602_O46-ATM.clean" ,"20170602_O11-B1.clean" ,"20170602_O07-wt.clean", "20170602_O45-B1.clean" ,"20170602_O44-B1p53.clean","20170602_O38-B1.clean","20170602_O37-palb1.clean","20160826_Organoids_BC1_O23_Clean_Debarcoded","20160826_Organoids_BC3_O25_Clean_Debarcoded")
maxR=data.frame()
for (culture in cultures){
	dataFromXshift_culture=dataFromXshift_org[dataFromXshift_org$File %in% culture,]
	print(culture)
	myData=rbind.data.frame(medianAbValueByCluster[,2:ncol(medianAbValueByCluster)],dataFromXshift_culture[,4:ncol(dataFromXshift_culture)])
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
maxR$sample2=maxR$sample
maxR$sample=gsub("20170602_","",maxR$sample )
maxR$sample=gsub("-.+.clean","",maxR$sample )
maxR$sample=gsub("_Clean_Debarcoded","",maxR$sample )
maxR$sample=gsub("20160826_Organoids_BC\\d_","",maxR$sample )
maxR$sample=gsub("O","Org",maxR$sample )

pdf("figure_5c.pdf")
stripchart(as.numeric(max)~factor(sample),data=maxR,
            vertical = TRUE, method = "jitter",pch=19,col='gray',group.names=NA,ylab="Max Correlation (Pearson's R)",main="",cex=.5,cex.axis=0.8)  
boxplot(as.numeric(max)~factor(sample),data=maxR,las=3,ylab="",xlab="",outpch = NA,lwd = 2,whisklty=1, outcex=0.1,boxwex=0.4,add=T,cex.axis=0.8)
dev.off()

png("figure_5c.png")
stripchart(as.numeric(max)~factor(sample),data=maxR,
            vertical = TRUE, method = "jitter",pch=19,col='gray',group.names=NA,ylab="Max Correlation (Pearson's R)",main="",cex=.5,cex.axis=0.8)  
boxplot(as.numeric(max)~factor(sample),data=maxR,las=3,ylab="",xlab="",outpch = NA,lwd = 2,whisklty=1, outcex=0.1,boxwex=0.4,add=T,cex.axis=0.8)
dev.off()

write.table(maxR,file='dataForFigure5c.txt',sep='\t',row.names=F)

