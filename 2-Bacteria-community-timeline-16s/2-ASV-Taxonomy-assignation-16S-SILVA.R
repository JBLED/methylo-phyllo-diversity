################################################################################################################################################################################
################
####ALTERNATIVE: TAXONOMIC ASSIGNATION USING SILVA

library(dada2)


#Work directory

#Work directories
pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-1"

pD2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-2"

####ASV tables

#Triming parameters
TFx=200
TRx=200
TLx=20

#sequences randomly sampled in each sample
NS=5000

Prare=0.003 #only keep ASV that have at least Prare% of relative abundance in at least one sample

############

setwd(pD1)

#Load corrected/rarefied ASV table
seqtab.nochim<-read.table(paste("RAREFIED_SeqTableNoChim-PL1_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

#ASV SEQUENCES
ASVs<-tolower(row.names(seqtab.nochim))

#Sample names
sample.names<-colnames(seqtab.nochim)


taxa <- assignTaxonomy(ASVs, paste(pD2,"/silva_nr_v138_train_set.fa.gz",sep=""), multithread=F)

setwd(pD2)

write.table(taxa,paste("SILVA-taxonomy_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,".txt",sep=""))

taxa<-as.matrix(read.table(paste("SILVA-taxonomy_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,".txt",sep=""),header=T))


rownames(taxa)<-NULL

taxa<-ifelse(is.na(taxa)==T,"unknown",taxa)


#Phyllum in Bacteria
PHYb<-as.vector(sort(unique(subset(taxa[,2],taxa[,1]=="Bacteria"))))


BAC.PHY.silva<-t(sapply(c(1:length(PHYb)),function(x){
	LX= paste("Bacteria", PHYb[x],sep="_")
	taX<-as.vector(paste(taxa[,1],taxa[,2],sep="_"))
	#Phyllo samples
	Fx<-subset(seqtab.nochim[,c(2:length(seqtab.nochim[1,]))], taX==LX)
	taxF<-subset(taX, taX==LX)
	taxF<-subset(taxF,as.vector(apply(Fx,1,max))!=0)
	Fx<-subset(Fx,as.vector(apply(Fx,1,max))!=0)
	#Meth community
	METHx<-as.vector(subset(seqtab.nochim[,1], taX==LX))
	taxM <-subset(taX, taX==LX)
	taxM<-subset(taxM, METHx!=0)
	METHx <-subset(METHx,METHx!=0)
	return(c(LX,length(Fx[,1]),sum(as.vector(apply(Fx,1,mean)))/NS,length(taxM),sum(METHx)/NS))
}))

colnames(BAC.PHY.silva)<-c("Taxa","ASV-Phyllosphere","Abundance-Phyllosphere","ASV-METHcomm.","Abundance-METHcomm.")


#Classes in Proteobacteria
CLAp<-as.vector(sort(unique(subset(taxa[,3],taxa[,2]=="Proteobacteria"))))


PROT.CLA.silva<-t(sapply(c(1:length(CLAp)),function(x){
	LX= paste("Proteobacteria", CLAp[x],sep="_")
	taX<-as.vector(paste(taxa[,2],taxa[,3],sep="_"))
	#Phyllo samples
	Fx<-subset(seqtab.nochim[,c(2:length(seqtab.nochim[1,]))], taX==LX)
	taxF<-subset(taX, taX==LX)
	taxF<-subset(taxF,as.vector(apply(Fx,1,max))!=0)
	Fx<-subset(Fx,as.vector(apply(Fx,1,max))!=0)
	#Meth community
	METHx<-as.vector(subset(seqtab.nochim[,1], taX==LX))
	taxM <-subset(taX, taX==LX)
	taxM<-subset(taxM, METHx!=0)
	METHx <-subset(METHx,METHx!=0)
	return(c(LX,length(Fx[,1]),sum(as.vector(apply(Fx,1,mean)))/NS,length(taxM),sum(METHx)/NS))
}))

colnames(PROT.CLA.silva)<-c("Taxa","ASV-Phyllosphere","Abundance-Phyllosphere","ASV-METHcomm.","Abundance-METHcomm.")


#Orders in AlphaProteobacteria
ORDa<-as.vector(sort(unique(subset(taxa[,4],taxa[,3]=="Alphaproteobacteria"))))


ALPH.ORD.silva<-t(sapply(c(1:length(ORDa)),function(x){
	LX= paste("Alphaproteobacteria", ORDa[x],sep="_")
	taX<-as.vector(paste(taxa[,3],taxa[,4],sep="_"))
	#Phyllo samples
	Fx<-subset(seqtab.nochim[,c(2:length(seqtab.nochim[1,]))], taX==LX)
	taxF<-subset(taX, taX==LX)
	taxF<-subset(taxF,as.vector(apply(Fx,1,max))!=0)
	Fx<-subset(Fx,as.vector(apply(Fx,1,max))!=0)
	#Meth community
	METHx<-as.vector(subset(seqtab.nochim[,1], taX==LX))
	taxM <-subset(taX, taX==LX)
	taxM<-subset(taxM, METHx!=0)
	METHx <-subset(METHx,METHx!=0)
	return(c(LX,length(Fx[,1]),sum(as.vector(apply(Fx,1,mean)))/NS,length(taxM),sum(METHx)/NS))
}))

colnames(ALPH.ORD.silva)<-c("Taxa","ASV-Phyllosphere","Abundance-Phyllosphere","ASV-METHcomm.","Abundance-METHcomm.")

#Famillies in Rhizobiales
FAMr<-as.vector(sort(unique(subset(taxa[,5],taxa[,4]=="Rhizobiales"))))

RHIZO.FAM.silva<-t(sapply(c(1:length(FAMr)),function(x){
	LX= paste("Rhizobiales", FAMr[x],sep="_")
	taX<-as.vector(paste(taxa[,4],taxa[,5],sep="_"))
	#Phyllo samples
	Fx<-subset(seqtab.nochim[,c(2:length(seqtab.nochim[1,]))], taX==LX)
	taxF<-subset(taX, taX==LX)
	taxF<-subset(taxF,as.vector(apply(Fx,1,max))!=0)
	Fx<-subset(Fx,as.vector(apply(Fx,1,max))!=0)
	#Meth community
	METHx<-as.vector(subset(seqtab.nochim[,1], taX==LX))
	taxM <-subset(taX, taX==LX)
	taxM<-subset(taxM, METHx!=0)
	METHx <-subset(METHx,METHx!=0)
	return(c(LX,length(Fx[,1]),sum(as.vector(apply(Fx,1,mean)))/NS,length(taxM),sum(METHx)/NS))
}))

colnames(RHIZO.FAM.silva)<-c("Taxa","ASV-Phyllosphere","Abundance-Phyllosphere","ASV-METHcomm.","Abundance-METHcomm.")

#Genera in Beijerinckiaceae
GENb<-as.vector(sort(unique(subset(taxa[,6],taxa[,5]=="Beijerinckiaceae"))))


BEIJ.GEN.silva<-t(sapply(c(1:length(GENb)),function(x){
	LX= paste("Beijerinckiaceae",GENb[x],sep="_")
	taX<-as.vector(paste(taxa[,5],taxa[,6],sep="_"))
	#Phyllo samples
	Fx<-subset(seqtab.nochim[,c(2:length(seqtab.nochim[1,]))], taX==LX)
	taxF<-subset(taX, taX==LX)
	taxF<-subset(taxF,as.vector(apply(Fx,1,max))!=0)
	Fx<-subset(Fx,as.vector(apply(Fx,1,max))!=0)
	#Meth community
	METHx<-as.vector(subset(seqtab.nochim[,1], taX==LX))
	taxM <-subset(taX, taX==LX)
	taxM<-subset(taxM, METHx!=0)
	METHx <-subset(METHx,METHx!=0)
	return(c(LX,length(Fx[,1]),sum(as.vector(apply(Fx,1,mean)))/NS,length(taxM),sum(METHx)/NS))
}))

colnames(BEIJ.GEN.silva)<-c("Taxa","ASV-Phyllosphere","Abundance-Phyllosphere","ASV-METHcomm.","Abundance-METHcomm.")

write.table(rbind(BAC.PHY.silva ,PROT.CLA.silva ,ALPH.ORD.silva ,RHIZO.FAM.silva ,BEIJ.GEN.silva),paste("Summary_tax-SILVA_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_Prare=", Prare,".txt",sep=""),row.names=F)

#Calculate ASV RELATIVE ABUNDANCE per sample

seqtab.nochimREL<-sapply(c(1:length(seqtab.nochim[1,])),function(x){
	return(seqtab.nochim[,x]/sum(seqtab.nochim[,x]))
	
})

rownames(seqtab.nochimREL)<-ASVs
colnames(seqtab.nochimREL)<-sample.names

write.table(seqtab.nochimREL,paste("RELATIVE_RAREFIED_SeqTableNoChim_PL1_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))



###################
###Summary of Taxonomy with METAdata

setwd(pD1)

####Load Metadata for each sample
METAsamp<-read.table("METADATA-PL1-wo-NEGcontrols.txt")

#Indexes for positive controls
vMETH<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]=="METH")

#Indexes for samples
vSAMP<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]!="METH")

#Taxa abundance per sample
taX<-as.vector(paste(taxa[,2],taxa[,3],taxa[,4],sep=";"))
CLAb<-sort(unique(taX))


FR<-t(sapply(c(1:length(CLAb)),function(x){
	LX= CLAb[x]
	#Phyllo samples
	Fx<-subset(seqtab.nochim, taX==LX)
	return(apply(Fx,2,sum)/NS)
}))

#Remove very rare taxa
CLAb<-subset(CLAb,apply(FR,1,max)>0.05)
FR<-subset(FR,apply(FR,1,max)> 0.05)

rownames(FR)<-CLAb

setwd(pD2)
COL<-read.table("COL.txt")

COLb<-round(seq(1,length(COL[,1]),length.out=length(CLAb)),0)
COLb <-COL[COLb,]
COLb <-rgb(COLb[,1], COLb[,2], COLb[,3])
COLb <-ifelse(CLAb =="unknown;unknown","grey", COLb)

#Description of main taxa relative abundance per factor:

pdf(paste("Summary_16s-ORDER-per-Factor_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_Prare=", Prare,"-wo-neg-controls-SILVA.pdf",sep=""))


par(mfrow=c(2,3),bty="n",mar=c(4,0,1,0))

###1 Taxonomic legend
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",main="Taxonomy",cex.main=1)

legend(x = "center",cex=0.8, horiz=F,legend= CLAb,fill= COLb,border=T,bg=rgb(1,1,1,0.9),box.col=NA)

#legend("center",legend=rownames(FRall),cex=0.6, fill= COLtax,bty="n",border=NA)

###2 METH communities and all samples

FRall<-cbind(FR[,vMETH] ,apply(FR[,vSAMP],1,mean))
 colnames(FRall)<-c("METH","PHYLLO")

par(mar=c(5,5,1,1))
barplot(FRall,las=1,ylab="Part of diversity",ylim=c(0,1.1),cex.axis=0.8,col= COLb,border=NA,main="All samples",cex.main=1)

###3 per SITE of sampling (excluding METH communities)

SITE<-as.vector(METAsamp[vSAMP,3])
SITEu<-sort(unique(SITE))

FRsite<-sapply(SITEu,function(x){
	return(apply(subset(t(FR[,vSAMP]), SITE==x),2,mean))
	
})

barplot(FRsite,las=1,ylab="Part of diversity",ylim=c(0,1.1),cex.axis=0.8,col= COLb,border=NA,main="Per site",cex.main=1)

###4 per SUBSITE of sampling (excluding METH communities)

SS<-paste(SITE,
	as.vector(METAsamp[vSAMP,4]),sep="-")
SSu<-sort(unique(SS))
#Add an empty class
SSu<-c(SSu,"")
#Manage subsite order (along geographic gradient)
SSo<-c(1,2,4,5,3)
SSu<-SSu[order(SSo)]

FRss<-sapply(SSu,function(x){
	return(apply(subset(t(FR[,vSAMP]), SS ==x),2,mean))
})

barplot(FRss,las=2,ylab="Part of diversity",ylim=c(0,1.1),cex.axis=0.8,col= COLb,border=NA,main="Per sub-site",cex.main=1,cex.names=0.8)

text((c(1.5,4.5)-0.5)*1.25,c(1.05,1.05),SITEu,cex=1)



###5 per HOST TREE SPECIES per site of sampling (excluding METH communities)

SPE<-as.vector(METAsamp[vSAMP,6])

#Replace tree species names by abbreviations
SPE <-ifelse(SPE=="Acer_saccharum","ACSA",
	ifelse(SPE=="Fagus_grandifolia","FAGR",
	ifelse(SPE=="Ostria_virginiana","OSVI",
	ifelse(SPE=="Abies_balsamea","ABBA",
	ifelse(SPE=="Acer_rubrum","ACRU",
	ifelse(SPE=="Acer_pennsylcanicum","ACPE",
	ifelse(SPE=="Querbus_rubra","QURU",SPE)))))))

#Add site	
SPEsite<-paste(SITE,SPE,sep="-")
SPEsiteu<-sort(unique(SPEsite))
#Add an empty class
SPEsiteu <-c(SPEsiteu," - ")
#Maname SPE order
SPEo<-c(1,2,3,4,6,7,8,9,10,5)
SPEsiteu <-SPEsiteu[order(SPEo)]
#Labels
SPEl<-matrix(unlist(strsplit(SPEsiteu,split="-")),ncol=2,byrow=T)[,2]

FRspe<-sapply(SPEsiteu,function(x){
	return(apply(subset(t(FR[,vSAMP]), SPEsite ==x),2,mean))
})

barplot(FRspe,las=2,ylab="Part of diversity",ylim=c(0,1.1),cex.axis=0.8,col= COLb,border=NA,main="Per tree species",cex.main=1,cex.names=0.8,names= SPEl)
text((c(2.5,8)-0.5)*1.25,c(1.05,1.05),SITEu,cex=1)

###6 per sampling time per site (excluding METH communities)

TIME<-as.numeric(METAsamp[vSAMP,8])
TIMEu<-sort(unique(TIME))

#Actual dates
TIMEru<-c("20 Jun.","27 Jun.","16 Jul.","6 Aug.","16 Aug.","7 Sep.","20 Sep.","18 Oct.")

TIMEr<-sapply(TIME,function(x){
	return(subset(TIMEru, TIMEu==x))
})

#Add an empty class
TIMEru <-c(TIMEru," ")
#Maname TIME order
TIMEo<-c(6,1,7,2,8,3,9,4,5)
TIMEru <-TIMEru[order(TIMEo)]

FRtime<-sapply(TIMEru,function(x){
	return(apply(subset(t(FR[,vSAMP]), TIMEr ==x),2,mean))
})

barplot(FRtime,las=2,ylab="Part of diversity",ylim=c(0,1.1),cex.axis=0.8,col= COLb,border=NA,main="Per sampling time",cex.main=1,cex.names=0.8)
text((c(2.5,7.5)-0.5)*1.25,c(1.05,1.05),SITEu,cex=1)

dev.off()



