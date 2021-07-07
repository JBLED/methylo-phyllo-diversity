library(seqinr)
library(dada2)

pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-1"

pD2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-2"

#Description of rpoB reference database built  for Alphaproteobacteria

#####1: retrieve rpoB reference database from https://link.springer.com/article/10.1186/s12866-019-1546-z#additional-information, downloaded from http://genoweb.toulouse.inra.fr/frogs_databanks/assignation/rpoB/

#####2: Recente Rhizobiales references were manually added:

## https://link.springer.com/article/10.1007/s10482-019-01357-6 >Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Lichenibacteriaceae;Lichenibacterium;Lichenibacterium_ramalinae;
## >Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Lichenibacteriaceae;Lichenibacterium;Lichenibacterium minor;
##>Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Beijerinckiaceae;GenusRH  #strains RH new genus AL1, AL8, CH11 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6912076/#!po=54.8780

#####3 problem with some ASVs assigned to Methylobacteriaceae while related to Beijerinckiaceae, and Bradyrhizobiaceae/Methylocystaceae being polyphyletic. Reformat Alphaproteobacteria database: remove all Methylobacteriaceae references and replace them by those used in first part (1-Phylo-of-plant-ass-Methylo-diversity), file rpoB-update-sept-2020.fas -  and ad strains from 2017 pilot survey for which rpoB was sequenced 

#####4 manually remove problematic sequences after quick and dirty phylogeny: 

# Salinarimonas rosea  (Methylobacteriaceae sequence classified as Bradyrhizobiaceae)
# Terasakiella pusilla and Beijerinckia sp. 28-YEA-48  (classified as Methylocystaceae and Beijerinckia but too divergent)
# Bosea sp. 117 , do not cluster with other Bosea
# Rhodoblastus acidophylus (Bradyrhizobiaceae but cluster with Beijerinckiaceae)

#NOTE THAT METHYLOCYSTACEAE, BEIJERINCKIACEAE AND BRADYRHIZOBIACEAE ARE STILL POLYPHYLETIC BUT GENERA ARE QUITE CONSISTENT WITHIN CLADES

#Step 1: initial taxonomic annotation
#Retrieve ASV table

#Trimming parameters
TFx=250
TRx=200
TLx=20

#sequences randomly sampled in each sample
NS=3000

Prare=0.005 #only keep ASV that have at least Prare% of relative abundance in at least one sample

setwd(pD1)

seqtab.nochimR<-read.table(paste("RAREFIED_SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)


#ASV SEQUENCES
ASVs<-row.names(seqtab.nochimR)

#Sample names
sample.names<-colnames(seqtab.nochimR)

###Taxonomic annotation

#Minimum bootstrapp support for taxonomic assignation
minBoot=50

taxa <- assignTaxonomy(ASVs, paste(pD2,"rpoB_Alphaproteobacteria-NEW.fasta",sep="/"), multithread=F,outputBootstraps=T,minBoot= minBoot)

boot<-taxa$boot #bootstrapp support for each taxonomic assignation
taxa <-taxa$tax # taxonomic assignation
rownames(taxa)<-NULL
rownames(boot)<-NULL

colnames(taxa)<-c("Kingdom","Phyllum","Class","Order","Family","Genus","Species")
colnames(boot)<-c("Kingdom","Phyllum","Class","Order","Family","Genus","Species")

setwd(pD2)

write.table(taxa,paste("Annotation-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

write.table(boot,paste("Boot-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

#######Step 2: Summary of initial ASV taxonomy

setwd(pD2)

taxa<-read.table(paste("Annotation-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

boot<-read.table(paste("Boot-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

###ASV taxonomy

seqASV<-sapply(ASVs,function(x){
	return(list(unlist(
		strsplit(x,split=""))))
}) #convert ASV in nucleotide sequence

tASVs<-sapply(c(1:length(ASVs)),function(x){
	tx<-as.vector(taxa[x,c(4:7)])
	tx<-subset(tx,is.na(tx)==F)
	return(paste(c(tx,"ASV",x),collapse="-"))
})	#ASV taxonomy

write.fasta(seqASV,names= tASVs,paste("Sequences-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.fas",sep=""))

###Number of ASVs and ASV relative abundance in METH communities and Phyllosphere samples

iMETH<-c(1:4)
aMETH<-as.vector(apply(seqtab.nochimR[, iMETH],1,max))
aMETH<-ifelse(aMETH==0,0,1)
fMETH<-as.vector(apply(seqtab.nochimR[, iMETH],1,mean)/NS)
iSAMP<-c(5:length(sample.names))
aSAMP<-as.vector(apply(seqtab.nochimR[, iSAMP],1,max))
aSAMP <-ifelse(aSAMP ==0,0,1)
fSAMP<-as.vector(apply(seqtab.nochimR[, iSAMP],1,mean)/NS)

ALL<-cbind(aMETH,fMETH,aSAMP,fSAMP)

COL<-read.table("COL.txt")

vTAX<-c("Alphaproteobacteria","Rhizobiales","Methylobacteriaceae","Methylobacterium")



par(mfrow=c(1,4))

SUMall<-matrix(unlist(sapply(c(4:7),function(x){
	LX=colnames(taxa)[x]
	VX= vTAX[x-3]
	VTX<-as.vector(taxa[,x-1])
	TX<-as.vector(taxa[,x])
	
	ALLx<-subset(ALL, VTX==VX)
	TX<-subset(TX, VTX==VX)
	
	TX<-ifelse(is.na(TX)==T,"unknown",TX)
	TU<-sort(unique(TX))
	COLtu<-round(seq(length(COL[,1]),1,length.out=length(TU)),0)
	COLtu<-COL[COLtu,]
	COLtu<-rgb(COLtu[,1], COLtu[,2], COLtu[,3])
	COLtu<-ifelse(TU =="unknown","grey", COLtu)
	SUM<-t(sapply(TU,function(x){
		return(apply(subset(ALLx,TX==x),2,sum))
	}))
	par(mar=c(4,4,2,1),bty="n")
plot(-10,-10,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
	legend(x = "right",cex=0.5, horiz=F,
		legend= TU,box.col="white",
		border="white",col= COLtu,text.font=1, fill=COLtu)
	par(mar=c(4,4,1,9),bty="n",new=T)
	barplot(SUM[,c(2,4)],las=1,ylab="Relative abundance",cex.axis=0.8,main= VX,col=COLtu,border=NA,names=c("METH","SAMP"),ylim=c(0,1))
	return(t(cbind(LX ,TU, SUM)))
})),ncol=6,byrow=T)	
	
colnames(SUMall)<-c("Level","Taxa","ASV-METH","R-METH","ASV-SAMP","F-SAMP")	

write.table(SUMall,paste("SUMMARY-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

# STEP 4
	
#Remove unaligned ASVs in MEGA and perform phylogeny --> Sequences-ref=Alphapro-minBoot=50_TF=250_TR=200_TL=20_NS=3000_prare=0.005-wo-NEGcontrols-wo-unaligned.fas

#Perform curation of taxonomy based on phylogeny: 
#Sequences-ref=Alphapro-minBoot=50_TF=250_TR=200_TL=20_NS=3000_prare=0.005-wo-NEGcontrols-wo-unaligned.mtsx
#Sequences-ref=Alphapro-minBoot=50_TF=250_TR=200_TL=20_NS=3000_prare=0.005-wo-NEGcontrols-wo-unaligned-names-for-TAXO.pdf

#In excel, generate Phylogeny-Taxo.txt by monifying manually tanoxomy of ASVs according to phylogenetic placement- Column 1: order of sequences in the phylogeny - Column 2: ASV old Taxonomy - Column 3: ASV name - Columns 5 to 8: manually currated taxonomy from phylogeny

#Retrieve Taxonomic assignation based on phylogeny
PhyloTaxo<-read.table("Phylogeny-Taxo.txt",header=T)

####Retrieve #Taxonomic assignation based on rpoB database
minBoot=50

TFx=250
TRx=200
TLx=20
#sequences randomly sampled in each sample
NS=3000
Prare=0.005 #only keep ASV that have at least Prare% of relative abundance in at least one sample


taxa <- read.table(paste("Annotation-ref=Alphapro-minBoot=", minBoot,"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

setwd(pD1)

seqtab.nochimR<-read.table(paste("RAREFIED_SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)
#ASV SEQUENCES
ASVs<-row.names(seqtab.nochimR)

#Sample names
sample.names<-colnames(seqtab.nochimR)

#ASVs that were removed from phylogeny (use assignation from initial rpoB database)
setwd(pD2)

iREM<-setdiff(c(1:length(ASVs)), PhyloTaxo[,3])

taxaREM<-cbind(iREM,taxa[iREM,c(1:6)])

colnames(taxaREM)<-c("ASV",colnames(taxaREM)[2:length(taxaREM[1,])])

#ASVs for which taxonomy was curated with phylogeny
taxaPHY <-cbind(PhyloTaxo[,3], "Bacteria","Proteobacteria",PhyloTaxo[,c(5,6,7,8)])

colnames(taxaPHY)<-colnames(taxaREM)

taxaALL<-rbind(taxaPHY, taxaREM)
taxaALL<-taxaALL[order(taxaALL[,1]),]

write.table(taxaALL ,paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))


###Number of ASVs and ASV relative abundance in METH communities and Phyllosphere samples

iMETH<-c(1:4)
aMETH<-as.vector(apply(seqtab.nochimR[, iMETH],1,max))
aMETH<-ifelse(aMETH==0,0,1)
fMETH<-as.vector(apply(seqtab.nochimR[, iMETH],1,mean)/NS)
iSAMP<-c(5:length(sample.names))
aSAMP<-as.vector(apply(seqtab.nochimR[, iSAMP],1,max))
aSAMP <-ifelse(aSAMP ==0,0,1)
fSAMP<-as.vector(apply(seqtab.nochimR[, iSAMP],1,mean)/NS)

ALL<-cbind(aMETH,fMETH,aSAMP,fSAMP)


COL<-read.table("COL.txt")

vTAX<-c("Alphaproteobacteria","Rhizobiales","Methylobacteriaceae","Methylobacterium")


par(mfrow=c(1,3))

SUMall<-matrix(unlist(sapply(c(5:7),function(x){
	LX=colnames(taxaALL)[x]
	VX= vTAX[x-4]
	VTX<-as.vector(taxaALL[,x-1])
	TX<-as.vector(taxaALL[,x])
	
	ALLx<-subset(ALL, VTX==VX)
	TX<-subset(TX, VTX==VX)
	
	TX<-ifelse(is.na(TX)==T,"unknown",TX)
	TU<-sort(unique(TX))
	COLtu<-round(seq(length(COL[,1]),1,length.out=length(TU)),0)
	COLtu<-COL[COLtu,]
	COLtu<-rgb(COLtu[,1], COLtu[,2], COLtu[,3])
	COLtu<-ifelse(TU =="unknown","grey", COLtu)
	SUM<-t(sapply(TU,function(x){
		return(apply(subset(ALLx,TX==x),2,sum))
	}))
	par(mar=c(4,4,2,1),bty="n")
plot(-10,-10,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
	legend(x = "right",cex=0.5, horiz=F,
		legend= TU,box.col="white",
		border="white",col= COLtu,text.font=1, fill=COLtu)
	par(mar=c(4,4,1,9),bty="n",new=T)
	barplot(SUM[,c(2,4)],las=1,ylab="Relative abundance",cex.axis=0.8,main= VX,col=COLtu,border=NA,names=c("METH","SAMP"),ylim=c(0,1))
	return(t(cbind(LX ,TU, SUM)))
})),ncol=6,byrow=T)	
	
colnames(SUMall)<-c("Level","Taxa","ASV-METH","R-METH","ASV-SAMP","F-SAMP")	

write.table(SUMall,paste("SUMMARY-PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))
