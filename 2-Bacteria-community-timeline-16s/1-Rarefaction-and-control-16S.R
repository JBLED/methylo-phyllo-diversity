library(seqinr)

#Work directory
pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-1"

####META DATA FOR SAMPLES
setwd(pD1)

META<-read.table("METADATA_16S.txt",header=T)

####ASV tables

p0 <-"DADA2-ON-PHYLLO-1-Analyze_16S_JB+ISA_Feb-2020-PSEUDO-bigdata"
TFx=200
TRx=200
TLx=20

seqtab.nochim<-read.table(paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""),header=T)

sample.names<-colnames(seqtab.nochim)


####STEP 1 : Merge ReadTrack files and add the number of reads after chimeras removing
#number of reads per sample
NR<-as.numeric(apply(seqtab.nochim,2,sum))

#Sequencing runs
MISEQr<-c("JB-PL0","JB-PL1","ISA-2014")

Reads_Track<-matrix(unlist(sapply(MISEQr,function(x){
	RT<-read.table(paste("Reads_Track_RUN=",x,
		"_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))
	sn<-row.names(RT)
	return(list(t(cbind(RT,sapply(sn,function(x){
		snx<-paste(unlist(strsplit(x,split="-")),collapse=".")
		return(subset(NR, sample.names == snx))
	})))))
})),ncol=6,byrow=T)

colnames(Reads_Track)<-c("input","filtered","denoisedF","denoisedR","merged","NoChim")
rownames(Reads_Track)<-sample.names

write.table(Reads_Track,paste("Reads_Track_16S-ALL_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

#########FROM THIS STEP, ONLY SAMPLES FROM BATCH "JB-PL1" WILL BE CONDIDERED  
####STEP 2: Extract info for ASVs from negative controls

Nneg<-as.vector(subset(META[,17],META[,4]=="NEG"))#names of negative control samples ordered by sequencing batch

Sneg<-sapply(Nneg ,function(x){
	return(subset(c(1:length(sample.names)), sample.names==x))	
})

POSneg<-sort(subset(c(1:length(seqtab.nochim[,1])),as.vector(seqtab.nochim[, Sneg])!=0))


NEGseq<-rownames(seqtab.nochim)[POSneg]
NSneg<-paste("ControlASV-", POSneg,"-NEG=",as.vector(seqtab.nochim[POSneg, Sneg]),sep="")

NEGseq<-sapply(NEGseq,function(x){
	return(list(unlist(strsplit(x,split=""))))
})

write.fasta(NEGseq,names= NSneg,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-NEGonly.fas",sep=""))

seqtab.nochimNEG<-t(t(seqtab.nochim[POSneg, Sneg]))

row.names(seqtab.nochimNEG)<-NSneg
colnames(seqtab.nochimNEG)<-Nneg
write.table(seqtab.nochimNEG,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-NEGonly.txt",sep=""))





####STEP 3: Extract info for ASVs from positive controls (Meth communities) 


Nmeth<-sort(as.vector(subset(META[,17],META[,1]=="METH"))) #names of METH community sample ordered by sequencing batch

Smeth<-sapply(Nmeth ,function(x){
	return(subset(c(1:length(sample.names)), sample.names==x))	
})

POSmeth<-sort(subset(c(1:length(seqtab.nochim[,1])),as.vector(seqtab.nochim[, Smeth])!=0))

POSseq<-rownames(seqtab.nochim)[POSmeth]

POSseq <-sapply(POSseq,function(x){
	return(list(unlist(strsplit(x,split=""))))
})

Npos<-paste("ControlASV-", POSmeth,"-METH=",as.vector(seqtab.nochim[POSmeth, Smeth]),sep="")

seqtab.nochimMETH<-t(t(seqtab.nochim[POSmeth, Smeth]))
colnames(seqtab.nochimMETH)<-Nmeth
row.names(seqtab.nochimMETH)<-Npos

write.fasta(POSseq,names= Npos,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-METHonly.fas",sep=""))


write.table(seqtab.nochimMETH,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-METHonly.txt",sep=""))

###STEP 4: (IN EXCELL) Determine a threshold for false ASVs. Be highly conservative and only keep ASVs that have at least 0.3% of relative abundance in at least one sample (Negative controls excluded), which allow to remove rare ASVs unrelated to the most abundant taxa found in METH community (details in Threshold-16S-METH.xlsx). 

#NOTE: METH COMMUNITY WAS CONTAMINATED WITH THERMOTOGALES (36.4% OF DIVERSITY, 3 ASVs)
#SPHINGOMONADALES REPRESENTS 5.4% (1 ASV)
#ENTEROBACTERALES (Escherichia)  REPRESENTS 18.5% (4 ASVs)
#METHYLOBACTERIUM REPRESENTS 38.3% (18 ASVs)

# STEP 5: rarefaction

#Samples from run PL1

Nsamp<-sort(as.vector(subset(META[,17],META[,19]=="PL1-16S")))

Ssamp<-sapply(Nsamp,function(x){ #index in SeqTab
	return(subset(c(1:length(sample.names)), sample.names==x))
})

#######Extract metadata for each sample

METAsamp<-t(sapply(Nsamp,function(x){
	Mx<-t(as.vector(subset(META, META[,17]==x)))
	return(Mx)
}))

# rarefaction
# Calculate ASV relative abundance for each sample
# Only keep ASVs with maximum frequency across samples >0.3% (conservative threshold; see Summary-16S-JB.xlsx/METH-Controls)
# Make rarefaction curves for samples, positive and negative controls (estimate the number of sequences to sample before reaching a plateau in diversity)
# Remove negative controls (less than 100 sequences, mostly cross contamination among samples)
# Perform rarefaction (sample randomly 3,000 sequences in each sample, conservative threshold corresponding to a bit less than the minimum number of sequences per sample)


seqtab.nochim <-seqtab.nochim[,c(Ssamp)]
sample.names<-Nsamp

##Calculate relative abundance in each sample
seqtab.nochimREL<-sapply(c(1:length(seqtab.nochim[1,])),function(x){
	return(as.numeric(seqtab.nochim[,x])/sum(as.numeric(seqtab.nochim[,x])))
})

##Calculate maximum relative abundance per ASV accross sample
maxREL<-apply(seqtab.nochimREL,1,max)

###Remove ASVs that do no reach the threshold calculated from = maximum frequency of false ASVd from METH positive controls  (in SUMMARY-community-analysis.xlsx/3-ASV-Meth-communities); conservative value

THx=0.003

seqtab.nochim<-subset(seqtab.nochim, maxREL> THx)

##RAREFACTION

#number of ASV and sequences in each sample

Nx<-cbind(t(sapply(c(1:length(seqtab.nochim[1,])),function(x){
	Fx<-as.vector(seqtab.nochim[,x])
	return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
})),
sapply(Ssamp,function(x){
	return(length(intersect(x,Smeth))+length(intersect(x,Sneg))*2)
})+1
)


pdf("Rarefaction-16S-PL1.pdf",width=4,height=4)

Kx=50 #Number of replicates resampling

Lx=20 #number of ressampling value to test in the range 1-total number of sequence in the sample

Cpx<-c("black","red","blue") #color code for actual values (samples, positive controls, negative controls)

Csx<-c(rgb(0.5,0.5,0.5,0.7),rgb(0.5,0,0,0.7),rgb(0,0,0.5,0.7)) #color code for estimated values after resampling (samples, positive controls, negative controls)

par(mar=c(4,4,1,1),bty="n")

plot(1,1,
	ylim=c(0,max(Nx[,1])),xlim=c(0,max(Nx[,2])),
	xlab="Sampled sequences",ylab="ASVs",cex.axis=0.8,las=1,log="",type="l")

sapply(c(1:length(seqtab.nochim[1,]))[order(Nx[,1],decreasing=T)],function(x){
	print(x)
	Fx<-as.vector(seqtab.nochim[,x])
	Fx<-subset(Fx,Fx!=0)
	Fx<-unlist(sapply(c(1:length(Fx)),function(x){
		return(seq(x,x,length.out=Fx[x]))
	}))
	Rx<-round(seq(1,Nx[x,2],length.out=Lx))
	RRk<-apply(sapply(Rx,function(x){
		#print(x)
		RRx=x
		return(sapply(c(1:Kx),function(x){
			return(length(unique(sample(Fx,
				RRx,replace=F))))
		}))
	}),2,mean)
	zc<-Nx[x,3]
	points(Rx,RRk,type="l",col=Csx[zc])
	points(Nx[x,2],Nx[x,1],col=Cpx[zc],pch=19,cex=0.5)
})	

dev.off()

#Remove negative controls

seqtab.nochim <-t(subset(t(seqtab.nochim),Nx[,3]!=3))
#seqtab.nochimREL <-t(subset(t(seqtab.nochimREL),Nx[,3]!=3))
METAsamp<-subset(METAsamp,Nx[,3]!=3)
sample.names<-subset(sample.names,Nx[,3]!=3)
#Nx<-subset(Nx,Nx[,3]!=3)

#Remove ASVs from negative controls
seqtab.nochim <-subset(seqtab.nochim,apply(seqtab.nochim,1,max)>=1)


#Sample randomly NS sequences in each sample (Check NS is lower than the minimum number of sequence per sample, second column of Nx)
NS=5000
NR=length(seqtab.nochim[,1])

seqtab.nochimR<-sapply(c(1:length(seqtab.nochim[1,])),function(x){
	#print(x)
	Fx<-as.vector(seqtab.nochim[,x])
	Ix=subset(c(1:length(Fx)),Fx!=0)
	Fx<-subset(Fx,Fx!=0)
	Fx<-unlist(sapply(c(1:length(Fx)),function(x){
		return(seq(Ix[x],Ix[x],length.out=Fx[x]))
	}))

	Fs<-sample(Fx, NS,replace=F)
	return(sapply(c(1: NR),function(x){return(length(subset(Fs,Fs==x)))}))
})
rownames(seqtab.nochimR)<-rownames(seqtab.nochim)

#Remove unsampled ASVs
seqtab.nochimR <-subset(seqtab.nochimR,apply(seqtab.nochimR,1,sum)>=1)



colnames(seqtab.nochimR)<-sample.names
rownames(METAsamp)<-sample.names
write.table(METAsamp,"METADATA-PL1-wo-NEGcontrols.txt")
write.table(seqtab.nochimR,paste("RAREFIED_SeqTableNoChim-PL1_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,"-wo-NEGcontrols.txt",sep=""))



#Compare ASV cumulated frequency before and after rarefaction

 #number of ASV and sequences in each sample after rarefaction
NRx<-t(sapply(c(1:length(seqtab.nochimR[1,])),function(x){
	Fx<-as.vector(seqtab.nochimR[,x])
	return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
}))


pdf(paste("16S-PL1-ASV-cumulated-frequency-per-sample-Prare=", THx,".pdf",sep=""),width=6,height=3)


par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))
plot(1,-1000,log="x",ylim=c(0.02,1),las=1,xlim=c(1,max(Nx[,1])),xlab="ASVs (log scale)",ylab="Cumulated frequency",cex.axis=0.8,main=paste("Before rarefaction (n= ",dim(seqtab.nochimREL)[1]," ASVs)",sep=""),cex.main=0.8)

sapply(c(1:length(seqtab.nochimREL[1,])),function(x){
	TYPEx= Nx[x,3]
	colx=ifelse(TYPEx==2,rgb(0.8,0,0,0.2),
		ifelse(TYPEx==3,rgb(0,0,0.8,0.2),
		rgb(0.4,0.4,0.4,0.2)))
	points(cumsum(sort(seqtab.nochimREL[,x],decreasing=T)),
		type="l",col= colx,lwd=2)
})

plot(1,-10000,log="x",ylim=c(0,NS),las=1,xlim=c(1,max(NRx[,1])),xlab="ASVs (log scale)",ylab="Cumulated abundance",cex.axis=0.8,main=paste("After rarefaction (n= ",dim(seqtab.nochimR)[1]," ASVs)",sep=""),cex.main=0.8)

sapply(c(1:length(seqtab.nochimR[1,])),function(x){
	TYPEx= METAsamp[x,1]
	colx=ifelse(TYPEx=="METH",rgb(0.8,0,0,0.2),
		rgb(0.4,0.4,0.4,0.2))
	points(cumsum(sort(seqtab.nochimR[,x],decreasing=T)),
		type="l",col= colx,lwd=2)
})

dev.off()





#STEP 6 Summarize the effect of rarefaction on the number of ASVs per sample

row.names(NRx)<-sample.names
NRx<-cbind(subset(Nx,Nx[,3]!=3),NRx)

NRneg<-cbind(subset(Nx,Nx[,3]==3),NA,NA)
row.names(NRneg)<-names(Sneg)

NRx<-rbind(NRneg, NRx)
NRx<-NRx[,c(3,1,2,4,5)]

seqtab.nochim<-read.table(paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""),header=T)

NRo<-t(sapply(c(1:length(seqtab.nochim[1,])),function(x){
	Fx<-as.vector(seqtab.nochim[,x])
	return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
}))

NRo<-t(sapply(row.names(NRx),function(x){
	return(subset(NRo,colnames(seqtab.nochim)==x))
}))

NRx<-cbind(NRx[,1], NRo, NRx[,c(2:5)])

NRx<-cbind(ifelse(NRx[,1]==3,'neg',ifelse(NRx[,1]==2,'meth','sample')),NRx[,c(2:7)])

colnames(NRx)<-c('Type',"ASVs-before-TH","Seq-before-TH","ASVs-before-raref","Seq-before-raref","ASVs-final","Seq-final")

write.table(NRx,paste("Summary-16S-PL1-raref_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,".txt",sep=""))