library(seqinr)

#Work directory
pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-1"

####META DATA FOR SAMPLES
setwd(pD1)

METArep<-read.table("METADATA_TO_REPLACE.txt",header=T)
META<-read.table("METADATA.txt",header=T) #some samples were inverted in plate pl2

####ASV tables

Mx="rpoB"
TFx=250
TRx=200
TLx=20

#MUST BE UNZIPPED (SeqTableNoChim_TF=250_TR=200_TL=20.txt.zip)
seqtab.nochim<-read.table(paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""),header=T)

sample.names<-colnames(seqtab.nochim)


####STEP 1 : Merge ReadTrack files and add the number of reads after chimeras removing
#number of reads per sample
NR<-as.numeric(apply(seqtab.nochim,2,sum))

#Sequencing runs
MISEQr<-paste("PL",c(0:4),sep="")

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

write.table(Reads_Track,paste("Reads_Track_ALL_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))


####STEP 2: Extract info for ASVs from negative controls

Nneg<-as.vector(subset(META[,17],META[,4]=="NEG"))[c(1,2,4,3)] #names of negative control samples ordered by sequencing batch

Sneg<-sapply(Nneg ,function(x){
	return(subset(c(1:length(sample.names)), sample.names==x))	
})

POSneg<-sort(subset(c(1:length(seqtab.nochim[,1])),as.vector(apply(seqtab.nochim[, Sneg],1,sum))!=0))

NEGseq<-rownames(seqtab.nochim)[POSneg]
NSneg<-paste("ControlASV-", POSneg,"-NEG=",as.vector(apply(seqtab.nochim[POSneg, Sneg],1,sum)),sep="")

seqtab.nochimNEG<-seqtab.nochim[POSneg, Sneg]
colnames(seqtab.nochimNEG)<-Nneg
write.table(seqtab.nochimNEG,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-NEGonly.txt",sep=""))


####STEP 3: Extract info for ASVs from positive controls (Meth communities) compare raw ASV relative abundance among samples


Nmeth<-sort(as.vector(subset(META[,17],META[,1]=="METH")))[c(2,3,4,1)] #names of METH community sample ordered by sequencing batch

Smeth<-sapply(Nmeth ,function(x){
	return(subset(c(1:length(sample.names)), sample.names==x))	
})

POSmeth<-sort(subset(c(1:length(seqtab.nochim[,1])),as.vector(apply(seqtab.nochim[, Smeth],1,sum))!=0))

POSseq<-rownames(seqtab.nochim)[POSmeth]
Npos<-paste("ControlASV-", POSmeth,"-METH=",as.vector(apply(seqtab.nochim[POSseq, Smeth],1,sum)),sep="")

seqtab.nochimMETH<-seqtab.nochim[POSseq, Smeth]
colnames(seqtab.nochimMETH)<-Nmeth

F1<-as.numeric(seqtab.nochimMETH[,1])

pdf("rpoB-METH-community-comparison.pdf",height=4,width=4)
par(bty="n",mar=c(4,4,1,1))
plot(-10,-10,xlim=c(0,max(F1)),ylim=c(0,max(seqtab.nochimMETH)),type="l",las=1,cex.axis=0.8,xlab="ASV abundance (batch 1)",ylab="ASV abundance (batches 2,3,4)", main="Positive controls (METH community)",log="")

sapply(c(2:length(seqtab.nochimMETH[1,])),function(x){
	Fx<-as.numeric(seqtab.nochimMETH[,x])
	points(F1, Fx,pch=19,cex=0.5,col=c("red","green","blue")[x-1])
	DF <- lm(Fx ~ F1)
	ld=sort(unique(F1))
	lines(ld, predict(DF,data.frame(F1 =ld)),
		col=c("red","green","blue")[x-1],lwd=1.5,lty=2)
		
	return(cor.test(Fx, F1)$estimate)	
})
dev.off()

write.table(seqtab.nochimMETH,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-METHonly.txt",sep=""))

###STEP 4: For positive controls, assign ASV taxonomy based on pairwise similarity (PS) between reads (reverse and forward separately) and reference sequences

###Sequences from actual strains that were used to generate the positive control + a couple reference sequences 

REFmeth<-read.fasta("REF-for-METH.fas")

REFmethN<-names(REFmeth)

REFmeth<-t(sapply(REFmeth,function(x){
return(unlist(strsplit(tolower(paste(as.vector(x),collapse="")),split="nnnnnnnnnn")))

}))


#For each ASV in the four meth communities, determine the closest reference sequence and PS - do separately the forward and the reverse to identify potential chimeras

SCORE<-t(sapply(c(1:length(seqtab.nochimMETH[,1])),function(x){
	X=x
	ASVx=unlist(strsplit(tolower(rownames(
		seqtab.nochimMETH))[x],split="nnnnnnnnnn"))
	ASVf<-unlist(strsplit(ASVx[1],split=""))
	ASVr<-unlist(strsplit(ASVx[2],split=""))
	PSx<-t(sapply(c(1:length(REFmeth[,1])),function(x){
		return(c(
		length(subset(ASVf, ASVf==
			unlist(strsplit(REFmeth[x,1],split="")))),
		length(subset(ASVr, ASVr==
			unlist(strsplit(REFmeth[x,2],split=""))))	
		))	
	}))
	return(c(X,paste(subset(c(1:length(REFmeth[,1])), 
		PSx[,1]==max(PSx[,1])),collapse="_"),
	max(PSx[,1]),	
	paste(subset(c(1:length(REFmeth[,1])), 
		PSx[,2]==max(PSx[,2])),collapse="_"),
	max(PSx[,2])))		
}))

colnames(SCORE)<-c("ASV","REFf","PSf","REFr","PSr")

write.table(SCORE ,paste("SCORE_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"-METHonly.txt",sep=""))


###Position of real samples
Nsamp<-setdiff(c(1:length(META[,1])),as.numeric(c(Sneg,Smeth)))

#####Fasta file only with sequences from positive controls and negative controls

SEQcont<-sapply(c(rownames(seqtab.nochimNEG),
	rownames(seqtab.nochimMETH)),function(x){
	return(strsplit(x,split=""))	
})

SEQnames<-c(
sapply(c(1:length(seqtab.nochimNEG[,1])),function(x){
	return(paste(
		paste(c("NEG",x,"F"),collapse="-"),
		paste(as.numeric(seqtab.nochimNEG[x,]),collapse="_"),
		sep="="))
}),
sapply(c(1:length(seqtab.nochimMETH[,1])),function(x){
	return(paste(
		paste(c("METH",x,"F"),collapse="-"),
		paste(as.numeric(seqtab.nochimMETH[x,]),collapse="_"),
		sep="="))
}))

write.fasta(SEQcont,names= SEQnames,paste("CONTROLsequences_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".fas",sep=""))

# STEP 5: In positive controls, calculate the proportion of sequences that are true, chimeras or contamination to determine a conservative threshold  (Rare-ASV-threshold-from-METH-community.xlsx from SCORE_TF=250_TR=200_TL=20-METHonly.txt) 

#######Extract metadata for each sample

#Because of an error in the PL2 sample sheet, replace METAdata from some samples by the right ones using METArep FILE 

METAsamp<-t(sapply(sample.names,function(x){
	Tx<-c(as.vector(subset(METArep[,2], METArep[,1]==x)),x)[1]

	Mx<-t(as.vector(subset(META, META[,17]==Tx)))
	return(Mx)
}))


# STEP 6: rarefaction
# Calculate ASV relative abundance for each sample
# Only keep ASVs with maximum frequency across samples >0.5% (conservative threshold; see SUMMARY-community-analysis.xlsx/sheet 3-ASV-Meth-communities)
# Make rarefaction curves for samples, positive and negative controls (estimate the number of sequences to sample before reaching a plateau in diversity)
# Remove negative controls (less than 100 sequences, mostly cross contamination among samples)
# Perform rarefaction (sample randomly 3,000 sequences in each sample, conservative threshold corresponding to a bit less than the minimum number of sequences per sample)


seqtab.nochim <-seqtab.nochim[,c(Smeth,Sneg,Nsamp)]
METAsamp<-METAsamp[c(Smeth,Sneg,Nsamp),]
sample.names<-sample.names[c(Smeth,Sneg,Nsamp)]

##Calculate relative abundance in each sample
seqtab.nochimREL<-sapply(c(1:length(seqtab.nochim[1,])),function(x){
	return(as.numeric(seqtab.nochim[,x])/sum(as.numeric(seqtab.nochim[,x])))
})

##Calculate maximum relative abundance per ASV accross sample
maxREL<-apply(seqtab.nochimREL,1,max)

###Remove ASVs that do no reach the threshold calculated from = maximum frequency of false ASVd from METH positive controls  (in SUMMARY-community-analysis.xlsx/3-ASV-Meth-communities); conservative value

THx=0.005

seqtab.nochim<-subset(seqtab.nochim, maxREL> THx)

##RAREFACTION

#number of ASV and sequences in each sample

Nx<-cbind(t(sapply(c(1:length(seqtab.nochim[1,])),function(x){
	Fx<-as.vector(seqtab.nochim[,x])
	return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
})),
c(seq(2,2,length.out=length(Smeth)),seq(3,3,length.out=length(Sneg)),seq(1,1,length.out=(length(seqtab.nochim[1,])-(length(Sneg)+length(Smeth))))
)) 


pdf("Rarefaction-rpoB.pdf",width=4,height=4)

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
NS=3000
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
write.table(METAsamp,"METADATA-wo-NEGcontrols.txt")
write.table(seqtab.nochimR,paste("RAREFIED_SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,"-wo-NEGcontrols.txt",sep=""))



#Compare ASV cumulated frequency before and after rarefaction

 #number of ASV and sequences in each sample after rarefaction
NRx<-t(sapply(c(1:length(seqtab.nochimR[1,])),function(x){
	Fx<-as.vector(seqtab.nochimR[,x])
	return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
}))


pdf(paste("rpoB-ASV-cumulated-frequency-per-sample-Prare=", THx,".pdf",sep=""),width=6,height=3)


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

#STEP 7: control rarefaction effect on data quality
#Estimate proportion of true ASV and contaminant after rarefaction from positive controls
#Perform phylogeny, assign ASV to reference sequences, identify contaminant (taxonomy from phylogeny) and estimate the proportion of true and false ASV (contaminant) IN MEGA
#Summarize the effect of rarefaction on the number of ASVs per sample



#Sequences form positive controls only
seqtab.nochimRmeth<-seqtab.nochimR[,c(1:4)]
seqtab.nochimRmeth<-subset(seqtab.nochimRmeth,as.vector(apply(seqtab.nochimRmeth,1,max))!=0)


pdf("rpoB-METH-community-comparison_after-rarefaction.pdf",height=4,width=4)
par(bty="n",mar=c(4,4,1,1))
plot(-10,-10,xlim=c(0,max(seqtab.nochimRmeth)),ylim=c(0,max(seqtab.nochimRmeth)),type="l",las=1,cex.axis=0.8,xlab="ASV abundance (batch 1)",ylab="ASV abundance (batches 2,3,4)", main="Positive controls (METH community)",log="")

F1=seqtab.nochimRmeth[,1]

sapply(c(2:length(seqtab.nochimRmeth[1,])),function(x){
	Fx<-as.numeric(seqtab.nochimRmeth[,x])
	points(F1, Fx,pch=19,cex=0.5,col=c("red","green","blue")[x-1])
	DF <- lm(Fx ~ F1)
	ld=sort(unique(F1))
	#lines(ld, predict(DF,data.frame(F1 =ld)),
		#col=c("red","green","blue")[x-1],lwd=c(7,3,1.5)[x-1],lty=2)
})

dev.off()

###Construct fasta file with sequences from positive contrals after rarefaction and reference requences
SEQcont<-sapply(rownames(seqtab.nochimRmeth),function(x){
	return(strsplit(x,split=""))	
})

SEQnames<-
sapply(c(1:length(seqtab.nochimRmeth[,1])),function(x){
	return(paste(
		paste(c("METH",x,"F"),collapse="-"),
		paste(as.numeric(seqtab.nochimRmeth[x,]),collapse="_"),
		sep="="))
})

SEQcont <-c(SEQcont ,read.fasta("REF-for-METH.fas"))

names(SEQcont)<-c(SEQnames,names(read.fasta("REF-for-METH.fas")))


write.fasta(SEQcont,names= names(SEQcont),paste("METH-sequences+Ref_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,"_after-rarefaction.fas",sep=""))

row.names(seqtab.nochimRmeth)<-SEQnames


write.table(seqtab.nochimRmeth,paste("METH-sequences+Ref_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,"_after-rarefaction.txt",sep=""))

##Perform phylogeny, assign ASV to reference sequences, identify contaminant (taxonomy from phylogeny) and estimate the proportion of true and false ASV (contaminant) IN MEGA

#Summarize the effect of rarefaction on the number of ASVs per sample

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

write.table(NRx,paste("Summary-raref_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=",THx,".txt",sep=""))