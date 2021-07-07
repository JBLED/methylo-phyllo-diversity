# Global meta analysis of diversity based on rpoB. In this script, all ASV are analyzed together (METH community discarded at this point)


pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-1"
pD2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-2"
pD3<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/2-Bacteria-community-timeline-16s/data-3"


########################################################
##STEP 1: Load ASV relative abundance table, Metadata for each sample and format them

setwd(pD1)

####Load Metadata for each sample
METAsamp<-read.table("METADATA-PL1-wo-NEGcontrols.txt")

#Indexes for positive controls
vMETH<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]=="METH")

#Indexes for samples
vSAMP<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]!="METH")

	
####Extract factors from metadata only for SAMPLES
#SITE of sampling
	SITE<-as.vector(METAsamp[vSAMP,3])
	SITEu<-sort(unique(SITE)) 
#SUBSITE of sampling
	SS<-paste(SITE,as.vector(METAsamp[vSAMP,4]),sep="-")
	SSu<-sort(unique(SS))
#HOST TREE SPECIES	
	SPE<-as.vector(METAsamp[vSAMP,6])
	#Replace tree species names by abbreviations
	SPE <-ifelse(SPE=="Acer_saccharum","ACSA",
		ifelse(SPE=="Fagus_grandifolia","FAGR",
		ifelse(SPE=="Ostria_virginiana","OSVI",
		ifelse(SPE=="Abies_balsamea","ABBA",
		ifelse(SPE=="Acer_rubrum","ACRU",
		ifelse(SPE=="Acer_pennsylcanicum","ACPE",
		ifelse(SPE=="Querbus_rubra","QURU",SPE)))))))
#Sampling TIME
	TIME<-as.numeric(METAsamp[vSAMP,8])
	TIMEu<-sort(unique(TIME))
	#Actual dates
	TIMEru<-c("20 Jun.","27 Jun.",
		"16 Jul.","6 Aug.",
		"16 Aug.","7 Sep.",
		"20 Sep.","18 Oct.")
	TIMEr<-sapply(TIME,function(x){
		return(subset(TIMEru, TIMEu==x))
	})
#Extraction batch
	EXT<-as.vector(METAsamp[vSAMP,16])
	EXTu<-sort(unique(EXT))
#PCR batch 
	PCR<-as.vector(METAsamp[vSAMP,18])
	PCRu<-sort(unique(PCR))
#Sequencing
	MISEQ <-as.vector(METAsamp[vSAMP,19])
	MISEQu<-sort(unique(MISEQ))	

#METAdata formated for permanova
FREQ.ENV<-as.data.frame(cbind(SITE,SS,SPE,TIME,EXT,PCR,MISEQ))	
	
##Load SeqTAB of relative abundances

TFx=200
TRx=200
TLx=20
NS=5000
Prare=0.003

seqtab.nochim<-read.table(paste("RAREFIED_SeqTableNoChim-PL1_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

###Only keep relative abundance for samples
FRall<-t(seqtab.nochim[,vSAMP])

###Apply Hellinger transformation on relative abundance to account for rare ASVs
library(vegan)

FRallH <-as.data.frame(decostand(FRall ,method="hellinger", MARGIN=1))

colnames(FRallH)<-paste("ASV", c(1:length(FRallH[1,])),sep="")
colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="")

##Load Taxonomic assignation for each ASV
setwd(pD2)
ASSsimp<-as.matrix(read.table(paste("SILVA-taxonomy_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_Prare=", Prare,".txt",sep=""),header=T))

ASSsimp<-ifelse(is.na(ASSsimp)==T,"unknown", ASSsimp)

taxA<-paste(ASSsimp[,2], ASSsimp[,3], ASSsimp[,4],sep=";")
TAXu<-sort(unique(taxA))


##Average relative abundance of each taxonomic level for each sample

FRsamp <-t(sapply(c(1:length(TAXu)),function(x){
	LX= TAXu[x]
	#Phyllo samples
	Fx<-subset(seqtab.nochim, taxA ==LX)
	return(apply(Fx,2,sum)/NS)
}))

FRall<-cbind(FRsamp[,vMETH] ,apply(FRsamp[,vSAMP],1,mean))
 colnames(FRall)<-c("METH","PHYLLO")

#Remove very rare taxa
FRall<-subset(FRall,apply(FRsamp,1,max)>0.05)
TAXu <-subset(TAXu,apply(FRsamp,1,max)>0.05)
FRsamp <-subset(FRsamp,apply(FRsamp,1,max)> 0.05)

rownames(FRall)<-TAXu
rownames(FRsamp)<-TAXu

setwd(pD2)
COL<-read.table("COL.txt")

COLb<-round(seq(1,length(COL[,1]),length.out=length(TAXu)),0)
COLb <-COL[COLb,]
COLb <-rgb(COLb[,1], COLb[,2], COLb[,3])
COLb <-ifelse(TAXu =="unknown;unknown","grey", COLb)


########################################################
##STEP 2: Perform  PCA and permanova on all samples together to weight the relative contribution of SITE, HOST TREE SPECIES and TIME


########################
######PCA ANALYSIS OF COMMUNITY

library(ade4)

#Display PCA
ACPh=dudi.pca(FRallH, scannf=F)
part<-round(100* ACPh $eig/sum(ACPh $eig),2)


x1=min(ACPh $li[,1])
x2=max(ACPh $li[,1])
y1=min(ACPh $li[,2])
y2=max(ACPh $li[,2])

#Factor to apply to Y coordinates in pie charts and contributions
	FACy=((y2-y1)/(x2-x1))

setwd(pD3)

pdf("PCA-per-sampling-site-Hellinger-transformation_MAIN-ORDER-frequency.pdf",height=6,width=6)
par(mar=c(4,4,1,1),bty="n")
plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="PCA per sampling site, Hellinger transformation",cex.main=1)


ordispider(ACPh $li[,c(1,2)],groups= SITE,col= "black",lty=1,label= F,cex=0.5)


#Pie Chart with relative abundance of each main tax per sample
sapply(c(1:length(ACPh $li[,1])),function(x){
	##Coordinates in the plot
	COORDc<-as.numeric(ACPh $li[x,c(1:2)])
	#relative abundance for each factor of the variable
	PCCc<-as.numeric(FRsamp[,x])
	#relative size of the pie chart (radius), proportional to the part of variance explained by the factor
	Sc=0.35
	CUMc<-c(0,cumsum(PCCc))/sum(PCCc)
	sapply(c(1:length(PCCc)),function(x){
		COLx<-COLb[x]
		z1<-CUMc[x]
		z2<-CUMc[x+1]
		z12<-seq(z1,z2,length.out=100)
		ang12<-2*pi*(z12)
		x12<-c(COORDc[1],Sc*sin(ang12)+COORDc[1])
		y12<-c(COORDc[2],Sc*cos(ang12)* FACy +COORDc[2])
		polygon(x12, y12,col= COLx,border=NA)	
	})	
})



ordiellipse(ACPh $li[,c(1,2)],groups= TIMEr,border="black",col= NA,lty=2,draw="polygon",label=T,cex=0.5,lwd=1,alpha=0)

ordiellipse(ACPh $li[,c(1,2)],groups= SPE,border=NA,col= rgb(0.5,0.5,0.5),lty=1,alpha=50,draw="polygon",label=T,cex=0.5)

#ordiellipse(ACPh $li[,c(1,2)],groups= TIME,col= NA,lty=2,draw="polygon",label=T,cex=0.5)
	

text(c(10,-10),c(20,20),SITEu,col="black",font=2,cex=1)

#Plot ASV contribution in variance (calculate mean per TAXA)

#Centre of the contribution plot in the main plot
Xc=4
Yc=15
nz=40 #number of discrete classes
Zc<-seq(0,1,length.out=(nz+1))
angc<-2*pi*(Zc)#convert in RAD

	#Factor for contribution
	Fz=0.25
	#"Scale" for strength
	#Size of the central circle (radius)
	Sct=0.5
	Xs= Sct*sin(2*pi*(seq(0,1,length.out=1000)))+Xc
	Ys= Sct*cos(2*pi*(seq(0,1,length.out=1000)))* FACy +Yc
	text(Xc, Yc-3,"Taxa contribution",cex=0.7,font=2)

#Convert contributions orientation in RAD
CONT<-ACPh$co
CONTl<-sqrt(CONT[,1]^2+CONT[,2]^2)
Cz<-ifelse(CONT[,1]>=0,
	acos(CONT[,2]/CONTl),
	(2*pi)-acos(CONT[,2]/CONTl))
Z<-Cz/(2*pi)
Z<-ifelse(is.nan(Z)==T,0,Z)
Z<-round(nz*Z)

#Display polygons of contribution in the decreasing order of ASV srelative abundance per Taxa
#Calculate total contribution per TAXA per RAD class
CONTtax<-t(sapply(c(1:length(TAXu)),function(x){
	CONTx<-subset(CONT, taxA == TAXu[x])
	COLx<-COLb[x]
	CONTlx<-subset(CONTl, taxA == TAXu[x])
	Zx<-subset(Z, taxA == TAXu[x])
	#Calculate total contribution per discret RAD class
	CONTz<-sapply(c(0:nz),function(x){
		return(sum(subset(CONTlx, Zx==x)))
	})
	return(CONTz)
}))

CONTtax <-sapply(c(1:length(CONTtax[1,])),function(x){
	return(cumsum(CONTtax[,x]))
})*Fz+ Sct


sapply(c(length(CONTtax[,1]):1),function(x){
	COLx= COLb[x]
	CONTz<-CONTtax[x,]
	#Convert in polygon coordinates
	Xx<-CONTz*sin(angc)+Xc
	Yx<-CONTz*cos(angc)* FACy +Yc
	#polygon(Xx,Yx,col=rgb(col2rgb(COLx)[1]/255,
	#		col2rgb(COLx)[2]/255,
	#		col2rgb(COLx)[3]/255,0.3),border= COLx,lwd=1)
	polygon(Xx,Yx,col=COLx,border= NA)	
})

	CONTz<-CONTtax[length(CONTtax[,1]),]
	#Convert in polygon coordinates
	Xx<-CONTz*sin(angc)+Xc
	Yx<-CONTz*cos(angc)* FACy +Yc
	polygon(Xx,Yx,col=NA,border= "black",lwd=0.5)


polygon(Xs,Ys,col="white",border="black",lwd=0.5)

dev.off()

###General permanova

np=10000 #number of permutations
nc= 2 #number of cores to use

print(date())	
ADONISallmiseq<-adonis(FRallH ~ SITE*SPE*TIME+SS+SS*SPE+SS*TIME, 
	data= FREQ.ENV,permutations=np,
	method="bray"	,parallel=nc)		
print(date())

write.table(as.data.frame(ADONISallmiseq$aov),paste("Global-PERMANOVA_",np,"_permutations_Hellinger-transformation.txt",sep=""))

###permanova on Proteobacteria

FRpH<-t(subset(t(FRallH),ASSsimp[,2]=="Proteobacteria"))

print(date())	
ADONISpmiseq<-adonis(FRpH ~ SITE*SPE*TIME+SS+SS*SPE+SS*TIME, 
	data= FREQ.ENV,permutations=np,
	method="bray"	,parallel=nc)	
print(date())

write.table(as.data.frame(ADONISpmiseq $aov),paste("Proteobacteria-PERMANOVA_",np,"_permutations_Hellinger-transformation.txt",sep=""))

###permanova on AlphaProteobacteria

FRaH<-t(subset(t(FRallH),ASSsimp[,3]=="Alphaproteobacteria"))

print(date())	
ADONISamiseq<-adonis(FRaH ~ SITE*SPE*TIME+SS+SS*SPE+SS*TIME, 
	data= FREQ.ENV,permutations=np,
	method="bray"	,parallel=nc)		
print(date())

write.table(as.data.frame(ADONISamiseq $aov),paste("Alphaproteobacteria-PERMANOVA_",np,"_permutations_Hellinger-transformation.txt",sep=""))

###permanova on Rhizobiales

FRrH<-t(subset(t(FRallH),ASSsimp[,4]=="Rhizobiales"))

print(date())	
ADONISrmiseq<-adonis(FRrH ~ SITE*SPE*TIME+SS+SS*SPE+SS*TIME, 
	data= FREQ.ENV,permutations=np,
	method="bray"	,parallel=nc)	
print(date())

write.table(as.data.frame(ADONISrmiseq $aov),paste("Rhizobiales-PERMANOVA_",np,"_permutations_Hellinger-transformation.txt",sep=""))

###permanova on Methylobacterium

FRmH<-t(subset(t(FRallH),ASSsimp[,6]=="Methylobacterium-Methylorubrum"))

print(date())	
ADONISmmiseq<-adonis(FRmH ~ SITE*SPE*TIME+SS+SS*SPE+SS*TIME, 
	data= FREQ.ENV,permutations=np,
	method="bray"	,parallel=nc)	
print(date())

write.table(as.data.frame(ADONISmmiseq $aov),paste("Methylobacterium-PERMANOVA_",np,"_permutations_Hellinger-transformation.txt",sep=""))


########################################################
##STEP 3: Test for significant increase or decrease of ASV abundance in function of time



	COLt<-as.vector(sapply(taxA,function(x){
		return(c(subset(COLb,TAXu==x),"white")[1])
	}))
pdf("cratere_plot_ASV-time-increase-SILVA-taxo.pdf",width=7,height=2.5)
par(mfrow=c(1,3),mar=c(4,4,1,1),bty="n")

COEF<-matrix(unlist(sapply(SITEu,function(x){
	X=x
	
	FRx<-subset(FRallH,SITE==X)
	SSx<-subset(SS,SITE==X)
	SPEx<-subset(SPE,SITE==X)
	TIMEx<-subset(TIME,SITE==X)
	
	#Remove ASVs absent from this site
	Ai<-subset(c(1:length(FRx[1,])),apply(FRx,2,max)!=0)

	COEFx<-t(sapply(Ai,function(x){
		mod<-lm(FRx[,x] ~ TIMEx*SPEx)
		return(c(x,mod$coefficients[c(1:2)],
			as.data.frame(anova(mod))[1,5]))
	}))
	
	Ao<-setdiff(c(1:length(FRx[1,])),Ai)
	COEFo<-cbind(Ao,0,0,NA)
	colnames(COEFo)<-colnames(COEFx)
	
	COEFx<-rbind(COEFo, COEFx)
	COEFx<-COEFx[order(COEFx[,1]),]
	
	COLx<-ifelse(COEFx[,4]<=0.05, COLt,rgb(0,0,0,0.1))
	CEXx<-sqrt(apply(FRx,2,mean)/pi)
	
	plot(COEFx[,3], COEFx[,4],
		cex.axis=0.8,las=1,
		xlab="estimate F~Date",ylab="p.value",
		log="y",pch=19,
		ylim=c(1,min(COEFx[,4],na.rm=T)),
		col=COLx,cex= CEXx*5,main=X)
	
	segments(0,1,0,1/10^10,col="black",lty=2)
		
	return(list(COEFx))
})),ncol=8)

COEF<-ifelse(is.nan(COEF)==T,1000, COEF)

CEX<-sqrt(apply(FRallH,2,mean)/pi)
pmin<-ifelse(COEF[,4]<=0.05,1,0)+ifelse(COEF[,8]<=0.05,1,0)

COLone<-sapply(COLt,function(x){
	return(rgb(col2rgb(x)[1]/255,
		col2rgb(x)[2]/255,
		col2rgb(x)[3]/255,0.3))
})

COLx<-ifelse(pmin==2,COLt,
	ifelse(pmin==1, COLone,rgb(0,0,0,0.1)))

plot(COEF[,3],COEF[,7],
	cex.axis=0.8,las=1,pch=19,
	xlab="Estimate in MSH",ylab="Estimate in SBL",
	cex= CEX*5,main="Site comparison",col=COLx)
	
segments(-1,-1,1,1,col="grey",lty=2)
segments(-1,0,1,0,col="black",lty=2)
segments(0,1,0,-1,col="black",lty=2)	
	
dev.off()

COEFint<-subset(COEF,COEF[,4]<=0.05)
COEFint<-subset(COEFint, COEFint[,8]<=0.05)

COEFmsh<-subset(COEF,COEF[,4]<=0.05)
COEFmsh <-subset(COEFmsh, COEFmsh[,8]>0.05)
COEFmsh <-subset(COEFmsh, COEFmsh[,3]>0)

COEFsbl<-subset(COEF,COEF[,4]>0.05)
COEFsbl <-subset(COEFsbl, COEFsbl[,8]<=0.05)
COEFsbl <-subset(COEFsbl, COEFsbl[,7]>0)

par(mar=c(1,4,1,1),bty='n',mfrow=c(5,5))
sapply(COEFint[,1],function(x){
	plot(as.factor(paste(SITE,TIME)),
	FRallH[,x],las=2,cex.lab=0.8,
	main=paste("ASV",x,sep=""),cex.main=1,
	col=COLone[x],border=COLt[x],xaxt="n",xlab="",ylab="")	
})


FRint<-sapply(COEFint[,1],function(x){
	print(x)
	FRx<-as.numeric(seqtab.nochim[x,vSAMP])/NS
	return(sapply(sort(unique(paste(SITE,TIME))),function(x){
		return(mean(subset(FRx,paste(SITE,TIME)==x)))
	}))
})	

FRsbl<-sapply(COEFsbl[,1],function(x){
	FRx<-as.numeric(seqtab.nochim[x,vSAMP])/NS
	return(sapply(sort(unique(paste(SITE,TIME))),function(x){
		return(mean(subset(FRx,paste(SITE,TIME)==x)))
	}))
})

FRmsh<-sapply(COEFmsh[,1],function(x){
	FRx<-as.numeric(seqtab.nochim[x,vSAMP])/NS
	return(sapply(sort(unique(paste(SITE,TIME))),function(x){
		return(mean(subset(FRx,paste(SITE,TIME)==x)))
	}))
})

iNONE<-setdiff(c(1:length(seqtab.nochim[,1])),
	c(COEFint[,1],COEFsbl[,1],COEFmsh[,1]))

FRnone<-sapply(iNONE,function(x){
	FRx<-as.numeric(seqtab.nochim[x,vSAMP])/NS
	return(sapply(sort(unique(paste(SITE,TIME))),function(x){
		return(mean(subset(FRx,paste(SITE,TIME)==x)))
	}))
})


###Summary of ASV relative abundance in fonction of TIME

pdf("ASV-time-frequency.pdf",width=5,height=5)

par(mar=c(2,4,2,1),mfrow=c(2,2),bty="n")

barplot(t(FRint)[order(COLt[COEFint[,1]]),],
	col=COLt[COEFint[,1]][order(COLt[COEFint[,1]])],
	main="Increase in MSH and SBL",xaxt="n",las=1,border=NA,
	ylab="Relative abundance",lwd=0.2)

barplot(t(FRsbl)[order(COLt[COEFsbl[,1]]),],
	col=COLt[COEFsbl[,1]][order(COLt[COEFsbl[,1]])],
	main="Increase in SBL only",xaxt="n",las=1,border= NA,
	ylab="Relative abundance",lwd=0.2)
	
barplot(t(FRmsh)[order(COLt[COEFmsh[,1]]),],
	col=COLt[COEFmsh[,1]][order(COLt[COEFmsh[,1]])],
	main="Increase in MSH only",xaxt="n",las=1,border= NA,
	ylab="Relative abundance",lwd=0.2)
	
barplot(t(FRnone)[order(COLt[iNONE]),],
	col=COLt[iNONE][order(COLt[iNONE])],
	main="no increase",xaxt="n",las=1,border= NA,
	ylab="Relative abundance",lwd=0.2)
dev.off()	
	
##############
###############
#############
