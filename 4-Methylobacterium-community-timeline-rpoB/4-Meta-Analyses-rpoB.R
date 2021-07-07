################################################
###################################################
###

#Work directories

pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-1"

pD2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-2"

pD3<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-3"

pD4<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-4"

pO<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"

###Phylogenetic tree
setwd(pD3)	
#minimum bootstrapp value
bx=0
#strain (tips of the tree) names	
namex <-as.vector(read.table("tree-names.txt")[,1])
	
#Reference sequence informations (names in second column)
setwd(pO)
CHAR<-read.table("sequence_characteristics.txt",header=T)
	
## Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs
setwd(pD1)

####Load Metadata for each sample
METAsamp<-read.table("METADATA-wo-NEGcontrols.txt")
#Indexes for positive controls
vMETH<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]=="METH")
#Indexes for samples
vSAMP<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]!="METH")

##EXCLUDE SAMPLES FROM SEQUENCING BATCH rpoB-1-step: WAS DONE WITH DIFFERENT EXTRACTION KIT AND SHOW A HIGHER PROPORTION OF UNASSIGNED ASV, CAULOBACTERALES AND BRADYRHIZOBIACEAE (SEE Summary_rpoB-tax-per-Factor_TF=250_TR=200_TL=20_NS=3000_Prare=0.005-wo-neg-controls-WITH-FIRST-BATCH.pdf)
vOUT<-subset(c(1:length(METAsamp[,1])),METAsamp[,19]=="rpoB-1-step")

#Remove outliers in PCA
EX<-c("T161.PL2.rpoB","T187"  ,        "T223.rpoB.PL3","T199.rpoB.PL3" ,"T152.PL2.rpoB")

sampx<-row.names(METAsamp)
sampx<-as.vector(sapply(sampx,function(x){
	return(length(intersect(x,EX)))
}))
sampx<-subset(c(1:length(METAsamp[,1])), sampx==1)

vSAMP<-setdiff(vSAMP,c(vOUT, sampx)	)

METAsamp<-METAsamp[vSAMP,]
	
##Load Rarefied SeqTAB of absolute abundances

TFx=250
TRx=200
TLx=20
NS=3000
Prare=0.005

seqtab.nochim<-read.table(paste("RAREFIED_SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

ASVs<-rownames(seqtab.nochim)

###Only keep abundances for samples
FRall<-t(seqtab.nochim[,vSAMP])
colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="-")

###Apply Hellinger transformation on relative abundance to account for rare ASVs
library(vegan)

FRallH <-as.data.frame(decostand(FRall ,method="hellinger", MARGIN=1))

colnames(FRallH)<-paste("ASV", c(1:length(FRallH[1,])),sep="-")

rm(seqtab.nochim)

##Load Taxonomic assignation for each ASV

setwd(pD2)

taxaALL <-read.table(paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

##Only keep Methylobacteriaceae

GENx<-as.vector(taxaALL[,7])
FAMx<-as.vector(taxaALL[,6])

rm(taxaALL)

#Methylobacteriaceae ASVs
ASVm<-subset(c(1:length(GENx)), FAMx =="Methylobacteriaceae")

#Genus name for these ASVs
GENx <-as.vector(GENx[ASVm])

#ASV names (genus names 4 first digit + ASV number)

ASVn<-paste(as.vector(sapply(GENx,function(x){
	return(paste(unlist(strsplit(x,split=""))[1:4],collapse=""))
})),"ASV", ASVm,sep="-")

ASVi<-cbind(NA,NA,"ASV",ASVn,GENx,"sp.","PHYLLOSPHERE",NA)
colnames(ASVi)<-colnames(CHAR)

CHAR<-rbind(CHAR,ASVi)

#tree information
setwd(pD3)
coordx<-read.table(paste("Tree-coordinates.txt",sep=""),
	header=T)
SCALE <-read.table(paste("SCALE.txt",sep=""),header=T)
PSn <-read.table(paste("nodes-attributes_MINBOOT=",
	bx,".txt",sep=""),header=T)

#Tip (strain) attributes
Tx<-c(1:length(namex))

#Reformat tree names (for references only)
#namer<-as.vector(sapply(namex,function(x){
#	nx<-unlist(strsplit(x,split="_"))
#	nnx<-paste(nx[3:(length(nx)-1)],
#		collapse="_")
#	return(ifelse(length(nx)==1,x, nnx))	
#}))
	
type<-t(sapply(Tx,function(x){
	#strain names
	stx=namex[x]
	
	CHARx<-as.vector(t(subset(CHAR,CHAR[,4]== stx )))
	typex<-ifelse(CHARx[3]=="ASV","ASV","ref")	
	return(c(typex ,CHARx[c(4:8)]))
}))
	
colnames(type)<-c("type","seq","genus",
	"species","env","clade")
	
#For each ASV, determine its consensus clade and supporting bootstrapp value
	
ASVi<-subset(c(1:length(type[,1])),type[,1]=="ASV")

	
LEVELx<-read.table(paste("Levels_MINBOOT=",
	bx,".txt",sep=""))
	
CLo<-cbind(NA,1,"un",0)
	
CLA<-cbind(type[ASVi,2],
	t(sapply(ASVi,function(x){
	print(x)
	Lx<-as.numeric(LEVELx[x,])
	Lvx<-setdiff(subset(Lx,Lx!=0),x)
	CLx<-t(sapply(Lvx,function(x){
		Ly<-subset(c(1:length(Lx)),Lx==x)
		Lsx<-subset(LEVELx[,1],LEVELx[, Ly]==x)
		Cx<-type[Lsx,6]
		Cx<-unique(setdiff(Cx,NA))
		return(c(x,length(Cx),Cx[1],
			subset(PSn[,8], PSn[,2]==x)))
	}))
	CLx<-rbind( CLx,CLo)
	CLx<-subset(CLx, CLx[,2]=="1")
	CLx<-subset(CLx, as.numeric(CLx[,4])==
		max(as.numeric(CLx[,4])))
	return(c(x,CLx[1,]))
})))
	
colnames(CLA)<-c("Name","Tree-order",
	"NodeAss","MatchAss","ClaAss","BootAss")

setwd(pD4)	
write.table(CLA,"ASV-Clade-assignment-WOout.txt")
	
#Only Methylobacterium ASV 
ASVo <-matrix(unlist(
	strsplit(CLA[,1],split="-")),ncol=3,byrow=T)

CLAo<-subset(CLA,ASVo[,1]=="Meth")
ASVo<-subset(as.numeric(ASVo[,3]),ASVo[,1]=="Meth")

	
#Sort samples by site, date, plot and tree
SAMPo<-order(paste(METAsamp[,3],
	METAsamp[,8],METAsamp[,4],METAsamp[,5]))
	
#Sort data accordingly
METAsamp<-METAsamp[SAMPo,]
FRall<-FRall[SAMPo, ASVo]
FRallH<-FRallH[SAMPo, ASVo]
ASVs<-ASVs[ASVo]

#Relative abundance per sample
FRallR<-t(sapply(c(1:length(FRall[,1])),function(x){
	return(FRall[x,]/sum(FRall[x,]))
}))
colnames(FRallR)<-colnames(FRall)
rownames(FRallR)<-rownames(FRall)


####Extract and format factors from sample metadata 
#SITE of sampling
	SITE<-as.vector(METAsamp[,3])
	SITEu<-sort(unique(SITE)) 
#SUBSITE of sampling
	SS<-paste(SITE,as.vector(METAsamp[,4]),sep="-")
	SSu<-sort(unique(SS))
#HOST TREE SPECIES	
	SPE<-as.vector(METAsamp[,6])
	#Replace tree species names by abbreviations
	SPE <-ifelse(SPE=="Acer_saccharum","ACSA",
		ifelse(SPE=="Fagus_grandifolia","FAGR",
		ifelse(SPE=="Ostria_virginiana","OSVI",
		ifelse(SPE=="Abies_balsamea","ABBA",
		ifelse(SPE=="Acer_rubrum","ACRU",
		ifelse(SPE=="Acer_pennsylcanicum","ACPE",
		ifelse(SPE=="Querbus_rubra","QURU",SPE)))))))
# TREE ID		
	ID<-paste(SS,as.vector(METAsamp[,5]),sep="-")	
#Sampling TIME
	TIME<-as.numeric(METAsamp[,8])
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
	EXT<-as.vector(METAsamp[,16])
	EXTu<-sort(unique(EXT))
#PCR batch 
	PCR<-as.vector(METAsamp[,18])
	PCRu<-sort(unique(PCR))
#Sequencing
	MISEQ <-as.vector(METAsamp[,19])
	MISEQu<-sort(unique(MISEQ))		

	#RELATIVE abundance of Clades per sample
	Fcla<-sapply(sort(unique(CLAo[,5])),function(x){	
		FRx<-subset(t(FRallR),CLAo[,5]==x)		
		return(apply(FRx,2,sum))
	})

	#Average relative abundance of Clades and number of ASVs
	FclaSum<-t(sapply(sort(unique(CLAo[,5])),function(x){	
		FRx<-subset(t(FRallR),CLAo[,5]==x)		
		return(c(length(FRx[,1]),mean(apply(FRx,2,sum))))
	}))
	

	#clade, color, font, ring (segment),cex
	COLcla<-rbind(
		cbind("Enterovirga",rgb(0.8,0.8,0.8),3,2,0.5),
		cbind("Microvirga",rgb(0.6,0.6,0.6),3,2,0.5),
		cbind("B","black",2,7,0.8),
		cbind("C",rgb(0.4,0.4,0.4),3,3,0.5),
		cbind("A9","red",2,7,0.8),
		cbind("A1","blue",2,7,0.8),
		cbind("A3","blue4",3,3,0.5),
		cbind("A6","orange",2,7,0.8),
		cbind("A2","cyan2",2,7,0.8),
		cbind("A4","brown",3,3,0.5),
		cbind("A5a","yellow3",3,3,0.5),
		cbind("A5b","yellow3",3,3,0.5),		
		cbind("A7","pink",3,3,0.5),
		cbind("A10","green3",2,7,0.8),
		cbind("A8","purple",3,3,0.5))
	
	COLclax<- sapply(sort(unique(CLAo[,5])),function(x){
		return(c(subset(COLcla[,2], COLcla[,1]==x),NA)[1])
	})
	
	par(mar=c(4,4,1,1),bty="n")
	FRsite<-sapply(SITEu,function(x){
		
		Fx<-subset(Fcla,SITE==x)
		return(apply(Fx,2,mean))
	})
	
	barplot(FRsite,las=2,col= COLclax,border=NA,las=1,cex.axis=0.8,ylab="Methylobacteriaceae ASV relative abundance")
	
	setwd(pD4)
	
	FRsum<-cbind(FclaSum, FRsite)
	colnames(FRsum)<-c("ASVs","Fr","MSH","SBL")
	
	write.table(FRsum,"Summary-Methylobacterium-ASV-WOout.txt")

###STEP 1: GLOBAL PERMANOVA ON METHYLOBACTERIUM DIVERSITY

#METAdata formated for permanova
FREQ.ENV<-as.data.frame(cbind(SITE,SS,SPE,TIME,EXT,PCR,MISEQ))

np=10000 #number of permutations
nc= 2 #number of cores to use

#Permanova (All samples)
print(date())	
ADONISallmiseq<-adonis(FRallH ~ SITE*SPE*TIME*SS, 
	data= FREQ.ENV, permutations=np,
	strata= MISEQ,method="bray"	,parallel=nc)		
print(date())

#Permanova (MSH samples)
FREQ.ENVm<-subset(FREQ.ENV,SITE=="MSH")
FRallHm<-subset(FRallH,SITE=="MSH")

print(date())	
ADONISallmiseqm<-adonis(FRallHm ~ SPE*TIME*SS, 
	data= FREQ.ENVm, permutations=np,
	strata= FREQ.ENVm$MISEQ,method="bray"	,parallel=nc)		
print(date())

#Permanova (SBL samples)
FREQ.ENVs<-subset(FREQ.ENV,SITE=="SBL")
FRallHs<-subset(FRallH,SITE=="SBL")

print(date())	
ADONISallmiseqs<-adonis(FRallHs ~ SPE*TIME*SS, 
	data= FREQ.ENVs, permutations=np,
	strata= FREQ.ENVs$MISEQ,method="bray"	,parallel=nc)		
print(date())

write.table(ADONISallmiseq$aov,"Permanova-Methylobacterium-All-samples-WOout.txt")
write.table(ADONISallmiseqm$aov,"Permanova-Methylobacterium-MSH-samples-WOout.txt")
write.table(ADONISallmiseqs$aov,"Permanova-Methylobacterium-SBL-samples-WOout.txt")

###Color codes for factors
COLfac<-rbind(
c("ABBA",rgb(0.35,0.35,0.35)),#Host tree species
c("ACRU","red3"),
c("ACSA","purple2"), 
c("OSVI","blue"), 
c("QURU","green"), 
c("FAGR","violet"), 
c("ACPE","orange3"),
c("SBL","orange"), #Sampling site
c("MSH","green3"), 
t(sapply(TIMEu,function(x){	#Sampling time
	z<-(x-min(TIMEu))/(max(TIMEu)-min(TIMEu))
	return(c(x,rgb(1-z,1-z,0)))
})),
t(sapply(c(1:length(EXTu)),function(x){	#Extraction batch (blue scale)
	z<-x/length(EXTu)
	return(c(EXTu[x],rgb(0,0,z)))
})),
t(sapply(c(1:length(PCRu)),function(x){	#PCR batch (red scale)
	z<-x/length(PCRu)
	return(c(PCRu[x],rgb(z,0,0)))
})),
t(sapply(c(1:length(MISEQu)),function(x){	#Sequencing batch (green scale)
	z<-x/length(MISEQu)
	return(c(MISEQu[x],rgb(0,z,0)))
})),
c("METH","red"), #Positive controls
c("0","red"))

########Spatial and temporal autocorrelation

##STEP 1: Calculate pairwise geographic distance among subsites and trees 

#Map within plots (relative position of trees within plots)
DISTss<-read.table("Map-per-subsite.txt",header=T)

#Calculate pairwise geographic distances between trees within subsites (euclidian distances)

DISTt<-matrix(unlist(sapply(as.vector(unique(DISTss[,2])),function(x){
	DISTx<-subset(DISTss,DISTss[,2]==x)
	return(t(matrix(unlist(sapply(c(1:(length(DISTx[,1])-1)),function(x){
		X=x
		return(sapply(c((X+1):(length(DISTx[,1]))),function(x){
			Y=x		
			return(c(
				t(DISTx[X,c(1:3)]),
				t(DISTx[Y,c(1:3)]),
				sqrt((DISTx[X,4]-DISTx[Y,4])^2
					+(DISTx[X,5]-DISTx[Y,5])^2)))
		}))
	})),ncol=7,byrow=T)))	
})),ncol=7,byrow=T)

##Calculate geographic distance among subsites from GPS coordinates
#Calculate a Distance Matrix for Geographic Points Using R (HOMEMADE SCRIPT)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
   # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
   # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
   if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
   if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
   else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
   else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
   m[tri] <- t(m)[tri]
   return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
   # Returns a matrix (M) of distances between geographic points.
   # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
   # (df.geopoints$lat[j], df.geopoints$lon[j]).
   # The row and column names are given by df.geopoints$name.
   GeoDistanceInMetres <- function(g1, g2){
      # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
      # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
      # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
      # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
      # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
      DistM <- function(g1, g2){
         require("Imap")
         return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
      }
      return(mapply(DistM, g1, g2))
   }
   n.geopoints <- nrow(df.geopoints)
   # The index column is used to ensure we only do calculations for the upper triangle of points
   df.geopoints$index <- 1:n.geopoints
   # Create a list of lists
   list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
   # Get a matrix of distances (in metres)
   mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
   # Set the row and column names
   rownames(mat.distances) <- df.geopoints$name
   colnames(mat.distances) <- df.geopoints$name
   return(mat.distances)
}

#GPS coordinates for each subsite
COORDS <- data.frame(name = c("MSH-06","MSH-01","SBL-03","MSH-03","SBL-02","SBL-04","SBL-01","MSH-04","MSH-H0","MSH-L0","MSH-02","MSH-05"),
                         lat  = c(45.53875,45.53964,45.99458,45.54111,45.99172,45.99268,45.98925,45.5418,45.5427,45.54146,45.541025,45.5433066666667),
                      lon  = c(73.154965,73.15775,73.98889,73.16395,73.99871,73.99133,74.0032,73.166335,73.156,73.16293,73.161085,73.1689466666667))

DISTS<-GeoDistanceInMetresMatrix(COORDS) 

DISTS<-matrix(sapply(c(1:length(DISTS[,1])),function(x){
	X=x
	return(sapply(c(1:length(DISTS[1,])),function(x){
		Sx<-rownames(DISTS)[X]
		Sy<-colnames(DISTS)[x]
		return(t(cbind(paste(sort(c(Sx,Sy)),collapse="_"),DISTS[X,x])))
	}))
}),ncol=2,byrow=T)

DISTS<-t(sapply(unique(DISTS[,1]),function(x){
	return(c(strsplit(x,split="_")[[1]],unique(subset(DISTS[,2], DISTS[,1]==x))))
}))

#Remove 2017 subsites
DISTS<-subset(DISTS,DISTS[,1]!="MSH-H0")
DISTS<-subset(DISTS,DISTS[,1]!="MSH-L0")
DISTS<-subset(DISTS,DISTS[,2]!="MSH-H0")
DISTS<-subset(DISTS,DISTS[,2]!="MSH-L0")

####Merge Distance among subsites and distance among trees within subsite: calculate all possible pairwise distances among trees

#All possible trees
TREE<-sort(unique(c(DISTt[,3], DISTt[,6])))

#Reformat the names so they match tree names in METAdata
TREEf<-as.vector(sapply(TREE,function(x){
	Tx<-unlist(strsplit(x,split="-"))
	return(paste(c(Tx[1],
		as.numeric(Tx[2:3])),collapse="-"))
}))

# calculate Pairwise geographic distance among all trees
DISTtree<-matrix(unlist(sapply(c(1:(length(TREE)-1)),function(x){
	X=x
	treeX= TREE[X]
	tx=TREEf[X]
	return(sapply(c((X+1):length(TREE)),function(x){
		Y=x
		treeY= TREE[Y]	
		ty=TREEf[Y]
		SSx<-unique(
			c(subset(DISTt[,2], DISTt[,3]== treeX),
			subset(DISTt[,5], DISTt[,6]== treeX)))
		SSy<-unique(
			c(subset(DISTt[,2], DISTt[,3]== treeY),
			subset(DISTt[,5], DISTt[,6]== treeY)))
		Dxy<-rbind(subset(DISTS, DISTS[,1]==SSx),
			subset(DISTS, DISTS[,2]==SSx))
		Dxy<-as.numeric(rbind(subset(Dxy, Dxy[,1]==SSy),
			subset(Dxy, Dxy[,2]==SSy))[,3])
		Txy<-rbind(subset(DISTt, DISTt[,3]== treeX),
			subset(DISTt, DISTt[,6]== treeX))
		Txy<-c(as.numeric(rbind(subset(Txy, Txy[,3]== treeY),
			subset(Txy, Txy[,6]== treeY))[,7]),NA)[1]	
		return(c(paste(sort(c(tx,ty)),collapse="_"),
			ifelse(is.na(Txy)==T, Dxy,Txy)))
	}))	
})),ncol=2,byrow=T)

#Add paiwise comparison between two samples from the same tree ( 0 m as a distance)
DISTtree <-rbind(DISTtree ,t(sapply(c(1:(length(TREE))),function(x){
	X=x
	treeX= TREE[X]
	tx=TREEf[X]
	SSx<-unique(c(subset(DISTt[,2], DISTt[,3]== treeX),
			subset(DISTt[,5], DISTt[,6]== treeX)))
	return(c(paste(c(tx,tx),collapse="_"),0))	
})))	


SITEu<-sort(unique(SITE))
CLAu<-setdiff(sort(unique(CLAo[,5])),"un")

#Step 2: Calculate pairwise geographic distance, time and Bray-Curtis dissimilarity
DISTtime<-matrix(unlist(sapply(SITEu,function(x){
	#print(x)
	SITEx=x
	Sx<-subset(c(1:length(SITE)), SITE == SITEx) #samples from this site
	zx<-t(FRallH[Sx,]) #corrected ASV relative abundance for this site
	DISTx<-matrix(unlist(sapply(c(1:(length(Sx)-1)),function(x){
		X=x
		SS1<-SS[Sx[X]] #Subsite for this sample
		T1= TIME[Sx[X]]#Sampling time for this sample
		S1<-as.numeric(zx[,X])#corrected ASV relative abundance for this sample
		SPE1<-SPE[Sx[X]] #host tree species for this sample
		TREE1<-ID[Sx[X]] #tree identity for this sample
		MISEQ1<-MISEQ[Sx[X]] #sequencing batch for this sample
		return(sapply(c((X+1):(length(Sx))),
			function(x){
			T2= TIME[Sx[x]]
			SS2<-SS[Sx[x]]
			S2<-as.numeric(zx[,x])
			SPE2<-SPE[Sx[x]]
			TREE2<-ID[Sx[x]]
			MISEQ2<-MISEQ[Sx[x]]	
			D12<-subset(DISTtree[,2],DISTtree[,1]
				==paste(sort(c(TREE1,TREE2)),
				collapse="_")) #Spatial distance
			B12=vegdist(t(cbind(S1,S2)),
				method="bray") #GLOBAL Bray dissimilarity
			#Bray dissimilarity for each clade	
			B12c<-sapply(CLAu,function(x){
				return(vegdist(t(
					subset(cbind(S1,S2), 
					CLAo[,5]==x)),
					method="bray"))
			})			
			return(c(SS1,SS2,TREE1,TREE2,
				SPE1,SPE2,MISEQ1,MISEQ2,
				D12,T1,T2,B12, B12c))		
		}))
	})),ncol=12+length(CLAu),byrow=T)
	return(t(cbind(SITEx,DISTx)))
})),ncol=13+length(CLAu),byrow=T)

colnames(DISTtime)<-c("Site","SS1","SS2","TREE1","TREE2","SPE1","SPE2","MISEQ1","MISEQ2","GEO-DIST","DAY1","DAY2","BD-All",paste("BD", CLAu,sep="-"))

####Models to test (each within sites):

# PC1 pairwise community similarity between trees within dates in function of geographic distance: 2 site, 4 dates (8 models)
DIST1<-subset(DISTtime,DISTtime[,11]==DISTtime[,12])

# PC2  pairwise community similarity in function of pairwise time : within sites (2 models), plots (2 models) and trees (2 models)
	
TIMEleg<-cbind(TIMEu ,TIMEru,
	sapply(TIMEu,function(x){
	return(c(subset(COLfac[,2],
		COLfac[,1]==x)))
}))

#For each clade, count the number of pairwise comparisons for which Bray-curtis distances are meaningless
nBD<-t(sapply(c(1:length(CLAu)),function(x){
	CATx= CLAu[x]
	BDi=subset(c(1:length(colnames(DIST1))),
		colnames(DIST1)==paste("BD",CATx,sep="-"))
	return(sapply(SITEu,function(x){
		D1<-subset(DIST1, DIST1[,1]== x)
		Bx<-as.numeric(D1[, BDi]) #Bray-Curtis 
		return(
			(length(subset(Bx,is.nan(Bx)==T))+
			length(subset(Bx,Bx==1)))/length(Bx))
	}))
}))
rownames(nBD)<-CLAu

#Category to test (remove clades with more than 99% of meaningless BD values in at least one site)
CATt<-c("All",subset(CLAu,apply(nBD,1,max)<0.99))

setwd(pD4)

pdf("Methylo-AUTOCORRELATION-WOout.pdf",height=8,width=8)

par(mfcol=c(5,6),mar=c(4,4,2,1),bty="n")


par(mar=c(0,0,0,0))
		plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")


#Plot legends
sapply(SITEu,function(x){
	SITEx=x
	Dx<-sort(unique(as.numeric(subset(DIST1[,11], DIST1[,1]== SITEx))))
	
	par(mar=c(4,0,1,0))
	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
		xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,ifelse(SITEx =="MSH","a","c"),cex=1.5,font=2)

	TIMElegx<- t(sapply(Dx,function(x){
		return(subset(TIMEleg,as.numeric(TIMEleg[,1])==x))
	}))
	legend("center",
	legend= TIMElegx[,2],cex=0.8, fill= TIMElegx[,3],bty="n",title=SITEx)

	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
		xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,ifelse(SITEx =="MSH","b","d"),cex=1.5,font=2)

	colx=subset(COLfac[,2], COLfac[,1]==SITEx)

	legend("center",
	legend= c("within sites","within plots","within trees"),cex=0.8, fill= c(
	rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255,1),
	rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400,1),
	rgb(col2rgb(colx)[1]/2000,
			col2rgb(colx)[2]/2000,
			col2rgb(colx)[3]/2000,1)				
	),bty="n",title=SITEx)	

})

#Plot real data and regressions

sapply(c(1:length(CATt)),function(x){

	CATx= CATt[x]
	BDi=subset(c(1:length(colnames(DIST1))),
		colnames(DIST1)==paste("BD",CATx,sep="-"))

	par(mar=c(0,3,0,0))
	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
		xaxt="n",yaxt="n",xlab="",ylab="")
		text(0.5,0.1, CATx,cex=0.8,font=2)


	SUMall<-sapply(SITEu,function(x){
	
		par(mar=c(4,4,2,1))
		print(c("SITE:",x))
		SITEx=x
	
		print("pairwise community similarity between trees within dates in function of geographic distance; all pairwise comparisons")
		D1<-subset(DIST1, DIST1[,1]== SITEx)
	
		Gx<-as.numeric(D1[,10]) #Geographic distance
		Bx<-as.numeric(D1[, BDi]) #Bray-Curtis distance
		Dx<-as.numeric(D1[,11]) #Sampling date
		Dmx<-Dx-min(Dx)
		M1<-lm(Bx ~Gx* Dmx)	
	
		Ax<-anova(M1)
		Sx<-summary(M1)
		Sx<-cbind(Sx$coefficients[,1],Sx$coefficients[,2],Sx$coefficients[,4])
		colnames(Sx)<-c("Est.","Est.SD","Est.p.value")
	#part of explained variance and significance (anova)
		M1<-cbind(round(Ax[,2]/sum(Ax[,2]),4),
			round(Ax[,5],4))
		row.names(M1)<-row.names(Ax)
		colnames(M1)<-c("Anova.R2","Anova.p.val.")
		M1<-M1[c(4,1,2,3),]
	
		M1<-cbind(M1,Sx)
		rownames(M1)<-c("Res./Int.",rownames(M1)[2:4])
		
		print("Detail per time point: linear model and ANOVA")
	
		Gmin=1
		Gmin2<- min(subset(Gx,Gx> Gmin))
	
		ld <- seq(Gmin,max(Gx), length.out=100)
		#COLt<-as.numeric(PC1x[,2])/max(TIME)
	
		la<-c(log10(Gmin),seq(log10(Gmin2),
			log10(max(Gx)), length.out=6))
		la<-round(10^la)

	#par(mar=c(0,0,0,0))
	#	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	#xaxt="n",yaxt="n",xlab="",ylab="")
	#text(0.05,0.95,ifelse(SITEx =="MSH","e","f"),cex=1.5,font=2)
	#par(new=T,mar=c(4,4,1,1))

	plot(10000,10000,
		log="x",xlab="pDistance (m)",
		ylab="BD",
		las=1,cex.axis=0.8,cex.lab=0.8,
		xlim=c(min(ld),max(ld)),
		ylim=c(min(as.numeric(DISTtime[,13])),
		max(as.numeric(DISTtime[,13]))),
		main="",
		cex.main=0.8,xaxt="n",font.lab=2)	
		
	axis(1, la,c(0,la[2:length(la)]),cex.axis=0.8)
	axis(1, c(Gmin,Gmin2),c(NA,NA),cex.axis=0.8,col="white",lty=3, lwd.ticks=0,lwd=2,lend=2)
	
	SUM1<-t(sapply(sort(unique(Dx)),function(x){
		colxx=subset(TIMEleg[,3], as.numeric(TIMEleg[,1])==x)
		colxx3=rgb(col2rgb(colxx)[1]/255,
			col2rgb(colxx)[2]/255,
			col2rgb(colxx)[3]/255,0.3)
		colxx6=rgb(col2rgb(colxx)[1]/255,
			col2rgb(colxx)[2]/255,
			col2rgb(colxx)[3]/255,0.1)	
		Gxx<-subset(Gx,Dx==x)
		Bxx<-subset(Bx,Dx==x)
		points(ifelse(Gxx==0,Gmin,Gxx),Bxx,pch=19,cex=0.3,
			col= colxx6)
		Mx<-lm(Bxx ~Gxx)
		px<-summary(Mx)$coefficient[2,4]
		COEF<-summary(Mx)$coefficient[,1]
		Axx<-c(round(anova(Mx)[,2]/
			sum(anova(Mx)[,2]),4),
			round(anova(Mx)[1,5],4))
		SD<-summary(Mx)$coefficient[,2]	
		#Intervals
		INTx<-sapply(ld,function(x){
			Lx=x
			VAR1x = (COEF[2]+SD[2])* Lx + (COEF[1]+SD[1])
			VAR2x = (COEF[2]-SD[2])* Lx + (COEF[1]-SD[1])
			VAR3x = (COEF[2]-SD[2])* Lx + (COEF[1]+SD[1])
			VAR4x = (COEF[2]+SD[2])* Lx + (COEF[1]-SD[1])
			return(c(VAR1x, VAR2x, VAR3x, VAR4x))
		})
		INTx<-cbind(apply(INTx,2,max),apply(INTx,2,min))	
		VARx = COEF[2]* ld + COEF[1]	
		polygon(c(ld,ld[length(ld):1]),
			c(INTx[,1],INTx[c(length(ld):1),2]),
			border="NA",col=colxx3)	
		lines(ld, VARx,
			col= colxx,
			lwd=1,lty=ifelse(px<=0.05,1,2),lend=2)		
		return(c(COEF,SD,px,Axx))
	}))
	TIMElegx<- t(sapply(sort(unique(Dx)),function(x){
		return(subset(TIMEleg,as.numeric(TIMEleg[,1])==x))
	}))
	
	colnames(SUM1)<-c("lm.Int.","lm.est.","lm.SD.Int.","lm.SD.est.","lm.p.val.est.","Anova.R2.est.","Anova.R2.res.","Anova.p.val.")
	rownames(SUM1)<-paste("Spatial autocorrelation; D", TIMElegx[,2],sep="=")
	
	SUM1 <-sapply(c(1:length(SUM1[,1])),function(x){	
		SUMx<-cbind(SUM1[x,c(7:6)],
			c(NA,SUM1[x,8]),
			SUM1[x,c(1:2)],
			SUM1[x,c(3:4)],
			c(NA,SUM1[x,5]))
		rownames(SUMx)<-c("Res./Int.","Gx")
		colnames(SUMx)<-colnames(M1)
		SUMx<-list(SUMx)
		names(SUMx)<-rownames(SUM1)[x]
		return(SUMx)
	})
	
	M1<-list(M1)
	names(M1)<-paste("Spatial autocorrelation; General model",sep=";")
	
	
	print("Pairwise community similarity in function of pairwise time; all pairwise comparisons within site")
	D2<-subset(DISTtime, DISTtime[,1]== SITEx)
	Bx<-as.numeric(D2[, BDi]) 
	Tx<-abs(as.numeric(D2[,11])-as.numeric(D2[,12]))#Pairwise time
	
	Tmin=5
	Tmin2<- min(subset(Tx,Tx>Tmin))
	
	ld <- seq(Tmin,max(Tx), length.out=100)
	#COLt<-as.numeric(PC1x[,2])/max(TIME)
	
	la<-c(log10(Tmin),seq(log10(Tmin2),log10(max(Tx)), length.out=6))
	la<-round(10^la)

	#par(mar=c(0,0,0,0))
	#	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	#xaxt="n",yaxt="n",xlab="",ylab="")
	#text(0.05,0.95,ifelse(SITEx =="MSH","g","h"),cex=1.5,font=2)
	#par(new=T,mar=c(4,4,1,1))

	plot(10000,10000,
		log="x",xlab="p-time (days)",
		ylab="BD",
		las=1,cex.axis=0.8,cex.lab=0.8,
		xlim=c(min(ld),max(ld)),
		ylim=c(min(as.numeric(DISTtime[,13])),
		max(as.numeric(DISTtime[,13]))),
		main="",xaxt="n",font.lab=2)	
		
	axis(1, la,c(0,la[2:length(la)]),cex.axis=0.8)
	axis(1, c(Tmin,Tmin2),c(NA,NA),cex.axis=0.8,col="white",lty=3, lwd.ticks=0,lwd=2,lend=2)	
	
	colx=subset(COLfac[,2], COLfac[,1]==SITEx)

	M2<-lm(Bx ~Tx)
	COEFx<-summary(M2)$coefficient[,1]
	SDx<-summary(M2)$coefficient[,2]
	px<-summary(M2)$coefficient[2,4]
	Ax<-anova(M2)
	Sx<-summary(M2)
	Sx<-cbind(Sx$coefficients[,1],Sx$coefficients[,2],Sx$coefficients[,4])
	colnames(Sx)<-c("Est.","Est.SD","Est.p.value")
	M2 <-cbind(round(Ax[,2]/sum(Ax[,2]),4),
		round(Ax[,5],4))
	row.names(M2)<-row.names(Ax)
	colnames(M2)<-c("Anova.R2","Anova.p.val.")
	
	M2 <-M2[c(2,1),]
	
	M2 <-cbind(M2,Sx)
	rownames(M2)<-c("Res./Int.",rownames(M2)[2])
	
	A= COEFx[2]
	B=COEFx[1]
	a=SDx[2]
	b=SDx[1]
	VARx = A*ld + B
	#Intervals
	INTx<-sapply(ld,function(x){
		Lx=x
		VAR1x = (A+a)* Lx + (B+b)
		VAR2x = (A-a)* Lx + (B-b)
		VAR3x = (A-a)* Lx + (B+b)
		VAR4x = (A+a)* Lx + (B-b)
		return(c(VAR1x, VAR2x, VAR3x, VAR4x))
	})
	INTx<-cbind(apply(INTx,2,max),apply(INTx,2,min))	
	points(ifelse(Tx==0,Tmin,Tx), Bx,pch=19,cex=0.3,col=
		rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255,0.01))
	polygon(c(ld,ld[length(ld):1]),
		c(INTx[,1],INTx[c(length(ld):1),2]),
		border="NA",
		col=rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255,0.3))	
	lines(ld, VARx,
		col= colx,
		lwd=1,lty=ifelse(px<=0.05,1,2),lend=2)	
	
	print("Pairwise community similarity in function of pairwise time; all pairwise comparisons within sub-site")
	D2a<-subset(D2, D2[,2]== D2[,3])
	Bx<-as.numeric(D2a[,BDi])
	Tx<-abs(as.numeric(D2a[,11])-as.numeric(D2a[,12]))
	M2a<-lm(Bx ~Tx)
	COEFx<-summary(M2a)$coefficient[,1]
	SDx<-summary(M2a)$coefficient[,2]
	px<-summary(M2a)$coefficient[2,4]
	Ax<-anova(M2a)
	Sx<-summary(M2a)
	Sx<-cbind(Sx$coefficients[,1],Sx$coefficients[,2],Sx$coefficients[,4])
	colnames(Sx)<-c("Est.","Est.SD","Est.p.value")
	M2a <-cbind(round(Ax[,2]/sum(Ax[,2]),4),
		round(Ax[,5],4))
	row.names(M2a)<-row.names(Ax)
	colnames(M2a)<-c("Anova.R2","Anova.p.val.")
	
	M2a <-M2a[c(2,1),]
	
	M2a <-cbind(M2a,Sx)
	rownames(M2a)<-c("Res./Int.",rownames(M2a)[2])
	
	A= COEFx[2]
	B=COEFx[1]
	a=SDx[2]
	b=SDx[1]
	VARx = A*ld + B
	#Intervals
	INTx<-sapply(ld,function(x){
		Lx=x
		VAR1x = (A+a)* Lx + (B+b)
		VAR2x = (A-a)* Lx + (B-b)
		VAR3x = (A-a)* Lx + (B+b)
		VAR4x = (A+a)* Lx + (B-b)
		return(c(VAR1x, VAR2x, VAR3x, VAR4x))
	})
	INTx<-cbind(apply(INTx,2,max),apply(INTx,2,min))	
	points(ifelse(Tx==0,Tmin,Tx), Bx,pch=19,cex=0.3,col=
		rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400,0.01))
	polygon(c(ld,ld[length(ld):1]),
		c(INTx[,1],INTx[c(length(ld):1),2]),
		border="NA",
		col=rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400,0.3))	
	lines(ld, VARx,
		col= rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400,1),
		lwd=1,lty=ifelse(px<=0.05,1,2),lend=2)	
	
	print("Pairwise community similarity in function of pairwise time; all pairwise comparisons within trees")
	D2b<-subset(D2a, D2a[,4]==D2a[,5])
	Bx<-as.numeric(D2b[,BDi])
	Tx<-abs(as.numeric(D2b[,11])-as.numeric(D2b[,12]))
	M2b<-lm(Bx ~Tx)
	COEFx<-summary(M2b)$coefficient[,1]
	SDx<-summary(M2b)$coefficient[,2]
	px<-summary(M2b)$coefficient[2,4]
	Ax<-anova(M2b)
	Sx<-summary(M2b)
	Sx<-cbind(Sx$coefficients[,1],Sx$coefficients[,2],Sx$coefficients[,4])
	colnames(Sx)<-c("Est.","Est.SD","Est.p.value")
	M2b <-cbind(round(Ax[,2]/sum(Ax[,2]),4),
		round(Ax[,5],4))
	row.names(M2b)<-row.names(Ax)
	colnames(M2b)<-c("Anova.R2","Anova.p.val.")
	
	M2b <-M2b[c(2,1),]
	
	M2b <-cbind(M2b,Sx)
	rownames(M2b)<-c("Res./Int.",rownames(M2b)[2])

	A= COEFx[2]
	B=COEFx[1]
	a=SDx[2]
	b=SDx[1]
	VARx = A*ld + B
	#Intervals
	INTx<-sapply(ld,function(x){
		Lx=x
		VAR1x = (A+a)* Lx + (B+b)
		VAR2x = (A-a)* Lx + (B-b)
		VAR3x = (A-a)* Lx + (B+b)
		VAR4x = (A+a)* Lx + (B-b)
		return(c(VAR1x, VAR2x, VAR3x, VAR4x))
	})
	INTx<-cbind(apply(INTx,2,max),apply(INTx,2,min))	
	points(ifelse(Tx==0,Tmin,Tx), Bx,pch=19,cex=0.3,col=
		rgb(col2rgb(colx)[1]/2000,
			col2rgb(colx)[2]/2000,
			col2rgb(colx)[3]/2000,0.01))
	polygon(c(ld,ld[length(ld):1]),
		c(INTx[,1],INTx[c(length(ld):1),2]),
		border="NA",
		col=rgb(col2rgb(colx)[1]/2000,
			col2rgb(colx)[2]/2000,
			col2rgb(colx)[3]/2000,0.3))	
	lines(ld, VARx,
		col= rgb(col2rgb(colx)[1]/2000,
			col2rgb(colx)[2]/2000,
			col2rgb(colx)[3]/2000,1),
		lwd=1,lty=ifelse(px<=0.05,1,2),lend=2)
		

	M2 <-list(M2)
	names(M2)<-paste("Temporal autocorrelation; within sites",sep=";")
	M2a <-list(M2a)
	names(M2a)<-paste("Temporal autocorrelation; within plots",sep=";")
	M2b <-list(M2b)
	names(M2b)<-paste("Temporal autocorrelation; within trees",sep=";")
	
	return(list(c(M1,SUM1,M2,M2a,M2b)))
	

})

SUMalln<-colnames(SUMall[[1]][[1]])

SUMall<-matrix(unlist(sapply(c(1:length(SUMall)),function(x){
	X=x
	SUMx<-SUMall[[X]]
	Nx=names(SUMx)
	
	SUMv<-matrix(unlist(sapply(c(1:length(SUMx)),function(x){
		SUMxx<-t(SUMx[[x]])
		colnames(SUMxx)<-paste(Nx[x],colnames(SUMxx),sep=";")
		return(SUMxx)
	})),ncol=5,byrow=T)
	SUMn<-unlist(sapply(c(1:length(SUMx)),function(x){
		SUMxx<-t(SUMx[[x]])
		return(paste(SITEu[X],Nx[x],colnames(SUMxx),sep=";"))
	}))
	return(list(t(cbind(SUMn, SUMv))))
	
})),ncol=6,byrow=T)


colnames(SUMall)<-c("Model/factor",SUMalln)

write.table(as.data.frame(SUMall),paste("Methylo-Spatial+Temporal-autocorrelation-Bray-", CATx ,"-ASV-WOout.txt",sep=""),row.names=F)

})

dev.off()


##############################
#####ANOVA PER ASV
library(ade4)
library(vegan)

setwd(pD4)

#Relative abundance after log transformation (for ANOVA per ASV)
FRallL<-sapply(
	c(1:length(FRall[1,])),function(x){
	fx<-as.vector(FRall[,x])
	fx<-fx/sum(fx)
	flx<-log10(fx)
	flx <-ifelse(fx==0,min(subset(flx, flx!="-Inf"))-0.5,flx)
	flx<-flx-min(flx)
	flx<-flx/sum(flx)
	return(flx)
})

FRallL <-ifelse(is.na(FRallL)==T,0, FRallL)


#Estimate part of variance per ASV for all factors and interactions: - Convert plots in numeric value as distance from plots MSH-6 and SBL-1 to test for the gradient effect, use also Time as numeric

#Plot in numeric values
SSd<-round(sapply(SS,function(x){
	DISTx<-rbind(subset(DISTtime, 
		DISTtime[,2]=="MSH-6"),
		subset(DISTtime, 
		DISTtime[,2]=="SBL-4"))
	DISTx<-subset(DISTx, DISTx[,3]==x)
	return(mean(as.numeric(DISTx[,10])))
}),0)

#Time in numeric values
TIMEn<-TIME-min(TIME)

#Factors
FAC<-c("SITE","TIMEn","SSd","SPE","SITE:TIMEn",
"SITE:SSd","TIMEn:SSd","SITE:SPE","TIMEn:SPE","SSd:SPE","SITE:TIMEn:SSd","SITE:TIMEn:SPE","SITE:SSd:SPE","TIMEn:SSd:SPE","SITE:TIMEn:SSd:SPE","Residuals")


ANOVA<-t(sapply(c(1:length(FRallL[1,])),
	function(x){
	
	fx<-FRallL[,x]	
	
	Mt<-lm(fx ~ SITE*TIMEn*SSd*SPE)	
	At<-anova(Mt)
	
	At <-t(sapply(FAC,function(x){
		return(as.numeric(subset(At,
			rownames(At)==x)[,c(2,5)]))
	}))
	
	return(c(At[,1]/sum(At[,1]),At[,2]))
}))


colnames(ANOVA)<-c(paste("VAR",FAC,sep="."),
	paste("pval",FAC,sep="."))

#Calculate adjusted pvalues

pADJ<-sapply(c(1:length(FAC)),function(x){
	pv<-as.vector(ANOVA[,c(length(FAC)+x)])
	return(p.adjust(pv,method="bonferroni"))
	
})

colnames(pADJ)<-paste("padj",FAC,sep=".")

write.table(cbind(CLAo, ANOVA, pADJ),"ANOVA-Association-per-Meth-ASV-WOout.txt",row.names=F)

###ANOVA per site (forest)

#Factors
FACs<-c("tn","ssd","spe","tn:ssd","tn:spe","ssd:spe",    "tn:ssd:spe","Residuals" )

sapply(c("MSH","SBL"),function(x){
	X=x
	tn<-subset(TIMEn,SITE==X)
	ssd<-subset(SSd,SITE==X)
	spe<-subset(SPE,SITE==X)
	anova<-t(sapply(c(1:length(FRallL[1,])),
		function(x){
		fx<-subset(FRallL[,x],SITE==X)	
		Mt<-lm(fx ~ tn*ssd*spe)	
		At<-anova(Mt)
		At <-t(sapply(FACs,function(x){
			return(as.numeric(subset(At,
			rownames(At)==x)[,c(2,5)]))
		}))
		return(c(At[,1]/sum(At[,1]),At[,2]))
	}))
	colnames(anova)<-
		c(paste("VAR",FACs,sep="."),
		paste("pval",FACs,sep="."))
	#Calculate adjusted pvalues
	padj<-sapply(c(1:length(FACs)),function(x){
		pv<-as.vector(anova[,
			c(length(FACs)+x)])
		return(p.adjust(pv,
			method="bonferroni"))
	})
	colnames(padj)<-paste("padj",FACs,sep=".")

	write.table(cbind(CLAo, anova, padj),
		paste(X,
		"ANOVA-Association-per-Meth-ASV-WOout.txt",
		sep="_"),row.names=F)	
})


#Average Relative abundance per ASV per factor
SUMasv<-cbind(sapply(sort(unique(SITE)),function(x){
	fx<-subset(FRall,SITE==x)
	fx<-t(sapply(c(1:length(fx[,1])),function(x){
		return(fx[x,]/sum(fx[x,]))
	}))
	return(apply(fx,2,mean))
}),
sapply(sort(unique(SS)),function(x){
	fx<-subset(FRall,SS==x)
	fx<-t(sapply(c(1:length(fx[,1])),function(x){
		return(fx[x,]/sum(fx[x,]))
	}))
	return(apply(fx,2,mean))
}),
sapply(sort(unique(paste(SITE,TIME))),function(x){
	fx<-subset(FRall,paste(SITE,TIME)==x)
	fx<-t(sapply(c(1:length(fx[,1])),function(x){
		return(fx[x,]/sum(fx[x,]))
	}))
	return(apply(fx,2,mean))
}),
sapply(sort(unique(paste(SITE,SPE))),function(x){
	fx<-subset(FRall,paste(SITE,SPE)==x)
	fx<-t(sapply(c(1:length(fx[,1])),function(x){
		return(fx[x,]/sum(fx[x,]))
	}))
	return(apply(fx,2,mean))
}))

write.table(cbind(CLAo, SUMasv),"Factor-abundance-per-Meth-ASV-WOout.txt",row.names=F)

#########################
#####SUP figure: PCA per site with emphasize on each main factor

#Contribution to PCA amplification factor
AC=5

#Size of factor contribution tested by anova
ANO=2

#Axis to display
A1=1
A2=2

pdf("ANOVA-per-Forest-Time-Plot-Host-WOout.pdf",width=5,height=7.5)

par(mfrow=c(3,2))


sapply(c("MSH","SBL"),function(x){

X=x

anova<-read.table(paste(X,
		"ANOVA-Association-per-Meth-ASV-WOout.txt",
		sep="_"),header=T)

FRallHs<-subset(FRallH,SITE==X)

ACPhx=dudi.pca(FRallHs, scannf=F,nf=10)
part<-round(100* ACPhx $eig/sum(ACPhx $eig),2)

par(mar=c(0,0,0,0),bty="n")
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,ifelse(X=="MSH","a","b"),cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

x1=min(ACPhx $li[, A1])
x2=max(ACPhx $li[, A1])
y1=min(ACPhx $li[, A2])
y2=max(ACPhx $li[,A2])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=1,
		xlab=paste("Axis",
		A1," (",part[1],"%)",sep=""),
		ylab=paste("Axis",
		A2," (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPhx $li[,c(A1,A2)],groups= subset(TIMEr ,SITE==X),col= rgb(0.5,0.5,0.5),lty=2,lwd=0.5)	

points(ACPhx $li[,c(A1,A2)],pch=19,cex=0.4,col=rgb(0,0,0,0.6))

ordiellipse(ACPhx $li[,c(A1,A2)],groups= subset(TIMEr ,SITE==X),border="grey",alpha=80,draw="polygon",col= rgb(1,1,1,0.01),lty=1,lwd=1,label=F,cex=0.8,font=1,kind="sd",conf=0.2)


sapply(unique(subset(TIMEr ,SITE==X)),function(x){
	coordx<-apply(subset(
		ACPhx $li[,c(A1,A2)],
		subset(TIMEr ,
		SITE==X) ==x),2,mean)
	text(coordx[1], coordx[2],x,cex=0.5,font=2)
})

#ASV contribution colored according to clades
CONT<-ACPhx$co

COLcont<-as.vector(sapply(CLA[,5],function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
}))

#Center of the contribution plot in the main plot
Xc=ifelse(X=="MSH",5,-6)
Yc=ifelse(X=="MSH",7,-5)

text(Xc,Yc+ifelse(X=="MSH",5,2),"ASV contribution & date association",cex=0.5,font=1)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,col= ifelse(anova[,23] <=0.1,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,pch=ifelse(anova[,23] <=0.05,19,21),cex= ifelse(anova[,23] <=0.1,sqrt(abs(anova[,7])/pi)*ANO,0),col= COLcont,lwd=0.75,bg="white")

return(subset(anova[,c(1,5)],anova[,23] <=0.05))

})

sapply(c("MSH","SBL"),function(x){

X=x

anova<-read.table(paste(X,
		"ANOVA-Association-per-Meth-ASV-WOout.txt",
		sep="_"),header=T)

FRallHs<-subset(FRallH,SITE==X)

ACPhx=dudi.pca(FRallHs, scannf=F,nf=10)
part<-round(100* ACPhx $eig/sum(ACPhx $eig),2)

par(mar=c(0,0,0,0),bty="n")
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,ifelse(X=="MSH","c","d"),cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

x1=min(ACPhx $li[, A1])
x2=max(ACPhx $li[, A1])
y1=min(ACPhx $li[, A2])
y2=max(ACPhx $li[,A2])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=1,
		xlab=paste("Axis",
		A1," (",part[1],"%)",sep=""),
		ylab=paste("Axis",
		A2," (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPhx $li[,c(A1,A2)],groups= subset(SS ,SITE==X),col= rgb(0.5,0.5,0.5),lty=2,lwd=0.5)	

points(ACPhx $li[,c(A1,A2)],pch=19,cex=0.4,col=rgb(0,0,0,0.6))

ordiellipse(ACPhx $li[,c(A1,A2)],groups= subset(SS ,SITE==X),border="grey",alpha=80,draw="polygon",col= rgb(1,1,1,0.01),lty=1,lwd=1,label=F,cex=0.8,font=1,kind="sd",conf=0.2)


sapply(unique(subset(SS ,SITE==X)),function(x){
	coordx<-apply(subset(
		ACPhx $li[,c(A1,A2)],
		subset(SS ,
		SITE==X) ==x),2,mean)
	text(coordx[1], coordx[2],x,cex=0.5,font=2)
})

#ASV contribution colored according to clades
CONT<-ACPhx$co

COLcont<-as.vector(sapply(CLA[,5],function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
}))

#Center of the contribution plot in the main plot
Xc=ifelse(X=="MSH",5,-6)
Yc=ifelse(X=="MSH",7,-5)

text(Xc,Yc+ifelse(X=="MSH",5,2),"ASV contribution & plot association",cex=0.5,font=1)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,col= ifelse(anova[,24] <=0.1,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,pch=ifelse(anova[,24] <=0.05,19,21),cex= ifelse(anova[,24] <=0.1,sqrt(abs(anova[,8])/pi)*ANO,0),col= COLcont,lwd=0.75,bg="white")

return(subset(anova[,c(1,5)],anova[,24] <=0.05))

})


sapply(c("MSH","SBL"),function(x){

X=x

anova<-read.table(paste(X,
		"ANOVA-Association-per-Meth-ASV-WOout.txt",
		sep="_"),header=T)

FRallHs<-subset(FRallH,SITE==X)

ACPhx=dudi.pca(FRallHs, scannf=F,nf=10)
part<-round(100* ACPhx $eig/sum(ACPhx $eig),2)

par(mar=c(0,0,0,0),bty="n")
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,ifelse(X=="MSH","e","f"),cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))


x1=min(ACPhx $li[, A1])
x2=max(ACPhx $li[, A1])
y1=min(ACPhx $li[, A2])
y2=max(ACPhx $li[,A2])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=1,
		xlab=paste("Axis",
		A1," (",part[1],"%)",sep=""),
		ylab=paste("Axis",
		A2," (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPhx $li[,c(A1,A2)],groups= subset(SPE ,SITE==X),col= rgb(0.5,0.5,0.5),lty=2,lwd=0.5)	

points(ACPhx $li[,c(A1,A2)],pch=19,cex=0.4,col=rgb(0,0,0,0.6))

ordiellipse(ACPhx $li[,c(A1,A2)],groups= subset(SPE ,SITE==X),border="grey",alpha=80,draw="polygon",col= rgb(1,1,1,0.01),lty=1,lwd=1,label=F,cex=0.8,font=1,kind="sd",conf=0.2)


sapply(unique(subset(SPE ,SITE==X)),function(x){
	coordx<-apply(subset(
		ACPhx $li[,c(A1,A2)],
		subset(SPE ,
		SITE==X) ==x),2,mean)
	text(coordx[1], coordx[2],x,cex=0.5,font=2)
})

#ASV contribution colored according to clades
CONT<-ACPhx$co

COLcont<-as.vector(sapply(CLA[,5],function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
}))

#Center of the contribution plot in the main plot
Xc=ifelse(X=="MSH",5,-6)
Yc=ifelse(X=="MSH",7,-5)

text(Xc,Yc+ifelse(X=="MSH",5,2),"ASV contribution & host association",cex=0.5,font=1)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,col= ifelse(anova[,25] <=0.1,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,A1]* AC,Yc+CONT[,A2]* AC,pch=ifelse(anova[,25] <=0.05,19,21),cex= ifelse(anova[,25] <=0.1,sqrt(abs(anova[,9])/pi)*ANO,0),col= COLcont,lwd=0.75,bg="white")

return(subset(anova[,c(1,5)],anova[,25] <=0.05))

})

dev.off()

##############################
#####SUMMARY FIGURE (MAIN)
setwd(pD4)

pdf("New-Main-Figure-WOout.pdf",height=6,width=7)

CEXleg=0.6

zones<-matrix(c(1,4,1,4,1,5,2,5,2,6,3,6,3,6),ncol=7)
layout(zones)
#layout.show(max((zones)))
par(bty="n")

#Color code for clades
ColCLAu<-sapply(c(CLAu,"un"),function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
})

#Color code for time
ColT<-sapply(TIME,function(x){
	return(subset(COLfac[,2],as.numeric(COLfac[,1])==x))
})

###Global PCA

ACPh=dudi.pca(FRallH, scannf=F)
part<-round(100* ACPh $eig/sum(ACPh $eig),2)

par(mar=c(0,0,0,0),bty="n")
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"a",cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

x1=min(ACPh $li[,1])
x2=max(ACPh $li[,1])+2
y1=min(ACPh $li[,2])
y2=max(ACPh $li[,2])+1


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=1,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(subset(ACPh $li[,c(1,2)],SITE=="SBL"),groups= subset(SITE,SITE=="SBL"),col= rgb(0.5,0.5,0.5),lty=1,lwd=0.5)

ordispider(subset(ACPh $li[,c(1,2)],SITE=="MSH"),groups= subset(TIMEr,SITE=="MSH"),col= rgb(0.8,0.8,0.8),lty=1,lwd=0.5)
	

points(ACPh $li[,c(1,2)],pch=ifelse(SITE=="SBL",17,24),cex=0.4,col="black",bg="white",lwd=0.5)

ordiellipse(subset(ACPh $li[,c(1,2)],SITE=="SBL"),groups= subset(SITE,SITE=="SBL"),border=NA,alpha=80,draw="polygon",col= "grey",lty=3,lwd=0.5,label=F,cex=0.8,font=1,kind="sd",conf=0.3)

ordiellipse(subset(ACPh $li[,c(1,2)],SITE=="MSH"),groups= subset(TIMEr,SITE=="MSH"),border="grey",alpha=130,draw="polygon",col= "white",lty=1,lwd=0.5,label=F,cex=0.8,font=1,kind="sd",conf=0.3)

sapply(sort(unique(SITE)),function(x){
	coordx<-apply(subset(ACPh $li[,c(1,2)],
		SITE==x),2,mean)
	text(coordx[1]+
		ifelse(x=="MSH",-7,4), 
		coordx[2],
		x,cex=1,font=2)
})

coordt<-sapply(unique(subset(TIMEr,SITE=="MSH")),function(x){
	coordx<-apply(subset(ACPh $li[,c(1,2)],
		TIMEr ==x),2,mean)
	text(coordx[1], 
		coordx[2],
		x,cex=0.6,font=1)
	return(coordx)	
})


#ASV contribution colored according to clades
CONT<-ACPh$co

COLcont<-as.vector(sapply(CLA[,5],function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
}))

#Center of the contribution plot in the main plot
Xc=4
Yc=-9.5

#Contribution amplification
AC=7

#Size of factor contribution tested by anova
ANO=4

text(Xc,Yc-5,"ASV contribution",cex=CEXleg,font=1)
text(Xc,Yc-6,"& Forest/Date association (ANOVA)",cex=CEXleg,font=1)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ[,1] <=0.05, rgb(0.5,0.5,0.5,0.5),NA),lwd=0.5,lend=2,lty=1)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ[,2] <=0.05, rgb(0.5,0.5,0.5,0.5),NA),lwd=0.5,lend=2)

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=21,cex= ifelse(pADJ[,1] <=0.05,sqrt(abs(apply(ANOVA[,c(1,2)],1,sum,na.rm=T))/pi)*ANO,0),col= COLcont,lwd=0.5,bg=rgb(1,1,1,0.9))

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=21,cex= ifelse(pADJ[,2] <=0.05,sqrt(abs(ANOVA[,2])/pi)*ANO,0),col= COLcont,lwd=0.5,bg= COLcont)

par(new=T,mar=c(4,4,2,2))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend(x = "topright",cex= CEXleg, horiz=F,legend= CLAu[c(1,3:8,2,9)],text.col= "black",box.col=NA,border=NA,pch=19,col= ColCLAu[c(1,3:8,2,9)], title="Clades",ncol=3,pt.lwd=0.5)

par(new=T,mar=c(5,5,1,1))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

VL<-c(0.04,0.1,0.2,0.4,0.3,0.3)

legend(x = "bottomleft",pt.cex =sqrt(VL/pi)*ANO,cex= CEXleg, horiz=F,legend= c(VL[1:4],"Forest","Date"),text.col= "black",box.col=NA,border=NA,pch=c(19,19,19,19,21,19),col= c(rgb(0.5,0.5,0.5),rgb(0.5,0.5,0.5),rgb(0.5,0.5,0.5),rgb(0.5,0.5,0.5),"black","black"), title="Part of var. / Factor",ncol=1,pt.lwd=0.5)

####Spatial autocorrelation

BDi=subset(c(1:length(colnames(DIST1))),
	colnames(DIST1)=="BD-All")

Gmin=1 #replace null distances by this value, for log scale
Gmax<-max(as.numeric(DISTtime[,10]))


par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,"b",
	cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

plot(10000,10000,
	xlab="pDistance (log scale, m)",
	ylab="Community dissimilarity (BC)",
	las=1,cex.axis=0.8,cex.lab=1,
	xlim=c(Gmin,Gmax),
	ylim=c(min(as.numeric(DISTtime[,13])),
	max(as.numeric(DISTtime[,13]))),
	log="x",xaxt="n",#Only for log scale on x-axis
	cex.main=0.8,font.lab=1)
	
	Gmin2<- min(subset(
		as.numeric(DISTtime[,10]),
		as.numeric(DISTtime[,10])> Gmin))
	la<-c(log10(Gmin),
		seq(log10(Gmin2),
		log10(Gmax), length.out=5))
	la<-round(10^la)
	
	axis(1, la,c(0,la[2:length(la)]),cex.axis=0.8)
	axis(1, c(Gmin,Gmin2),c(NA,NA),cex.axis=0.8,col="white",lty=3, lwd.ticks=0,lwd=2,lend=2)		

COLs<-sapply(SITEu,function(x){
	SITEx=x
	D1<-subset(DIST1, DIST1[,1]== SITEx)
	Gx<-as.numeric(D1[,10]) #Geographic distance
	Gx=ifelse(Gx==0,Gmin,Gx)
	Bx<-as.numeric(D1[, BDi]) #Bray-Curtis dissimilarity
	Dx<-as.numeric(D1[,11]) #Sampling date
	Dmx<-Dx-min(Dx)
	#pairwise community similarity between trees within dates in function of geographic distance; all pairwise comparisons
	M1<-lm(Bx ~Gx)	
	Sx<-summary(M1)
	ld <- seq(Gmin,max(Gx), length.out=100)
	colx<-subset(COLfac[,2],COLfac[,1]==SITEx)
	COLint<-rgb(col2rgb(colx)[1]/255,
		col2rgb(colx)[2]/255,
		col2rgb(colx)[3]/255,0.1)
	points(Gx,
		Bx,pch=19,cex=0.3,
		col= COLint)	
	px<-Sx$coefficients[2,4]
	COEF<-Sx$coefficients[,1]
	SD<-Sx$coefficients[,2]
	
	#Intervals
	sapply(c(1:5),function(x){
		X=x
		INTx<-sapply(ld,function(x){
			Lx=x
			VAR1x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR2x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			VAR3x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR4x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			return(c(VAR1x, VAR2x, VAR3x, VAR4x))
		})
		INTx<-cbind(apply(INTx,2,max),
			apply(INTx,2,min))	
		polygon(c(ld,ld[length(ld):1]),
			c(INTx[,1],INTx[c(length(ld):1),2]),
			border="NA",col=COLint)
			
	})	
	VARx = COEF[2]* ld + COEF[1]		
	lines(ld, VARx,col= colx,
		lwd=1.5,lty=ifelse(px<=0.05,1,2),lend=2)		
	return(colx)
})

legend(x = "topright",cex= CEXleg, horiz=F,legend= SITEu,text.col= "black",box.col=NA,border=NA,pch=19,col= COLs, title="Forest of origin")


#Temporal autocorrelation 

Dmin=5 #replace null time by this value, for log scale
Dmax<-max(abs(as.numeric(DISTtime[,11])-as.numeric(DISTtime[,12])))
Dmin2<- min(subset(abs(as.numeric(DISTtime[,11])-as.numeric(DISTtime[,12])), abs(as.numeric(DISTtime[,11])-as.numeric(DISTtime[,12]))> Dmin))

la<-c(log10(Dmin),seq(log10(Dmin2),log10(Dmax), length.out=5))
la<-round(10^la)
#la<-ifelse(la<Dmin2,NA,la)

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,"c",
	cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

plot(10000,10000,
	xlab="pTime (log scale, days)",
	ylab="Community dissimilarity (BC)",
	las=1,cex.axis=0.8,cex.lab=1,
	xlim=c(Dmin-1, Dmax),
	ylim=c(min(as.numeric(DISTtime[,13])),
	max(as.numeric(DISTtime[,13]))),
	log="x",xaxt="n",#Only for log scale on x-axis
	cex.main=0.8,font.lab=1)	

	axis(1, la,c(0,la[2:length(la)]),cex.axis=0.8)
	axis(1, c(Dmin,Dmin2),c(NA,NA),cex.axis=0.8,col="white",lty=3, lwd.ticks=0,lwd=2,lend=2)

COLs<-sapply(SITEu,function(x){
	SITEx=x
	D1<-subset(DISTtime, DISTtime[,1]== SITEx)
	Bx<-as.numeric(D1[, BDi]) #Bray-Curtis dissimilarity
	Tx<-abs(as.numeric(D1[,11])-as.numeric(D1[,12])) #Time distance
	
	
	Tx =ifelse(Tx ==0,Dmin, Tx)
	
	#pairwise community similarity between trees within dates in function of geographic distance; all pairwise comparisons
	M1<-lm(Bx ~ Tx)	
	Sx<-summary(M1)
	ld <- seq(Dmin,Dmax, length.out=100)
	colx<-subset(COLfac[,2],COLfac[,1]==SITEx)
	COLint<-rgb(col2rgb(colx)[1]/255,
		col2rgb(colx)[2]/255,
		col2rgb(colx)[3]/255,0.1)
	COLp<-rgb(col2rgb(colx)[1]/255,
		col2rgb(colx)[2]/255,
		col2rgb(colx)[3]/255,0.05)
	aTx<-sample(seq(-0.5,0.5,length.out=101),
		size=length(Tx),replace=T)	
	points(Tx+ aTx,
		Bx,pch=19,cex=0.3,
		col= COLp)	
	px<-Sx$coefficients[2,4]
	COEF<-Sx$coefficients[,1]
	SD<-Sx$coefficients[,2]
	
	#Intervals
	sapply(c(1:5),function(x){
		X=x
		INTx<-sapply(ld,function(x){
			Lx=x
			VAR1x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR2x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			VAR3x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR4x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			return(c(VAR1x, VAR2x, VAR3x, VAR4x))
		})
		INTx<-cbind(apply(INTx,2,max),
			apply(INTx,2,min))	
		polygon(c(ld,ld[length(ld):1]),
			c(INTx[,1],INTx[c(length(ld):1),2]),
			border="NA",col=COLint)
			
	})	
	VARx = COEF[2]* ld + COEF[1]		
	lines(ld, VARx,col= colx,
		lwd=1.5,lty=ifelse(px<=0.05,1,2),lend=2)		
	return(colx)
})

legend(x = "topright",cex= CEXleg, horiz=F,legend= SITEu,text.col= "black",box.col=NA,border=NA,pch=19,col= COLs, title="Forest of origin")

#BC dissimilarity in function of site and time


par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"d",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

D<-as.numeric(DIST1[,11])
SITE1<-as.vector(DIST1[,1])
colx<-sapply(SITE1,function(x){
	return(subset(COLfac[,2],COLfac[,1]==x))
})
BC<-as.numeric(DIST1[,BDi])


COLp<-sapply(colx ,function(x){
	return(rgb(col2rgb(x)[1]/255,
		col2rgb(x)[2]/255,
		col2rgb(x)[3]/255,0.2))
})

plot(D+sample(seq(-1.5,1.5,length.out=length(D))), 
	BC,col= COLp,pch=19,cex=0.3,
	lwd=1,log="",las=1,
	ylab="Community dissimilarity (BC)",
	xlab="Sampling date",cex.axis=0.8,xaxt="n",
	xlim=c(min(D)-5,max(D)+5))

#Month categories
DM<-c("June","June","July","Aug.","Aug.","Sept.","Sept.","Oct.")
Du<-sort(unique(D))
M<-as.vector(sapply(D,function(x){
	return(subset(DM,Du==x))
}))
Mu<-unique(DM)
colmu<-sapply(Mu,function(x){
	mean(unique(subset(D, M==x)))
})
colmu <-(colmu-min(colmu))
colmu <-colmu/max(colmu)
colmu <-rgb(1-colmu,1-colmu,0)

segM<-sapply(Mu,function(x){
	return(mean(subset(D,M==x)))
})
par(cex.axis=0.8)
axis(1,at= segM,labels=Mu)

#density (violin plot) parameters
DENS=100
SVP=1

colsu<-sapply(SITEu,function(x){
	SITEx=x
	Dx<-subset(D, SITE1== SITEx)
	BCx<-subset(BC, SITE1== SITEx)
	colsx=unique(subset(colx, SITE1 ==SITEx))
	COLpx=unique(subset(COLp, SITE1 ==SITEx))
	
	sapply(sort(unique(Dx)),function(x){
		X=x
		BCxx<-subset(BCx,Dx==X)
		XX<-density(BCxx,n= DENS,na.rm=T)$y/SVP
		YY<-density(BCxx,n= DENS,na.rm=T)$x
		xx=subset(XX,abs(median(BCxx)-YY)==min(abs(median(BCxx)-YY)))
		yy=subset(YY,abs(median(BCxx)-YY)==min(abs(median(BCxx)-YY)))
		polygon(
			c(XX+X,(X-XX)[DENS:1]),
			c(YY,YY[DENS:1]),
			col= COLpx,
			border=unique(colsx),lwd=0.5)
			
		segments(xx+X,yy,X-xx,yy,col= colsx,lwd=2.5,lend=2)	
				
	})
	return(unique(colsx))
})

legend(x = "topright",cex= CEXleg, horiz=F,
	legend= unique(SITE),
	box.col=NA,pch=19,
	border="white",col= colsu,text.font=1,
	text.col="black",title="Forest of origin")

 


#Detail spatial autocorrelation per date in MSH only

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"e",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

x="MSH"
SITEx=x
	
D1<-subset(DIST1, DIST1[,1]== SITEx)
Gx<-as.numeric(D1[,10]) #Geographic distance
Bx<-as.numeric(D1[, BDi]) #Bray-Curtis distance
Dx<-as.numeric(D1[,11]) #Sampling date
Dmx<-Dx-min(Dx)

TIMElegx<- t(sapply(sort(unique(Dx)),function(x){
		return(subset(TIMEleg,as.numeric(TIMEleg[,1])==x))
	}))

		
Gmin=1
Gmin2<- min(subset(Gx,Gx> Gmin))
	
ld <- seq(Gmin,max(Gx), length.out=100)
	
la<-c(log10(Gmin),
	seq(log10(Gmin2),
	log10(max(Gx)), length.out=5))
la<-round(10^la)


plot(10000,10000,
	log="x",xlab="pDistance (log scale, m)",
	ylab="Community dissimilarity (BC)",
	las=1,cex.axis=0.8,cex.lab=1,
	xlim=c(min(ld),max(ld)),
	ylim=c(min(Bx),
	max(Bx)),
	main="",
	cex.main=1,xaxt="n",font.lab=1)	
		
	axis(1, la,c(0,la[2:length(la)]),cex.axis=0.8)
	axis(1, c(Gmin,Gmin2),c(NA,NA),cex.axis=0.8,col="white",lty=3, lwd.ticks=0,lwd=2,lend=2)
	
	SUM1<-t(sapply(sort(unique(Dx)),function(x){
		colxx=subset(COLfac[,2], COLfac[,1]==x)
		colxx3=rgb(col2rgb(colxx)[1]/255,
			col2rgb(colxx)[2]/255,
			col2rgb(colxx)[3]/255,0.3)
		colxx6=rgb(col2rgb(colxx)[1]/255,
			col2rgb(colxx)[2]/255,
			col2rgb(colxx)[3]/255,0.05)	
		Gxx<-subset(Gx,Dx==x)
		Bxx<-subset(Bx,Dx==x)
		points(ifelse(Gxx==0,Gmin,Gxx),Bxx,pch=19,cex=0.3,
			col= colxx3)
		Mx<-lm(Bxx ~Gxx)
		px<-summary(Mx)$coefficient[2,4]
		COEF<-summary(Mx)$coefficient[,1]
		Axx<-c(round(anova(Mx)[,2]/
			sum(anova(Mx)[,2]),4),
			round(anova(Mx)[1,5],4))
		SD<-summary(Mx)$coefficient[,2]	
	#Intervals
	sapply(c(1:4),function(x){
		X=x
		INTx<-sapply(ld,function(x){
			Lx=x
			VAR1x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR2x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			VAR3x = (COEF[2]-SD[2]*X)* Lx + 
				(COEF[1]+SD[1]*X)
			VAR4x = (COEF[2]+SD[2]*X)* Lx + 
				(COEF[1]-SD[1]*X)
			return(c(VAR1x, VAR2x, VAR3x, VAR4x))
		})
		INTx<-cbind(apply(INTx,2,max),
			apply(INTx,2,min))	
		polygon(c(ld,ld[length(ld):1]),
			c(INTx[,1],INTx[c(length(ld):1),2]),
			border="NA",col= colxx6)
			
	})
		VARx = COEF[2]* ld + COEF[1]		
		lines(ld, VARx,
			col= colxx,
			lwd=1.5,lty=ifelse(px<=0.05,1,2),lend=2)		
	}))
	
	legend("topleft",
	legend= TIMElegx[,2],cex= CEXleg, col= TIMElegx[,3],bty="n",title="MSH",pch=19)


####Test for node significance: BC dissimilarity in function of pairwise time


###rTREE with nodes labelled with PERMANOVA part of variance

#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=1.25

#Compression of tree from the root (>1) 
COMP=1.1

#Margin size
MAR=1.1

#tip attribute coordinates (concentric circles)
Ai=1.01+seq(0.05,1,length.out=10)
Ae=Ai+0.8*(Ai[2]-Ai[1])

#Number of graduations in the PS scale
Nscale=6

#Branch color
COLbr<-rgb(0.9,0.9,0.9)

Lmin=min(c(coordx[,5], coordx[,7]))
Lmax=max(coordx[,5])*COMP

#Limit of the graphic
rLIM=MAR* Lmax
z=max(c(coordx[,1], coordx[,2],coordx[,3]))


#PS Scale
	
#Level in the scale
Lscale=round(seq(Lmin,Lmax,length.out=Nscale),0)
#Corresponding PS values
PSscale<-unlist(sapply(Lscale,function(x){
	return(subset(as.vector(SCALE[,4]), as.vector(SCALE[,2])==x))
}))
Xscale<-sapply(c(1:(Nscale-1)),function(x){
	X=Lmax-Lscale[x]
	segments(-0.01*rLIM,X,-0.025*rLIM,X)
	return(X)
})
maxTIP=max(Xscale)

PSscale<-PSscale[c(1:length(Xscale))]
segments(-0.01*rLIM,min(Xscale),
	-0.01*rLIM,max(Xscale))
text(-0.005*rLIM, Xscale, round(PSscale,2)
	,cex=0.6,pos=2,adj=1,font=1)
text(-0.005*rLIM, Ae[1]* maxTIP, "PS",
	cex=0.8,font=1,pos=2,adj=1,)

#1 ASVs
	
Tmx<-intersect(subset(Tx, type[,1]=="ASV"),subset(Tx, type[,3]=="Methylobacterium"))

#ASV clade assignation
CA<-sapply(Tmx,function(x){
	nx<-type[x,2]
	cx=subset(CLAo[,5],CLAo[,1]== nx)
	return(cx)
})


par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"f",cex=1.5,font=2)

par(mar=c(1,1,1,1),bty="n",new=T)

#Frame
plot(-10000,-10000,xlim=c(-rLIM, rLIM),
	ylim=c(-rLIM, rLIM),xlab="",
	ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n")


#Remove nodes after strain N228 (last Methylobacterium in the list of tree tip labels)

nmx<- subset(c(1:length(type[,1])),type[,2]=="NS228")


nmcx<-subset(c(1:length(coordx[,1])),coordx[,2]<= nmx)
	
MOx<-subset(coordx[nmcx,6],coordx[nmcx,5]==max(coordx[nmcx,5]))	
	
nmcx<-c(nmcx,subset(c(1:length(coordx[,1])), coordx[,4]== MOx))	

COLbr=rgb(0.6,0.6,0.6)
	
#TREE
COORD<-t(sapply(nmcx,function(x){
	#Nodes
	Yb=coordx[x,1]
	Xb1=Lmax-coordx[x,5]
	angb<-2*pi*(Yb/z)*OPEN
	XXb1<-Xb1*sin(angb)
	YYb1<-Xb1*cos(angb)
	#Branches
	Xb2=Lmax-coordx[x,7]
	Xb2<-ifelse(Xb2<0,0,Xb2)
	XXb2<-Xb2*sin(angb)
	YYb2<-Xb2*cos(angb)
	##rake
	Yn=seq(coordx[x,2],coordx[x,3],length.out=100)
	angn=2*pi*(Yn/z)*OPEN
	XXn<-Xb1*sin(angn)
	YYn<-Xb1*cos(angn)
	points(XXn, YYn, type="l",
		col=COLbr,lwd=0.75)
	segments(XXb1, YYb1,XXb2, YYb2,
		col=COLbr,lwd=0.75)
	#points(XXb1, YYb1,pch=19,col=rgb(0,0,0,1),cex=0.5)
		return(c(x,XXb1, YYb1))
}))


#Nodes (supported by at least 30% of bootstrapps)

setwd(pD3)
#Tree topology in level format
LEVELS<-read.table(paste("Levels_MINBOOT=",
	0,".txt",sep=""))
colnames(LEVELS)<-paste("L",
	c(1:(length(LEVELS[1,]))),sep="-")
	
#remove empty levels
LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))

#Only nodes with at least 30% of boot support
Nbx<-subset(PSn[,2],PSn[,8]>=0.3)



#Number of ASV per nodes: only keep nodes with at least two ASVs
Nbx<-subset(Nbx ,sapply(Nbx,function(x){
	#Node
	NX=x
	#Level
	LEVELx <-subset(PSn[,1],PSn[,2]==NX)	
	#Column index in the level file
	Lx<-subset(c(1:length(LEVELS[1,])),paste("L", LEVELx,sep="-")==colnames(LEVELS))
	#Nodes in this level
	ND<-as.vector(LEVELS[, Lx])
	#ASV line indexes in the level file for this node
	STx<-subset(c(1:length(LEVELS[,1])),ND== NX)
	#Corresponding ASV names
	nx= type[STx,2]	
	return(length(intersect(nx,
		paste("Meth",colnames(FRallH),sep="-"))))	
})>1)

#Calculate average absolute abundance per node and per site
Fbx<-t(sapply(Nbx,function(x){
	#Node
	NX=x
	#Level
	LEVELx <-subset(PSn[,1],PSn[,2]==NX)	
	#Column index in the level file
	Lx<-subset(c(1:length(LEVELS[1,])),paste("L", LEVELx,sep="-")==colnames(LEVELS))
	#Nodes in this level
	ND<-as.vector(LEVELS[, Lx])
	#ASV line indexes in the level file for this node
	STx<-subset(c(1:length(LEVELS[,1])),ND== NX)
	#Corresponding ASV names
	nx= type[STx,2]
	
	FRx<-sapply(colnames(FRall),function(x){
		z<-length(intersect(nx ,
			paste("Meth",x,sep="-")))
		FRxx<-as.vector(t(subset(t(FRall),
			colnames(FRall)==x)))
		return(z* FRxx)
	})
	
	FRx<-apply(FRx,1,sum)
		
	return(sapply(SITEu,function(x){
		return(mean(subset(FRx,SITE==x)))
	}))
}))	

####Spatial and temporal autocorrelation test for each node in each site

#Threshold for frequency per node per site
FT=0
#Relative size of variance explained by pairwise time
ANO=4

setwd(pD4)
sapply(SITEu,function(x){
	SITEx=x
	#Node that can be tested in this site
	Nbs=subset(Nbx,
		as.vector(subset(t(Fbx),
		SITEu==SITEx))> FT)
	#SAMPLE index for this site
	Sx<-subset(c(1:length(SITE)),SITE==SITEx)
	#Shape of the node (square= MSH; circle=SBL)
	PCH<-ifelse(SITEx=="MSH",0,1)
	
	TEMPx<-t(sapply(Nbs,function(x){
		print(x)
		#Node
		NX=x
		#Level
		LEVELx <-subset(PSn[,1],PSn[,2]==NX)
		#Column index in the level file
		Lx<-subset(c(1:length(LEVELS[1,])),
			paste("L", LEVELx,sep="-")
			==colnames(LEVELS))
		#Nodes in this level
		ND<-as.vector(LEVELS[, Lx])
		#ASV line indexes in the level file for this node
		STx<-intersect(Tmx,
			subset(c(1:length(LEVELS[,1])),
				ND== NX))
		#Corresponding ASV names
		nx= intersect(type[STx,2], ASVn)
		#Abundance (Hellinger transformation)
		fx<-as.data.frame(sapply(nx,function(x){
			return(subset(t(FRallH),
				paste("Meth",
				colnames(FRallH),
				sep="-")==x))
		}))[Sx,]
		
		#Pairwise index	
		DISTx<-matrix(unlist(
			sapply(c(1:(length(Sx)-1)),
			function(x){
			X=x
			SS1<-SS[Sx[X]] #Subsite for this sample
			T1= TIME[Sx[X]]#Sampling time for this sample
			S1<-as.numeric(fx[X,])#corrected ASV relative abundance for this sample
			SPE1<-SPE[Sx[X]] #host tree species for this sample
			TREE1<-ID[Sx[X]] #tree identity for this sample
			MISEQ1<-MISEQ[Sx[X]] #sequencing batch for this sample
			return(sapply(c((X+1):(length(Sx))),
				function(x){
				T2= TIME[Sx[x]]
				SS2<-SS[Sx[x]]
				S2<-as.numeric(fx[x,])
				SPE2<-SPE[Sx[x]]
				TREE2<-ID[Sx[x]]
				MISEQ2<-MISEQ[Sx[x]]	
				D12<-subset(DISTtree[,2],DISTtree[,1]
					==paste(sort(c(TREE1,TREE2)),
					collapse="_")) #Spatial distance
				B12=vegdist(t(cbind(S1,S2)),
					method="bray") #Bray dissimilarity
				return(c(SS1,SS2,TREE1,TREE2,
					SPE1,SPE2,MISEQ1,MISEQ2,
					D12,T1,T2,B12))		
			}))
		})),ncol=12,byrow=T)
		
		#Temporal correlation
		Dx<-abs(as.numeric(DISTx[,10])-as.numeric(DISTx[,11]))
		BCx<-as.numeric(DISTx[,12])
		#Linear model
		Mx<-lm(BCx ~Dx)
		Ax<-as.data.frame(anova(Mx))
		SLx<-summary(Mx)$coefficients
		pvar=Ax[1,2]/sum(Ax[,2]) #part of variance explained by pairwise time
		svar=Ax[1,5] #significance of this variance
		BCm<-SLx[1,c(1:2)] #average and sd BC dissimilarity index
		Sm<-SLx[2,c(1:2)] #average and sd slope of BC in function of pairwise time
		
		#Coordinates of the node in the graphic
		COORDxx<-subset(COORD,coordx[nmcx,4]==NX)
		#COORDxx<-subset(COORD,COORD[,1]==NX)
		XX=COORDxx[1,2]
		YY=COORDxx[1,3]	
		
		
		return(c(NX,LEVELx,XX,YY, pvar,svar, BCm, Sm))
		
	}))
	
	colnames(TEMPx)<-c("node","level","X","Y","pVAR","sVAR","meanBC","sdBC","meanSlope","sdSlope")
	setwd(pD4)
	write.table(TEMPx,paste("Temp-autocorrelation-per-node-Site=",SITEx,"-WOout.txt",sep=""))
})	
	
	#Plot part of variance in BC explained by pairwise time. If not significant, no plot. If significant but the slope is negative (unexpected), empty point. If significant and the slope is positive (expected), full point. Size proportional to the part of variance. Shade proportional to the strenght of the (positive) slope
	setwd(pD4)
	TEMPm<-read.table("Temp-autocorrelation-per-node-Site=MSH-WOout.txt",header=T)
	TEMPs<-read.table("Temp-autocorrelation-per-node-Site=SBL-WOout.txt",header=T)

	
	#relative size of points proportional to slope
	MAXslope<-max(c(TEMPm[,9],TEMPs[,9]))
	Mslope<-TEMPm[,9]/MAXslope
	Sslope<-TEMPs[,9]/MAXslope
	RS=5
	
	#Bonferonni correction of pvalue
	pvalM<-p.adjust(TEMPm[,6])
	pvalS<-p.adjust(TEMPs[,6])
	
	sapply(c(1:length(TEMPm[,1])),function(x){
		XX<-TEMPm[x,3]
		YY<-TEMPm[x,4]
		slopex<-abs(Mslope[x])
		pvalx<-pvalM[x]
		slx<-slopex	
		#slx<-sqrt(slopex/pi)	
		pvalx<-ifelse(pvalx<=0.001,1,
			ifelse(pvalx<=0.01,0.75,
			ifelse(pvalx<=0.05,0.5,0)))
		colx<-"green2"
		colxx<-rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255, pvalx)
		colx<-rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400, pvalx)			
		points(XX,YY,pch=21,
			bg= ifelse(Mslope[x]<0,NA, colxx),
			col= ifelse(Mslope[x]<0,NA, colx),
			lwd=0.5,cex=slx* RS)
	})
	
	sapply(c(1:length(TEMPs[,1])),function(x){
		XX<-TEMPs[x,3]
		YY<-TEMPs[x,4]
		slopex<-abs(Sslope[x])
		pvalx<-pvalS[x]
		slx<-slopex
		#slx<-sqrt(slopex/pi)	
		pvalx<-ifelse(pvalx<=0.001,1,
			ifelse(pvalx<=0.01,0.75,
			ifelse(pvalx<=0.05,0.5,0)))
		colx<-"orange"
		colxx<-rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255, pvalx)
		colx<-rgb(col2rgb(colx)[1]/400,
			col2rgb(colx)[2]/400,
			col2rgb(colx)[3]/400, pvalx)			
		points(XX,YY,pch=21,
			bg= ifelse(Sslope[x]<0,NA, colxx),
			col= ifelse(Sslope[x]<0,NA, colx),
			lwd=0.5,cex=slx* RS)
	})


#Consensus clade label

Aei2 <-mean(c(Ae[2],Ai[2]))

sapply(setdiff(COLcla[,1],c("Enterovirga","Microvirga")),
	function(x){
	ST1<-min(c(subset(Tmx,CA==x),subset(Tx, type[,6]==x)))
	ST2<-max(c(subset(Tmx,CA==x),subset(Tx, type[,6]==x)))
	nax<-length(subset(Tmx,CA==x))
	angx=2*pi*(mean(c(ST1, ST2))/z)*OPEN
	XX1i<-Ae[3]* maxTIP*sin(angx)
	YY1i<-Ae[3]* maxTIP*cos(angx)						
	text(XX1i, YY1i,x,cex=ifelse(nax==0,0.4,0.8),col= "black",font=1)
	
	ST12<-seq(ST1, ST2,length.out=200)
	angx=2*pi*(ST12/z)*OPEN
	ang1=2*pi*(ST1/z)*OPEN
	ang2=2*pi*(ST2/z)*OPEN
	XX1i<-c(Aei2* maxTIP*sin(ang1),
		Ae[2]* maxTIP*sin(angx),
		Aei2* maxTIP*sin(ang2))
	YY1i<-c(Aei2* maxTIP*cos(ang1),
		Ae[2]* maxTIP* cos(angx),
		Aei2* maxTIP* cos(ang2))
	points(XX1i, YY1i,type="l",lwd=0.5)	
})

	#Strain label
	
	ang.lab<-seq(360,(1-OPEN)*360,length.out=length(Tx))-270
	adjx<-ifelse(ang.lab>-90,0, 1)
	ang.lab<-ifelse(ang.lab>-90,ang.lab, ang.lab+180)
	
	ATx=1
	y1=Ai[ATx]
	y2=(Ae-((Ae-Ai)/2))[ATx]
	
	
	sapply(Tx[1: nmx],function(x){
		print(x)
		cx<-type[x,6]
		nx<-type[x,2]
		cx<-ifelse(is.na(cx)==T,subset(CLA[,5],CLA[,1]== nx),cx)
		cx<-c(subset(COLcla[,2],COLcla[,1]==cx),"grey")[1]
		angx=2*pi*(x/z)*OPEN
		XX1i<-y1* maxTIP*sin(angx)
		YY1i<-y1* maxTIP*cos(angx)						
		
		points(1.02*maxTIP*sin(angx), 1.02*maxTIP*cos(angx),pch=1,cex=0.2,col=ifelse(type[x,1]=="ASV", "black",NA),lwd=0.5)
		
		text(XX1i, YY1i,nx,cex=0.15,srt= ang.lab[x],adj= adjx[x],col=ifelse(type[x,1]=="ASV", "black","grey"),font=1)	
			
	})
	
	#Legend for slope
	
	TEMPscale<-seq(1,0,length.out=5)[1:4]
	TEMPn<-round(seq(MAXslope,0,length.out=5)[1:4],3)
	
	legend("topright",pt.cex = TEMPscale* RS,
		cex= CEXleg, horiz=F,legend= TEMPn,
		text.col="black",box.col=NA,border=NA,
		pch=21,col= "white",pt.bg="grey",
		title= expression("slope lm(BC~pTime)"), 
		pt.lwd=0.5)
		
	legend(x = "bottomright",cex= CEXleg, horiz=F,
		legend= unique(SITE),
		box.col=NA,pch=19,
		border="white",col= colsu,text.font=1,
		text.col="black",title="Forest of origin")	
	
dev.off()




	
################
################
################
################
################

