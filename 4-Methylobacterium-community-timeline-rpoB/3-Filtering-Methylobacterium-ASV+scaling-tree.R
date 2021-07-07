#Work directories

pD1<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-1"

pD2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-2"

pD3<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/4-Methylobacterium-community-timeline-rpoB/data-3"

########################################################
##STEP 1: Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs

setwd(pD1)

####Load Metadata for each sample
METAsamp<-read.table("METADATA-wo-NEGcontrols.txt")

#Indexes for positive controls
vMETH<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]=="METH")

#Indexes for samples
vSAMP<-subset(c(1:length(METAsamp[,1])),METAsamp[,1]!="METH")

##EXCLUDE SAMPLES FROM SEQUENCING BATCH rpoB-1-step: WAS DONE WITH DIFFERENT EXTRACTION KIT AND SHOW A HIGHER PROPORTION OF UNASSIGNED ASV, CAULOBACTERALES AND BRADYRHIZOBIACEAE (SEE Summary_rpoB-tax-per-Factor_TF=250_TR=200_TL=20_NS=3000_Prare=0.005-wo-neg-controls-WITH-FIRST-BATCH.pdf)

vOUT<-subset(c(1:length(METAsamp[,1])),METAsamp[,19]=="rpoB-1-step")
vSAMP<-setdiff(vSAMP,vOUT)	
	
##Load Rarefied SeqTAB of absolute abundances

setwd(pD1)

TFx=250
TRx=200
TLx=20
NS=3000
Prare=0.005

seqtab.nochim<-read.table(paste("RAREFIED_SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""))

ASVs<-rownames(seqtab.nochim)

###Only keep abundances for samples
FRall<-t(seqtab.nochim[,vSAMP])

###Apply Hellinger transformation on relative abundance to account for rare ASVs
library(vegan)

FRallH <-as.data.frame(decostand(FRall ,method="hellinger", MARGIN=1))

colnames(FRallH)<-paste("ASV", c(1:length(FRallH[1,])),sep="")
colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="")

##Load Taxonomic assignation for each ASV
setwd(pD2)

taxaALL <-read.table(paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

##Simplify taxonomy (Emphasize on Order, and Familly within Rhizobiales)

ASSsimp<-ifelse(is.na(as.vector(taxaALL[,5]))==T,"unknown",ifelse(is.na(as.vector(taxaALL[,6]))==T,as.vector(taxaALL[,5]),ifelse(as.vector(taxaALL[,5])=="Rhizobiales",paste(taxaALL[,5], taxaALL[,6],sep="_"),as.vector(taxaALL[,5]))))

#Methylobacteriaceae ASVs
ASVm<-subset(c(1:length(ASSsimp)), ASSsimp=="Rhizobiales_Methylobacteriaceae")

#Export Methylobacteriaceae ASVs in fasta format

SEQm<-sapply(ASVm,function(x){
	return(list(unlist(strsplit(ASVs[x],split=""))))
})

library(seqinr)

setwd(pD3)

write.fasta(SEQm ,paste("Methylobacteriaceae-ASVs_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,".fas",sep=""),names=paste(taxaALL[ASVm,7],"ASV", ASVm,sep="-"))

#Manually combine and align ASVs sequences with reference sequences from complete genomes (1-Phylo-of-plant-ass-Methylo-diversity/data/rpoB-update-sept-2020.fas) and 2017 strains used for rpoB marker 

#Simplify reference names and retrieve characteristics

p2<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"

#Reference sequence informations (names in second column)
setwd(p2)
CHAR<-read.table("sequence_characteristics.txt",header=T)

setwd(pD3)
MethASV<-read.fasta("ASV+REF-FULL.fas")

#rename reference sequences according to strain names before performing phylogeny (for ASVs, manually simplify Methylobacterium by Meth, Microvirga by Micr and Enterovirga by Ente)

NEWnames<-as.vector(sapply(names(MethASV),function(x){
	Nx=x
	Ny<-c(subset(
	
		paste(as.vector(CHAR[,4]),as.vector(CHAR[,5]),as.vector(CHAR[,3]),as.vector(CHAR[,7]))
		,
		as.vector(CHAR[,2]== Nx)),Nx)[1]
	return(Ny)
}))

write.fasta(MethASV,names= NEWnames,"ASV+REF-FULL.fas")

########################################################
##STEP 2: Perform phylogeny in MEGA (ML tree, 200 permutations, pairwise deletion) of Methylobacteriaceae ASVs + references : ASV+REF-FULL-from-nwk.mtsx

########################################################
##STEP 3: Retrieve and format the ML tree (no branch length): scale the tree on pairwise nucleotide similarity

###FUNCTIONS
num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

empty.plot<-function(x1,x2,y1,y2,Nx,Ny){
	plot(-1000,-1000,
	xlim=c(x1,x2),ylim=c(y1,y2),
	xaxt="n",yaxt="n",xlab=Nx,ylab=Ny)
}

library(ape)
library(seqinr)

setwd(pD3)

treex<-read.tree("ASV+REF-from-nwk-boot-only.nwk")

###Corresponding Fasta file
fastax<-read.fasta("ASV+REF.fas")


#Extract node coordinates, sequence names and bootstrapp values from the newick tree
coordx<-treex$edge #tree topology: nodes (first column) and their daugther nodes. Daughther node include tips (tips have indexes 1 to n and nodes, n+1 to n+N where n is the number of tips and N the number of nodes, tips excluded)
namex<-treex$tip.label #sequence (tip) indexes
bootx<-treex$node.label #bootstrapp value for each node
nodex<-sort(unique(coordx[,1])) #internal nodes indexes

#read p-distance matrix among sequences (all sites) Was calculated in MEGA7 (all positions included)

pdist<-as.matrix(read.csv("ASV+REF-pdist-formated.csv",header=T,sep=";"))


###Each node will be associated to a level, according to its hierachical position in the tree 
#number of levels (must be high enough to cover all the nodes)
NL=50
#number of levels to add between each level (will be used to adjust levels when searching for the best correlation between pdistance and levels)
NLex=10

#Minimum bootstrapp value to keep nodes
bx=0


#Sequence indexes
tips<-c(1:length(namex))

###Build a matrix with each sequence in rows and each level in column; fill the matrix iteratively for each level with nodes names (empty positions as 0)


#initial files (sequence index only)
write.table(cbind(tips),"treetemp.txt",row.names=F,col.names=F)
write.table(cbind(tips),"boottemp.txt",row.names=F,col.names=F)


#iteratively fill the matrix (stacked from the tips) with node names
sapply(c(1:NL),function(x){
	print(x)
	tempx<-read.table("treetemp.txt")
	tempb<-read.table("boottemp.txt")
	#nodes in this level
	Lx<-as.vector(tempx[,x])
	#retrieve mother nodes
	Ln<-sapply(Lx,function(x){
		return(c(subset(coordx[,1],coordx[,2]==x),0)[1])
	})
	#retrieve bootstrapp value for each mother node
	Bn<-sapply(Ln,function(x){
		return(c(subset(bootx,nodex==x),0)[1])
	})
	Bn<-ifelse(Bn=="",0,Bn)
	write.table(cbind(tempx, Ln),
		"treetemp.txt",row.names=F,col.names=F)
	write.table(cbind(tempb, Bn),
		"boottemp.txt",row.names=F,col.names=F)		
})


tempx<-as.matrix(read.table("treetemp.txt"))
tempb<-as.matrix(read.table("boottemp.txt"))

#Remove nodes with bootstrapp value lower than bx but keet the root node
tempx<-ifelse(tempx ==(length(namex)+1),tempx,ifelse(tempb<bx,0,tempx))


#stack the matrix to the root, so each node is associated to a single level 
tempx<-sapply(c(1:length(tempx[,1])),function(x){
	tx<-as.numeric(tempx[x,])
	return(c(tx[1],
		subset(tx[2:length(tx)],tx[2:length(tx)]==0),
		subset(tx[2:length(tx)],tx[2:length(tx)]!=0)))
})
#Remove empty levels
tempx<-t(subset(tempx,apply(tempx,1,max)!=0))

#Sort sequences according to new tree topology

NTx<-length(tempx[1,])

tempx <-sort(sapply(c(1:length(tempx[,1])),function(x){
	return(paste(tempx[x,c(NTx:1)],collapse="_"))
}))

tempx<-as.matrix(t(sapply(tempx,function(x){
	return(as.numeric(unlist(strsplit(x,split="_"))[NTx:1]))
})))


###New tree topology mother and daughter nodes
coordnx <-matrix(as.numeric(unlist(strsplit(unique(unlist(sapply(c(1:length(tempx[,1])),function(x){
	tx<-as.vector(tempx[x,])
	tx<-subset(tx,tx!=0)
	return(sapply(c(2:length(tx)),function(x){
		return(paste(tx[x],tx[x-1],sep="_"))
	}))
}))),split="_"))),ncol=2,byrow=T)

#New order of sequence names
namex<-sapply(as.vector(tempx[,1]),function(x){
	return(namex[x])
})

#re-index tips in tempx, accordingly
tempx<-cbind(tips,tempx[,c(2:length(tempx[1,]))])

write.table(tempx[,1],"tree-topo.txt",row.names=F,col.names=F)

###Insert NLex empty levels between each levels with nodes information
tempex<-sapply(c(1:NLex),function(x){
	return(seq(0,0,length.out=length(tempx[,1])))
})

sapply(c(2:length(tempx[1,])),function(x){
	tx<-as.matrix(read.table("tree-topo.txt"))
	write.table(cbind(tx,tempex,tempx[,x]),"tree-topo.txt",row.names=F,col.names=F)	
})

###Ad several NLex empty levels at the root of the tree (add more if necessary)
tempx <-cbind(as.matrix(read.table("tree-topo.txt")), tempex,tempex,tempex, tempex, tempex
, tempex,tempex,tempex, tempex, tempex, tempex,tempex,tempex, tempex, tempex
, tempex,tempex,tempex, tempex, tempex, tempex,tempex,tempex, tempex, tempex
, tempex,tempex,tempex, tempex, tempex, tempex,tempex,tempex, tempex, tempex
, tempex,tempex,tempex, tempex, tempex)

NL=length(tempx[1,]) #number of levels
write.table(tempx,"tree-topo.txt",row.names=F,col.names=F)

write.table(namex,"tree-names.txt",row.names=F,col.names=F)



###Generate a file with the following information for each node (tips excluded)
# "Level": the initial level for this node
# "node": the name of the node
# "mother_node": the name of the mother node in the tree
# "sequences": the number of sequences in this node
# "PS_median": the median of pairwise similarity (PS) among sequences from this node
# "PS_max": the maximum PS
# "PS_min": the minimum PS
# "bootstrapp": the bootstrapp value for this node



PS<-matrix(unlist(sapply(c(2:length(tempx[1,])),function(x){
	print(x)
	X=x
	vx<-setdiff(unique(tempx[,X]),0)	
	sapply(vx,function(x){
		#print(x)
		stx<-subset(tempx[,1], tempx[,X]==x)
		#Names in nwk file (strain names)
		nx<-namex[stx]
		
		#retrieve pdistance among these strains	
		IDx<-as.vector(sapply(nx ,function(x){
			as.numeric(subset(pdist[,2],pdist[,1]==x))
		}))
		PSx <-as.vector(as.numeric(pdist[IDx,(IDx+2)]))
		
		PSx<-subset(PSx,is.na(PSx)==F)
		
		#CONVERT pdistance in psimilarity (PS)
		PSx=1-PSx
				
		return(c(X,x,length(nx),median(PSx),
			max(PSx),min(PSx)))
	})	
})),ncol=6,byrow=T)

emb<-sapply(PS[,2],function(x){
	return(c(subset(coordnx[,1],coordnx[,2]==x),0)[1])
}) #Mother nodes

dau<-sapply(PS[,2],function(x){
	return(paste(subset(coordnx[,2],coordnx[,1]==x),collapse="_"))
}) #daughter nodes

PS<-cbind(PS[,c(1:2)],emb,PS[,c(3:6)])

boot<-as.numeric(sapply(PS[,2],function(x){
	return(subset(bootx,nodex==x))
})) #bootstrapp value for each node

PS<-cbind(PS,boot)

colnames(PS)<-c("Level","node","mother_node","sequences","PS_median","PS_max","PS_min","bootstrapp")

write.table(PS,paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),row.names=F,col.names=F)


###This step will allow to optimize match between level and PS. WARNING: only use if the evolutionnary rate seems quite similar in each branch. Each node being associated to a level and to a median PS value, levels will be iteratively moved for each node (while respecting the tree topology) until a maximal correlation is found between PS and Level (then one can assume that node from the same level have very similar PS)

#Convert PS value in negative log scale to zoom at the tips of the tree 
CORx=1.005 #to avoid -inf value when conveerting in inverse log (must be >1)
PSc<-(-log10(CORx-PS[,5]))

##Initial correlation between PS and level
#global correlation between PCs and log(level)
CORi=cor.test(log(PS[,1]), PSc)$estimate

pdf(paste("Tree-iterations_MINBOOT=",bx,".pdf",sep=""))

par(mfrow=c(2,2),mar=c(4,4,2,1),bty="n")
plot(log(PS[,1]), PSc,log="",
	xlab = "Level (log)",
	ylab="modified PS",las=1,cex.axis=0.8,pch=19,cex=0.5,main="before iterations")

par(new=T)
empty.plot(0,1,0,1,"","")	
text(0.25,0.25,paste("r2 =",round(CORi,4)),font=2)	
	
##Number of iterations to find the best correlation
K=10000

CORall<-t(sapply(c(1:K),function(x){
	PSk<-as.matrix(read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep="")))
	tempk<-as.matrix(read.table("tree-topo.txt"))	

	#plot(PSk[,1], PSc,log="x",
	#	xlab = "Level (log scale)",
	#	ylab="modified PS",las=1,cex.axis=0.8,pch=19,cex=0.5)

	#global correlation between PCs and log(level)
	CORi=cor.test(log(PSk[,1]), PSc)$estimate
	
	#par(new=T)
	#empty.plot(0,1,0,1,"","")	
	#text(0.25,0.25,paste("r2 =",round(CORi,4)),font=2)	

	
	#Nodes that can be moved
	MOV<-t(sapply(c(1:length(PSk[,1])),function(x){
		#print(x)
		Gx=PSk[x,2]
		Nx=PSk[x,1]
		Go=subset(tempk[,Nx-1],tempk[,Nx]==Gx)
		Gi=subset(tempk[,Nx+1],tempk[,Nx]==Gx)
		Ge=length(subset(Go,Go==0))
		Gf=length(subset(Gi,Gi==0))
		#return(c(length(Go),Ge))
		return(c(length(Go),Ge,length(Gi),Gf))
	}))
	#Sample one node to move
	#MOV<-sample(subset(PSk[,2],MOV[,1]==MOV[,2]),1)
	MOV<-sample(unique(c(
		subset(PSk[,2],MOV[,1]==MOV[,2]),
		subset(PSk[,2],MOV[,3]==MOV[,4])))
		,1)
	#Initial level of this node
	Ni=subset(PSk[,1],PSk[,2]==MOV)

	#Possible new levels for this node
	tempi<-subset(tempk,tempk[,Ni]==MOV)
	tempi<-ifelse(tempi==MOV,0,tempi)
	Np<-subset(c(1:NL),apply(tempi,2,max)==0)
	Np<-Np[order(Np,decreasing=T)]
	Np<-cbind(Np,c(0,cumsum(Np[1:(length(Np)-1)]-Np[2:(length(Np))]-1)))
	Gp<-subset(Np[,2],Np[,1]==Ni)
	Np<-sort(subset(Np[,1],Np[,2]==Gp))
	
	#Test all possible level for this node and calculate global correlation between PCs and log(level); keep the level with the higest coefficient correlation
	CORp<-sapply(Np,function(x){
		Nnew=ifelse(PSk[,2]==MOV,x,PSk[,1])
		return(cor.test(log(Nnew), PSc)$estimate)
	})
	Nmax<-min(subset(Np,CORp==min(CORp)))
	#print(c(x,CORi,min(CORp)))
	#construct a new tempx file

	tempn<-t(sapply(c(1:length(tempi[,1])),function(x){
		return(c(tempi[x,subset(c(1:NL),c(1:NL)<Nmax)],
		MOV,
		tempi[x,subset(c(1:NL),c(1:NL)>Nmax)]))
	}))

	tempn<-rbind(tempn,tempk[setdiff(c(1:length(tempk[,1])),tempn[,1]),])

	tempn<-tempn[order(tempn[,1]),]

	#construct a new PS file

	PSn<-cbind(ifelse(PSk[,2]==MOV,Nmax,PSk[,1]),
		PSk[,c(2:length(PSk[1,]))])

	write.table(tempn,"tree-topo.txt",row.names=F,col.names=F)
	write.table(PSn,paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),row.names=F,col.names=F)

	print(c(x,round(as.vector(CORi),4)))
	return(c(x,MOV,CORi))
}))

plot(CORall[,1],CORall[,3],type="l",xlab="Iteration",ylab="r2",cex.axis=0.8,main="iterations",las=1)

#point of stabilization
points(min(subset(CORall[,1],CORall[,3]==min(CORall[,3]))),min(CORall[,3]),pch=19,col="red")

PSn<-as.matrix(read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep="")))
tempc<-as.matrix(read.table("tree-topo.txt"))	

colnames(PSn)<-colnames(PS)

write.table(PSn,paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),row.names=F)
write.table(tempc,paste("Levels_MINBOOT=",bx,".txt",sep=""),row.names=F,col.names=F)
	
####Predict PSc for each level using a linear model

#Transform level in log
Lc=log(PSn[,1])

plot(Lc, PSc,
	xlab = "Level (log)",
	ylab="modified PS",las=1,cex.axis=0.8,pch=19,cex=0.5,main="after iterations")

modPS <- glm(PSc ~ Lc)
summary(modPS)
ld <- log(seq(1,NL, length.out=NL))
PSp<-predict(modPS, data.frame(Lc = ld),type = "response")#Predicted modified PS from log level
lines(ld, PSp,col="red",lwd=2)

par(new=T)
empty.plot(0,1,0,1,"","")	
text(0.25,0.25,paste("r2 =",round(cor.test(Lc, PSc)$estimate,4)),font=2,col="red")

plot(PSn[,1], PSn[,5],log="",
	xlab = "Level",
	ylab="PS",las=1,cex.axis=0.8,pch=19,cex=0.5,main="real values")

dev.off()


##########Overview of the tree topology

###Generate a file with nodes and tips coordinate (X=levels; Y=sequence coordinates)
#Coordinates of tips branches
coord0<-t(sapply(c(1:length(tempc[,1])),function(x){
	Sx=tempc[x,1]#tips index
	Nx=subset(tempc[x,],tempc[x,]!=0)[2]#node name
	Lx=subset(c(1:NL),tempc[x,]!=0)[2]#lower level
	#segments(0,x,Lx,x)
	return(c(Sx, Sx, Sx, Sx,min(PSn[,1]),Nx,Lx))
}))	

colnames(coord0)<-c("mean_seq","max_seq","min_seq","seq_or_node","level","mother_node","mother_level")

write.table(coord0,paste("coord_MINBOOT=",bx,".txt",sep=""),row.names=F)
	
#Coordinates of nodes
NLx<-setdiff(subset(c(1:NL),apply(tempc,2,max)!=0),1)

sapply(setdiff(NLx,max(NLx)),function(x){
#sapply(NLx,function(x){
	print(x)
	X=x
	coordx<-as.matrix(read.table(
		paste("coord_MINBOOT=",bx,".txt",sep=""),
		header=T))
	coordn<-subset(coordx,coordx[,7]==X)
	Nm<-unique(coordn[,6])#node names
	coordn <-t(sapply(Nm,function(x){
		mx=mean(subset(coordn[,1],coordn[,6]==x))#mean_seq for this node	
		minx<-min(subset(coordn[,1],coordn[,6]==x))#min_seq for this node	
		maxx<-max(subset(coordn[,1],coordn[,6]==x))#max_seq for this node	
		mo<-subset(PSn[,3],PSn[,2]==x)#mother node
		lo<-c(subset(PSn[,1],PSn[,2]==mo),max(NLx))[1]#level of mother node; if no mother node, root on a artificial level
		return(c(mx, maxx,minx,x,X,mo,lo))
	}))
	colnames(coordn)<-colnames(coordx)
	coordn<-rbind(coordx,coordn)
	write.table(coordn,
		paste("coord_MINBOOT=",bx,".txt",sep=""),
		row.names=F)
})

coordx<-as.matrix(read.table(
	paste("coord_MINBOOT=",bx,".txt",sep=""),header=T))

#Coordinate for the root node
coordn<-subset(coordx,coordx[,7]==max(coordx[,7]))
#coordn<-subset(coordx,coordx[,7]==0)

#coordx<-subset(coordx,coordx[,7]!=max(coordx[,7]))

Nm<-unique(coordn[,6])#node names

coordn<-t(sapply(Nm,function(x){
	mx=mean(subset(coordn[,1],coordn[,6]== x))#mean_seq for this node	
	minx<-min(subset(coordn[,1],coordn[,6]== x))#min_seq for this node	
	maxx<-max(subset(coordn[,1],coordn[,6]== x))#max_seq for this node	
	mo<-0#mother node
	lo=max(coordn[,7])#level of mother node
	return(c(mx, maxx,minx,x,lo,mo,NL))
}))	
	
colnames(coordn)<-colnames(coordx)
coordx<-rbind(coordx,coordn)

write.table(coordx,paste("Tree-coordinates.txt",sep=""),row.names=F)


####Plot the tree

PSr = CORx-(10^(-PSp))#real PS value

#Level scale
SCALE<-as.matrix(cbind(round(ld,4),round(exp(ld),0),round(PSp,4),round(PSr,4)))
colnames(SCALE)<-c("log(Level)","Level","(-log10(1.005-PS))","PS")

write.table(SCALE,paste("SCALE.txt",sep=""),row.names=F)

Lmin=min(PSn[,1])
Lmax=max(PSn[,1])
#Number of graduations in the scale
Nscale=6
Lscale=round(seq(Lmin,Lmax,length.out=Nscale),0)
#Corresponding PS values
PSscale<-sapply(Lscale,function(x){
	return(subset(as.vector(SCALE[,4]), as.vector(SCALE[,2])==x))
})

pdf(paste("Tree.pdf",sep=""),height=12,width=6)

#Start of the plot
par(mar=c(4,1,1,1),bty="n",mfrow=c(1,1))

#Frame
plot(-100,-100,xlim=c(Lmax, Lmin-0.25*(Lmax-Lmin)),ylim=c(0,length(namex)),xlab="PS",ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n")

#PS scale
axis(1, Lscale, round(PSscale,2),cex.axis=0.8)

#Plot the tree
sapply(c(1:length(coordx[,1])),function(x){
	#branch coordinates
	Yb=coordx[x,1]
	Xb1=coordx[x,5]
	Xb2=coordx[x,7]
	segments(Xb1,Yb,Xb2,Yb,col="grey")
	#node coordinates
	Yn1=coordx[x,2]
	Yn2=coordx[x,3]
	segments(Xb1, Yn1, Xb1, Yn2,col="grey")
	#bootstrapp value (if aplicable)
	Nx=coordx[x,4]#node/sequence
	bx=subset(PSn[,8],PSn[,2]==Nx)
	bx<-ifelse(length(bx)==0,0,ifelse(is.na(bx)==T,0,bx))
	points(Xb1,Yb,pch=19,col=rgb(bx,0,0,1),cex=bx/2)
})

#Sequence labels

text(seq(Lmin,Lmin,length.out=length(namex)),
	c(1:length(namex)),
	labels=namex,cex=0.2,pos=4,adj=1)	
	

dev.off()


#########################################
#########################################
## Step 4 - summary of Methylobacterium ASVs and comparison of diversity with isolated strains
 #############################
#########################################
#########################################



pO="/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"


library(parallel)

###FUNCTIONS
num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

reduce.group<-function(groups,ncor) {
	reduce.gx<-function(i){
		gx<-num.strsplit.u(groups[i])
		z<-sapply(groups,function(x){
			return(length(intersect(
				gx,num.strsplit.u(x))))
		})
		return(paste(sort(unique(
			num.strsplit.u(names(subset(z,z!=0))))),
			collapse="_"))
	}
	return(unlist(mclapply(c(1:length(groups)), 
		reduce.gx,mc.cores=ncor)))
}


##### STEP 6a
# LOAD and format meta data and ASV table only for Methylobacteriaceae
###Load data

	#strain (tips of the tree) names
	setwd(pD3)
	
	#minimum bootstrapp value
	bx=0
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
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
vSAMP<-setdiff(vSAMP,vOUT)	

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

##Only keep genus

GENx<-as.vector(taxaALL[,7])

rm(taxaALL)

#Methylobacterium ASVs
ASVm<-subset(c(1:length(GENx)), GENx =="Methylobacterium")

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
	
CLA<-cbind(type[ASVi,2],t(sapply(ASVi,function(x){
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
	
#Sort Methylobacteriaceae ASV data according to tree order
ASVo<-as.numeric(matrix(unlist(
	strsplit(CLA[,1],split="-")),ncol=3,byrow=T)[,3])
	
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
	Fcla<-sapply(sort(unique(CLA[,5])),function(x){	
		FRx<-subset(t(FRallR),CLA[,5]==x)		
		return(apply(FRx,2,sum))
	})

	#Average relative abundance of Clades and number of ASVs
	FclaSum<-t(sapply(sort(unique(CLA[,5])),function(x){	
		FRx<-subset(t(FRallR),CLA[,5]==x)		
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
	
	COLclax<- sapply(sort(unique(CLA[,5])),function(x){
		return(c(subset(COLcla[,2], COLcla[,1]==x),NA)[1])
	})
	
	par(mar=c(4,4,1,1),bty="n")
	FRsite<-sapply(SITEu,function(x){
		
		Fx<-subset(Fcla,SITE==x)
		return(apply(Fx,2,mean))
	})
	
	barplot(FRsite,las=2,col= COLclax,border=NA,las=1,cex.axis=0.8,ylab="Methylobacteriaceae ASV relative abundance")
	
	setwd(pD3)
	
	FRsum<-cbind(FclaSum, FRsite)
	colnames(FRsum)<-c("ASVs","Fr","MSH","SBL")
	
	write.table(FRsum,"Summary-Methylobacteriaceae-ASV.txt")
	
#Fasta file combining aligned strains rpoB sequences and ASVs (nnn removed)

setwd(pO)
strains<-read.table("strain-metadata.txt",header=T)

library(seqinr)
setwd(pD3)
AsvSt<-read.fasta("ASV+2018-Strains-short.fas")

AsvStn<-names(AsvSt)
#strain names
STn<-as.vector(strains[,1])

#Calculate p-distance between sequences in MEGA

pDIST<-as.matrix(read.table("ASV+2018-Strains-short-pdist.txt",header=T))

#Remove non-methylo ASVs

Ix<-as.numeric(pDIST[,1])

Ix<-subset(Ix,sapply(as.vector(pDIST[,2]),function(x){
	nx<-unlist(strsplit(x,split="-"))
	vx=length(nx)
	gx<-nx[1]
	return(ifelse(vx==2,1,ifelse(gx=="Meth",1,0)))
})==1)

#Perform comparison isolation/ASV for different levels of pDIST

sapply(c(0:7)*0.003,function(x){

#Group sequences that have les than pn nucleotide divergence
pn=x
print(pn)

GROUP<-unique(sapply(Ix,function(x){
	px<-as.numeric(pDIST[x,Ix+2])
	px<-ifelse(is.na(px)==T,as.numeric(pDIST[Ix,(x+2)]),px)
	Iy<-subset(Ix, px<=pn)
	return(paste(sort(c(x, Iy)),collapse="_"))
}))


print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))
GROUP =unique(reduce.group(GROUP,ncor=2))
print(length(GROUP))

GROUP<-as.vector(sapply(GROUP,function(x){
	Gx<-num.strsplit.u(x)
	return(paste(as.vector(pDIST[Gx,2]),collapse="_"))
}))

#For each group, count the number of strains and ASVs

COUNT<-t(sapply(GROUP,function(x){
	Gx<-unlist(strsplit(x,split="_"))
	return(c(length(intersect(ASVn, Gx)),length(intersect(STn, Gx))))
	
}))

#Group with only ASVs or strains
GROUPa<-subset(GROUP,as.vector(apply(COUNT,1,min))==0)

# group with at least one strain and one ASV
GROUP<-subset(GROUP,as.vector(apply(COUNT,1,min))!=0)

#For each comparable group, count the number of strains and ASVs abundance

COMP<-t(sapply(GROUP,function(x){
	#print(x)
	Gx <-unlist(strsplit(x,split="_"))
	Sx<-intersect(STn, Gx)
	ST<-length(Sx)
	Ax<-intersect(ASVn, Gx)
	Fx<-sapply(Ax,function(x){
		return(as.vector(subset(t(FRallR),
			paste("Meth",colnames(FRallR),sep="-")==x)))
	})
	
	#Frequency in site where strains come from
	Fsite<-sapply(Sx,function(x){
		X<-as.vector(subset(strains[,3],strains[,1]==x))
		return(sum(apply(subset(Fx, SITE== X),2,mean)))
	})
	#Consensus clade
	CX<-unique(sapply(Ax,function(x){
		return(subset(CLA[,5],CLA[,1]==x))
	}))
	COLx<-subset(COLcla[,2], COLcla[,1]==CX)
	
	return(as.vector(c(ST,length(Ax), 
		mean(apply(Fx,1,sum)),
		mean(Fsite),CX,COLx)))
}))
row.names(COMP)<-NULL


#Unmatched diversity
COUNTa<-t(sapply(GROUPa,function(x){
	#print(x)
	Gx <-unlist(strsplit(x,split="_"))
	Sx<-intersect(STn, Gx)
	Ax<-intersect(ASVn, Gx)
	return(c(length(Sx),length(Ax)))
}))

GROUPs<-as.vector(subset(GROUPa, COUNTa[,1]!=0))
SUMs<-t(sapply(GROUPs,function(x){
	#print(x)
	Gx <-unlist(strsplit(x,split="_"))
	Sx<-intersect(STn, Gx)
	ST<-length(Sx)
	#Consensus clade
	CX<-as.vector(unique(sapply(Sx,function(x){
		return(subset(strains[,2], strains[,1]==x))
	})))
	COLx<-subset(COLcla[,2], COLcla[,1]==CX)	
	return(c(ST,CX,COLx))
}))

GROUPa<-as.vector(subset(GROUPa, COUNTa[,2]!=0))
SUMa<-t(sapply(GROUPa,function(x){
	#print(x)
	Gx <-unlist(strsplit(x,split="_"))
	Ax<-intersect(ASVn, Gx)
	Fx<-sapply(Ax,function(x){
		return(as.vector(subset(t(FRallR),
			paste("Meth",colnames(FRallR),sep="-")==x)))
	})
	#Consensus clade
	CX<-unique(sapply(Ax,function(x){
		return(subset(CLA[,5],CLA[,1]==x))
	}))
	COLx<-c(subset(COLcla[,2], COLcla[,1]==CX),"white")[1]
	return(c(length(Ax),mean(apply(Fx,1,sum)),CX,COLx))
}))
row.names(SUMa)<-NULL

SUMa<-t(sapply(sort(unique(SUMa[,3])),function(x){
	return(c(sum(subset(as.numeric(SUMa[,1]),SUMa[,3]==x)),
		sum(subset(as.numeric(SUMa[,2]),SUMa[,3]==x)),
		unique(subset(SUMa[,4],SUMa[,3]==x))))
}))

SUMa<-rbind(SUMa,cbind(sum(as.numeric(COMP[,2])),sum(as.numeric(COMP[,3])),"grey"))

SUMs <-t(sapply(unique(SUMs[,2]),function(x){
	return(c(sum(subset(as.numeric(SUMs[,1]), SUMs[,2]==x)),unique(subset(SUMs[,3], SUMs[,2]==x))))
}))


SUMs<-rbind(SUMs,cbind(sum(as.numeric(COMP[,1])),"grey"))

CLAu<-setdiff(sort(unique(c(COMP[,5],rownames(SUMa),rownames(SUMs)))),"")
CLAu<-sapply(CLAu,function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
})

#Summary figure

setwd(pD3)

pdf(paste("Comparison-ISO-ASV-pn=",pn,".pdf",sep=""),width=5,height=5)

par(mar=c(4,4,5,5),bty="n")

plot(as.numeric(COMP[,1]),as.numeric(COMP[,3]),log="xy",las=1,cex.axis=0.8,xlab="Methylobacterium strains",ylab="Methylobacterium diversity (barcoding)",pch=19,col=COMP[,6],lwd=1.5,cex=0.8)

par(mar=c(4,3,1,1),new=T)
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")

legend(x = "topright",cex=0.8, horiz=F,legend= names(CLAu),text.col= "black",box.col=NA,border=NA,pch=19,col= CLAu)


#Pie chart for unmatched strains	
	xs=0.85
	ys=0.1
	
	PCCc<-as.numeric(SUMs[,1])
	COLc<-as.vector(SUMs[,2])
	CUMc<-c(0,cumsum(PCCc))/sum(PCCc)
	sapply(c(1:length(PCCc)),function(x){
		COLx<-COLc[x]
		z1<-CUMc[x]
		z2<-CUMc[x+1]
		z12<-seq(z1,z2,length.out=100)
		ang12<-2*pi*(z12)
		x12<-c(xs,sin(ang12)*0.1+ xs)
		y12<-c(ys,cos(ang12)*0.1+ys)
		polygon(x12, y12,
			col= ifelse(COLx=="grey","white",COLx),
			border= ifelse(COLx=="grey","black",NA))	
	})	
	
	text(0.85,0.275,"Unmatched strains",cex=0.8)
	text(0.85,0.225,paste("(",round(100*(1-as.numeric(SUMs[length(SUMs[,1]),1])/sum(as.numeric(SUMs[,1]))),2),"%)",sep=""),cex=0.8)

#Pie chart for unmatched ASVs	
	xs=0.075
	ys=0.925
	
	PCCc<-as.numeric(SUMa[,2])
	COLc<-as.vector(SUMa[,3])
	CUMc<-c(0,cumsum(PCCc))/sum(PCCc)
	sapply(c(1:length(PCCc)),function(x){
		COLx<-COLc[x]
		z1<-CUMc[x]
		z2<-CUMc[x+1]
		z12<-seq(z1,z2,length.out=100)
		ang12<-2*pi*(z12)
		x12<-c(xs,sin(ang12)*0.1+ xs)
		y12<-c(ys,cos(ang12)*0.1+ys)
		polygon(x12, y12,
			col= ifelse(COLx=="grey","white",COLx),
			border= ifelse(COLx=="grey","black",NA))	
	})	

	text(0.375,0.95,"Unmatched diversity",cex=0.8)
	text(0.375,0.9,paste("(",round(100*(1-as.numeric(SUMa[length(SUMa[,1]),2])/sum(as.numeric(SUMa[,2]))),2),"%)",sep=""),cex=0.8)



dev.off()

SUMas<-t(sapply(names(CLAu),function(x){
	return(c(
	
	sum(as.numeric(subset(COMP[,1],COMP[,5]==x))),
	sum(as.numeric(subset(COMP[,2],COMP[,5]==x))),
	sum(as.numeric(subset(COMP[,3],COMP[,5]==x))),
	sum(as.numeric(subset(SUMs[,1],rownames(SUMs)==x))),
	sum(as.numeric(subset(SUMa[,1],rownames(SUMa)==x))),
	sum(as.numeric(subset(SUMa[,2], rownames(SUMa)==x)))))
}))

colnames(SUMas)<-c("Matched-Strains","Matched-ASVs",
	"Matched-Diversity","Unmatched-Strains",
	"Unmatched-ASVs","Unmatched-Diversity")

COMP<-cbind(COMP[,c(1,2,3,5)],GROUP)
colnames(COMP)<-c("Strains","ASV","Diversity","CLade","SeqNames")

write.table(SUMas,paste("Comparison-ISO-ASV-pn=",pn,".txt",sep=""))
write.table(COMP,paste("MATCH-ISO-ASV-pn=",pn,".txt",sep=""))


})
