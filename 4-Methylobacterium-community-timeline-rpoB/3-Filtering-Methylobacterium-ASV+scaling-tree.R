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
## Step 4 - Assign ASV to clades based on phylogeny and format metadata for next analyses
 #############################
#########################################
#########################################


library(ape)
library(seqinr)

pM="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis"

pO="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/6-Strains-2018"

pT<-"/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis/rTREE"

###Load data

	#strain (tips of the tree) names
	setwd(pT)
	
	#minimum bootstrapp value
	bx=0
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
	#Reference sequence informations (names in second column)
	setwd(pO)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	
pR="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/2-rarefaction&control"

pC="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/5-Comparative-community-analysis"

## Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs

setwd(pR)

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

rm(seqtab.nochim)

###Apply Hellinger transformation on relative abundance to account for rare ASVs
library(vegan)

FRallH <-as.data.frame(decostand(FRall ,method="hellinger", MARGIN=1))

colnames(FRallH)<-paste("ASV", c(1:length(FRallH[1,])),sep="")
colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="")

##Load Taxonomic assignation for each ASV
pA="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/3-Taxonomic-assignation"

setwd(pA)

taxaALL <-read.table(paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

##Simplify taxonomy (Emphasize on Order, and Familly within Rhizobiales)

ASSsimp<-ifelse(is.na(as.vector(taxaALL[,5]))==T,"unknown",ifelse(is.na(as.vector(taxaALL[,6]))==T,as.vector(taxaALL[,5]),ifelse(as.vector(taxaALL[,5])=="Rhizobiales",paste(taxaALL[,5], taxaALL[,6],sep="_"),as.vector(taxaALL[,5]))))

#Genus
GENx<-taxaALL[,7]

rm(taxaALL)

#Methylobacteriaceae ASVs
ASVm<-subset(c(1:length(ASSsimp)), ASSsimp=="Rhizobiales_Methylobacteriaceae")

rm(ASSsimp)

GENx <-as.vector(GENx[ASVm])

#ASV names (genus names 4 first digit + ASV number)

ASVn<-paste(as.vector(sapply(GENx,function(x){
	return(paste(unlist(strsplit(x,split=""))[1:4],collapse=""))
})),"ASV", ASVm,sep="-")

ASVi<-cbind(NA,"ASV",ASVn,GENx,"sp.","PHYLLOSPHERE",NA)
colnames(ASVi)<-colnames(CHAR)

CHAR<-rbind(CHAR,ASVi)

#tree information
setwd(pT)
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
	CHARx<-as.vector(t(subset(CHAR,CHAR[,3]== stx )))
	typex<-ifelse(CHARx[2]=="ASV","ASV","ref")	
	return(c(typex ,CHARx[c(3:7)]))
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
FRallH<-FRallH[SAMPo, ASVo]
FRall<-FRall[SAMPo, ASVo]
ASVs<-ASVs[ASVo]

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

#METAdata formated for permanova
FREQ.ENV<-as.data.frame(cbind(SITE,SS,SPE,TIME,EXT,PCR,MISEQ))

##For each ASV, calculate average relative abundance 
AAll<-apply(FRall,2,mean)/NS

##For each ASV, calculate relative abundance per SITE
ASite<-sapply(SITEu,function(x){
	ix<-subset(c(1:length(FRall[,1])),SITE==x)
	return(apply(FRall[ix,],2,mean)/NS)
})

##For each ASV, calculate relative abundance per PLOT
APlot<-sapply(SSu[c(6,1,2,3,4,5,7,8,10,9)],function(x){
	ix<-subset(c(1:length(FRall[,1])),SS==x)
	return(apply(FRall[ix,],2,mean)/NS)
})

##For each ASV, calculate relative abundance per Host tree specie per SITE
SiteSPEu<-sort(unique(paste(SITE,SPE,sep="-")))

ASpe<-sapply(SiteSPEu,function(x){
	ix<-subset(c(1:length(FRall[,1])),
		paste(SITE,SPE,sep="-")== x)
	return(apply(FRall[ix,],2,mean)/NS)
})	

##For each ASV, calculate relative abundance per Time point per SITE
SiteTIMEu<-sort(unique(paste(SITE,TIME,sep="-")))

ATime<-sapply(SiteTIMEu,function(x){
	ix<-subset(c(1:length(FRall[,1])),
		paste(SITE,TIME,sep="-")== x)
	return(apply(FRall[ix,],2,mean)/NS)
})

ACladex<-as.numeric(matrix(unlist(strsplit(CLA[,1],split="-")),ncol=3,byrow=T)[,3])

AClade<-as.vector(sapply(colnames(FRall),function(x){
	Ax<-as.numeric(unlist(strsplit(x,split="ASV")))[2]
	return(subset(CLA[,5], ACladex==Ax))
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
		cbind("A5","yellow3",3,3,0.5),
		cbind("A7","pink",3,3,0.5),
		cbind("A10","green3",2,7,0.8),
		cbind("A8","purple",3,3,0.5))
		
COLcla<-sapply(AClade,function(x){
	return(c(subset(COLcla[,2], COLcla[,1]==x),"white")[1])
})		

barplot(t(t(ASite)/apply(ASite,2,sum)),las=2,col= COLcla,border=NA)
barplot(t(t(APlot)/apply(APlot,2,sum)),las=2,col= COLcla,border=NA)
barplot(t(t(ATime)/apply(ATime,2,sum)),las=2,col= COLcla,border=NA)
barplot(t(t(ASpe)/apply(ASpe,2,sum)),las=2,col= COLcla,border=NA)

#########################################
## STEP 5 - for each level in the tree redefine taxonomic divisions according to significant bootstrapp in the tree and perform PERMANOVA on metadata,   ################



#Number of permutations for PERMANOVA
K=10000
	
#Tree topology in level format
LEVELS<-read.table(paste("Levels_MINBOOT=",
	bx,".txt",sep=""))
colnames(LEVELS)<-paste("L",
	c(1:(length(LEVELS[1,]))),sep="-")
	
#remove empty levels
LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))
	
#for each sequence, fill empty levels with node from the closest level
	
LEVELf<-t(sapply(c(1:length(LEVELS[,1])),function(x){
	Sx=x
	Lx<-as.vector(LEVELS[Sx,])
	lx<-c(1:length(Lx))
	Lx<-sapply(lx,function(x){
		return(max(subset(Lx,lx>=x)))
	})
	return(c(Sx,Lx[2:length(Lx)]))
}))

colnames(LEVELf)<-colnames(LEVELS)

#Only Keep levels for ASVs
LEVELf<-LEVELf[ASVi,]
	
###For each level, perform a PERMANOVA

###Permanova on diversity per level
LV=c(2:length(LEVELf[1,]))

#Number of ASVs per level 
STn<-sapply(LV,function(x){
	#print(x)
	Lx=x
	return(length(unique(LEVELf[, Lx])))
})	

#Remove levels with less than 2 ASVs
LV2<-subset(LV, STn>=2)

#####################
### PERMANOVA PER LEVEL


AOV<-sapply(LV2,function(x){
	print(x)
	Lx=x
	
	#Level
	LEVELx<-as.numeric(unlist(strsplit(colnames(LEVELf)[x],split="-"))[2])
	
	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	NDu<-setdiff(unique(ND),0)
	
	####For each node, retrieve ASV abundance (From FRall)
	N1<-sapply(NDu,function(x){
		NDx=x
		#ASV index in this node
		STx<-subset(c(1:length(LEVELf[,1])),ND== NDx)
		#total abundance per sample for these ASVs
		FRx<-cbind(FRall[,STx],0)
		return(apply(FRx,1,sum))
	})
	colnames(N1)<-paste("node", NDu,sep="-")
	
	#remove nodes without observation
	N1<-as.data.frame(t(subset(t(N1),apply(N1,2,max)!=0)))
	
	#Hellinger transformation
	N1 <-as.data.frame(decostand(N1 ,method="hellinger", MARGIN=1))
	
				
	MOD<-adonis(N1 ~ SITE*SS*SPE*TIME, data= FREQ.ENV, permutations=K,strata= FREQ.ENV $MISEQ,method="bray")$aov
	
	return(c(dim(N1), MOD[,5], MOD[,6]))
	
})

DESC<-rownames(MOD)


rownames(AOV)<-c("Fac","Nodes",paste("R2", DESC),paste("Pr(>F)", DESC))
colnames(AOV)<-colnames(LEVELf)[LV2]

#Remove levels with less than two nodes
AOV<-t(subset(t(AOV),AOV[2,]>1))

setwd(pM)

write.table(AOV,paste("Permanova-tree-min-bootstrap=",bx,".txt",sep=""))



####Display result

#Levels
LV<-as.numeric(matrix(unlist(strsplit(colnames(AOV),split="-")),ncol=2,byrow=T)[,2])

#Levels for x-scale
LVscale<-round(seq(min(LV),max(LV),length.out=6))

#Corresponding PS values to show on the scale
PSscale<-round(sapply(LVscale,function(x){
	return(subset(SCALE[,4],SCALE[,2]==x))
}),2)

#Color code for factors
COLfac<-c("blue","grey","red","yellow2",
	"cyan2","purple","red3","green",
	"yellow4","orange","purple3","green3",
	"brown","orange3","black","grey",
	NA)


pdf(paste("Permanova-tree-min-bootstrap=",bx,".pdf",sep=""),width=5,height=4)

par(mar=c(4,4,1,1),bty="n")

plot(-10,-10,xlim=c(max(LV),min(LV)),ylim=c(0,1.25*max(AOV[c(3:11),])),xlab="PS",ylab="Part of variance",las=1,cex.axis=0.8,xaxt="n")

axis(1, LVscale, PSscale,cex.axis=0.8)

PMx<-sapply(c(1:(length(DESC)-1))+2,function(x){
	
	colx<-COLfac[x-2]
	
	DESCx<-rownames(AOV)[x]
	
	R2x<-as.vector(AOV[x,])
	px=as.vector(AOV[x+length(DESC),])
	pm=min(px)
	colm<-ifelse(DESCx=="R2 Residuals","grey",ifelse(pm>0.05,NA, colx))
	
	points(LV,R2x,type="l",col= colm,lty=ifelse(DESCx=="R2 Residuals",3,3))
	
	#colma =ifelse(is.na(colm)==T,NA,rgb(col2rgb(colm)[1]/255,
	#		col2rgb(colm)[2]/255,
	#		col2rgb(colm)[3]/255,0.3))
	
	points(LV,R2x,pch=19,
		cex=ifelse(px<=0.001,1,
		ifelse(px<=0.01,0.75,
		ifelse(px<=0.05,0.5,
		ifelse(px<=0.1,0,0)))),bg= ifelse(is.na(colm)==T,NA,"white"),col= colm)	
	#points(LV,R2x,pch=21,
	#	cex=ifelse(px<=0.001,1,
	#	ifelse(px<=0.01,0.75,
	#	ifelse(px<=0.05,0.5,
	#	ifelse(px<=0.1,0.25,0)))),bg= colma,col= colm)
	
	return(pm)		
})

#Legend only for significant factors
DESCleg<-c(subset(DESC[c(1:(length(DESC)-1))],PMx<=0.05),"residuals")
COLleg<-c(subset(COLfac[c(1:(length(DESC)-1))],PMx<=0.05),"grey")


legend(x = "topleft",cex=0.4, horiz=F,legend= DESCleg,text.col= COLleg,box.col=NA,pch=0,border=NA)

legend(x = "topright",cex=0.6, horiz=F,legend= c(
	expression("p"<= "0.001"),
	expression("p"<= "0.01"),
	expression("p"<= "0.05")),
	#expression("p"<= "0.1")),
	border=T,box.col=NA,pch=21,pt.cex=c(1,0.75,0.5,0.25))


dev.off()


################################################
###################################################
### STEP 6: summary of Methylobacterium ASVs and comparison with isolated strains

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

pM="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis"

pO="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/6-Strains-2018"

pT<-"/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis/rTREE"

###Load data

	#strain (tips of the tree) names
	setwd(pT)
	
	#minimum bootstrapp value
	bx=0
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
	#Reference sequence informations (names in second column)
	setwd(pO)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	
pR="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/2-rarefaction&control"

pC="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/5-Comparative-community-analysis"

## Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs

setwd(pR)

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
pA="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/3-Taxonomic-assignation"

setwd(pA)

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

ASVi<-cbind(NA,"ASV",ASVn,GENx,"sp.","PHYLLOSPHERE",NA)
colnames(ASVi)<-colnames(CHAR)

CHAR<-rbind(CHAR,ASVi)

#tree information
setwd(pT)
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
	CHARx<-as.vector(t(subset(CHAR,CHAR[,3]== stx )))
	typex<-ifelse(CHARx[2]=="ASV","ASV","ref")	
	return(c(typex ,CHARx[c(3:7)]))
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
		cbind("A5","yellow3",3,3,0.5),
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
	
	setwd(pM)
	
	FRsum<-cbind(FclaSum, FRsite)
	colnames(FRsum)<-c("ASVs","Fr","MSH","SBL")
	
	write.table(FRsum,"Summary-Methylobacteriaceae-ASV.txt")
	
#Fasta file combining aligned strains rpoB sequences and ASVs (nnn removed)

pO="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/6-Strains-2018"

setwd(pO)
strains<-read.table("strain-metadata.txt",header=T)

library(seqinr)
setwd(pM)
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
	
	return(c(ST,length(Ax), mean(apply(Fx,1,sum)),
		mean(Fsite),CX,COLx))
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

setwd(pM)

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
##############################
#####Meta analysis of Methylobacterium ASVs with emphazise on time


library(ade4)
library(vegan)

setwd(pM)

pdf("PCA-Methylo-SUMMARY-TAXA.pdf",height=7.5,width=7.5)

zones<-matrix(c(1,3,2,4),ncol=2)
layout(zones)
#layout.show(max((zones)))
par(bty="n")

CLAu<-sapply(sort(unique(CLA[,5])),function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
})

###rTREE with nodes labelled with PERMANOVA part of variance

#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=1.2

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

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"a",cex=1.5,font=2)

par(mar=c(1,1,1,1),bty="n",new=T)

#Frame
plot(-10000,-10000,xlim=c(-rLIM, rLIM),
	ylim=c(-rLIM, rLIM),xlab="",
	ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n")


#Remove nodes after strain N228 (last Methylobacterium)

nmx<- subset(c(1:length(type[,1])),type[,2]=="NS228")


nmcx<-subset(c(1:length(coordx[,1])),coordx[,2]<= nmx)
	
MOx<-subset(coordx[nmcx,6],coordx[nmcx,5]==max(coordx[nmcx,5]))	
	
nmcx<-c(nmcx,subset(c(1:length(coordx[,1])), coordx[,4]== MOx))	
	
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
PSscale<-PSscale[c(1:length(Xscale))]
segments(-0.01*rLIM,min(Xscale),
	-0.01*rLIM,max(Xscale))
text(-0.005*rLIM, Xscale, round(PSscale,2)
	,cex=0.6,pos=2,adj=1,font=1)
text(-0.25*rLIM, mean(Xscale), "PS",cex=0.8,font=1)
	
maxTIP=max(Xscale)

#PART of variance explained by time in an anova for each node with at least 50% of support and ASVs

#1 ASVs
	
Tmx<-intersect(subset(Tx, type[,1]=="ASV"),subset(Tx, type[,3]=="Methylobacterium"))

CA<-sapply(Tmx,function(x){
	nx<-type[x,2]
	cx=subset(CLA[,5],CLA[,1]== nx)
	return(cx)
})


#Nodes (supported by at least 30% of bootstrapps)

setwd(pT)
#Tree topology in level format
LEVELS<-read.table(paste("Levels_MINBOOT=",
	0,".txt",sep=""))
colnames(LEVELS)<-paste("L",
	c(1:(length(LEVELS[1,]))),sep="-")
	
#remove empty levels
LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))

setwd(pM)

#Only nodes with at least 30% of boot support
Nbx<-subset(PSn[,2],PSn[,8]>=0.3)

#format metadata for permanova
METAp<-as.data.frame(cbind(SITE,SS,SPE,TIME,MISEQ))
	
#Number of ASV per nodes: only keep nodes with at least two ASVs (for PERMANOVA)
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

#Count the number of samples and variable for each nodes

VRcount<-t(sapply(Nbx,function(x){
	#print(x)
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
	nx= intersect(type[STx,2], ASVn)
	fx<-as.data.frame(sapply(nx,function(x){
		return(subset(t(FRallH),
			paste("Meth",colnames(FRallH),sep="-")==x))
	}))
	METAx<-subset(METAp,apply(fx,1,sum)!=0)
	fxx<-subset(fx,apply(fx,1,sum)!=0)
	return(c(NX ,length(fxx[,1]),
		sapply(c(1:length(METAx[1,])),function(x){
			return(length(unique(METAx[,x])))
		})))
}))

#Remove nodes that have only one variable in at least one factor
VRcount<-subset(VRcount,apply(VRcount,1,min)>1)

#Remove nodes that are present in less than nsamp samples
nsamp=11
VRcount<-subset(VRcount, VRcount[,2]>= nsamp)
Nbx <-VRcount[,1]

ANO=4

#For permanova, perform K permutations
K=10000

FACad<-c("SITE","SS","SPE","TIME","SITE:SPE","SS:SPE","SS:TIME","SPE:TIME","SS:SPE:TIME","Residuals","Total")

ANOnode<-t(sapply(Nbx,function(x){
	print(x)
	#Node
	NX=x
	#Level
	LEVELx <-subset(PSn[,1],PSn[,2]==NX)		
	#Column index in the level file
	Lx<-subset(c(1:length(LEVELS[1,])),paste("L", LEVELx,sep="-")==colnames(LEVELS))
	#Nodes in this level
	ND<-as.vector(LEVELS[, Lx])
	#ASV line indexes in the level file for this node
	STx<-intersect(Tmx,
		subset(c(1:length(LEVELS[,1])),ND== NX))
	#Corresponding ASV names
	nx= intersect(type[STx,2], ASVn)

	fx<-as.data.frame(sapply(nx,function(x){
		return(subset(t(FRallH),
			paste("Meth",colnames(FRallH),sep="-")==x))
	}))
	
	METAx<-subset(METAp,apply(fx,1,sum)!=0)
	fxx<-subset(fx,apply(fx,1,sum)!=0)

	Mt<-adonis(fxx ~ SITE*SS*SPE*TIME, data= METAx, permutations=K,strata= METAx$MISEQ,method="bray")$aov

	#write.table(as.data.frame(Mt),paste("ANOVA-node-",NX,".txt",sep=""))

	VARx<-Mt[,5]
	Px<-Mt[,6]
	
	#Mt<-lm(fx ~ SITE*TIME*SS*SPE)	
	#At<-anova(Mt)
	VARtx<-subset(VARx,row.names(Mt)=="TIME")
	Ptx<-subset(Px,row.names(Mt)=="TIME")
	VARsx<-subset(VARx,row.names(Mt)=="SITE")
	Psx<-subset(Px,row.names(Mt)=="SITE")
	VARhx<-subset(VARx,row.names(Mt)=="SPE")
	Phx<-subset(Px,row.names(Mt)=="SPE")
	VARssx<-subset(VARx,row.names(Mt)=="SS")
	Pssx<-subset(Px,row.names(Mt)=="SS")		
	
	VARtx<-ifelse(Ptx<=0.05, VARtx,0)
	VARsx<-ifelse(Psx<=0.05, VARsx,0)
	VARhx<-ifelse(Phx<=0.05, VARhx,0)
	VARssx<-ifelse(Psx<=0.05, VARssx,0)
	
	varx<-cumsum(c(VARtx, VARhx , VARssx ,VARsx))
	VARtx<-varx[1]
	VARhx <-varx[2]	
	VARssx <-varx[3]
	VARsx <-varx[4]
	
	#Consensus clade
	ncx<-unique(sapply(nx,function(x){
		return(subset(CLA[,5],CLA[,1]==x))
	}))
	cx<-ncx[1]
	ncx=length(ncx)
	cx=ifelse(ncx==1,cx,"un")
	colx<-subset(CLAu,names(CLAu)== cx)	

	colH<-rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255,0.3)
	colSS<-rgb(col2rgb(colx)[1]/255,
			col2rgb(colx)[2]/255,
			col2rgb(colx)[3]/255,0.1)	

	COORDxx<-subset(COORD,coordx[nmcx,4]==NX)
	XX=COORDxx[1,2]
	YY=COORDxx[1,3]	
	
	points(XX, YY,pch=1,cex= sqrt(abs(VARsx)/pi)* ANO,
		col= colx,lwd=0.25,bg="white")
	points(XX, YY,pch=19,cex=sqrt(abs(VARssx)/pi)* ANO,
		col= colSS,lwd=0.25,bg="white")	
	points(XX, YY,pch=19,cex=sqrt(abs(VARhx)/pi)* ANO,
		col= colH,lwd=0.25,bg="white")			
	points(XX, YY,pch=19,cex=sqrt(abs(VARtx)/pi)* ANO,
		col= colx,lwd=0.25,bg="white")		
	
	VARx<-sapply(FACad,function(x){
		return(c(subset(Mt[,5],rownames(Mt)==x),NA)[1])
	})
	Px<-sapply(FACad,function(x){
		return(c(subset(Mt[,6],rownames(Mt)==x),NA)[1])
	})
	
	return(c(NX,LEVELx,length(nx), cx,min(STx),max(STx), VARx,Px))
	
}))

colnames(ANOnode)<-c("Node","Level","ASVs","Clade","FistASV","LastASV",
paste("VAR", FACad,sep="."),paste("pval", FACad,sep="."))

write.table(ANOnode,"PERMANOVA-per-node.txt",row.names=F)
ANOnode<-read.table("PERMANOVA-per-node.txt",header=T)

#Consensus clade label

sapply(setdiff(COLcla[,1],c("A5","Enterovirga","Microvirga")),
	function(x){
	colx=subset(COLcla[,2], COLcla[,1]==x)
	ST1<-min(c(subset(Tmx,CA==x),subset(Tx, type[,6]==x)))
	ST2<-max(c(subset(Tmx,CA==x),subset(Tx, type[,6]==x)))
	nax<-length(subset(Tmx,CA==x))
	angx=2*pi*(mean(c(ST1, ST2))/z)*OPEN
	XX1i<-Ae[3]* maxTIP*sin(angx)
	YY1i<-Ae[3]* maxTIP*cos(angx)						
	text(XX1i, YY1i,x,cex=ifelse(nax==0,0.4,0.8),col= colx,font=1)
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
		
		points(1.02*maxTIP*sin(angx), 1.02*maxTIP*cos(angx),pch=19,cex=0.2,col=ifelse(type[x,1]=="ASV", cx,NA))
		
		text(XX1i, YY1i,nx,cex=0.2,srt= ang.lab[x],adj= adjx[x],col=cx,font=1)	
			
	})

par(new=T,mar=c(0,0,0.2,0.5))

VP<-c(0.04,0.1,0.2,0.4)

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend("topright",pt.cex =sqrt(VP/pi)*ANO,cex=0.5, horiz=T,legend= VP,text.col= "black",box.col=NA,border=NA,pch=19,col= "black", title="Cumulative part of var. (ANOVA)", pt.lwd=0.25)

par(new=T,mar=c(0,0,1.5,0.5))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend("topright",pt.cex =sqrt(max(VP)/pi)*ANO,cex=0.5, horiz=T,legend= c("Site","Plot","Host","Time"),text.col= "black",box.col=NA,border=NA,pch=c(1,19,19,19),col= c("black",rgb(0,0,0,0.1),rgb(0,0,0,0.3),"black"), title="Factors", pt.lwd=0.25)

###Global PCA

ACPh=dudi.pca(FRallH, scannf=F)
part<-round(100* ACPh $eig/sum(ACPh $eig),2)

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"b",cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

x1=min(ACPh $li[,1])
x2=max(ACPh $li[,1])
y1=min(ACPh $li[,2])
y2=max(ACPh $li[,2])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPh $li[,c(1,2)],groups= SITE,col= rgb(0.9,0.9,0.9),lty=1,lwd=0.5)	

points(ACPh $li[,c(1,2)],pch=19,cex=0.4,col=rgb(0,0,0,0.6))

ordiellipse(ACPh $li[,c(1,2)],groups= SITE,border="NA",alpha=40,draw="polygon",col= "grey",lty=1,label=T,cex=0.5,font=1)	

#ASV contribution colored according to clades
CONT<-ACPh$co

COLcont<-as.vector(sapply(CLA[,5],function(x){
	return(c(subset(COLcla[,2],COLcla[,1]==x),"grey")[1])
}))

#Estimate amd part of variance Sie/Fr per ASV
COR<-t(sapply(c(1:length(FRallH[1,])),function(x){
	fx<-FRallH[,x]
	Mt<-lm(fx ~ SITE*TIME*SS*SPE)	
	At<-anova(Mt)
	return(c(
		subset(At[,2]/sum(At[,2]),row.names(At)=="SITE"),
		subset(At[,5],row.names(At)=="SITE")
	))
	
}))

colnames(COR)<-c("AN.var","AN.pval")

pADJ<-p.adjust(COR[,2])


#Center of the contribution plot in the main plot
Xc=-2.5
Yc=-12.5

#Contribution amplification
AC=7

#Size of factor contribution tested by anova
ANO=2

text(Xc,Yc+4,"ASV contribution & site association",cex=0.5,font=2)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ <=0.05,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=19,cex= ifelse(pADJ <=0.05,sqrt(abs(COR[,1])/pi)*ANO,0),col= COLcont,lwd=0.75,bg="white")

par(new=T,mar=c(4,5,1,1))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend(x = "bottomleft",cex=0.5, horiz=F,legend= names(CLAu),text.col= "black",box.col=NA,border=NA,pch=19,col= CLAu, title="Clades")

par(new=T,mar=c(4,7,1,1))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

VL<-c(0.04,0.1,0.2,0.4)

legend(x = "bottomleft",pt.cex =sqrt(VL/pi)*ANO,cex=0.5, horiz=F,legend= VL,text.col= "black",box.col=NA,border=NA,pch=19,col= "black", title="Part of var. Site (ANOVA)")




#PCA per site (return ASVs significantly associated with time)

SUMtime<-matrix(unlist(sapply(SITEu,function(x){
	X=x
	
	FRx<-subset(FRallH,SITE==X)
	SSx<-subset(SS,SITE==X)
	IDx<-subset(ID,SITE==X)
	SSx=matrix(unlist(strsplit(SSx,split="-")),
		ncol=2,byrow=T)[,2]
		
	SPEx<-subset(SPE,SITE==X)
	TIMEx<-subset(TIME,SITE==X)
	TIMErx<-subset(TIMEr,SITE==X)
	
	
	#Display PCA
	ACPx=dudi.pca(FRx, scannf=F)
	part<-round(100* ACPx $eig/sum(ACPx $eig),2)
	#Graphic limit
	x1=min(ACPx $li[,1])
	x2=max(ACPx $li[,1])
	y1=min(ACPx $li[,2])
	y2=max(ACPx $li[,2])
	
	par(mar=c(0,0,0,0))
		plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,ifelse(X=="MSH","c","d"),cex=1.5,font=2)
	par(new=T,mar=c(4,4,1,1))
	
	plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=0.8)
		
	ordispider(ACPx $li[,c(1,2)],groups= TIMErx,
		col= "grey",lty=1,lwd=0.5,label=F)
	points(ACPx $li[,c(1,2)],pch=19,cex=0.4,col="black")
	ordiellipse(ACPx $li[,c(1,2)],groups= TIMErx,border="NA",alpha=40,draw="polygon",col="black",lty=1,label=T,cex=0.5,font=1)

#Estimate amd part of variance Time/Fr per ASV
COR<-t(sapply(c(1:length(FRx[1,])),function(x){
	fx<-FRx[,x]
	Mt<-lm(fx ~ TIMEx*SSx*SPEx)	
	At<-anova(Mt)
	Ct<-cor.test (TIMEx, fx)
	return(c(
		subset(At[,2]/sum(At[,2]),row.names(At)=="TIMEx"),
		subset(At[,5],row.names(At)=="TIMEx"),
		Ct$estimate
	))
	
}))

colnames(COR)<-c("AN.var","AN.pval","Est")

pADJ<-p.adjust(COR[,2])
	
#ASV contribution colored according to clades
CONT<-ACPx$co

#Center of the contribution plot in the main plot
Xc=ifelse(X=="SBL",-15,5)
Yc=ifelse(X=="SBL",-10,7.5)

#Contribution amplification
AC=7

ANO=3

text(Xc,Yc+ifelse(X=="SBL",3,-3),"ASV contribution & time shift",cex=0.5,font=2)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ <=0.05,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=ifelse(COR[,3]<0,21,19),cex= ifelse(pADJ <=0.05,sqrt(abs(COR[,1])/pi)* ANO,0),col= COLcont,lwd=0.75,bg="white")
	
VL<-c(0.02,0.04,0.1,0.2)

text(ifelse(X=="SBL",-20,-7.5),ifelse(X=="SBL",4,11),X,font=2,cex=0.8)

par(new=T,mar=c(
ifelse(X=="SBL",5,4),
ifelse(X=="SBL",7,4),
1,
ifelse(X=="SBL",1,2)
))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")


legend(x = ifelse(X=="SBL","bottomleft","topright"),pt.cex =sqrt(VL/pi)*ANO,cex=0.5, horiz=T,legend= VL,text.col= "black",box.col=NA,border=NA,pch=19,col= "black", title="Part of var. Time (ANOVA)")

return(t(cbind(x,subset(cbind(CLA,pADJ),pADJ<=0.05))))
			
})),ncol=8,byrow=T)


dev.off()

#####ASV significantly associated with time: display time dynamics separatly in MSH and SBL

SUMtime<-SUMtime[order(SUMtime[,6]),]

ASVtu<-unique(SUMtime[,2])

pdf("ASV-Sign-Increase-Time.pdf",height=5,width=5)

par(mar=c(1,3,1,1),mfrow=c(5,5),bty="n")

FRt<-apply(sapply(ASVtu,function(x){
	
	CLAx<-subset(CLA[,5],CLA[,1]==x)
	COLx<-subset(COLcla[,2],COLcla[,1]== CLAx)
	COLx2<-rgb(col2rgb(COLx)[1]/255,
			col2rgb(COLx)[2]/255,
			col2rgb(COLx)[3]/255,0.2)
	
	
	FRx<-as.numeric(subset(t(FRallH),
		paste("Meth",colnames(FRallH),sep="-")==x))
	
	plot(-10,-10,xlim=c(min(TIME),max(TIME)),
		ylim=c(0,max(FRx)),xaxt="n",xlab="time",
		ylab="Rel. ab.",las=1,cex.axis=0.7,log="")
	
	sapply(SITEu,function(x){
		Tx<-subset(TIME, SITE==x)
		Fx<-subset(FRx, SITE==x)
		Fm<-sapply(sort(unique(Tx)),function(x){
			return(mean(subset(Fx, Tx==x)))
		})
		
		points(Tx,Fx,pch=ifelse(x=="MSH",19,1),cex=0.1,col= COLx2,lwd=0)
		points(sort(unique(Tx)), Fm,type="l",lty=ifelse(x=="MSH",1,2),col=COLx,lwd=0.5)
	})

	SUMx<-subset(SUMtime, SUMtime[,2]==x)
	legend(x = "top",pt.cex =0,cex=0.4, 
		horiz=F,legend= paste(SUMx[,1],"; p=",
		round(as.numeric(SUMx[,8]),4)),
		text.col= "black",box.col=NA,border=NA,title=x)	
	
	return(as.numeric(subset(t(FRallR),
		paste("Meth",colnames(FRallR),sep="-")==x)))
		
		
}),1,sum)

dev.off()

par(mfrow=c(1,2))
sapply(SITEu,function(x){
	ASVt<-subset(SUMtime[,2],SUMtime[,1]==x)
	SITEx=x
	
	FRt<-t(sapply(ASVt,function(x){
		FRx<-as.vector(subset(t(FRallR),
			paste("Meth",colnames(FRallR),sep="-")==x))
		Tx<-subset(TIME,SITE==SITEx)	
		FRx <-subset(FRx,SITE==SITEx)		
		return(sapply(sort(unique(Tx)),function(x){
			return(mean(subset(FRx, Tx ==x)))
		}))
	}))

	CLAt<-t(sapply(ASVt,function(x){
		CLAx<-subset(CLA[,5],CLA[,1]==x)
		COLx<-subset(COLcla[,2],COLcla[,1]== CLAx)
		COLx2<-rgb(col2rgb(COLx)[1]/255,
			col2rgb(COLx)[2]/255,
			col2rgb(COLx)[3]/255,0.2)
		return(c(CLAx, COLx, COLx2))
	}))
	
	barplot(FRt,col= CLAt[,3],border=CLAt[,2])
})	

#Supplementary figure: ASV associated with plots and host tree species

pdf("Sup-PCA-Methylo-SUMMARY-TAXA.pdf",height=7.5,width=7.5)

zones<-matrix(c(1,3,2,4),ncol=2)
layout(zones)
#layout.show(max((zones)))
par(bty="n")



#PCA per Host tree species

sapply(SITEu,function(x){
	X=x
	
	FRx<-subset(FRallH,SITE==X)
	SSx<-subset(SS,SITE==X)
	IDx<-subset(ID,SITE==X)
	SSx=matrix(unlist(strsplit(SSx,split="-")),
		ncol=2,byrow=T)[,2]
		
	SPEx<-subset(SPE,SITE==X)
	TIMEx<-subset(TIME,SITE==X)
	TIMErx<-subset(TIMEr,SITE==X)
	
	#Display PCA
	ACPx=dudi.pca(FRx, scannf=F)
	part<-round(100* ACPx $eig/sum(ACPx $eig),2)
	#Graphic limit
	x1=min(ACPx $li[,1])
	x2=max(ACPx $li[,1])
	y1=min(ACPx $li[,2])
	y2=max(ACPx $li[,2])
	
	par(mar=c(0,0,0,0))
		plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,ifelse(X=="MSH","a","b"),cex=1.5,font=2)
	par(new=T,mar=c(4,4,1,1))
	
	plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=0.8)
		
	ordispider(ACPx $li[,c(1,2)],groups= SPEx,
		col= "grey",lty=1,lwd=0.5,label=F)
	points(ACPx $li[,c(1,2)],pch=19,cex=0.4,col="black")
	ordiellipse(ACPx $li[,c(1,2)],groups= SPEx,border="NA",alpha=40,draw="polygon",col="black",lty=1,label=T,cex=0.5,font=1)

#Estimate amd part of variance Spe/Fr per ASV
COR<-t(sapply(c(1:length(FRx[1,])),function(x){
	fx<-FRx[,x]
	Mt<-lm(fx ~ TIMEx*SSx*SPEx)	
	At<-anova(Mt)
	return(c(
		subset(At[,2]/sum(At[,2]),row.names(At)=="SPEx"),
		subset(At[,5],row.names(At)=="SPEx")	))
	
}))

colnames(COR)<-c("AN.var","AN.pval")

pADJ<-p.adjust(COR[,2])
	
#ASV contribution colored according to clades
CONT<-ACPx$co

#Center of the contribution plot in the main plot
Xc=ifelse(X=="SBL",-15,5)
Yc=ifelse(X=="SBL",-10,7.5)

#Contribution amplification
AC=7

ANO=3

text(Xc,Yc+ifelse(X=="SBL",3,-3),"ASV contribution & tree species",cex=0.5,font=2)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ <=0.05,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=19,cex= ifelse(pADJ <=0.05,sqrt(abs(COR[,1])/pi)* ANO,0),col= COLcont,lwd=0.75,bg="white")
	
VL<-c(0.02,0.04,0.1,0.2)

text(ifelse(X=="SBL",-20,-7.5),ifelse(X=="SBL",4,11),X,font=2,cex=0.8)

par(new=T,mar=c(
ifelse(X=="SBL",5,4),
ifelse(X=="SBL",7,4),
1,
ifelse(X=="SBL",1,2)
))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend(x = ifelse(X=="SBL","bottomleft","topright"),pt.cex =sqrt(VL/pi)*ANO,cex=0.5, horiz=T,legend= VL,text.col= "black",box.col=NA,border=NA,pch=19,col= "black", title="Part of var. Time (ANOVA)")

})


#PCA per Plot

sapply(SITEu,function(x){
	X=x
	
	FRx<-subset(FRallH,SITE==X)
	SSx<-subset(SS,SITE==X)
	IDx<-subset(ID,SITE==X)
	SSx=matrix(unlist(strsplit(SSx,split="-")),
		ncol=2,byrow=T)[,2]
		
	SPEx<-subset(SPE,SITE==X)
	TIMEx<-subset(TIME,SITE==X)
	TIMErx<-subset(TIMEr,SITE==X)
	
	#Display PCA
	ACPx=dudi.pca(FRx, scannf=F)
	part<-round(100* ACPx $eig/sum(ACPx $eig),2)
	#Graphic limit
	x1=min(ACPx $li[,1])
	x2=max(ACPx $li[,1])
	y1=min(ACPx $li[,2])
	y2=max(ACPx $li[,2])
	
	par(mar=c(0,0,0,0))
		plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
	text(0.05,0.95,ifelse(X=="MSH","c","d"),cex=1.5,font=2)
	par(new=T,mar=c(4,4,1,1))
	
	plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=0.8)
		
	ordispider(ACPx $li[,c(1,2)],groups= SSx,
		col= "grey",lty=1,lwd=0.5,label=F)
	points(ACPx $li[,c(1,2)],pch=19,cex=0.4,col="black")
	ordiellipse(ACPx $li[,c(1,2)],groups= SSx,border="NA",alpha=40,draw="polygon",col="black",lty=1,label=T,cex=0.5,font=1)

#Estimate amd part of variance Spe/Fr per ASV
COR<-t(sapply(c(1:length(FRx[1,])),function(x){
	fx<-FRx[,x]
	Mt<-lm(fx ~ TIMEx*SSx*SPEx)	
	At<-anova(Mt)
	return(c(
		subset(At[,2]/sum(At[,2]),row.names(At)=="SSx"),
		subset(At[,5],row.names(At)=="SSx")	))
	
}))

colnames(COR)<-c("AN.var","AN.pval")

pADJ<-p.adjust(COR[,2])
	
#ASV contribution colored according to clades
CONT<-ACPx$co

#Center of the contribution plot in the main plot
Xc=ifelse(X=="SBL",-15,5)
Yc=ifelse(X=="SBL",-10,7.5)

#Contribution amplification
AC=7

ANO=3

text(Xc,Yc+ifelse(X=="SBL",3,-3),"ASV contribution & plot association",cex=0.5,font=2)

points(Xc,Yc,pch=10,cex= AC,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,col= ifelse(pADJ <=0.05,"grey",NA),lwd=0.75,lend=2)

points(Xc+ CONT[,1]* AC,Yc+CONT[,2]* AC,pch=19,cex= ifelse(pADJ <=0.05,sqrt(abs(COR[,1])/pi)* ANO,0),col= COLcont,lwd=0.75,bg="white")
	
VL<-c(0.02,0.04,0.1,0.2)

text(ifelse(X=="SBL",-20,-7.5),ifelse(X=="SBL",4,11),X,font=2,cex=0.8)

par(new=T,mar=c(
ifelse(X=="SBL",5,4),
ifelse(X=="SBL",7,4),
1,
ifelse(X=="SBL",1,2)
))

plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")

legend(x = ifelse(X=="SBL","bottomleft","topright"),pt.cex =sqrt(VL/pi)*ANO,cex=0.5, horiz=T,legend= VL,text.col= "black",box.col=NA,border=NA,pch=19,col= "black", title="Part of var. Time (ANOVA)")

})


dev.off()


	
###########
###########
###########
###########
###########
###########OLD

ACPh $co




setwd(pM)

pn=0.018

COMP<-read.table(paste("MATCH-ISO-ASV-pn=",pn,".txt",sep=""),header=T)


#Comparison of MILAG values for strains and time dynamics for equivalent ASV

Pt="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/8-temperature-screen"

setwd(pM)

pn=0.018

COMP<-read.table(paste("MATCH-ISO-ASV-pn=",pn,".txt",sep=""),header=T)

setwd(Pt)

MiLag<-read.table("MiLag.txt",header=T)
ST<-rownames(MiLag)
ST<-ifelse(ST=="E-061-B","E-061",ifelse(ST=="E-061-R","E-061",ST))

#Only keep Groups with data in both conparable ASV and temperature screen

GROUP<-as.vector(COMP[,5])

nST<-as.vector(sapply(GROUP,function(x){
	return(length(intersect(unlist(strsplit(x,split="_")), ST)))
	
}))

GROUP<-subset(GROUP, nST!=0)

#Metadata for each strain

STa<-intersect(ST,unlist(strsplit(GROUP,split="_")))

MetaS<-t(sapply(STa,function(x){
	return(t(subset(strains, strains[,1]==x))[,1])
}))

SSt<-as.vector(sapply(MetaS[,4],function(x){
	SSx<-unlist(strsplit(x,split="-"))
	return(paste(SSx[1],as.numeric(SSx[2]),sep="-"))
}))

TIMEt<-as.numeric(MetaS[,10])
TIMEt<-TIMEt-min(TIMEt)

COLt<-rgb(1-TIMEt/max(TIMEt),1-TIMEt/max(TIMEt),0)

#For each strain, return average relative abundance of its corresponding ASV(s) in samples with same subsite and same sampling time

Fst<-sapply(STa,function(x){
	X=x
	Gx<-subset(GROUP,sapply(GROUP,function(x){
		return(length(intersect(X,
			unlist(strsplit(x,split="_")))))
	})==1)
	Gx<-unlist(strsplit(Gx,split="_"))
	Ax<-intersect(ASVn, Gx)
	#ASV relative abundance per sample
	FRx<-apply(sapply(Ax,function(x){
		return(as.vector(subset(t(FRallR),
			paste("Meth",colnames(FRallR),sep="-")==x)))
	}),1,sum)
	TIMEx<-subset(TIMEt, STa ==x)
	SStx<-subset(SSt, STa ==x)
	FR0<-mean(subset(FRx ,paste(SS,TIME-min(TIME))==
		paste(SStx, TIMEx)))
	#Relative abundance at the same subsite a the previous time poin
	TIMEp<-subset(TIMEt, SSt ==SStx)
	TIMEp<-max(subset(TIMEp, TIMEp<TIMEx))
	FRp<-mean(subset(FRx ,paste(SS,TIME-min(TIME))==
		paste(SStx, TIMEp)))	
		
	return(c(FRp, FR0))
})

STmilag<-t(sapply(STa,function(x){
	Sx<-ifelse(x=="E-061","E-061-B", x)
	return(as.numeric(subset(MiLag,rownames(MiLag)==Sx)))
}))
colnames(STmilag)<-colnames(MiLag)

library(ade4)

ACPmilag=dudi.pca(STmilag, scannf=T)
4
part<-round(100* ACPh $eig/sum(ACPh $eig),2)

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"a",cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

X=1
Y=2

x1=min(ACPmilag $li[, X])
x2=max(ACPmilag $li[, X])
y1=min(ACPmilag $li[,Y])
y2=max(ACPmilag $li[,Y])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis ",X," (",part[X],"%)",sep=""),
		ylab=paste("Axis ",Y," (",part[Y],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPmilag $li[,c(X,Y)],groups= TIMEt,col= rgb(0.9,0.9,0.9),lty=1,lwd=0.5)	
ordiellipse(ACPmilag $li[,c(X,Y)],groups= TIMEt,border="NA",alpha=40,draw="polygon",col= "grey",lty=1,label=T,cex=0.5,font=1)	

points(ACPmilag $li[,c(X,Y)],pch=19,cex=10*Fst,col=COLt)






#Average Milag values per group
ASVmilag<-sapply(GROUP,function(x){
	#print(x)
	Gx<-unlist(strsplit(x,split="_"))
	Sx<-intersect(Gx,ST)
	#MiLag values
	Sx<-ifelse(Sx=="E-061","E-061-B_E-061-R", Sx)
	Sx<-unlist(strsplit(Sx,split="_"))
	MiLagx<-apply(sapply(Sx,function(x){
		return(as.numeric(subset(MiLag,rownames(MiLag)==x)))
	}),1,mean)
	return(MiLagx)
})
colnames(ASVmilag)<-NULL
rownames(ASVmilag)<-colnames(MiLag)

#Total ASV relative abundance per group

ASVfr<-sapply(GROUP,function(x){
	#print(x)
	Gx<-unlist(strsplit(x,split="_"))
	Ax<-intersect(ASVn, Gx)
	#ASV relative abundance per sample
	FRx<-apply(sapply(Ax,function(x){
		return(as.vector(subset(t(FRallR),
			paste("Meth",colnames(FRallR),sep="-")==x)))
	}),1,sum)	
	return(FRx)
})
colnames(ASVfr)<-NULL
rownames(ASVfr)<-rownames(FRallR)

FRst<-sapply(sort(unique(paste(SITE,TIME))),function(x){
	return(
		apply(
		subset(ASVfr,paste(SITE,TIME)==x),2,mean))
	})



#Consensus clade per group
CX<-as.vector(sapply(GROUP,function(x){
	#print(x)
	Gx<-unlist(strsplit(x,split="_"))
	Ax<-intersect(ASVn, Gx)
	return(unique(sapply(Ax,function(x){
		return(subset(CLA[,5],CLA[,1]==x))
	})))
}))

COLx<-sapply(CX,function(x){
	return(subset(COLcla[,2], COLcla[,1]==x))
})


	
#Total ASV absolute abundance per group

ASVfa<-sapply(GROUP,function(x){
	#print(x)
	Gx<-unlist(strsplit(x,split="_"))
	Ax<-intersect(ASVn, Gx)
	#ASV relative abundance per sample
	FRx<-apply(sapply(Ax,function(x){
		return(as.vector(subset(t(FRall),
			paste("Meth",colnames(FRall),sep="-")==x)))
	}),1,sum)	
	return(FRx)
})
colnames(ASVfa)<-NULL
rownames(ASVfa)<-rownames(FRallR)

library(ade4)
library(vegan)

ASVfaH <-as.data.frame(decostand(ASVfr ,method="hellinger", MARGIN=1))

ACPh=dudi.pca(ASVfaH, scannf=F)
part<-round(100* ACPh $eig/sum(ACPh $eig),2)

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"a",cex=1.5,font=2)
par(new=T,mar=c(4,4,1,1))

x1=min(ACPh $li[,1])
x2=max(ACPh $li[,1])
y1=min(ACPh $li[,2])
y2=max(ACPh $li[,2])


plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1)

ordispider(ACPh $li[,c(1,2)],groups= TIME,col= rgb(0.9,0.9,0.9),lty=1,lwd=0.5)	

ordiellipse(ACPh $li[,c(1,2)],groups= TIME,border="NA",alpha=40,draw="polygon",col= "grey",lty=1,label=T,cex=0.5,font=1)	

points(ACPh $li[,c(1,2)],pch=19,cex=0.4,col=rgb(0,0,0,0.6))


#Group contribution
CONT<-ACPh$co

#Center of the contribution plot in the main plot
Xc=0
Yc=-3

points(Xc,Yc,pch=10,cex=5,lwd=0.5,col="black",lend=2)

segments(Xc,Yc,Xc+ CONT[,1]*1,Yc+CONT[,2]*1,col= COLx,lwd=1.5,lend=2)

points(Xc+ CONT[,1],Yc+CONT[,2],pch=21,cex=1.4,col=MAINtaxCOL,lwd=0.75,bg="white")

text(Xc+ CONT[,1],Yc+CONT[,2], MAINtaxCOD,col= "black",cex=0.3)
	




################################################
###################################################
### OLD: perform spatial and temporal autocorrelation analysis on diversity within each significant node of the tree


##### STEP 6a
# Calculate pairwise geographic distance among subsites and trees 

p0="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/4-Global-community-analysis"

setwd(p0)

#Map within subsite #in Echantillonnage_2018-2019.xlsx
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

##Calculate geographic distance among subsites
#Calculate a Distance Matrix for Geographic Points Using R (HOMEMADE SCRIPT FOUND IN INTERNET)

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

#Pairwise geographic distance among all trees
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

##### STEP 6b
# LOAD and format meta data and ASV table only for Methylobacteriaceae

pM="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis"

pO="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/6-Strains-2018"

pT<-"/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis/rTREE"

###Load data

	#strain (tips of the tree) names
	setwd(pT)
	
	#minimum bootstrapp value
	bx=0
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
	#Reference sequence informations (names in second column)
	setwd(pO)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	
pR="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/2-rarefaction&control"

pC="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/5-Comparative-community-analysis"

## Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs

setwd(pR)

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
colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="")

###Apply Hellinger transformation on relative abundance to account for rare ASVs
library(vegan)

FRallH <-as.data.frame(decostand(FRall ,method="hellinger", MARGIN=1))

colnames(FRallH)<-paste("ASV", c(1:length(FRallH[1,])),sep="")

rm(seqtab.nochim)

##Load Taxonomic assignation for each ASV
pA="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/3-Taxonomic-assignation"

setwd(pA)

taxaALL <-read.table(paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

##Simplify taxonomy (Emphasize on Order, and Familly within Rhizobiales)

ASSsimp<-ifelse(is.na(as.vector(taxaALL[,5]))==T,"unknown",ifelse(is.na(as.vector(taxaALL[,6]))==T,as.vector(taxaALL[,5]),ifelse(as.vector(taxaALL[,5])=="Rhizobiales",paste(taxaALL[,5], taxaALL[,6],sep="_"),as.vector(taxaALL[,5]))))

#Genus
GENx<-taxaALL[,7]

rm(taxaALL)

#Methylobacteriaceae ASVs
ASVm<-subset(c(1:length(ASSsimp)), ASSsimp=="Rhizobiales_Methylobacteriaceae")

rm(ASSsimp)

GENx <-as.vector(GENx[ASVm])

#ASV names (genus names 4 first digit + ASV number)

ASVn<-paste(as.vector(sapply(GENx,function(x){
	return(paste(unlist(strsplit(x,split=""))[1:4],collapse=""))
})),"ASV", ASVm,sep="-")

ASVi<-cbind(NA,"ASV",ASVn,GENx,"sp.","PHYLLOSPHERE",NA)
colnames(ASVi)<-colnames(CHAR)

CHAR<-rbind(CHAR,ASVi)

#tree information
setwd(pT)
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
	CHARx<-as.vector(t(subset(CHAR,CHAR[,3]== stx )))
	typex<-ifelse(CHARx[2]=="ASV","ASV","ref")	
	return(c(typex ,CHARx[c(3:7)]))
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

#METAdata formated for permanova
FREQ.ENV<-as.data.frame(cbind(SITE,SS,SPE,TIME,EXT,PCR,MISEQ))

#Tree topology in level format
LEVELS<-read.table(paste("Levels_MINBOOT=",
	bx,".txt",sep=""))
colnames(LEVELS)<-paste("L",
	c(1:(length(LEVELS[1,]))),sep="-")
	
#remove empty levels
LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))

#Only Keep levels for ASVs
LEVELf<-LEVELS[ASVi,]

###Pairwise distance (geo and time among all samples)

SITEu<-sort(unique(SITE))

DISTtime<-matrix(unlist(sapply(SITEu,function(x){
	#print(x)
	SITEx=x
	Sx<-subset(c(1:length(SITE)), SITE == SITEx) #samples from this site
	
	DISTx<-matrix(unlist(sapply(c(1:(length(Sx)-1)),function(x){
		X=x
		SS1<-SS[Sx[X]] #Subsite for this sample
		T1= TIME[Sx[X]]#Sampling time for this sample
		SPE1<-SPE[Sx[X]] #host tree species for this sample
		TREE1<-ID[Sx[X]] #tree identity for this sample
		MISEQ1<-MISEQ[Sx[X]] #sequencing batch for this sample
		return(sapply(c((X+1):(length(Sx))),function(x){
			T2= TIME[Sx[x]]
			SS2<-SS[Sx[x]]
			SPE2<-SPE[Sx[x]]
			TREE2<-ID[Sx[x]]
			MISEQ2<-MISEQ[Sx[x]]	
			D12<-subset(DISTtree[,2],DISTtree[,1]
				==paste(sort(c(TREE1,TREE2)),collapse="_")) #Spatial distance			
			return(c(X,x,SS1,SS2,TREE1,TREE2,
				SPE1,SPE2,MISEQ1,MISEQ2,
				D12,T1,T2))		
		}))
	})),ncol=13,byrow=T)
	return(t(cbind(SITEx,DISTx)))
})),ncol=14,byrow=T)

colnames(DISTtime)<-c("Site","I1","I2","SS1","SS2","TREE1","TREE2","SPE1","SPE2","MISEQ1","MISEQ2","GEO-DIST","DAY1","DAY2")


#STEP 6c: perform spatial autocorrelation analyzes for each level of the tree

#Test only nodes with at least BP bootstrap
BP=0.3

node<-subset(PSn[,2],PSn[,8]>=BP)

#For each node, and each site, count the number of sample with enough observations and number of ASV
OBS<-t(sapply(node,function(x){
	NX=x
	#Level
	LEVELx <-subset(PSn[,1], PSn[,2]== NX)
	#Column index in the level file
	Lx<-subset(c(1:length(LEVELf[1,])),paste("L", LEVELx,sep="-")==colnames(LEVELf))
	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	#ASV line indexes in the level file for this node
	STx<-subset(c(1:length(LEVELf[,1])),ND== NX)
	#ASV abundance
	Fx<-cbind(FRall[,STx],0)
	
	Fsite<-sapply(SITEu,function(x){
		SITEx=x
		Fxx<-subset(Fx, SITE== SITEx)
		return(length(subset(as.vector(apply(Fxx,1,sum)),
			as.vector(apply(Fxx,1,sum))!=0)))
	})		
	
	return(c(x,length(STx), Fsite))
})	)	

#Only keep nodes with at least 2 ASVs
OBS <-subset(OBS, OBS[,2]>=2)

#Only keep nodes with at least Nsamp observation in at least MSH or SBL
Nsamp=20

setwd(paste(pM,"AUTOCORRELATION-TREE",sep="/"))

#Perform autocorrelation analysis for each node

SITEx="SBL"


nodex<-subset(OBS[,1],OBS[,4]>= Nsamp)

SUMcor<-t(sapply(nodex,function(x){
	NX=x
	print(NX)
	#Level
	LEVELx <-subset(PSn[,1], PSn[,2]== NX)
	#Column index in the level file
	Lx<-subset(c(1:length(LEVELf[1,])),paste("L", LEVELx,sep="-")==colnames(LEVELf))
	
	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	
	#ASV line indexes in the level file for this node
	STx<-subset(c(1:length(LEVELf[,1])),ND== NX)
	
	#ASV abundance
	Fx<-FRall[,STx]
	
	#ASV abundance (Hellinger)
	FHx<-FRallH[,STx]
	
	#Node PS median value
	PSx<-subset(PSn[,5],PSn[,2]==NX)
	
	#Node clade Assignation
	CLAx<-paste(unique(subset(CLA[,5],ND==NX)),collapse="_")
	
	#Calculate pairwise Bray-curtis dissimilarities between samples for these ASVs (Hellinger)	
	BDc<-as.numeric(unlist(sapply(SITEu,function(x){
		#print(x)
		SITEx=x
		Sx<-subset(c(1:length(SITE)), SITE == SITEx) #samples from this site
		FS<-subset(FHx,SITE == SITEx) #abundances for this site	
		return(unlist(sapply(c(1:(length(Sx)-1)),function(x){
			X=x
			S1<-as.vector(t(FS[X,]))
			return(sapply(c((X+1):(length(Sx))),function(x){
				S2<-as.vector(t(FS[x,]))
				vegdist(t(cbind(S1,S2)),method="bray")			
			}))
		})))
	})))
		
	#Remove pairwise comparisons with no observation
	DISTx<-subset(DISTtime,is.nan(BDc)==F)
	BDc<-subset(BDc,is.nan(BDc)==F)
	#Remove pairwise comparisons with observation in only one pair
	DISTx<-subset(DISTx,BDc!=1)
	BDc<-subset(BDc,BDc!=1)	
		
	
		D1<-subset(DISTx, DISTx[,1]== SITEx)
		Bx<-subset(BDc, DISTx[,1]== SITEx)
		
		Gx<-as.numeric(D1[,12]) #Geographic distance
		Dx<-as.numeric(D1[,13]) #Sampling date
		Tx<-abs(as.numeric(D1[,13])-as.numeric(D1[,14]))#Pairwise time
		Dmx<-Dx-min(Dx)
		
		Mg<-lm(Bx ~Gx* Dmx*Tx)	
		Ag<-anova(Mg)
		Cg<-coefficients(Mg)
		
	return(c(NX,PSx,CLAx,Ag[,2]/sum(Ag[,2]),Ag[,5],Cg))
}))





################################################
###################################################
################################################
###################################################
################################################
###################################################
################################################
###################################################
###OLD:: perform permutation for each node to test for significant association with factors. Permutations will be done on absolute abundance. Because of the high number of ASVs, permutations will be donne on random subsamples of KS sequences

library(ape)
library(seqinr)

pM="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis"

pO="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/6-Strains-2018"

pT<-"/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/6-Methylobacteriaceae-community-analysis/rTREE"

###Load data

	#strain (tips of the tree) names
	setwd(pT)
	
	#minimum bootstrapp value
	bx=0
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
	#Reference sequence informations (names in second column)
	setwd(pO)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	
pR="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/2-rarefaction&control"

pC="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/5-Comparative-community-analysis"

## Load ASV relative abundance table, Metadata for each sample and format them - Extract Methylobacteriaceae ASVs

setwd(pR)

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

rm(seqtab.nochim)

colnames(FRall)<-paste("ASV", c(1:length(FRall[1,])),sep="")

##Load Taxonomic assignation for each ASV
pA="/Users/jean-baptisteleducq/Dropbox/Methylobacterium/Article/7-rpoB-community-analysis/3-Taxonomic-assignation"

setwd(pA)

taxaALL <-read.table(paste("PHYLLOAnnotation_TF=",TFx,"_TR=",TRx,"_TL=",TLx,"_NS=",NS,"_prare=", Prare,"-wo-NEGcontrols.txt",sep=""),header=T)

##Simplify taxonomy (Emphasize on Order, and Familly within Rhizobiales)

ASSsimp<-ifelse(is.na(as.vector(taxaALL[,5]))==T,"unknown",ifelse(is.na(as.vector(taxaALL[,6]))==T,as.vector(taxaALL[,5]),ifelse(as.vector(taxaALL[,5])=="Rhizobiales",paste(taxaALL[,5], taxaALL[,6],sep="_"),as.vector(taxaALL[,5]))))

#Genus
GENx<-taxaALL[,7]

rm(taxaALL)

#Methylobacteriaceae ASVs
ASVm<-subset(c(1:length(ASSsimp)), ASSsimp=="Rhizobiales_Methylobacteriaceae")

rm(ASSsimp)

GENx <-as.vector(GENx[ASVm])

#ASV names (genus names 4 first digit + ASV number)

ASVn<-paste(as.vector(sapply(GENx,function(x){
	return(paste(unlist(strsplit(x,split=""))[1:4],collapse=""))
})),"ASV", ASVm,sep="-")

ASVi<-cbind(NA,"ASV",ASVn,GENx,"sp.","PHYLLOSPHERE",NA)
colnames(ASVi)<-colnames(CHAR)

CHAR<-rbind(CHAR,ASVi)

#tree information
setwd(pT)
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
	CHARx<-as.vector(t(subset(CHAR,CHAR[,3]== stx )))
	typex<-ifelse(CHARx[2]=="ASV","ASV","ref")	
	return(c(typex ,CHARx[c(3:7)]))
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

#METAdata formated for permanova
FREQ.ENV<-as.data.frame(cbind(SITE,SS,SPE,TIME,EXT,PCR,MISEQ))

#Tree topology in level format
LEVELS<-read.table(paste("Levels_MINBOOT=",
	bx,".txt",sep=""))
colnames(LEVELS)<-paste("L",
	c(1:(length(LEVELS[1,]))),sep="-")
	
#remove empty levels
LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))
	
#for each sequence, fill empty levels with node from the closest level
	
LEVELf<-t(sapply(c(1:length(LEVELS[,1])),function(x){
	Sx=x
	Lx<-as.vector(LEVELS[Sx,])
	lx<-c(1:length(Lx))
	Lx<-sapply(lx,function(x){
		return(max(subset(Lx,lx>=x)))
	})
	return(c(Sx,Lx[2:length(Lx)]))
}))

colnames(LEVELf)<-colnames(LEVELS)

#Only Keep levels for ASVs
LEVELf<-LEVELf[ASVi,]

#Level column index
LV=c(2:length(LEVELf[1,]))

#Number of ASVs per level 
STn<-sapply(LV,function(x){
	#print(x)
	Lx=x
	return(length(unique(LEVELf[, Lx])))
})	

#Remove levels with less than 2 ASVs
LV2<-subset(LV, STn>=2)


#Permutations tests performed separately in MSH and SBL

#Column index in FREQ.ENV, for factors to test (Site, Plot, Spe, Time)
	FACi<-c(2,3,4)

##Column index in FREQ.ENV for categories of permutation (Plot, Spe, Time, sequencing batch). For every factor but site, all categories different than the tested factor must be took into account as permutation categories, assuming that factors are well balanced within site. For site, only consider permutations among sequencing batches, as nested factors (BUT FAGR and ACSA in SPE) are completely independant. We assume that the difference between sites is so big that the effect of shared species between them is neglectible 
	CATi<-c(2,3,4)

#Number of permutations
K=10000

#Number of sequences to subsample or permutations (redo sampling for every permutation)
KS=1000

#Test only nodes with at least BP bootstrap
BP=0.5

#Perform permutation for each level

setwd(paste(pM,"PERMUTATIONS-TREE",sep="/"))

sapply(LV2,function(x){
	Lx=x
	#Level
	LEVELx<-as.numeric(unlist(strsplit(colnames(LEVELf)[x],split="-"))[2])
	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	#Bootstrap values
	BD<-sapply(ND,function(x){
		return(subset(PSn[,8], PSn[,2]==x))
	})
	#Nodes with less than BP bootstrap support
	NDno<-subset(ND ,BD<BP)

	NDu<-setdiff(unique(ND),NDno)
	####For each supported node, retrieve ASV abundance (From FRall, wil be used for permutations)
	N1<-sapply(NDu,function(x){
		NDx=x
		#ASV index in this node
		STx<-subset(c(1:length(LEVELf[,1])),ND== NDx)
		#total abundance per sample for these ASVs
		FRx<-cbind(FRall[,STx],0)
		return(apply(FRx,1,sum))
	})
	#Ad total abundance of non supported nodes
	N0<-as.vector(apply(FRall,1,sum)-apply(N1,1,sum))
	N1<-cbind(N1, N0)
	
	colnames(N1)<-paste("node", c(NDu,"un"),sep="-")
	
	#For each node, retrieve ASV relative abundance (From FRallR, will only be used to estimate real values)
	N1R<-t(sapply(c(1:length(N1[,1])),function(x){
		return(N1[x,]/sum(N1[x,]))
	}))
	colnames(N1R)<-colnames(N1)	
	
	#remove nodes without observation
	N1<-as.data.frame(t(subset(t(N1),apply(N1,2,max)!=0)))
	N1R<-as.data.frame(t(subset(t(N1R),apply(N1R,2,max)!=0)))	
	#Nodes 
	MB1<-colnames(N1)
	
	#Probability of sampling corrected per leaf sample (so that the number of sequences per leaf samples will be equal)
	PR<-matrix(sapply(c(1:length(N1[,1])),
		function(x){
		Fx<-as.vector(t(N1)[,x])
		Tx<-sum(Fx)
		Sx=1/length(N1[,1])
		Fx<-Sx* Fx/Tx
		return(t(cbind(x,c(1:length(Fx)),Fx)))
	}),ncol=3,byrow=T) #First column: sample; second column: node; third column: probability of sampling
	
	#Remove nodes/sample with null probability
	PR<-subset(PR,PR[,3]!=0)
	
	
	####1 : Test for association with sites
	print(date())
	print(paste("Level=", LEVELx,"; Fac=",colnames(FREQ.ENV)[1]))
	#Unique factors
	FACu<-sort(unique(as.vector(FREQ.ENV[, 1])))

	#Perform K permutation of nodes within categories
	PERM<-t(sapply(c(1:K),function(x){
		#print(x)
		X=x	
		####Rarefy the node table (N1) to have KS sequences
		#1 sample KS sequences per sample
		vs<-sample(c(1:length(PR[,1])),size=KS,prob= PR[,3],replace=T)
		#2 attribute to nodes
		SRx<-PR[vs,2]
		#3 attribute to samples
		LRx<-PR[vs,1]
		#Sampled Nodes
		MBx<-MB1[SRx]
		MB<-sort(unique(MBx))
		#Factor to test
		FAC<-as.vector(FREQ.ENV[LRx,1])
		#category within wich permutations have to be done to avoid interaction among factors
		CAT<-as.vector(FREQ.ENV[LRx, 7])
		CATperm<-sort(unique(CAT)) 
		# Permute nodes within CAT 
		MBp<-MB1[unlist(sapply(CATperm,function(x){
			return(sample(
				subset(SRx, CAT==x),
				size=length(subset(c(1:length(MBx)), 
				CAT==x))))		
		}))]
		# Retrieve unpermuted factors	
		FACp<-FAC[unlist(sapply(CATperm,function(x){
			return(subset(c(1:length(FAC)), CAT==x))		
		}))]
		# Return expected frequency per node per factor	
		return(sapply(MB1,function(x){
			FACx<-subset(FACp, MBp==x)
			return(sapply(FACu,function(x){
				return(length(subset(FACx, FACx==x)))
			}))
		}))
	}))/KS	
	#Real values (relative abundance per factor)
	FAC<-as.vector(FREQ.ENV[, 1])
	REAL<-as.vector(t(sapply(FACu,function(x){
		Nx<-subset(N1R, FAC ==x)
		return(apply(Nx,2,mean))
	})))/length(FACu)
	#REAL<-KS*REAL/sum(REAL)
	#Expected values (permutations)
	EXP<-paste(round(apply(PERM,2,mean),2),
		round(apply(PERM,2,sd),2),sep="±")
	#Factor_node association
	DESC<-matrix(sapply(MB1,function(x){
		X=x
		return(sapply(FACu,function(x){
			return(paste(x,X,sep="_"))
		}))
	}),nrow=1)
	#P value calculated from permutations for each association
	pv<-(sapply(c(1:length(DESC)),function(x){
		return(	length(subset(PERM[,x], 
		PERM[,x]>=REAL[x])))
	})+1)/(K+1)
	#Adjusted pvalue
	pu<-p.adjust(pv, method="bonferroni")

	#summary file
	SUMsite<-cbind("SITE",as.vector(DESC) ,
		as.vector(REAL), EXP,round(pv,5),round(pu,5))
	colnames(SUMsite)<-c("fac","ass.","obs",
		"exp","p.obs>exp","padj.obs>exp")

	####2 : Test for association with other factors, within sites
	SUMnode<-matrix(unlist(sapply(SITEu,function(x){
		SITEx=x
		FREQ.ENVx<-subset(FREQ.ENV, FREQ.ENV[,1]== SITEx)
		N1x<-subset(N1, FREQ.ENV[,1]== SITEx)
		N1Rx<-subset(N1R, FREQ.ENV[,1]== SITEx)
		#remove nodes without observation
		N1x <-as.data.frame(t(subset(t(N1x),apply(N1x,2,max)!=0)))
		N1Rx <-as.data.frame(t(subset(t(N1Rx),apply(N1Rx,2,max)!=0)))
		MB1x<-colnames(N1x)	
		return(t(matrix(unlist(sapply(FACi,function(x){
			print(paste("SITE=",SITEx,"; Level=", LEVELx,"; Fac=",colnames(FREQ.ENV)[x]))
			FACx=x
			#Column index in data for categories of permutation, excluding the tested factor
			CATx<-setdiff(CATi,FACx)
			#Unique factors
			FACu<-sort(unique(as.vector(FREQ.ENVx[, FACx])))

			#Probability of sampling corrected per leaf sample (so that the number of sequences per leaf samples will be equal)
			PR<-matrix(sapply(c(1:length(N1Rx[,1])),
				function(x){
				Fx<-as.vector(t(N1Rx)[,x])
				Tx<-sum(Fx)
				Sx=1/length(N1Rx[,1])
				Fx<-Sx* Fx/Tx
				return(t(cbind(x,c(1:length(Fx)),Fx)))
			}),ncol=3,byrow=T) #First column: sample; second column: node; third column: probability of sampling
	
			#Remove nodes/sample with null probability
			PR<-subset(PR,PR[,3]!=0)

			#Perform K permutation of nodes within categories
			PERM<-t(sapply(c(1:K),function(x){
				#print(x)
				X=x			
				####Rarefy the node table (N1) to have KS sequences
				#1 sample KS sequences per sample
				vs<-sample(c(1:length(PR[,1])),size=KS,prob= PR[,3],replace=T)
				#2 attribute to nodes
				SRx<-PR[vs,2]
				#3 attribute to samples
				LRx<-PR[vs,1]				
				#Nodes
				MBx<-MB1x[SRx]
				MB<-sort(unique(MBx))
				#Factor to test
				FAC<-as.vector(FREQ.ENVx[LRx,FACx])
				#category within wich permutations have to be done to avoid interaction among factors
				CAT<-as.vector(paste(FREQ.ENVx[LRx, CATx[1]],
					FREQ.ENVx[LRx, CATx[2]],
					FREQ.ENVx[LRx, CATx[3]],sep="_"))
				CATperm<-sort(unique(CAT)) 
				MBp<-MBx[unlist(sapply(CATperm,function(x){
					return(sample(subset(c(1:length(MBx)), CAT==x),
						size=length(subset(c(1:length(MBx)), 
							CAT==x))))		
				}))]	
				FACp<-FAC[unlist(sapply(CATperm,function(x){
					return(subset(c(1:length(FAC)), CAT==x))		
				}))]	
				return(sapply(MB1x,function(x){
					FACx<-subset(FACp, MBp==x)
					return(sapply(FACu,function(x){
						return(length(subset(FACx, FACx==x)))
					}))
				}))
			}))	/KS	
			#Real values (relative abundance per factor)
			FAC<-as.vector(FREQ.ENVx[, FACx])
			REAL<-as.vector(t(sapply(FACu,function(x){
				Nx<-subset(N1Rx, FAC ==x)
				return(apply(Nx,2,mean))
			})))/length(FACu)
			#Expected values (permutations)
			EXP<-paste(round(apply(PERM,2,mean),2),
				round(apply(PERM,2,sd),2),sep="±")

			#Factor_node association
			DESC<-matrix(sapply(MB1x,function(x){
				X=x
				return(sapply(FACu,function(x){
					return(paste(x,X,sep="_"))
				}))
			}),nrow=1)
			#P value calculated from permutations for each association
			pv<-(sapply(c(1:length(DESC)),function(x){
				return(	length(subset(PERM[,x], PERM[,x]>=REAL[x])))
			})+1)/(K+1)
			#Adjusted pvalue
			pu<-p.adjust(pv, method="bonferroni")
			#summary file
			SUMxx<-cbind(as.vector(DESC) ,as.vector(REAL), 
				EXP,round(pv,5),round(pu,5))
			colnames(SUMxx)<-c("ass.","obs","exp",
				"p.obs>exp","padj.obs>exp")
			return(t(cbind(paste(SITEx,colnames(FREQ.ENVx)[FACx],sep="_"), SUMxx)))
		})),ncol=6,byrow=T)))
	})),ncol=6,byrow=T)
	
	SUMnode<-rbind(SUMnode, SUMsite)
	
	print(subset(SUMnode,as.numeric(SUMnode[,6])<=0.05))
	
	write.table(SUMnode,
		paste("PERM_LEVEL=", LEVELx,"_K=",K,"_boot=",bx,".txt",sep=""),row.names=F)

})
	
#arrivé ici

############################################################################################################################################################




#Display metadata on a circular tree
#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=0.85

#Compression of tree from the root (>1) 
COMP=1.1

#Margin size
MAR=1.5

#tip attribute coordinates (concentric circles)
Ai=1.01+seq(0.05,1,length.out=10)
Ae=Ai+0.8*(Ai[2]-Ai[1])

#Number of graduations in the PS scale
Nscale=6

#Branch color
COLbr<-rgb(0.7,0.7,0.7)

Lmin=min(c(coordx[,5], coordx[,7]))
Lmax=max(coordx[,5])*COMP

#Limit of the graphic
rLIM=MAR* Lmax
z=max(c(coordx[,1], coordx[,2],coordx[,3]))

pdf(paste("Circular-tree_MINBOOT=",bx,".pdf",sep=""))
par(mar=c(1,1,1,1),bty="n",mfrow=c(1,1))

#Frame
plot(-10000,-10000,xlim=c(-rLIM, rLIM),
	ylim=c(-rLIM, rLIM),xlab="",
	ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n")
	
#TREE
COORD<-t(sapply(c(1:length(coordx[,1])),function(x){
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
		col=COLbr<-rgb(0.7,0.7,0.7),lwd=1)
	segments(XXb1, YYb1,XXb2, YYb2,
		col=COLbr<-rgb(0.7,0.7,0.7),lwd=1)
	#points(XXb1, YYb1,pch=19,col=rgb(0,0,0,1),cex=0.5)
		return(c(x,XXb1, YYb1))
}))
	
#Node attributes (bootstrapp)
Nb=11
Nbl=3
	
Nx<-c(1:length(PSn[,1]))
colx<-seq(0.8,0,length.out=Nb)
colx<-rgb(colx, colx, colx)
	
sapply(Nx,function(x){
	bx=PSn[x,8]
	bx<-ifelse(is.na(bx)==T,0,bx)
	cx<-subset(colx,seq(0,1,length.out=Nb)==round(bx,1))
	nx=PSn[x,2]
	COORDxx<-subset(COORD,coordx[,4]==nx)
	points(COORDxx[1,2], 
		COORDxx[1,3],pch=19,col= cx,cex=bx/3)
})
	
#Bootstrapp legend
	
sapply(c(1:Nb),function(x){
	X=seq(0.7*rLIM,0.85*rLIM,length.out=Nb)[x]
	bbx=seq(bx,1,length.out=Nb)[x]
	points(X,0.86*rLIM,pch=19,col=colx[x],cex=bbx/3)		
})
text(0.775*rLIM,0.9*rLIM,"Node support",cex=0.6,font=1)
sapply(c(1: Nbl),function(x){
	X=seq(0.7*rLIM,0.85*rLIM,length.out= Nbl)[x]
	text(X,0.86*rLIM,
		round(seq(bx*100,100,length.out=Nbl)[x],0),
		cex=0.5,pos=1,adj=1)			
})	

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
PSscale<-PSscale[c(1:length(Xscale))]
segments(-0.01*rLIM,min(Xscale),
	-0.01*rLIM,max(Xscale))
text(-0.005*rLIM, Xscale, round(PSscale,2)
	,cex=0.6,pos=2,adj=1,font=1)
text(-0.17*rLIM, mean(Xscale), "PS",cex=0.8,font=1)
	
maxTIP=max(Xscale)

#####PLOT Tips attributes
Nlab=12 #maximum number of labels in legend
l1=-0.05 #maximum x coord for legend
l2=-0.75 #minimum x coord for legend

#######environment of sampling
	
ATx=2
	
#Display biome of origin for each sequence
colECO<-cbind(
	c("PLANT","PHYLLOSPHERE","ENDOPHYTE",
	"SEED","RHIZHOSPHERE","SOIL","LICHEN-FUNGI",
	"WATER","OCEAN","MICROBIOME"),
	c("green3","green","yellow3",
	"yellow","brown","red","orange",
	"cyan","blue","purple"),
	c("Plant","Phyllo.","Endop.",
	"Seed","Rhizo.","Soil","Fung.",
	"Water","Ocean","Gut"))	
	
sapply(Tx,function(x){
	typex<-unlist(strsplit(type[x,5],split="/"))
	antx<-length(subset(typex, typex=="ANTHROPO"))
	typex<-setdiff(typex,"ANTHROPO")
	colECOx<-subset(colECO[,2], colECO[,1]== typex)
	ang1=2*pi*((x-0.5)/z)*OPEN
	ang2=2*pi*((x+0.5)/z)*OPEN
	XX1i<-Ai[ATx]* maxTIP*sin(ang1)
	YY1i<-Ai[ATx]* maxTIP*cos(ang1)
	XX2i<-Ai[ATx]* maxTIP*sin(ang2)
	YY2i<-Ai[ATx]* maxTIP*cos(ang2)
	XX1e<-Ae[ATx]* maxTIP*sin(ang1)
	YY1e<-Ae[ATx]* maxTIP*cos(ang1)
	XX2e<-Ae[ATx]* maxTIP*sin(ang2)
	YY2e<-Ae[ATx]* maxTIP*cos(ang2)		
	polygon(c(XX1i,XX2i,XX2e,XX1e),
		c(YY1i,YY2i,YY2e,YY1e),
		col= colECOx,border=ifelse(antx==1,
		"black",NA),lwd=0.5)
})
	
#Scale 
Ns=length(colECO[,1])+1	
colx=c(colECO[,2],"white")
bordx<-c(sapply(c(1:(Ns-1)),function(x){
	return(NA)
}),"black")
	
y1=Ae[ATx]
y2=(Ae-((Ae-Ai)/2))[ATx]

Nmin=min(sapply(c(1:Ns),function(x){
	x1=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x]
	x2=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x+1]
	polygon(c(x1,x1,x2,x2),
		c(y1, y2, y2, y1)*maxTIP,
		col= colx[x],border= bordx[x],lwd=0.5)
	text((x1+x2)/2,(y1 + y2)*maxTIP/2,
		c(colECO[,3],"Anthr.")[x],cex=0.35)	
	return(min(c(x1,x2)))				
}))

text(Nmin, ((y1+y2)/2)*maxTIP, "Environment"
	,cex=0.6,pos=2,adj=1,font=1)	

#Second tips attribute (ASV relative abundance, Hellinger transformation, average accross samples)
	
Cx3=log(as.numeric(type[,3]))
	
	Cx2=ifelse(as.numeric(type[,3])==0,min(subset(Cx3,as.numeric(type[,3])!=0)),Cx3)	
	
	Cx1=Cx2+abs(min(Cx2))
	Cx<-Cx1/max(Cx1)
	
	sapply(Tx,function(x){
		ang1=2*pi*((x-0.5)/z)*OPEN
		ang2=2*pi*((x+0.5)/z)*OPEN
		XX1i<-Ai[1]* maxTIP*sin(ang1)
		YY1i<-Ai[1]* maxTIP*cos(ang1)
		XX2i<-Ai[1]* maxTIP*sin(ang2)
		YY2i<-Ai[1]* maxTIP*cos(ang2)
		XX1e<-Ae[1]* maxTIP*sin(ang1)
		YY1e<-Ae[1]* maxTIP*cos(ang1)
		XX2e<-Ae[1]* maxTIP*sin(ang2)
		YY2e<-Ae[1]* maxTIP*cos(ang2)		
		polygon(c(XX1i,XX2i,XX2e,XX1e),
			c(YY1i,YY2i,YY2e,YY1e),
			col=rgb(0,0,1,Cx[x]),border=NA)
	})
	
	text(-0.005*rLIM, ((Ai[1]+Ae[1])/2)*maxTIP, "ASV average abundance"
		,cex=0.6,pos=2,adj=1,font=1)
	
	#Scale first tip attribute
	Ns=100	
	Nl=4
	
	colx=seq(min(Cx),max(Cx),length.out=Ns)
	
	sapply(c(1:Ns),function(x){
		x1=seq(-0.9*rLIM,-0.5*rLIM,length.out=(Ns+1))[x]
		x2=seq(-0.9*rLIM,-0.5*rLIM,length.out=(Ns+1))[x+1]
		polygon(c(x1,x1,x2,x2),
			c(Ai[1],Ae[1],Ae[1],Ai[1])*maxTIP,
			col=rgb(0,0,1, colx[x]),border=NA)			
	})
	
	colx=seq(0,1,length.out=Nl)
	val=round(exp((colx*max(Cx1))-abs(min(Cx2))),2)
	
	sapply(c(1:Nl),function(x){
		X=seq(-0.9*rLIM,-0.5*rLIM,length.out=(Nl))[x]
		Y=Ai[1]*maxTIP
		valx =paste(round(val[x],3),"%")
		text(X,Y, valx,cex=0.5,pos=1,adj=1)			
	})	
	
	
	VAR<-c("SITE","SPE","TEMP","DAY14")
	
	##Color codes
	Tu<-sort(setdiff(as.numeric(type[,c(subset(c(1:length(type[1,])),colnames(type)=="DAY14"))]),c(NA)),decreasing=T)
		
	COLtu<-(Tu-min(Tu))/(max(Tu)-min(Tu))
	COLtu<-rgb(1-COLtu,1-COLtu,0)
	
	VARc<-t(rbind(
	cbind(1,"Sampling site",
		c("MSH","SBL"),
		c("green3","orange"),
		c("black","black"),
		c("MSH","SBL")),
	
	cbind(2,"Tree species",
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU"),
		c("purple2","violet","blue","green",rgb(0.35,0.35,0.35),
		"orange3","red3"),
		c("black","black","black","black","black","black","black"),
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU")),
	
	cbind(3,"Isolation temperature",
		c(20,30),
		c("cyan2","red2"),
		c("black","black"),
		c("20°C","30°C")),
	
	cbind(4,"Sampling time",Tu,COLtu,
		c("white",NA,NA,NA,"black"),
		c("Oct.","","","","June"))		
	))

	row.names(VARc)<-c("ring","fac","var","col","label.col","label")
	
	
	sapply(c(1:length(VAR)),function(x){
		X=x
		Aix<-Ai[X+1]
		Aex<-Ae[X+1]
		varx<-VAR[X]
		
		VARcx<-subset(t(VARc),as.numeric(VARc[1,])==X)
		typex<-as.vector(subset(t(type),colnames(type)==varx))
		
		colx<-as.vector(VARcx[,4])
		
		colx<-as.vector(sapply(typex,function(x){
			return(c(subset(as.vector(VARcx[,4]),
				as.vector(VARcx[,3])==x),NA)[1])
		}))
		ATx=X+2
		
		#Strain label
		sapply(Tx,function(x){
			cx<-colx[x]
			ang1=2*pi*((x-0.5)/z)*OPEN
			ang2=2*pi*((x+0.5)/z)*OPEN
			XX1i<-Ai[ATx]* maxTIP*sin(ang1)
			YY1i<-Ai[ATx]* maxTIP*cos(ang1)
			XX2i<-Ai[ATx]* maxTIP*sin(ang2)
			YY2i<-Ai[ATx]* maxTIP*cos(ang2)
			XX1e<-Ae[ATx]* maxTIP*sin(ang1)
			YY1e<-Ae[ATx]* maxTIP*cos(ang1)
			XX2e<-Ae[ATx]* maxTIP*sin(ang2)
			YY2e<-Ae[ATx]* maxTIP*cos(ang2)		
			polygon(c(XX1i,XX2i,XX2e,XX1e),
				c(YY1i,YY2i,YY2e,YY1e),
				col= cx,border=NA,lwd=0.5)
		})
		
		#Legend
		VARx<-unique(as.vector(VARcx[,2]))
		
		labx<-as.vector(VARcx[,6])
		colx<-as.vector(VARcx[,4])
		lab.colx<-as.vector(VARcx[,5])
			
		Ns=length(labx)
	
		y1=Ae[ATx]
		y2=(Ae-((Ae-Ai)/2))[ATx]
	
		Nmin=min(sapply(c(1:Ns),function(x){
			x1=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x]
			x2=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x+1]	
			polygon(c(x1,x1,x2,x2),
				c(y1,y2,y2,y1)*maxTIP,
				col=colx[x],border=NA)
			text(mean(c(x1,x2)),
				mean(c(y1,y2))*maxTIP, 
				labx[x],col=lab.colx[x],cex= 0.35)
			return(min(c(x1,x2)))				
		}))
		text(Nmin, ((y1+y2)/2)*maxTIP, VARx
		,cex=0.6,pos=2,adj=1,font=1)			
	})

	
	
	#Absolute abundance of Clades per sample
	Fcla<-sapply(sort(unique(CLA[,5])),function(x){	
		FRx<-subset(t(FRall),CLA[,5]==x)		
		return(apply(FRx,2,sum))
	})

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
		cbind("A5","yellow3",3,3,0.5),
		cbind("A7","pink",3,3,0.5),
		cbind("A10","green3",2,7,0.8),
		cbind("A8","purple",3,3,0.5))
	
	COLclax<- sapply(sort(unique(CLA[,5])),function(x){
		return(c(subset(COLcla[,2], COLcla[,1]==x),NA)[1])
	})
	
	barplot(sapply(sort(unique(paste(METAsamp[,3],METAsamp[,8]))),function(x){
		
		Fx<-subset(Fcla,paste(METAsamp[,3],METAsamp[,8])==x)/NS
		return(apply(Fx,2,mean))
	}),las=2,col= COLclax,border=NA)
	
	
	Tmsh<-sort(unique(subset(METAsamp[,8],METAsamp[,3]=="MSH")))
	
	plot(-10,1000,xlim=c(min(Tmsh), max(Tmsh)),ylim=c(0,2),log="",xlab="time",ylab='ab.',las=2)
	
	sapply(sort(unique(CLA[,5])),function(x){
		FRx<-t(subset(t(FRall),CLA[,5]==x))
		COLclax<-c(subset(COLcla[,2], COLcla[,1]==x),NA)[1]
		FRx <-rbind(sapply(Tmsh,function(x){
			FRxx<-subset(FRx,METAsamp[,8]==x)/NS
			return(apply(FRxx,2,mean))
		}),cbind(10,10,10,10))
		FRx<-subset(FRx ,apply(FRx,1,max)>0.01)
		
		sapply(c(1:length(FRx[,1])),function(x){
			points(Tmsh, FRx[x,]/mean(FRx[x,]),type="l",col= COLclax)
		})
	})	
	
#Arrivé ici

	write.table(type,"Metadata-tree.txt")
	
#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=0.88

#Compression of tree from the root (>1) 
COMP=1.1

#Margin size
MAR=1.5

#tip attribute coordinates (concentric circles)
Ai=1.01+seq(0.05,1,length.out=10)
Ae=Ai+0.8*(Ai[2]-Ai[1])

#Number of graduations in the PS scale
Nscale=6

	#Branch color
	COLbr<-rgb(0.7,0.7,0.7)
	

	Lmin=min(c(coordx[,5], coordx[,7]))
	Lmax=max(coordx[,5])*COMP

	#Limit of the graphic
	rLIM=MAR* Lmax
	z=max(c(coordx[,1], coordx[,2],coordx[,3]))
	
	#Start of the plot
	pdf(paste("Circular-tree_MINBOOT=",bx,".pdf",sep=""))
	par(mar=c(1,1,1,1),bty="n",mfrow=c(1,1))

	#Frame
	plot(-10000,-10000,xlim=c(-rLIM, rLIM),ylim=c(-rLIM, rLIM),xlab="",ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n")
	
	#TREE
	COORD<-t(sapply(c(1:length(coordx[,1])),function(x){
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
		points(XXn, YYn, type="l",col=COLbr<-rgb(0.7,0.7,0.7),lwd=1)
		segments(XXb1, YYb1,XXb2, YYb2,col=COLbr<-rgb(0.7,0.7,0.7),lwd=1)
		#points(XXb1, YYb1,pch=19,col=rgb(0,0,0,1),cex=0.5)
		return(c(x,XXb1, YYb1))
	}))
	
	
	#Node attributes (bootstrapp)
	Nb=11
	Nbl=3
	
	Nx<-c(1:length(PSn[,1]))
	colx<-seq(0.8,0,length.out=Nb)
	colx<-rgb(colx, colx, colx)
	
	sapply(Nx,function(x){
		bx=PSn[x,8]
		bx<-ifelse(is.na(bx)==T,0,bx)
		
		cx<-subset(colx,seq(0,1,length.out=Nb)==round(bx,1))
		
		nx=PSn[x,2]
		COORDxx<-subset(COORD,coordx[,4]==nx)
		
		points(COORDxx[1,2], COORDxx[1,3],pch=19,col= cx,cex=bx/3)
	})
	
	#Bootstrapp legend

	
	sapply(c(1:Nb),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out=Nb)[x]
		
		bbx=seq(bx,1,length.out=Nb)[x]
		points(X,0.86*rLIM,pch=19,col=colx[x],cex=bbx/3)			
	})
	text(0.775*rLIM,0.9*rLIM,"Node support",cex=0.6,font=1)
	sapply(c(1: Nbl),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out= Nbl)[x]
		text(X,0.86*rLIM,round(seq(bx*100,100,length.out=Nbl)[x],0),cex=0.5,pos=1,adj=1)			
	})	
	
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
	PSscale<-PSscale[c(1:length(Xscale))]
	
	segments(-0.01*rLIM,min(Xscale),
		-0.01*rLIM,max(Xscale))
	text(-0.005*rLIM, Xscale, round(PSscale,2)
		,cex=0.6,pos=2,adj=1,font=1)
	text(-0.17*rLIM, mean(Xscale), "PS",cex=0.8,font=1)
	
	maxTIP=max(Xscale)
	
	
	#####PLOT Tips attributes
	
	Nlab=12 #maximum number of labels in legend
	l1=-0.05 #maximum x coord for legend
	l2=-0.75 #minimum x coord for legend
		
	#######environment of sampling
	
	ATx=2
	
	#Display biome of origin for each strain
	colECO<-cbind(c("PLANT","PHYLLOSPHERE","ENDOPHYTE","SEED","RHIZHOSPHERE","SOIL","LICHEN-FUNGI","WATER","OCEAN","MICROBIOME"),
	c("green3","green","yellow3","yellow","brown","red","orange","cyan","blue","purple"),
	c("Plant","Phyllo.","Endop.","Seed","Rhizo.","Soil","Fung.","Water","Ocean","Gut"))	
	
	
	sapply(Tx,function(x){
		typex<-unlist(strsplit(type[x,4],split="/"))
		antx<-length(subset(typex, typex=="ANTHROPO"))
		typex<-setdiff(typex,"ANTHROPO")
		colECOx<-subset(colECO[,2], colECO[,1]== typex)
		ang1=2*pi*((x-0.5)/z)*OPEN
		ang2=2*pi*((x+0.5)/z)*OPEN
		XX1i<-Ai[ATx]* maxTIP*sin(ang1)
		YY1i<-Ai[ATx]* maxTIP*cos(ang1)
		XX2i<-Ai[ATx]* maxTIP*sin(ang2)
		YY2i<-Ai[ATx]* maxTIP*cos(ang2)
		XX1e<-Ae[ATx]* maxTIP*sin(ang1)
		YY1e<-Ae[ATx]* maxTIP*cos(ang1)
		XX2e<-Ae[ATx]* maxTIP*sin(ang2)
		YY2e<-Ae[ATx]* maxTIP*cos(ang2)		
		polygon(c(XX1i,XX2i,XX2e,XX1e),
			c(YY1i,YY2i,YY2e,YY1e),
			col= colECOx,border=ifelse(antx==1,"black",NA),lwd=0.5)
	})
	
	#Scale first tip attribute
	Ns=length(colECO[,1])+1	
	colx=c(colECO[,2],"white")
	bordx<-c(sapply(c(1:(Ns-1)),function(x){
		return(NA)
	}),"black")
	
	y1=Ae[ATx]
	y2=(Ae-((Ae-Ai)/2))[ATx]

	Nmin=min(sapply(c(1:Ns),function(x){
		x1=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x]
		x2=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x+1]
		polygon(c(x1,x1,x2,x2),
			c(y1, y2, y2, y1)*maxTIP,
			col= colx[x],border= bordx[x],lwd=0.5)
		text((x1+x2)/2,(y1 + y2)*maxTIP/2,
			c(colECO[,3],"Anthr.")[x],cex=0.35)	
		return(min(c(x1,x2)))				
	}))
	
	text(Nmin, ((y1+y2)/2)*maxTIP, "Environment"
		,cex=0.6,pos=2,adj=1,font=1)	
	
	##Variables for strains only: SITE (site of sampling), SPE (host tree species), TEMP (temperature of isolation) and TIME (day of sampling) 
	
	VAR<-c("SITE","SPE","TEMP","DAY14")
	
	##Color codes
	Tu<-sort(setdiff(as.numeric(type[,c(subset(c(1:length(type[1,])),colnames(type)=="DAY14"))]),c(NA)),decreasing=T)
		
	COLtu<-(Tu-min(Tu))/(max(Tu)-min(Tu))
	COLtu<-rgb(1-COLtu,1-COLtu,0)
	
	VARc<-t(rbind(
	cbind(1,"Sampling site",
		c("MSH","SBL"),
		c("green3","orange"),
		c("black","black"),
		c("MSH","SBL")),
	
	cbind(2,"Tree species",
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU"),
		c("purple2","violet","blue","green",rgb(0.35,0.35,0.35),
		"orange3","red3"),
		c("black","black","black","black","black","black","black"),
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU")),
	
	cbind(3,"Isolation temperature",
		c(20,30),
		c("cyan2","red2"),
		c("black","black"),
		c("20°C","30°C")),
	
	cbind(4,"Sampling time",Tu,COLtu,
		c("white",NA,NA,NA,"black"),
		c("Oct.","","","","June"))		
	))

	row.names(VARc)<-c("ring","fac","var","col","label.col","label")
	
	
	sapply(c(1:length(VAR)),function(x){
		X=x
		Aix<-Ai[X+1]
		Aex<-Ae[X+1]
		varx<-VAR[X]
		
		VARcx<-subset(t(VARc),as.numeric(VARc[1,])==X)
		typex<-as.vector(subset(t(type),colnames(type)==varx))
		
		colx<-as.vector(VARcx[,4])
		
		colx<-as.vector(sapply(typex,function(x){
			return(c(subset(as.vector(VARcx[,4]),
				as.vector(VARcx[,3])==x),NA)[1])
		}))
		ATx=X+2
		
		#Strain label
		sapply(Tx,function(x){
			cx<-colx[x]
			ang1=2*pi*((x-0.5)/z)*OPEN
			ang2=2*pi*((x+0.5)/z)*OPEN
			XX1i<-Ai[ATx]* maxTIP*sin(ang1)
			YY1i<-Ai[ATx]* maxTIP*cos(ang1)
			XX2i<-Ai[ATx]* maxTIP*sin(ang2)
			YY2i<-Ai[ATx]* maxTIP*cos(ang2)
			XX1e<-Ae[ATx]* maxTIP*sin(ang1)
			YY1e<-Ae[ATx]* maxTIP*cos(ang1)
			XX2e<-Ae[ATx]* maxTIP*sin(ang2)
			YY2e<-Ae[ATx]* maxTIP*cos(ang2)		
			polygon(c(XX1i,XX2i,XX2e,XX1e),
				c(YY1i,YY2i,YY2e,YY1e),
				col= cx,border=NA,lwd=0.5)
		})
		
		#Legend
		VARx<-unique(as.vector(VARcx[,2]))
		
		labx<-as.vector(VARcx[,6])
		colx<-as.vector(VARcx[,4])
		lab.colx<-as.vector(VARcx[,5])
			
		Ns=length(labx)
	
		y1=Ae[ATx]
		y2=(Ae-((Ae-Ai)/2))[ATx]
	
		Nmin=min(sapply(c(1:Ns),function(x){
			x1=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x]
			x2=seq(l1*rLIM,l2*rLIM,length.out= Nlab)[x+1]	
			polygon(c(x1,x1,x2,x2),
				c(y1,y2,y2,y1)*maxTIP,
				col=colx[x],border=NA)
			text(mean(c(x1,x2)),
				mean(c(y1,y2))*maxTIP, 
				labx[x],col=lab.colx[x],cex= 0.35)
			return(min(c(x1,x2)))				
		}))
		text(Nmin, ((y1+y2)/2)*maxTIP, VARx
		,cex=0.6,pos=2,adj=1,font=1)			
	})

	
	#Clade
	
	ATx=7
	
	CLA<-type[,c(subset(c(1:length(type[1,])),colnames(type)=="clade"))]
	
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
		cbind("A5","yellow3",3,3,0.5),
		cbind("A7","pink",3,3,0.5),
		cbind("A10","green3",2,7,0.8),
		cbind("A8","purple",3,3,0.5))
	
	sapply(unique(CLA),function(x){
		cx<-subset(COLcla[,2], COLcla[,1]==x)
		fontx<-as.numeric(subset(COLcla[,3], COLcla[,1]==x))
		ATx <-as.numeric(subset(COLcla[,4], COLcla[,1]==x))
		cexx<-as.numeric(subset(COLcla[,5], COLcla[,1]==x))
		zx<-subset(Tx, CLA==x)
		
		#remove isolated tips
		zx<-subset(zx,c(sapply(c(1:(length(zx)-1)),function(x){
			return(zx[x+1]-zx[x])
		}),1)==1)
		
		z1<-min(zx)
		z2<-max(zx)
		zx<-seq(z1,z2,length.out=1000)
		zm=mean(zx)
		
		yz<-((Ai+Ae)[ATx])/2
		ym<-Ae[ATx+1]
		
		angz=2*pi*(zx/z)*OPEN
		angm=2*pi*(zm/z)*OPEN
			
		XX1i<-yz* maxTIP*sin(angz)
		YY1i<-yz* maxTIP*cos(angz)
		XXm<-ym* maxTIP*sin(angm)
		YYm<-ym* maxTIP*cos(angm)		
	
		points(XX1i,
				YY1i,
				col= cx,lwd=3,type="l")
				
		text(XXm,YYm,x,font=fontx,col=cx,cex=cexx)
		
	})
	
	#Title
	text(-0.01*rLIM,Ae[7]*maxTIP,
		"Clades",font=1,cex=0.6,pos=2,adj=1)
	
	
	######Sampling year and strain name/clade
	
	batchx<-as.vector(subset(t(type),colnames(type)=="BATCH"))
	batchx<-ifelse(is.na(batchx)==T,"",ifelse(batchx=="Iso-2018",2017,2018))
	pchx<-ifelse(is.na(batchx)==T,0,
		ifelse(batchx==2018,19,
		ifelse(batchx==2017,21,0)))
	fontx<-ifelse(batchx=="",3,
		ifelse(batchx==2018,2,
		ifelse(batchx==2017,2,NA)))
	
	y1=(Ai[1]+1)/2
	
	#Strain label
	
	ang.lab<-seq(360,(1-OPEN)*360,length.out=length(Tx))-270
	adjx<-ifelse(ang.lab>-90,0, 1)
	ang.lab<-ifelse(ang.lab>-90,ang.lab, ang.lab+180)
	
	sapply(Tx,function(x){
		cx<-subset(COLcla[,2], COLcla[,1]== CLA[x])
		nx<-ifelse(batchx[x]=="",namex[x],paste(namex[x]," (",batchx[x],")",sep=""))
		angx=2*pi*(x/z)*OPEN
		XX1i<-y1* maxTIP*sin(angx)
		YY1i<-y1* maxTIP*cos(angx)						
		text(XX1i, YY1i,nx,cex=0.2,srt= ang.lab[x],adj= adjx[x],col=cx,font=fontx[x])	
			
	})
	
	#Title
	text(-0.01*rLIM,Ae[1]*maxTIP,
		"Isolated strains (year)",font=2,cex=0.6,pos=2,adj=1)
	text(-0.01*rLIM,Ai[1]*maxTIP,
		"Reference strains",font=3,cex=0.6,pos=2,adj=1)	
	
	
	dev.off()

