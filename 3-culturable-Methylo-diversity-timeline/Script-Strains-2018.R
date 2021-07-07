#Draw a circular phylogenetic tree from newick format and test for node association with environmental variables 

library(ape)
library(seqinr)
library(TreeTools)
library(ips)

###FUNCTIONS
num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

empty.plot<-function(x1,x2,y1,y2,Nx,Ny){
	plot(-1000,-1000,
	xlim=c(x1,x2),ylim=c(y1,y2),
	xaxt="n",yaxt="n",xlab=Nx,ylab=Ny)
}


#########################################
#########################################
## A - convert newick tree and find ##### 
###### best correlation between tree ####
###### topology and pairwise similarity # 
###### (PS) #############################
#########################################
#########################################

######################
###work directories

#Directory for data
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"

##load tree in newick format and corresponding fasta file
setwd(pD)

###Newick file (with Nodal support values, no branch length) 
treex<-read.tree("rpoB_mrbayes_consensus.nwk")

#round Nodal support
treex$node.label<-round(as.numeric(treex$node.label))

#Collapse nodes with less of bx support
bx=30

#node names
treex<-collapseUnsupportedEdges(treex,value="node.label",cutoff= bx)


#Root on Microvirga/Enterovirga (outgroup)

ox<-c(grep("Microvirga",treex$tip.label),
	grep("Enterovirga",treex$tip.label))
ox<-c(min(ox):max(ox))

treex<-root(treex,outgroup=treex$tip.label[ox])

#add branch length (Grafen method; https://pubmed.ncbi.nlm.nih.gov/2575770/)
treex <-compute.brlen(treex,power=0.5)

###Corresponding Fasta file
fastax<-read.fasta("rpoB_Complete-ref+Partial-strains.fas")

#Extract node coordinates, sequence names and bootstrapp values from the newick tree
coordx<-treex$edge #tree topology: nodes (first column) and their daugther nodes. Daughther node include tips (tips have indexes 1 to n and nodes, n+1 to n+N where n is the number of tips and N the number of nodes, tips excluded)
namex<-treex$tip.label #sequence (tip) indexes
bootx<-treex$node.label #bootstrapp value for each node
nodex<-sort(unique(coordx[,1])) #nodes indexes

#read p-distance matrix among sequences (all sites) Was calculated in MEGA7 (all positions included)

pdist<-read.table("rpoB_Complete-ref+Partial-strains-pdistance-formated.txt",header=T)


#Sequences characteristics (in column Name_in_pdist, matching names in the pdist file, spaces were replaced by underscores). Column fasta_names matches the fasta file names. Column Name_in_phyllo matches the tree file (namex)

CHAR<-read.table("sequence_characteristics.txt",header=T)

LAB<-as.vector(sapply(treex$tip.label,function(x){
	CHARx<-subset(CHAR,CHAR[,1]==x)
	return(paste(CHARx[,4], CHARx[,5], CHARx[,6], ifelse(is.na(CHARx[,8])==T,"",CHARx[,8])))
}))

#clade color
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
		cbind("A5b","yellow3",3,7,0.5),
		cbind("A7","pink",3,3,0.5),
		cbind("A10","green3",2,7,0.8),
		cbind("A8","purple",3,3,0.5))

#color and font of the tip label
COL.FONT<-t(sapply(treex$tip.label,function(x){
	CHARx<-subset(CHAR,CHAR[,1]==x)
	
	lab1x<-ifelse(is.na(CHARx[,8])==T,
		1,3)
	lab2x<- ifelse(is.na(CHARx[,8])==T,
		CHARx[,4],
		paste(CHARx[,4], 
		CHARx[,5], CHARx[,6]))
	
	colx<-ifelse(is.na(CHARx[,8])==T,"black",
		subset(COLcla[,2],
		COLcla[,1]==CHARx[,8]))
	
	return(c(lab1x, lab2x,colx))
}))



treex$tip.label<-as.vector(COL.FONT[,2])
treex2<-treex
treex2$tip.label<-ifelse(as.vector(COL.FONT[,1])=="1","O","")

pdf("tree-Fig.pdf")

par(mar=c(0,0,0,0))
plot.phylo(treex,"fan" ,show.tip.label=T,col="white",show.node.label = F, edge.color = "white",cex=0.25,tip.color=COL.FONT[,3],label.offset=0.05,x.lim=c(-1.3,1.3),y.lim=c(-1.3,1.3))
par(new=T)
plot.phylo(treex2,"fan" ,show.tip.label=T,col="grey",show.node.label = T, edge.color = "grey",cex=0.25,tip.color="black",font=1,x.lim=c(-1.3,1.3),y.lim=c(-1.3,1.3))


dev.off()

###Each node will be associated to a level, according to its hierachical position in the tree 
#number of levels (must be high enough to cover all the nodes)
NL=50
#number of levels to add between each level (will be used to adjust levels when searching for the best correlation between pdistance and levels)
NLex=10


#Sequence indexes
tips<-c(1:length(namex))

###Build a matrix with each sequence in rows and each level in column; fill the matrix iteratively for each level with nodes names (emprty positions as 0)

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

#stack the matrix to the root, so each node is associated to a single level 
tempx<-sapply(c(1:length(tempx[,1])),function(x){
	tx<-as.numeric(tempx[x,])
	return(c(tx[1],
		subset(tx[2:length(tx)],
			tx[2:length(tx)]==0),
			subset(tx[2:length(tx)],
			tx[2:length(tx)]!=0)))
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
		print(x)
		stx<-subset(tempx[,1], tempx[,X]==x)
		#Names in nwk file (strain names)
		nx<-namex[stx]
		#Corresponding names in pdist file
		npx<-as.vector(sapply(nx,function(x){
			return(subset(CHAR[,2],CHAR[,1]==x))
		}))
		#retrieve pdistance among these strains
		pdistI<-as.vector(sapply(npx,function(x){
			return(subset(pdist[,1],pdist[,2]==x))
		}))
		
		PSx<-as.vector(unlist(matrix(pdist[pdistI,(pdistI+2)],ncol=1)))
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

##Number of iterations to find the best correlation
K=3000

pdf(paste("Tree-iterations_MINBOOT=",bx,".pdf",sep=""))

par(mfrow=c(2,2),mar=c(4,4,2,1),bty="n")
plot(log(PS[,1]), PSc,log="",
	xlab = "Level (log)",
	ylab="modified PS",las=1,cex.axis=0.8,pch=19,cex=0.5,main="before iterations")

par(new=T)
empty.plot(0,1,0,1,"","")	
text(0.25,0.25,paste("r2 =",round(CORi,4)),font=2)	
	


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

pdf(paste("Tree2.pdf",sep=""),height=12,width=6)

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
	bx<-ifelse(length(bx)==0,0,ifelse(is.na(bx)==T,0,bx/100))
	points(Xb1,Yb,pch=19,col=rgb(bx,0,0,1),cex=bx/2)
})

LABn<-as.vector(sapply(namex,function(x){
	CHARx<-subset(CHAR,CHAR[,1]==x)
	return(paste(CHARx[,4], CHARx[,5], CHARx[,6], ifelse(is.na(CHARx[,8])==T,"",CHARx[,8])))
}))


#Sequence labels

text(seq(Lmin,Lmin,length.out=length(namex)),
	c(1:length(namex)),
	labels= LABn,cex=0.2,pos=4,adj=1)	
	

dev.off()


#########################################
#########################################
## B - draw a circular tree with metadata
 #############################
#########################################
#########################################


library(ape)
library(seqinr)

###Load data
#Directory for data
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"
	#strain (tips of the tree) names
	setwd(pD)
	
	#minimum bootstrapp value
	bx=30
	
	namex <-as.vector(read.table("tree-names.txt")[,1])
	#namex <-read.tree(paste("rpoB_Complete-ref+Partial-strains_ML100-2.nwk",sep=""))[2]$tip.label
	
	#strains information
	strains<-read.table("strain-metadata.txt",header=T)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	#add references
	refn<-setdiff(namex,gsub("-","_",strains[,1]))
	refn<-cbind(refn,"ref",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
	colnames(refn)<-colnames(strains)
	
	strains<-cbind(gsub("-","_",strains[,1]),strains[,c(2:length(strains[1,]))])
	
	colnames(strains)<-colnames(refn)
	
	strains<-rbind(strains,refn)
	rm(refn)
	
	#tree information
	coordx<-read.table(paste("Tree-coordinates.txt",sep=""),header=T)
	SCALE <-read.table(paste("SCALE.txt",sep=""),header=T)
	PSn <-read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),header=T)


	#Tip (strain) attributes
	Tx<-c(1:length(namex))
	
	type<-t(sapply(Tx,function(x){
		
		#strain names
		stx=namex[x]
		
		CHARx<-subset(CHAR,CHAR[,1]== stx )
		
		STRAINx<-subset(strains,strains[,1]== stx)
		
		##Characteristics for strains and references
		#clade et reference (if not, strain)
		clx=c(as.vector(CHARx[,8]),as.vector(STRAINx[,2]))
		typex<-ifelse(clx[2]=="ref","ref","strain")
		clx<-setdiff(clx,c("ref",NA))
		#environment
		envx<-as.vector(CHARx[,7])
		#genus
		genx<-as.vector(CHARx[,5])
		#species
		spex<-as.vector(CHARx[,6])
		
		###sampling and isolation characteristics for straina only
		#SITE (site of sampling)
		#PLOT (plot of sampling within site)	
		#TREE (sampled tree)
		#SPE (host tree species)
		#BATCH of isolation
		#SAMPLE
		#Temperature of isolation
		#day of sampling (on a subsective scale)
		#day of sampling (rounded for 2 weeks)
		#date of sampling (actual)
		STRAINx<-as.vector(t(STRAINx))[3:12]
		
		return(as.vector(c(stx,typex,clx,envx,genx,spex,STRAINx)))
	}))
	colnames(type)<-c("Strain","type","clade","env","genus","species","SITE","PLOT","TREE","SPE","BATCH","SAMPLE","TEMP","DAY","DAY14","DATE")

	write.table(type,"Metadata-tree.txt")
	
#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=0.88

#Compression of tree from the root (>1) 
COMP=1.1

#Margin size
MAR=1.55

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
		bx=PSn[x,8]/100
		bx<-ifelse(is.na(bx)==T,0,bx)
		
		cx<-subset(colx,seq(0,1,length.out=Nb)==round(bx,1))
		
		nx=PSn[x,2]
		COORDxx<-subset(COORD,coordx[,4]==nx)
		
		points(COORDxx[1,2], COORDxx[1,3],pch=19,col= cx,cex=bx/3)
	})
	
	#Bootstrapp legend

	
	sapply(c(1:Nb),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out=Nb)[x]
		
		bbx=seq(bx/100,1,length.out=Nb)[x]
		points(X,0.86*rLIM,pch=19,col=colx[x],cex=bbx/3)			
	})
	text(0.775*rLIM,0.9*rLIM,"Nodal support",cex=0.6,font=1)
	sapply(c(1: Nbl),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out= Nbl)[x]
		text(X,0.86*rLIM,round(seq(bx,100,length.out=Nbl)[x],0),cex=0.5,pos=1,adj=1)			
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
	c("green3","green","yellow3","yellow2","brown","red","orange","cyan","blue","purple"),
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
	
	#clade color, 	COLcla<-rbind(
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
		cbind("A5b","yellow3",3,7,0.5),
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
	LABn<-as.vector(sapply(namex,function(x){
		CHARx<-subset(CHAR,CHAR[,1]==x)
		return(CHARx[,4])
	}))
	
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
		nx<-ifelse(batchx[x]=="", LABn[x],paste(LABn[x]," (",batchx[x],")",sep=""))
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

#########################################
#########################################
## C - for each level redefine taxonomic divisions according to significant bootstrapp in the tree and perform PERMANOVA on metadata,   ################
#########################################
#########################################

library(vegan)
library(ape)
library(seqinr)

#Directory for data
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"

#Minimum bootstrapp value to keep nodes
bx=30

#Number of permutations for PERMANOVA
K=10000
	
	#tree information
	setwd(pD)
	#coordx<-read.table(paste("Tree-coordinates.txt",sep=""),header=T)
	SCALE <-read.table(paste("SCALE.txt",sep=""),header=T)
	namex <-as.vector(read.table("tree-names.txt")[,1]) #strain (tips of the tree) names
	
	PSn <-read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),header=T)
	type<-read.table("Metadata-tree.txt",header=T)
	
	#Tree topology in level format
	LEVELS<-read.table(paste("Levels_MINBOOT=",bx,".txt",sep=""))
	colnames(LEVELS)<-paste("L",c(1:(length(LEVELS[1,]))),sep="-")
	
	#remove empty levels
	LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))
	
	#for each node in the level table, retrieve the bootstrapp value and replace node by 0 if bootstrapp is lower than bx
	#LEVELb<-cbind(LEVELS[,1],sapply(c(2:length(LEVELS[1,])),function(x){
	#	Lx<-as.vector(LEVELS[,x])
	#	return(sapply(Lx,function(x){
	#		return(x*length(intersect(x,PSb[,2])))
	#	}))	
	#}))
	#colnames(LEVELb)<-colnames(LEVELS)
	
	#remove empty levels
	#LEVELb <-t(subset(t(LEVELb),apply(LEVELb,2,max)!=0))
	
	#for each sequence, fill empty levels with node from the closest level
	
	LEVELf<-t(sapply(c(1:length(LEVELS[,1])),function(x){
		Sx=x
		Lx<-as.vector(LEVELS[Sx,])
		lx<-c(1:length(Lx))
		Lx<-sapply(lx,function(x){
			return(setdiff(
				subset(Lx,lx>=x),0)[1])
		})
		return(c(Sx,Lx[2:length(Lx)]))
	}))
	
	colnames(LEVELf)<-colnames(LEVELS)
	 
	
###For each level, perform a PERMANOVA
#use year of isolation as strata and test different models with combinations of following variables:

# SITE (site of sampling)
# HOST (host tree species)
# TEMP (temperature of isolation)
#DAY14 (day of sampling, simplified over 14 days)
# BATCH OF isolation (strata)


####Categories (considered as samples)
#CAT<-paste(type[,7],type[,8],type[,9],type[,10],type[,13],type[,15],ifelse(is.na(type[,11])==T,NA,ifelse(type[,11]=="Iso-2018",2018,2019)))

CAT<-paste(type[,7],type[,8],type[,9],type[,10],type[,13],type[,15],ifelse(is.na(type[,11])==T,NA,ifelse(type[,11]=="Iso-2018",2018,ifelse(type[,11]=="Iso-2019-1",2019,2020))))


CAT1<-setdiff(sort(unique(CAT)),"NA NA NA NA NA NA NA")

ENV1<-as.data.frame(t(sapply(CAT1,function(x){
	return(unlist(strsplit(x,split=" ")))
})))

colnames(ENV1)<-c("SITE","PLOT","TREE","HOST","TEMP","DAY","YEAR")

###Permanova with clades
CLAu<-as.vector(sort(unique(type[,3])))

####For each clade, retrieve number of strains per category
	NC<-sapply(CLAu,function(x){
		Cx=x
		#Strain index in this node
		STx<-subset(c(1:length(type[,1])),type[,3]== Cx)
		#Metadata for these strains
		typex<-type[STx,]
		type1<-CAT[STx]
		#retrieve number of strains for each category within this node
		return(sapply(CAT1,function(x){
			return(length(subset(type1, type1==x)))
		}))	
	})
	#remove nodes without observation
	NC <-as.data.frame(t(subset(t(NC),apply(NC,2,max)!=0)))
	
	write.table(NC,"exemple-frequencies.txt")
	write.table(ENV1,"exemple-factors.txt")
	#Number of strains per node per combination of factors
	NC<-as.data.frame(read.table("exemple-frequencies.txt"),header=T)
	#Factors
	ENV1 <-as.data.frame(read.table("exemple-factors.txt"),header=T)	

	#Complete model	
	MOD<-adonis(NC ~ SITE*HOST*TEMP*DAY, data= ENV1, permutations=K,strata= ENV1$YEAR,method="bray")$aov.tab

	write.table(MOD,"Strains-PERMANOVA-per-clade.txt")

###Permanova on diversity per level
LV=c(2:length(LEVELf[1,]))


#Number of strains per node per level (references excluded)
STn<-sapply(LV,function(x){
	#print(x)
	Lx=x
	ND<-as.vector(LEVELf[, Lx])
	ND<-subset(ND,type[,2]!="ref")
	#strains per node 	
	stn<-summary(as.factor(ND),
		max.sum=length(unique(ND)))
	paste(stn,collapse="_")	
})

#remove redondant levels
LV<-as.vector(sapply(unique(STn),function(x){
	LVx<-subset(LV,STn==x)
	return(LVx[length(LVx)])
}))	

#####################
### COMPLETE MODEL PER LEVEL("SITE","HOST","TEMP","DAY","YEAR")
AOV<-sapply(LV,function(x){
	print(x)
	Lx=x
	
	#Level
	LEVELx<-as.numeric(unlist(strsplit(colnames(LEVELf)[x],split="-"))[2])
	
	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	NDu<-setdiff(unique(ND),0)
	
	####For each node, retrieve number of strains per category
	N1<-sapply(NDu,function(x){
		NDx=x
		#Strain index in this node
		STx<-subset(LEVELf[,1],ND== NDx)
		#Metadata for these strains
		typex<-type[STx,]
		type1<-CAT[STx]
		#retrieve number of strains for each category within this node
		return(sapply(CAT1,function(x){
			return(length(subset(type1, type1==x)))
		}))	
	})
	colnames(N1)<-paste("node", NDu,sep="-")
	#remove nodes without observation
	N1<-as.data.frame(t(subset(t(N1),apply(N1,2,max)!=0)))
			
	write.table(N1,"exemple-frequencies.txt")
	write.table(ENV1,"exemple-factors.txt")
	#Number of strains per node per combination of factors
	N1<-as.data.frame(read.table("exemple-frequencies.txt"),header=T)
	#Factors
	ENV1 <-as.data.frame(read.table("exemple-factors.txt"),header=T)	
		
	MOD<-adonis(N1 ~ SITE*HOST*TEMP*DAY, data= ENV1, permutations=K,strata= ENV1$YEAR,method="bray")$aov
	
	return(c(dim(N1), MOD[,5], MOD[,6]))
	
})

DESC<-rownames(MOD)


rownames(AOV)<-c("Fac","Nodes",paste("R2", DESC),paste("Pr(>F)", DESC))
colnames(AOV)<-colnames(LEVELf)[LV]


write.table(AOV,paste("Permanova-tree-min-bootstrap=",bx,".txt",sep=""))

##########################
##########################
### D # Test for signifant association between diversity and factors for each level of the tree (permutations): Site, Host tree species, Temperature of isolation, sampling time
##########################
##########################

library(ape)

######################
###work directories
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"


	#Minimum bootstrapp value
	bx=30
	
	#tree information
	setwd(pD)
	#coordx<-read.table(paste("Tree-coordinates.txt",sep=""),header=T)
	SCALE <-read.table(paste("SCALE.txt",sep=""),header=T)

	namex <-as.vector(read.table("tree-names.txt")[,1]) #strain (tips of the tree) names
	
	PSn <-read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),header=T)
	type<-read.table("Metadata-tree.txt",header=T)

	#only nodes supported by at least bx bootstrapps
	#PSb<-subset(PSn,PSn[,8]>=bx)
	
	#Tree topology in level format
	LEVELS<-read.table(paste("Levels_MINBOOT=",bx,".txt",sep=""))
	colnames(LEVELS)<-paste("L",c(1:(length(LEVELS[1,]))),sep="-")
	
	#remove empty levels
	LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))
	
	
	#for each sequence, fill empty levels with node from the closest level
	
	LEVELf<-t(sapply(c(1:length(LEVELS[,1])),function(x){
		Sx=x
		Lx<-as.vector(LEVELS[Sx,])
		lx<-c(1:length(Lx))
		Lx<-sapply(lx,function(x){
			return(setdiff(
				subset(Lx,lx>=x),0)[1])
		})
		return(c(Sx,Lx[2:length(Lx)]))
	}))
	
	colnames(LEVELf)<-colnames(LEVELS)
	 
	
###Factors to test for association by permutation
# SITE (site of sampling)
# HOST (host tree species)
# TEMP (temperature of isolation)
# DAY14 (day of sampling, simplified over 14 days)
# BATCH OF isolation (strata)

#CATo<-paste(type[,7],type[,8],type[,9],type[,10],type[,13],type[,15],ifelse(is.na(type[,11])==T,NA,ifelse(type[,11]=="Iso-2018",2018,2019)))

CATo<-paste(type[,7],type[,8],type[,9],type[,10],type[,13],type[,15],ifelse(is.na(type[,11])==T,NA,ifelse(type[,11]=="Iso-2018",2018,ifelse(type[,11]=="Iso-2019-1",2019,2020))))



CAT1<-setdiff(sort(unique(CATo)),"NA NA NA NA NA NA NA")

ENV1<-as.data.frame(t(sapply(CAT1,function(x){
	return(unlist(strsplit(x,split=" ")))
})))

colnames(ENV1)<-c("SITE","PLOT","TREE","HOST","TEMP","DAY","YEAR")

###levels
LV=c(2:length(LEVELf[1,]))


#Number of strains per node per level (references excluded)
STn<-sapply(LV,function(x){
	#print(x)
	Lx=x
	ND<-as.vector(LEVELf[, Lx])
	ND<-subset(ND,type[,2]!="ref")
	#strains per node 	
	stn<-summary(as.factor(ND),
		max.sum=length(unique(ND)))
	paste(stn,collapse="_")	
})

#remove redondant levels
LV<-as.vector(sapply(unique(STn),function(x){
	LVx<-subset(LV,STn==x)
	return(LVx[length(LVx)])
}))	

#Number of permutations
K=100000

#For each level, test for significant association between strains distribution among taxa and factors. For each tested factor, perform permutation of strains within combined categories of over factors


sapply(LV,function(x){
	
	Lx=x
	
	#Level
	LEVELx<-as.numeric(unlist(strsplit(colnames(LEVELf)[x],split="-"))[2])
	
	

	#Nodes in this level
	ND<-as.vector(LEVELf[, Lx])
	NDu<-setdiff(unique(ND),0)
	
	####For each node, retrieve number of strains per category
	
	N1<-sapply(NDu,function(x){
		NDx=x
		#Strain index in this node
		STx<-subset(LEVELf[,1],ND== NDx)
		#Metadata for these strains
		typex<-type[STx,]
		type1<-CATo[STx]
		#retrieve number of strains for each category within this node
		return(sapply(CAT1,function(x){
			return(length(subset(type1, type1==x)))
		}))	
	})	
	
	
	colnames(N1)<-paste("node", NDu,sep="-")
	#remove nodes without observation
	N1<-as.data.frame(t(subset(t(N1),apply(N1,2,max)!=0)))
		
	write.table(N1,"exemple-frequencies.txt")
	write.table(ENV1,"exemple-factors.txt")
	#Number of strains per node per combination of factors
	N1<-as.data.frame(read.table("exemple-frequencies.txt"),header=T)
	#Factors
	ENV1 <-as.data.frame(read.table("exemple-factors.txt"),header=T)	
	
	#number of nodes (taxa)
	NN=length(N1[1,])
	
	#Format data for permutations
	data<-matrix(unlist(sapply(c(1:NN),function(x){
		X=x
		NDx<-colnames(N1)[X]
		Nx<-N1[,X]
		ENVx <-subset(ENV1,Nx!=0)
		Nx<-subset(Nx,Nx!=0)
		return(t(matrix(unlist(sapply(c(1:length(Nx)),function(x){
			Y=x
			nx<-Nx[Y]
			ex<-as.vector(t(ENVx[Y,]))
			sapply(c(1:nx),function(x){
				return(c(NDx, ex))
			})
		})),ncol=8,byrow=T)))
	})),ncol=8,byrow=T)
	
	#Nodes
	MBx<-as.vector(data[,1])
	MB<-sort(unique(MBx))
		
###For node association, permutations should be done in order to remove data structure. For instance, when testing for association with sampling site, permutations as done within temp/year categories. If for instance we have 8 strains from nodes A and B: A1, A2, A3, A4, B1, B2, B3 and B4. Strains A1, A2, B3 and B4 come from site X, the others from site Y. Strains A1, A2, A3, A4 were isolated at 20°C, the others at 30°C. Strains A1 and B1 were sampled in 2017, the others in 2018, we obtain the following data table

# Strain	Node	Site	Temp	Year
# A1		A		X		20		2017
# A2		A		X		20		2018	
# A3		A		Y		20		2018
# A4		A		Y		20		2018
# B1		B		Y		30		2017
# B2		B		Y		30		2018
# B3		B		X		30		2018
# B4		B		X		30		2018

###If we want to test the association of Node with site, permutations have to we done within Temp/year categories, resulting in 4 categories

# Strain	Node	Site
#CAT1: 	20/2017
# A1		A		X
#CAT2: 	20/2018
# A2		A		X
# A3		A		Y
# A4		A		Y
#CAT3: 	30/2017
# B1		B		Y
#CAT4: 	30/2018
# B2		B		Y
# B3		B		X
# B4		B		X

##In categories CAT1 (strain A1) and CAT3 (strain B1), there is only one strain, so permutations won't affect expectation for these strains: they will alway be associated with the same site in permutations 



#####Test for association with SITE: permutations in TEMP/YEAR categories

	#Factor to test
	FAC<-as.vector(data[,2])

	#category within wich permutations have to be done to avoid interaction among factors
	CAT<-as.vector(paste(data[,6],data[,8],sep="_"))
	CATperm<-sort(unique(CAT)) 

	#Perform K permutation of nodes within categories
	PERM<-t(sapply(c(1:K),function(x){
		X=x	
		MBp<-MBx[unlist(sapply(CATperm,function(x){
			return(sample(subset(c(1:length(MBx)), CAT==x),
			size=length(subset(c(1:length(MBx)), CAT==x))))		
		}))]	
		FACp<-FAC[unlist(sapply(CATperm,function(x){
			return(subset(c(1:length(FAC)), CAT==x))		
		}))]	
		return(sapply(MB,function(x){
			FACx<-subset(FACp, MBp==x)
			return(sapply(unique(FAC),function(x){
				return(length(subset(FACx, FACx==x)))
			}))
		}))
	}))	
	
	#Real values
	REAL<-matrix(sapply(MB,function(x){
		FACx<-subset(FAC, MBx==x)
		return(sapply(unique(FAC),function(x){
			return(length(subset(FACx, FACx==x)))
		}))
	}),nrow=1)	
	
	#Expected values (permutations)
	EXP<-paste(round(apply(PERM,2,mean),2),round(apply(PERM,2,sd),2),sep="±")	
	
	#Factor_node association
	DESC<-matrix(sapply(MB,function(x){
		X=x
		return(sapply(unique(FAC),function(x){
			return(paste(x,X,sep="_"))
		}))
	}),nrow=1)
	
	#P value calculated from permutations for each association
	pu<-p.adjust((sapply(c(1:length(DESC)),function(x){
		return(	length(subset(PERM[,x], PERM[,x]>=REAL[x])))
	})+1)/(K+1), method="bonferroni")
	
	#summary file
	SUMsite<-cbind(as.vector(DESC) ,as.vector(REAL), EXP,round(pu,5))
	colnames(SUMsite)<-c("ass.","obs","exp","p.obs>exp")		
	


#####Test for association with TEMP: permutations in SITE/YEAR categories

	#Factor to test
	FAC<-as.vector(data[,6])

	#category within wich permutations have to be done to avoid interaction among factors
	CAT<-as.vector(paste(data[,2],data[,8],sep="_"))
	CATperm<-sort(unique(CAT)) 

	#Perform K permutation of nodes within categories
	PERM<-t(sapply(c(1:K),function(x){
		X=x	
		MBp<-MBx[unlist(sapply(CATperm,function(x){
			return(sample(subset(c(1:length(MBx)), CAT==x),
			size=length(subset(c(1:length(MBx)), CAT==x))))		
		}))]	
		FACp<-FAC[unlist(sapply(CATperm,function(x){
			return(subset(c(1:length(FAC)), CAT==x))		
		}))]	
		return(sapply(MB,function(x){
			FACx<-subset(FACp, MBp==x)
			return(sapply(unique(FAC),function(x){
				return(length(subset(FACx, FACx==x)))
			}))
		}))
	}))	
	
	#Real values
	REAL<-matrix(sapply(MB,function(x){
		FACx<-subset(FAC, MBx==x)
		return(sapply(unique(FAC),function(x){
			return(length(subset(FACx, FACx==x)))
		}))
	}),nrow=1)	
	
	#Expected values (permutations)
	EXP<-paste(round(apply(PERM,2,mean),2),round(apply(PERM,2,sd),2),sep="±")	
	
	#Factor_node association
	DESC<-matrix(sapply(MB,function(x){
		X=x
		return(sapply(unique(FAC),function(x){
			return(paste(x,X,sep="_"))
		}))
	}),nrow=1)
	
	#P value calculated from permutations for each association
	pu<-p.adjust((sapply(c(1:length(DESC)),function(x){
		return(	length(subset(PERM[,x], PERM[,x]>=REAL[x])))
	})+1)/(K+1), method="bonferroni")
	
	#summary file
	SUMtemp<-cbind(as.vector(DESC) ,as.vector(REAL), EXP,round(pu,5))
	colnames(SUMtemp)<-c("ass.","obs","exp","p.obs>exp")			
	PSx<-subset(SCALE[,4],SCALE[,2]== LEVELx)
	
	SUMall<-rbind(
		cbind("site",SUMsite),
		cbind("temp",SUMtemp))
		
	colnames(SUMall)<-c("test",colnames(SUMsite))	
	
	print(c("Level:",LEVELx,"significant tests:",length(subset(SUMall[,1], as.numeric(SUMall[,5])<=0.05))))
	
	
	write.table(SUMall,
		paste("PERM2018_LEVEL=", LEVELx,"_PS=", PSx,"_K=",K,"_boot=",bx,".txt",sep=""))
	
})


#Summary only with significant tests

SUMo<-cbind("none","none",0,0,0)

SUM<-matrix(unlist(sapply(LV,function(x){
	
	Lx=x
	LEVELx<-as.numeric(unlist(strsplit(colnames(LEVELf)[x],split="-"))[2])
	PSx<-subset(SCALE[,4],SCALE[,2]== LEVELx)
	
	SUMx<-read.table(paste("PERM2018_LEVEL=", LEVELx,"_PS=", PSx,"_K=",K,"_boot=",bx,".txt",sep=""),header=T)
	
	colnames(SUMo)<-colnames(SUMx)
	SUMx<-cbind(LEVELx ,PSx,rbind(SUMx,SUMo))
	colnames(SUMx)<-c("Level","PS",colnames(SUMx[3:length(SUMx[1,])]))
	
	return(t(subset(SUMx,SUMx[,7]<=0.05)))

})),ncol=7,byrow=T)

SUM<-subset(SUM,SUM[,3]!="none")


colnames(SUM)<-c("Level","PS","test","ass.","obs","exp", "p.obs.exp")

write.table(SUM,
		paste("PERM2018_summary_K=",K,"_boot=",bx,".txt",sep=""))

############################
############################
### E # Summary figure of permutation tests and permanova. For permutation, only display results for significant tests: site2 amd temp
############################
############################

library(ape)

#Directory for data
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"
	
	#Minimum bootstrapp value
	bx=30
	
	
	#tree information
	setwd(pD)
	coordx<-read.table(paste("Tree-coordinates.txt",sep=""),header=T)
	SCALE <-read.table(paste("SCALE.txt",sep=""),header=T)
	namex <-as.vector(read.table("tree-names.txt")[,1]) #strain (tips of the tree) names
	
	PSn <-read.table(paste("nodes-attributes_MINBOOT=",bx,".txt",sep=""),header=T)
	Tx<-c(1:length(namex))
	type<-read.table("Metadata-tree.txt",header=T)
	CHAR<-read.table("sequence_characteristics.txt",header=T)
	
	#Number of permutation (permutation test for node association)
	K=100000

	#Tree topology in level format
	LEVELS<-read.table(paste("Levels_MINBOOT=",bx,".txt",sep=""))
	colnames(LEVELS)<-paste("L",c(1:(length(LEVELS[1,]))),sep="-")
	
	#remove empty levels
	LEVELS<-t(subset(t(LEVELS),apply(LEVELS,2,max)!=0))
	
	#for each sequence, fill empty levels with node from the closest level
	
	LEVELf<-t(sapply(c(1:length(LEVELS[,1])),function(x){
		Sx=x
		Lx<-as.vector(LEVELS[Sx,])
		lx<-c(1:length(Lx))
		Lx<-sapply(lx,function(x){
			return(setdiff(
				subset(Lx,lx>=x),0)[1])
		})
		return(c(Sx,Lx[2:length(Lx)]))
	}))
	
	colnames(LEVELf)<-colnames(LEVELS)

###levels
LV=c(2:length(LEVELf[1,]))

#Number of strains per node per level (references excluded)
STn<-sapply(LV,function(x){
	#print(x)
	Lx=x
	ND<-as.vector(LEVELf[, Lx])
	ND<-subset(ND,type[,2]!="ref")
	#strains per node 	
	stn<-summary(as.factor(ND),
		max.sum=length(unique(ND)))
	paste(stn,collapse="_")	
})

STu<-unique(STn)

#Identify level with exactly the same nodes
LevND<-as.vector(sapply(STu,function(x){
	return(paste(subset(LV,STn==x),collapse="_"))
}))

LevNDx <-as.vector(sapply(LevND ,function(x){
	return(paste(as.numeric(
		setdiff(unlist(strsplit(
			colnames(LEVELf)[
			as.numeric(unlist(
			strsplit(x,split="_")))],
			split="-")),"L")),collapse="_"))
}))




#remove redundant levels
LV<-as.vector(sapply(STu,function(x){
	LVx<-subset(LV,STn==x)
	return(LVx[length(LVx)])
}))	

	
	#Significant associated nodes 
	SUM<-read.table(paste("PERM2018_summary_K=",K,"_boot=",bx,".txt",sep=""),header=T)

#########Graphical parameters
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=0.75

#Compression of tree from the root (>1) 
COMP=1.1

#Margin size
MAR=1.3

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

#####Parameters for Tips attributes
	
	Nlab=12 #maximum number of labels in legend
	l1=-0.05 #maximum x coord for legend
	l2=-0.75 #minimum x coord for legend	
	
	##Variables for strains only and significant factors in permutations tests only: SITE (site of sampling) and TEMP (temperature of isolation)
		
	TEST<-sort(as.vector(unique(SUM[,3])),decreasing=T) 
		
	nFAC=length(TEST)
		
	##Color codes
	Tu<-sort(setdiff(as.numeric(type[,c(subset(c(1:length(type[1,])),colnames(type)=="DAY14"))]),c(NA)),decreasing=T)
		
	COLtu<-(Tu-min(Tu))/(max(Tu)-min(Tu))
	COLtu<-rgb(1-COLtu,1-COLtu,0)
	
	VARc<-t(rbind(
	
	cbind("SITE","site","Forest of origin",
		c("MSH","SBL"),
		c("green3","orange"),
		c("black","black"),
		c("MSH","SBL")),	
	
	cbind("SPE","host","Tree species",
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU"),
		c("purple2","violet","blue","green",rgb(0.35,0.35,0.35),
		"orange3","red3"),
		c("black","black","black","black","black","black","black"),
		c("ACSA","FAGR","OSVI","QURU","ABBA","ACPE","ACRU")),
	
	cbind("TEMP","temp","Isolation temperature",
		c(20,30),
		c("cyan2","red2"),
		c("black","black"),
		c("20°C","30°C")),
	
	cbind("DAY14","day","Sampling time",Tu,COLtu,
		c("white",NA,NA,NA,"black"),
		c("Oct.","","","","June"))		
	))

	row.names(VARc)<-c("FAC","test","fac","var","col","label.col","label")

###Load PERMANOVA results

AOV<-read.table(paste("Permanova-tree-min-bootstrap=",bx,".txt",sep=""),header=T)

#Remove tests for Levels with less than 7 nodes (meaningless)

AOV<-t(subset(t(AOV),AOV[2,]>6))

#Levels
LV<-as.numeric(matrix(unlist(strsplit(colnames(AOV),split="[.]")),ncol=2,byrow=T)[,2])

#Graphic just with node labels
	pdf(paste("Tree-NodeLabels_boot=",bx,".pdf",sep=""))
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
	#Node attributes (name)
	Nx<-c(1:length(PSn[,1]))
	sapply(Nx,function(x){
		nx=PSn[x,2]
		COORDxx<-subset(COORD,coordx[,4]==nx)
		points(COORDxx[1,2], COORDxx[1,3],
			pch=19,col= rgb(0,0,0,0.3),cex=0.5)
		text(COORDxx[1,2], COORDxx[1,3],nx,
			cex=0.2,col="white",font=2)
	})
					
dev.off()


	
######Graphic

	#Start of the plot
	pdf(paste("Summary-tree_Permutations_and_AOV_boot=",bx,"_K=",K,".pdf",sep=""))
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
		bx=PSn[x,8]/100
		bx<-ifelse(is.na(bx)==T,0,bx)
		cx<-subset(colx,seq(0,1,length.out=Nb)==round(bx,1))
		nx=PSn[x,2]
		COORDxx<-subset(COORD,coordx[,4]==nx)
		points(COORDxx[1,2], COORDxx[1,3],pch=19,col= cx,cex=bx/3)
	})
	sapply(c(1:Nb),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out=Nb)[x]
		bbx=seq(bx/100,1,length.out=Nb)[x]
		points(X,0.86*rLIM,pch=19,col=colx[x],cex=bbx/3)			
	})
	text(0.775*rLIM,0.9*rLIM,"Nodal support",cex=0.6,font=1)
	sapply(c(1: Nbl),function(x){
		X=seq(0.7*rLIM,0.85*rLIM,length.out= Nbl)[x]
		text(X,0.86*rLIM,round(seq(bx,100,length.out=Nbl)[x],0),cex=0.5,pos=1,adj=1)			
	})		
	#PS Scale (shared between the tree and that the AOV graph)
	#Level in the scale
	#Lscale=round(seq(Lmin,Lmax,length.out=Nscale),0)
	Lscale=round(seq(Lmin,max(LV),length.out=Nscale),0)
	#Corresponding PS values
	PSscale<-unlist(sapply(Lscale,function(x){
		return(subset(as.vector(SCALE[,4]), as.vector(SCALE[,2])==x))
	}))
	xs1=-0.01
	xs2=-0.2
	Xscale<-sapply(c(1:Nscale),function(x){
		X=Lmax-Lscale[x]
		segments(xs1*rLIM,X,(xs1-0.04)*rLIM,X)
		segments(xs2*rLIM,X,(xs2+0.04)*rLIM,X)
		return(X)
	})
	maxTIP=max(Xscale)
	PSscale<-PSscale[c(1:length(Xscale))]
	segments(xs2*rLIM,min(Xscale),
		xs2*rLIM,max(Xscale))
	segments(xs1*rLIM,min(Xscale),
		xs1*rLIM,max(Xscale))	
	text(0.5*(xs1+xs2)*rLIM, Xscale, round(PSscale,2)
		,cex=0.5,font=1)
	text(0.5*(xs1+xs2)*rLIM, maxTIP*(Ai[1]+Ae[1])/2, "PS",cex=0.8,font=1)
	#Display Clades	
	CLA<-type[,c(subset(c(1:length(type[1,])),colnames(type)=="clade"))]
	#clade color
	COLcla<-rbind(
		cbind("Enterovirga",rgb(0.8,0.8,0.8),3,3,0.5),
		cbind("Microvirga",rgb(0.6,0.6,0.6),3,3,0.5),
		cbind("B","black",1,nFAC+2,0.8),
		cbind("C",rgb(0.4,0.4,0.4),1,3,0.5),
		cbind("A9","black",1,nFAC+2,0.8),
		cbind("A1","black",1,nFAC+2,0.8),
		cbind("A3",rgb(0.4,0.4,0.4),1,3,0.5),
		cbind("A6","black",1,nFAC+2,0.8),
		cbind("A2","black",1,nFAC+2,0.8),
		cbind("A4",rgb(0.4,0.4,0.4),1,3,0.5),
		cbind("A5a",rgb(0.4,0.4,0.4),1,3,0.5),
		cbind("A5b","black",1,nFAC+2,0.8),
		cbind("A7",rgb(0.4,0.4,0.4),1,3,0.5),
		cbind("A10","black",1,nFAC+2,0.8),
		cbind("A8",rgb(0.4,0.4,0.4),1,3,0.5))
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
		z1<-min(zx)+0.25
		z2<-max(zx)-0.25
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
				col= cx,lwd=2,type="l")
		text(XXm,YYm,x,font=fontx,col=cx,cex=cexx)
	})
	######Sampling year and strain name/clade

	LABn<-as.vector(sapply(namex,function(x){
		CHARx<-subset(CHAR,CHAR[,1]==x)
		return(CHARx[,4])
	}))
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
		nx<-ifelse(batchx[x]=="", LABn[x],paste(LABn[x]," (",batchx[x],")",sep=""))
		angx=2*pi*(x/z)*OPEN
		XX1i<-y1* maxTIP*sin(angx)
		YY1i<-y1* maxTIP*cos(angx)						
		text(XX1i, YY1i,nx,cex=0.2,srt= ang.lab[x],adj= adjx[x],col=cx,font=fontx[x])	
	})	
	#Display species name for reference genomes, when available
	SPE<-as.vector(type[,c(subset(c(1:length(type[1,])),colnames(type)=="species"))])
	SPE<-ifelse(CLA=="Microvirga","",ifelse(CLA=="Enterovirga","",
		ifelse(SPE =="sp.","",paste("M.",SPE))))
	y1=Ai[2]
	sapply(Tx,function(x){
		cx<-subset(COLcla[,2], COLcla[,1]== CLA[x])
		nx<-SPE[x]
		angx=2*pi*(x/z)*OPEN
		XX1i<-y1* maxTIP*sin(angx)
		YY1i<-y1* maxTIP*cos(angx)						
		text(XX1i, YY1i,nx,cex=0.2,srt= ang.lab[x],adj= adjx[x],col=cx,font=3)	
	})
#Display factors and permutation test results

COORD<-sapply(c(length(TEST):1),function(x){
	print(x)	
	FACx=TEST[x]
		X=x
		Aix<-Ai[X+1]
		Aex<-Ae[X+1]
		VARcx<-subset(t(VARc),VARc[2,]==FACx)
		varx<-unique(VARcx[,1])
		typex<-as.vector(subset(t(type),colnames(type)==varx))
		colx<-as.vector(VARcx[,5])
		colx<-as.vector(sapply(typex,function(x){
			return(c(subset(as.vector(VARcx[,5]),
				as.vector(VARcx[,4])==x),NA)[1])
		}))
		ATx=X+1
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
		VARx<-unique(as.vector(VARcx[,3]))
		labx<-as.vector(VARcx[,7])
		colx<-as.vector(VARcx[,5])
		lab.colx<-as.vector(VARcx[,6])
		Ns=length(labx)
		y1=Ae[ATx]
		y2=Ai[ATx]
		Nmin=min(sapply(c(1:Ns),function(x){
			x1=seq(l1*rLIM,l2*rLIM,
				length.out= Nlab)[x]
			x2=seq(l1*rLIM,l2*rLIM,
				length.out= Nlab)[x+1]	
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
#####Signficantly associated nodes == ONLY DISPLAY MONOPHYLETIC GROUPS (REMOVE POLYPHYLETIC GROUPS)
	SUMx<-subset(SUM,SUM[,3]== FACx)
	#Remove test that were done on polyphyletic groups
	POL<-sapply(c(1:length(SUMx[,1])),function(x){
		nodex<-as.numeric(unlist(strsplit(unlist(strsplit(as.vector(SUMx[x,4]),split="_"))[2],split="[.]"))[2])
		#Level at which it was tested
		levelx<-as.numeric(SUMx[x,1])
		#Strains that were tested
		STx<-subset(LEVELf[,1],
			as.vector(subset(t(LEVELf),
			colnames(LEVELf)==paste("L", levelx,sep="-")))== nodex)		
		#Real level 
		levelt<-subset(PSn[,1],PSn[,2]==nodex)
		#Real Strains
		STt<-subset(LEVELS[,1],
			as.vector(subset(t(LEVELS),
			colnames(LEVELS)==paste("L", levelt,sep="-")))== nodex)
		#IF NOT REAL STRAINS OF THE NODES WERE INCLUDED IN THE TEST, CONSIDER AS POLYPHYLETIC AND REJECT THE TEST
		return(ifelse(length(intersect(STx, 
			STt))<length(STt),
			1,0))
	})
	SUMx<-subset(SUMx,POL==0)
	#Variable and node
	ND<-matrix(unlist(strsplit(as.vector(SUMx[,4]),split="_")),ncol=2,byrow=T)
	vx<-ND[,1]
	ND<-as.numeric(matrix(unlist(
		strsplit(ND[,2],split="[.]")),
		ncol=2,byrow=T)[,2])
	NDu<-unique(ND)
	#P. value
	px<-as.numeric(SUMx[,7])
	#Levels
	Lx<-as.numeric(SUMx[,1])
	##ARRIVÉ ICITTE
	COORDx<-t(sapply(NDu[length(NDu):1],
		function(x){
		print(x)
		NDx=x	
		#Levels
		Lxx<-subset(Lx,ND== NDx)
		#Retrieve untested levels (same node composition)
		Lxx<-sort(as.numeric(unlist(
			strsplit(LevNDx[
			sapply(Lxx,function(x){
			grep(paste("_",x,"_",sep=""),
				paste("_", LevNDx,
				"_",sep=""))
		})],split="_"))))
		L0<-min(as.numeric(unlist(
			strsplit(LevNDx,split="_"))))
		Lxx<-c(Lxx,max(Lxx)+1)
		#Pvalues
		#pxx<-subset(px,ND== NDx)
		pxx<-max(subset(px,ND== NDx))
		axx<-ifelse(pxx <=0.001,0.25,ifelse(pxx <=0.01,0.1,0.05))
		sxx<-ifelse(pxx<=0.001,"***",
			ifelse(pxx<=0.01,"**","*"))
		#Variable
		#vxx<-subset(vx,ND== NDx)
		vxx<-unique(subset(vx,ND== NDx))
		colvx<-sapply(c(1:length(vxx)),
			function(x){
			colxx<-subset(VARcx[,5],
				VARcx[,4]== vxx[x])
			return(rgb(col2rgb(colxx)[1]/255,
				col2rgb(colxx)[2]/255,
				col2rgb(colxx)[3]/255, axx[x]))	
		})
		colbx<-sapply(c(1:length(vxx)),
			function(x){
			return(subset(VARcx[,5],
				VARcx[,4]== vxx[x]))	
		})
		colpx<-sapply(c(1:length(vxx)),
			function(x){
			colxx<-subset(VARcx[,5],
				VARcx[,4]== vxx[x])
			return(rgb(col2rgb(colxx)[1]/450,
				col2rgb(colxx)[2]/450,
				col2rgb(colxx)[3]/450))	
		})		
		Stx<-as.vector(t(subset(t(LEVELf),
			colnames(LEVELf)==
			paste("L", Lxx[1],sep="-"))))
		Stx<-subset(LEVELf[,1],Stx== NDx)
		Yn=seq(min(Stx),
			max(Stx),length.out=100)
		angn=2*pi*(Yn/z)*OPEN
		#Xb1=Lmax-Lxx[1]
		Xb1=Lmax-L0
		Xb2=Lmax-Lxx[length(Lxx)]
		XX1n<-Xb1*sin(angn)
		YY1n<-Xb1*cos(angn)
		XX2n<-Xb2*sin(angn)
		YY2n<-Xb2*cos(angn)
		polygon(c(XX1n,XX2n[100:1]),
			c(YY1n,YY2n[100:1]),
			border=colbx,
			col= ifelse(X==2,colvx,NA),
			lwd=1.2)
		Xbm=mean(c(Xb1, Xb2))
		XXm<-Xbm*sin(mean(angn))
		YYm<-Xbm*cos(mean(angn))
		#text(XXm, YYm, sxx,font=2,cex=1.2,col=colpx)
		return(as.vector(c(FACx, NDx,
			XXm, YYm, sxx,vxx,colpx,pxx)))
	}))
	return(t(COORDx))
})		

	COORD<-matrix(unlist(COORD),ncol=8,byrow=T)
	
	colnames(COORD)<-c("fac","node","X","Y","sig","var","color","pval")
	
	#Display p-values (add manually)
	#text(as.numeric(COORD[,3]),as.numeric(COORD[,4])+ifelse(COORD[,1]=="temp",-1.5*xs1*maxTIP,+1.5*xs1*maxTIP),paste(ifelse(COORD[,6]=="20","20°C",ifelse(COORD[,6]=="30","30°C",COORD[,6])), COORD[,5],sep=""),cex=0.5,font=2,col=COORD[,7])


###### Add PERMANOVA on the top left corner

#Levels for x-scale
LscaleAOV<-round(seq(min(LV),max(LV),length.out=6))

#Corresponding PS values to show on the scale
PSscaleAOV<-round(sapply(LscaleAOV,function(x){
	return(subset(SCALE[,4],SCALE[,2]==x))
}),2)

#Color code for factors
COLfac<-c("blue","grey","red","yellow2",
	"cyan2","purple","red3","green","yellow4","orange",
	"purple3","green3","brown","orange3",
	"black","grey",NA)


####Set graphical limits
#minimal Part of variance
AOVx1=(xs2-0.01)*rLIM
#maximal Part of variance
AOVx2=-maxTIP*Ae[2]
#minimal PS
AOVy1=min(Xscale)
#minimal PS
AOVy2=max(Xscale)

#polygon(c(AOVx1, AOVx1, AOVx2, AOVx2),
	#c(AOVy1, AOVy2, AOVy2, AOVy1),col="white",border="grey")

#Maximal part of variance
MaxVAR=round(1.1*max(AOV[c(3:17),]),2)

DESC<-setdiff(unique(unlist(strsplit(rownames(AOV[c(3:length(AOV[,1])),]),split=" "))),c("R2","Pr(>F)"))

#Convert LV in graphique coordinates
#LVreal<-c(min(Lscale):max(Lscale))
#LVgraph<-seq(max(Xscale),min(Xscale),length.out=(max(Lscale)-min(Lscale)+1))

#LVgraph <-sapply(LV,function(x){
#	return(subset(LVgraph, LVreal==x))
#})

LVgraph <-(Lmax-LV)

##Plot Variance axis
		
	VARscale<-sapply(c(1:(Nscale-1)),function(x){
		X=seq(AOVx2, AOVx1,length.out=(Nscale-1))[x]
		N=seq(MaxVAR,0,length.out=(Nscale-1))[x]
		segments(X, min(LVgraph)-(0.015*rLIM), X, min(LVgraph)-2*(0.015*rLIM))
		text(X, min(LVgraph)-(0.015*rLIM), round(N,2),cex=0.5,pos=1)
		return(X)
	})
	
	segments(min(VARscale), min(LVgraph)-(0.015*rLIM),max(VARscale), min(LVgraph)-(0.015*rLIM))
	
	text(mean(VARscale), (min(LVgraph)-2*(0.015*rLIM))/2,"Part of variance",cex=0.8,font=1)



#Plot variance explained by each factor in function of PS
PMx<-sapply(c(1:(length(DESC)-1))+2,function(x){
	#color
	colx<-COLfac[x-2]
	#factor
	DESCx<-rownames(AOV)[x]
	#R2 (part of explained variance)
	R2x<-as.numeric(AOV[x,])
	#Convert R2 in graphic coordinates
	
	R2graph<-R2x/MaxVAR
	R2graph<-AOVx1+(AOVx2-AOVx1)* R2graph
	
	#pvalue
	px=as.vector(AOV[x+length(DESC),])
	pm=min(px)
	colm<-ifelse(DESCx=="R2 Residuals","grey",ifelse(pm>0.05,NA, colx))
	
	points(R2graph, LVgraph ,type="l",col= colm,lty=ifelse(DESCx=="R2 Residuals",3,2),lwd=0.5)
	
	#colma =ifelse(is.na(colm)==T,NA,rgb(col2rgb(colm)[1]/255,
	#		col2rgb(colm)[2]/255,
	#		col2rgb(colm)[3]/255,0.3))
	
	points(R2graph, LVgraph,pch=19,
		cex=ifelse(px<=0.001,0.5,
		ifelse(px<=0.01,0.3,
		ifelse(px<=0.05,0.15,
		ifelse(px<=0.1,0,0)))),bg= ifelse(is.na(colm)==T,NA,"white"),col= colm)	
	#points(LV,R2x,pch=21,
	#	cex=ifelse(px<=0.001,1,
	#	ifelse(px<=0.01,0.75,
	#	ifelse(px<=0.05,0.5,
	#	ifelse(px<=0.1,0.25,0)))),bg= colma,col= colm)
	
	return(pm)		
})

	legend(x = -Ae[3]*maxTIP,y=0.65*Lmax,cex=0.5, horiz=F,legend= c(
	expression("p"<= "0.001"),
	expression("p"<= "0.01"),
	expression("p"<= "0.05")),
	border=T,box.col=NA,pch=19,pt.cex=c(0.5,0.3,0.15))

dev.off()

#Factor that are significant in permanoval for at least one level with color code for the Venn diagram to add manually in Figure 3a (Figure3.pptx; see also file Permanova-tree-min-bootstrap=30.xlsx for details)
subset(cbind(DESC, COLfac),PMx<=0.05)

#Node and corresponding clades that are significantly associated with at least one factor, to add manually in figure 3b (Figure3.pptx)

COORD #node positions in the tree can be retrieved in file Tree-NodeLabels_boot=30.pdf 

