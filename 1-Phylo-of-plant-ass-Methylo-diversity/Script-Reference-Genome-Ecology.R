#load libraries
library(ape)
library(seqinr)

#Directory for data
pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/1-Phylo-of-plant-ass-Methylo-diversity/data"

setwd(pD)

#rpoB sequences
SEQ<-read.fasta("rpoB-update-sept-2020.fas")

#sequence names
nsx<-names(SEQ)

#Consensus rpoB phylogenetic tree (newick format)
tx<-read.tree("consensus-rpob-Jack.nwk")

#convert node support
tx $node.label<-round(100*as.numeric(tx $node.label))

#retrieve tip labels in the tree
ntx<-tx$tip.label

plot(tx,"fan",cex=0.2)

#retrieve genome information matching tip labels
orderP<-read.table("order_phylo-update.txt",header=T)

orderP<-orderP[order(orderP[,1]),]

#Only keep Methylobacterium tips

tlx<-keep.tip(tx,subset(orderP[,1],orderP[,5]=='Methylobacterium'))

orderP <-subset(orderP,orderP[,5]=='Methylobacterium')

#add branch length (Grafen method; https://pubmed.ncbi.nlm.nih.gov/2575770/)
tlx<-compute.brlen(tlx,power=0.5)
plot(tlx,"fan",cex=0.2)


######FIGURE 1
 
pdf("Figure-1-TreeMethylo-ref-rpoB.pdf",width=7,height=7)

#Frame for the tree (adjust marx in function of the size of the graphic)

marx=8.5

par(bty="n",mar=c(marx, marx, marx, marx))

plot.phylo(tlx,"fan",cex=0.2,show.tip.label = F,rotate.tree=90,show.node.label=T,adj=0.5,font=2, edge.color ="grey")

#########Graphical parameters for tips information (added separately)
#Openness of the circle (value between 0 and 1; 1 = full circle; 0 = fully compressed)
OPEN=1
#Margin size
MAR=2
#tip attribute coordinates (concentric circles)
Ai=1.01+seq(0.05,1,length.out=24)
Ae=Ai+0.8*(Ai[2]-Ai[1])

#Limit of the graphic
rLIM=MAR

#Start of the plot
par(mar=c(1,1,1,1),bty="n",mfrow=c(1,1),new=T)

#Frame
plot(-10000,-10000,xlim=c(-rLIM, rLIM),ylim=c(-rLIM, rLIM),xlab="",ylab="",cex.axis=0.8,las=1,xaxt="n",yaxt="n",new=T)	

maxTIP=1

	#######environment of sampling
	LABx=-11
	ATx=-12
	LEGx=6
		

#Genome order
z<-orderP[,1]

z<-z[length(z):1]

#Display biome of origin for each strain
	colECO<-cbind(c("PLANT","PHYLLOSPHERE","ENDOPHYTE","SEED","RHIZHOSPHERE","SOIL","LICHEN-FUNGI","WATER","OCEAN","MICROBIOME"),
	c("green3","green","yellow3","yellow2","brown","red","orange","cyan","blue","purple"),
	c("Plant (unknown part)","Phyllosphere","Endophyte","Seed","Rhizosphere","Soil","Fungi or lichen","Water","Ocean","Gut"))	

#Display clades and species

CLu<-as.vector(sort(unique(orderP[,8])))

#CLu<-setdiff(as.vector(sort(unique(orderP[,8]))),c("Enterovirga","Microvirga"))

CLx<-ifelse(CLu=="Microvirga",10,
	ifelse(CLu=="Enterovirga", 10, 15))

SPu<-setdiff(as.vector(sort(unique(subset(orderP[,6], orderP[,5]=="Methylobacterium")))),
	c("sp.",
	"thiocyanatum",#same than M. populi
	"zatmanii")#same than M. extorquens
	) 

	ang.lab<-seq(360,(1-OPEN)*360,length.out=length(z))-270
	adjx<-ifelse(ang.lab>-90,1, 0)
	adjy<-ifelse(ang.lab>-90,0, 1)
	ang.lab<-ifelse(ang.lab>-90,ang.lab, ang.lab+180)

sapply(SPu,function(x){	
	clx<-unique(subset(orderP[,8],orderP[,6]==x))
	ATxx <-subset(CLx, CLu ==clx)
	zx<-sort(subset(orderP[,1], orderP[,6]==x))	
	atu<-ifelse(x=="platani",2,
		ifelse(x=="tarhaniae",2,3))		
	z1<-subset(c(1:length(z)),z==min(zx))
	z2<-subset(c(1:length(z)),z==max(zx))
	
	zx<-seq(z1,z2,length.out=1000)
	zm=mean(zx)	
	yz<-((Ai+Ae)[ATxx-atu])/2
	ym<-Ae[ATxx-atu]	
	angz=2*pi*(zx/length(z))*OPEN
	angm=2*pi*(zm/length(z))*OPEN			
	XX1i<-yz* maxTIP*sin(angz)
	YY1i<-yz* maxTIP*cos(angz)
	XXm<-ym* maxTIP*sin(angm)
	YYm<-ym* maxTIP*cos(angm)			
	points(XX1i,YY1i,
		lwd=1.5,type="l")		
	text(XXm,YYm,paste("M.",x),font=4,cex=0.3,srt=ang.lab[round(mean(c(z1:z2)))],adj= adjy[round(mean(c(z1:z2)))])
})



sapply(CLu,function(x){	
	ATxx <-subset(CLx, CLu ==x)
	zx<-sort(subset(orderP[,1], orderP[,8]==x))	
	#remove isolated tips
	zx<-subset(zx,c(sapply(c(1:(length(zx)-1)),
		function(x){
			return(zx[x+1]-zx[x])
	}),1)==1)		
	z1<-subset(c(1:length(z)),z==min(zx))
	z2<-subset(c(1:length(z)),z==max(zx))
	zx<-seq(z1,z2,length.out=1000)
	zm=mean(zx)	
	yz<-((Ai+Ae)[ATxx+6])/2
	ym<-Ae[ATxx +8]	
	angz=2*pi*(zx/length(z))*OPEN
	angm=2*pi*(zm/length(z))*OPEN			
	XX1i<-yz* maxTIP*sin(angz)
	YY1i<-yz* maxTIP*cos(angz)
	XXm<-ym* maxTIP*sin(angm)
	YYm<-ym* maxTIP*cos(angm)			
	points(XX1i,YY1i,
		lwd=2,type="l")		
	text(XXm,YYm,x,font=2,cex=0.8)	
})

#Strain label

	
sapply(c(1:length(z)),function(x){
	nx=as.vector(subset(orderP[,3],orderP[,1]==z[x]))
	cx<-as.vector(subset(orderP[,8],orderP[,1]==z[x]))
	ATxx<-subset(CLx,CLu==cx)
	typex<-unlist(strsplit(as.vector(
		subset(orderP[,7],orderP[,1]==z[x])),
		split="/"))
	#print(c(x,typex))	
	antx<-length(subset(typex, typex=="ANTHROPO"))
	typex<-setdiff(typex,"ANTHROPO")
	colECOx<-subset(colECO[,2], colECO[,1]== typex)
	ang1=2*pi*((x-0.4)/length(z))*OPEN
	ang2=2*pi*((x+0.4)/length(z))*OPEN
	XX1i<-Ai[ATx+ ATxx]* maxTIP*sin(ang1)
	YY1i<-Ai[ATx+ ATxx]* maxTIP*cos(ang1)
	XX2i<-Ai[ATx+ ATxx]* maxTIP*sin(ang2)
	YY2i<-Ai[ATx+ ATxx]* maxTIP*cos(ang2)
	XX1e<-Ae[ATx+ ATxx]* maxTIP*sin(ang1)
	YY1e<-Ae[ATx+ ATxx]* maxTIP*cos(ang1)
	XX2e<-Ae[ATx+ ATxx]* maxTIP*sin(ang2)
	YY2e<-Ae[ATx+ ATxx]* maxTIP*cos(ang2)		
	polygon(c(XX1i,XX2i,XX2e,XX1e),
		c(YY1i,YY2i,YY2e,YY1e),
		col= colECOx,
		border=ifelse(antx==1,"black",NA),
		lwd=1)
	#construction lines (optional)	
	#segments(0,0,(XX1i+ XX2i)/2,(YY1i+ YY2i)/2,col="grey",lwd=0.5)	
	#Names	
	angx=2*pi*(x/length(z))*OPEN
	XX1i<-Ae[LABx+ ATxx]* maxTIP*sin(angx)
	YY1i<-Ae[LABx+ ATxx]* maxTIP*cos(angx)						
	text(XX1i, YY1i,nx,cex=0.3,srt= ang.lab[x],adj= ifelse(adjx[x]==1,0,1),col=ifelse(length(colECOx) ==0,"grey", colECOx),font=2)	
})

legend("topright",cex=0.5, horiz=F,
	legend= c(colECO[,3],"Anthopogenic","unknown"),
	text.col= c(colECO[,2],"black","grey"),box.col=NA,pch=22,pt.bg=c(colECO[,2],"white","white"),
	border=NA,col= c(colECO[,2],"black","white"),title="Environmental sources",title.col="black",ncol=1)

	dev.off()






######################
######################
######################
######################
###Pairwise similarity among clades and within clades (calculated in MEGA from rpoB-update-sept-2020.meg)


DIST<-read.table("Distance-matrix.txt",header=T)

MINx=1-max(as.matrix(DIST[,c(4:length(DIST[1,]))]),na.rm=T)
MAXx=1-min(as.matrix(DIST[,c(4:length(DIST[1,]))]),na.rm=T)

CLAu<-summary(as.factor((DIST[,3])))
CLAu<-subset(names(CLAu), as.numeric(CLAu)>2)
CLAu<-setdiff(CLAu,"Microvirga")

#density (violin plot) parameters
DENS=150
SVP=30


pdf("rpoB-PS-clades.pdf",height=5,width=6)


par(mar=c(4,1,1,1),bty="n")

plot(-10,-10,xlab="PS",ylim=c(length(CLAu)+6,0),xlim=c(MINx-0.11,MAXx),yaxt="n",ylab="",las=1,cex.axis=0.8)

#PS within Methylobacterium clades
sapply(c(1:length(CLAu)),function(x){
	X=x
	CLAx=CLAu[X]
	ix=as.vector(subset(DIST[,1],DIST[,3]== CLAx))
	DISTx<-as.vector(1-as.matrix(DIST[ix, ix+3]))
	DISTx<-subset(DISTx,is.na(DISTx)==F)
	z<-density(DISTx,n= DENS,na.rm=T)
	Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
	Xvp <-Xvp/(3*max(Xvp))
	Yvp<-c(z$x, z$x[DENS:1])
	polygon(Yvp,Xvp+X,col=rgb(1,0,0,0.3),border=NA)
	segments(median(DISTx),min(Xvp)+X,median(DISTx),max(Xvp)+X,lwd=1.5,col="red")
	text(MINx,X,  CLAx,cex=0.6,pos=2,col="red")
})

text(MINx-0.03,mean(c(1:length(CLAu))),  "Within clades",cex=0.6,col="red",srt=90,font=2)

segments(MINx-0.02,1,MINx-0.02,length(CLAu),lwd=1.5,col="red")



#PS among clades from clade A

claA<-as.vector(unique(DIST[,3]))
claA<-setdiff(claA,c("B","C","Microvirga","Enterovirga"))

PSa<-as.vector(unlist(sapply(c(1:(length(claA)-1)),function(x){
	X=x
	CLAx= claA[X]
	ix=as.vector(subset(DIST[,1],DIST[,3]== CLAx))
	return(sapply(c((X+1):(length(claA))),function(x){
		CLAy= claA[x]
		iy=as.vector(subset(DIST[,1],DIST[,3]== CLAy))	
		DISTx<-1-as.matrix(DIST[ix, iy+3])
		DISTy<-1-as.matrix(DIST[iy, ix+3])	
		return(c(as.vector(DISTx),as.vector(DISTy)))
	}))
})))
PSa <-subset(PSa,is.na(PSa)==F)

X=length(CLAu)+1
denx= PSa
NX="among A1-A9"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(1,0,1,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(1,0,1,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(1,0,1,1))


#PS between A and B
PSab<-as.vector(unlist(sapply(c(1:(length(claA))),function(x){
	X=x
	CLAx= claA[X]
	ix=as.vector(subset(DIST[,1],DIST[,3]== CLAx))
	iy=as.vector(subset(DIST[,1],DIST[,3]== "B"))	
	DISTx<-1-as.matrix(DIST[ix, iy+3])
	DISTy<-1-as.matrix(DIST[iy, ix+3])	
	return(c(as.vector(DISTx),as.vector(DISTy)))
})))
PSab <-subset(PSab,is.na(PSab)==F)

X=length(CLAu)+2
denx= PSab
NX="between A and B"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(1,0,1,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(1,0,1,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(1,0,1,1))


#PS between A and C
PSac<-as.vector(unlist(sapply(c(1:(length(claA))),function(x){
	X=x
	CLAx= claA[X]
	ix=as.vector(subset(DIST[,1],DIST[,3]== CLAx))
	iy=as.vector(subset(DIST[,1],DIST[,3]== "C"))	
	DISTx<-1-as.matrix(DIST[ix, iy+3])
	DISTy<-1-as.matrix(DIST[iy, ix+3])	
	return(c(as.vector(DISTx),as.vector(DISTy)))
})))
PSac <-subset(PSac,is.na(PSac)==F)

X=length(CLAu)+3
denx= PSac
NX="between A and C"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(1,0,1,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(1,0,1,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(1,0,1,1))

#PS between B and C

	ix=as.vector(subset(DIST[,1],DIST[,3]== "B"))
	iy=as.vector(subset(DIST[,1],DIST[,3]== "C"))	
	DISTx<-1-as.matrix(DIST[ix, iy+3])
	DISTy<-1-as.matrix(DIST[iy, ix+3])	
	PSbc<-c(as.vector(DISTx),as.vector(DISTy))
	PSbc <-subset(PSbc,is.na(PSbc)==F)

X=length(CLAu)+4
denx= PSbc
NX="between B and C"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(1,0,1,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(1,0,1,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(1,0,1,1))

text(MINx-0.08,mean(c(1:(4+length(CLAu)))),  "Within Methylobacterium",cex=0.6,col=rgb(1,0,1,1),srt=90,font=2)

segments(MINx-0.07,1,MINx-0.07,length(CLAu)+4,lwd=1.5,col=rgb(1,0,1,1))



#Within Microvirga

ix=as.vector(subset(DIST[,1],DIST[,3]== "Microvirga"))
DISTx<-as.vector(1-as.matrix(DIST[ix, ix+3]))
DISTx<-subset(DISTx,is.na(DISTx)==F)

X=length(CLAu)+5
denx= DISTx
NX="within Microvirga"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(0.5,0.5,0.5,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(0.5,0.5,0.5,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(0.5,0.5,0.5,1),font=2)

# Between Methylobacterium and Microvirga

claM<-as.vector(unique(DIST[,3]))
claM <-setdiff(claM,c("Microvirga","Enterovirga"))

PSm<-as.vector(unlist(sapply(c(1:(length(claM))),function(x){
	X=x
	CLAx= claM[X]
	ix=as.vector(subset(DIST[,1],DIST[,3]== CLAx))
	iy=as.vector(subset(DIST[,1],DIST[,3]== "Microvirga"))	
	DISTx<-1-as.matrix(DIST[ix, iy+3])
	DISTy<-1-as.matrix(DIST[iy, ix+3])	
	return(c(as.vector(DISTx),as.vector(DISTy)))
})))
PSm <-subset(PSm,is.na(PSm)==F)

X=length(CLAu)+6
denx= PSm
NX="between Methylobacterium and Microvirga"

z<-density(denx,n= DENS,na.rm=T)
Xvp<-c(-(z$y)/SVP,(z$y[DENS:1])/SVP)
Xvp <-Xvp/(3*max(Xvp))
Yvp<-c(z$x, z$x[DENS:1])
polygon(Yvp,Xvp+X,col=rgb(0,0,0,0.3),border=NA)
segments(median(denx),min(Xvp)+X,median(denx),max(Xvp)+X,lwd=1.5,col=rgb(0,0,0,1))
text(MINx,X,  NX,cex=0.6,pos=2,col=rgb(0,0,0,1),font=2)

dev.off()