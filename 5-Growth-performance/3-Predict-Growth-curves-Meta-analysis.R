#############################
#### TEMPERATURE SCREEN
#############################

pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/5-Growth-performance/data"

library(gplots)
library(ade4)
library(vegan)

###Retrieve strain characteristics
pO<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/3-culturable-Methylo-diversity-timeline/data"

setwd(pO)

GROUP<-read.table("strain-metadata.txt",header=T)

#clade color
COLcla<-rbind(
	cbind("Enterovirga",rgb(0.8,0.8,0.8)),
	cbind("Microvirga",rgb(0.6,0.6,0.6)),
	cbind("B","black"),
	cbind("C",rgb(0.4,0.4,0.4)),
	cbind("A9","red"),
	cbind("A1","blue"),
	cbind("A3","blue4"),
	cbind("A6","orange"),
	cbind("A2","cyan2"),
	cbind("A4","brown"),
	cbind("A5","yellow3"),
	cbind("A7","pink"),
	cbind("A10","green3"),
	cbind("A8","purple"),
	cbind("NA",NA),
	cbind("Neg",NA))

###Retrieve metadata for temperature screen

setwd(pD)

#Petri plan
PLAN<-read.table("plan_petris.txt",header=T)
#Strain list
STRAIN<-read.table("strain.txt",header=T)
#Filter to remove problematic plates
FILT<-read.table("filter.txt",header=T)

#Clade names
NG<-sapply(as.vector(STRAIN[,1]),function(x){
	Sx<-paste(unlist(
		strsplit(x,split="-"))[c(1:2)],
		collapse="-")
	Gx<-subset(as.vector(GROUP[,2]),
		as.vector(GROUP[,1])== Sx)
	return(ifelse(length(Gx)==0,"NA",Gx))
})

STRAIN<-cbind(STRAIN[,c(1:2)],
	as.vector(NG),
	STRAIN[,c(4:length(STRAIN[1,]))])
	
###### For each time point, each strain and each condition, return median intensities corrected by background (one file per time point)

#Plate names
NAME<-c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24") 

#Tested temperatures
TP<-c("T20","T30")

#Date of pictures (3 tim points; inoculation on 15-02-2019)
TM<-c("22-02-2019","28-02-2019","11-03-2019") 

#Time (days) between inoculation and pictures
ND<-c(7,13,24)

#Temperature treatments
TT<-c("P20M20","P20M30","P30M20","P30M30")

#Pool raw data from 3 time points and return name of strains 
ST<-unique(as.vector(sapply(TM,function(x){
	TMx<-x
	NDx<-subset(ND,TM==TMx)
	
	#RAW DATA FOR THIS TIME POINT 
	#(retrieved for each plate)		
	INTall<-matrix(sapply(NAME,function(x){
		NAMEx<-x
		
		#RAW DATA FOR MONITORING TEMPERATURE = 20°C
		TX="T20" 
		P1<-paste(pD,"/Pic_",TMx,
			"/",TX,"/",NAMEx,sep="")
		setwd(P1)
		INT<-read.table("spot_intensities.txt",
			header=T)	
		FILTx<-t(subset(t(FILT),
			colnames(FILT)== 
			paste(TX,NAMEx,sep="_")))[,1]
		PLx<-t(subset(t(PLAN),
			colnames(PLAN)== NAMEx))[,1]
			
		mINT2020<-cbind(PLx,20,20,NDx,
			NAMEx ,INT[PLAN[,1],], 
			FILTx[PLAN[,1]])
		colnames(mINT2020)<-
			c("strain","T1","T2",
			"day","Petri",
			colnames(INT),"filter")
		mINT3020<-cbind(PLx,30,20,NDx, 
			NAMEx ,INT[PLAN[,2],], 
			FILTx[PLAN[,2]])
		colnames(mINT3020)<-
			c("strain","T1","T2",
			"day","Petri",
			colnames(INT),"filter")
		
		#RAW DATA FOR MONITORING TEMPERATURE = 30°C
		TX="T30" 
		P1<-paste(pD,"/Pic_",TMx,
			"/",TX,"/",NAMEx,sep="")
		setwd(P1)
		INT<-read.table("spot_intensities.txt",
			header=T)	
		FILTx<-t(subset(t(FILT),
			colnames(FILT)== 
			paste(TX,NAMEx,sep="_")))[,1]
		PLx<-t(subset(t(PLAN),
			colnames(PLAN)== NAMEx))[,1]
			
		mINT2030<-cbind(PLx,20,30,NDx,
			NAMEx ,INT[PLAN[,1],], 
			FILTx[PLAN[,1]])
		colnames(mINT2030)<-
			c("strain","T1","T2",
			"day","Petri",
			colnames(INT),"filter")
		mINT3030<-cbind(PLx,30,30,NDx, 
			NAMEx ,INT[PLAN[,2],], 
			FILTx[PLAN[,2]])
		colnames(mINT3030)<-
			c("strain","T1","T2",
			"day","Petri",
			colnames(INT),"filter")		

		return(t(rbind(mINT2020, 
			mINT3020, mINT2030, mINT3030)))
	}),ncol=14,byrow=T)
		
	colnames(INTall)<-c("strain","T1","T2",
		"day","Petri","spot",
		"X_position","Y_position","radius",
		"av_intensity","sd_intensity",
		"av_intensity_back","sd_intensity_back",
		"filter")

	#Strains
	STx<-INTall[,1]
		
	#Retrieve strain clade
	CLx<-as.vector(sapply(STx,function(x){
		return(c(subset(as.vector(STRAIN[,3]),
			as.vector(STRAIN[,1])==x),"Neg")[1])
	}))
		
	print(length(CLx))
	setwd(pD)
	write.table(cbind(CLx, INTall),
		paste("Detail-Growth-DAY=",NDx,
		".txt",sep=""))
	return(STx)
		
})))	


setwd(pD)
	
#Corrected intensity per strain, per replicate, per treatment, per date
STTx<-sapply(ND,function(x){
	NDx=x
	STTm<-read.table(paste(
		"Detail-Growth-DAY=",NDx,
		".txt",sep=""),header=T)
	return(as.numeric(STTm[,11])-
		as.numeric(STTm[,13]))
})

#Filter per strain, per replicate, per treatment (determine if replicates should be removed; e.g.: because spot touch each over in the plate)
Fx<-sapply(ND,function(x){
	NDx=x
	STTm<-read.table(paste(
		"Detail-Growth-DAY=",NDx,
		".txt",sep=""),header=T)
	return(as.numeric(STTm[,15]))
})

#Replicat characteristics: Clade, strain name, treatment, day, petri and spot number in the petri
CHARx<-read.table(paste("Detail-Growth-DAY=",
	ND[1],".txt",sep=""),
	header=T)[,c(1,2,3,4,6,7,8,9)]
	
#Negative controls mean and sd intensity
mean(rbind(subset(STTx,CHARx[,2]=="CN2"),
	subset(STTx,CHARx[,2]=="CN1")))
sd(rbind(subset(STTx,CHARx[,2]=="CN2"),
	subset(STTx,CHARx[,2]=="CN1")))
	
#Replicate mean and sd intensity (negative controls excluded)
mean(subset(STTx,CHARx[,1]!="Neg"))
sd(subset(STTx,CHARx[,1]!="Neg"))


#####Correct intensity by position in petri dish: predict intensity according to position on petri dish (36 position per petri) and make correction accordingly (because of more competition in the center, less growth)
	
CENTER<-c(15,16,21,22) #Spots located in the middle
MED<-c(8:11,14,17,20,23,26:29) #Spots located at intermediate distance between center and border
BORDER<-setdiff(c(1:36),c(MED,CENTER)) #Spots located at the border

#Values in the center	
vC<-(unlist(sapply(CENTER,function(x){
	return(subset(STTx,CHARx[,6]==x))
})))
#Values in the intermediate zone	
vM<-(unlist(sapply(MED,function(x){
	return(subset(STTx,CHARx[,6]==x))
})))
#Values at the border
vB<-(unlist(sapply(BORDER,function(x){
	return(subset(STTx,CHARx[,6]==x))
})))
	
#Mean and sd value in each zone
vC<-c(mean(vC),sd(vC))	
vM <-c(mean(vM),sd(vM))	
vB <-c(mean(vB),sd(vB))		

#x and y position of each spot in each petri	
xx<-as.vector(cbind(as.numeric(CHARx[,7]),
	as.numeric(CHARx[,7]),
	as.numeric(CHARx[,7])))
yy<-as.vector(cbind(as.numeric(CHARx[,8]),
	as.numeric(CHARx[,8]),
	as.numeric(CHARx[,8])))
#Intensity per spot	
zz<-as.vector(STTx)

cx=ifelse(CHARx[,2]=="0",0,   #remove empty positions
	ifelse(apply(Fx,1,mean)!=1,0, #remove spots that touch each other
	ifelse(CHARx[,2]=="J-028",0,
	ifelse(CHARx[,2]=="E-024",0,1)))) #remove contaminated strains
cx<-c(cx,cx,cx)

#####Display spot intensities before and after correction by position
pdf("prediction.pdf")
	
	#All spots before correction (size proportional to growth)
	plot(xx,yy,cex= cx* zz/5,pch=19,
		col=rgb(0,0,0,0.1),lwd=0,las=2,
		cex.axis=0.8,main="Observed I",
		xlab="X",ylab="Y")
	
	#Model predicting growth in function of spot position
	model<-lm(zz ~
		I((xx^2)*(yy^2))
		+I((xx^2)*(yy))
		+I((xx)*(yy^2))
		+I((xx)*(yy)))
	COEF<-model $coefficients
	xp<-c(min(xx):max(xx))
	yp<-c(min(yy):max(yy))
	
	#Predicted values from the model for each posible position in a petri and display
	PRED<-sapply(xp,function(x){
		return(
			COEF[2]*x^2*yp^2
			+COEF[3]*x^2*yp
			+COEF[4]*x*yp^2
			+COEF[5]*x*yp
			+COEF[1])
	})
	rownames(PRED)<-yp
	colnames(PRED)<-xp
	heatmap(PRED,Rowv=NA,Colv=NA,scale="none",
		main="Predicted I",xlab="X",ylab="Y")
	
	#Predicted values for each spot
	PRED<-as.vector(
		sapply(c(1:length(zz)),function(x){
		xo= xx[x]
		yo=yy[x]
		return(
			COEF[2]* xo ^2* yo ^2
			+COEF[3]* xo ^2* yo
			+COEF[4]* xo* yo ^2
			+COEF[5]* xo* yo
			+COEF[1])
	}))
	#Correct growth intensity by predicted values
	STTn<-lm(zz ~ PRED)$residuals
	STTn<-STTn-min(STTn)
	
	#All spots after correction (size proportional to growth)
	plot(xx,yy,cex=cx* STTn/5,
		pch=19,col=rgb(0,0,0,0.1),
		lwd=0,las=2,cex.axis=0.8,
		main="Corrected I",xlab="X",ylab="Y")

dev.off()	

#Format corrected values
STTn<-matrix(STTn,ncol=3)
	
#remove positions were spots touch each other
CHARx<-subset(CHARx,apply(Fx,1,mean)==1)
STTx <-subset(STTx,apply(Fx,1,mean)==1)
STTn<-subset(STTn,apply(Fx,1,mean)==1)
	
#Remove empty positions
STTx <-subset(STTx,CHARx[,2]!="0")
STTn<-subset(STTn,CHARx[,2]!="0")	
CHARx<-subset(CHARx,CHARx[,2]!="0")

#Remove strains with very weird outcome
STTx <-subset(STTx,CHARx[,2]!="E-024")
STTn<-subset(STTn,CHARx[,2]!="E-024")	
CHARx<-subset(CHARx,CHARx[,2]!="E-024")	
STTx <-subset(STTx,CHARx[,2]!="J-028")
STTn<-subset(STTn,CHARx[,2]!="J-028")	
CHARx<-subset(CHARx,CHARx[,2]!="J-028")

#Correct by negative control intensities	
NEGm<-mean(subset(STTn,CHARx[,1]=='Neg'))
STTn<-STTn-NEGm	
STTn<-ifelse(STTn<=0,0, STTn)

#Keep only spot with growth in at least one time point after corrections	
CHARx<-subset(CHARx,apply(STTn,1,max)>0)	
STTx <-subset(STTx,apply(STTn,1,max)>0)
STTn <-subset(STTn ,apply(STTn,1,max)>0)	
	
#Average and sd per time point
apply(STTn,2,mean)
apply(STTn,2,sd)
	
#For each replicate, predict spot intensity over time using logistic growth, return predicted max intensities and Lag time

NDo<-c(0,ND) #ARTIFICIALY ADD TIME 0 (GROWTH = 0)

fdx=1.5 #Factor to extend the prediction after day 24 (tested in the range 1-2; the chosen value is the one with less outliers and best correlation between estimated and real yields)
TDo<-seq(min(NDo),max(NDo),length.out=10001)
TD1<-seq(min(NDo),max(NDo)* fdx,length.out=40001)

iTD<-sapply(ND,function(x){
	return(subset(c(1:length(TD1)),
		abs(x-TD1)==min(abs(x-TD1)))[1])
})

#Assume that maximum intensity cannot be reached within Mmin days
Nmin=6


Io=0.1 #Replace 0 by Io
zz=50

#range Number of days after inoculations assuming growth is neglectible 
NDm<-(7*seq(0.1,1,length.out=zz)^3)[1:zz]

RateYieldNew<-t(sapply(c(1:length(CHARx[,1])),
	function(x){

setwd(paste(pD,"temp",sep="/"))

jpeg(paste("Pred-",x,".jpeg",sep=""),width=500,height=500,quality=100)

par(mfrow=c(2,2),bty="n",mar=c(4,4,1,1),cex.axis=0.8,las=1)



#Real values for this replicate and this treatment (ad 0 at time 0)
Ix= as.numeric(STTn[x,])
Ix<-ifelse(Ix<Io,Io,Ix)
Ix= c(0, Ix)


plot(NDo,c(Ix),pch=c(1,19,19,19),col="black",
	ylim=c(0,max(Ix)*1.5),xlab="time (day)",
	xlim=c(0,max(TD1)),
	ylab="intensity",main=x)
	

legend(x = "topright",cex=0.8, horiz=F,
	legend= c("observed","assumed"),
	text.col= "black",box.col=NA,pch=c(19,1),
	border=NA,col= "black",title="real values")	
	
#1) Polynomial model fitting these values
model<-lm(Ix ~I(NDo[1:4]^2)+I(NDo[1:4]))
COEF<-model $coefficients
#Predicted values from the polynomial model
y=  COEF[2]*(TDo ^2) + COEF[3]* TDo + COEF[1]
df<-data.frame(x= TDo,y=y)

plot(NDo,c(Ix),pch=c(1,19,19,19),col="black",
	ylim=c(0,max(Ix)*1.5),xlab="time (day)",
	xlim=c(0,max(TD1)),
	ylab="intensity")
points(df,type="l",col="purple",lwd=1.5,lty=3)

legend(x = "topright",cex=0.8, horiz=F,
	legend= c("range 0-24"),
	text.col= "black",box.col=NA,lwd=1.5,lty=3,
	border=NA,col= "purple",title="quadratic fit on real")

#2) Fit a log normal model on the polynomial residual (will be use only to adjust predicted values <NDm days)
opt <- optim(c(1, 1, 1), function(p) 
	sum((dlnorm(df$x, p[1], p[2]) * p[3] - df$y)^2))
zY<-opt$par[3] * dlnorm(TDo, opt$par[1], opt$par[2])

plot(-10,-10,
	ylim=c(0,max(Ix)*1.5),xlab="time (day)",
	xlim=c(0,max(TD1)),
	ylab="intensity")
points(df,type="l",col="purple",lwd=1.5,lty=3)
	
points(TDo, zY,type="l",col="blue",lwd=1.5)

legend(x = "topright",cex=0.8, horiz=F,
	legend= c("range 0-24"),
	text.col= "black",box.col=NA,lwd=1.5,
	border=NA,col= "blue",title="log normal fit on quatradic")

	#3) Keep only predicted values in the range 0-(<NDm) days and combine with real values for > NDm days

plot(NDo,Ix,pch=c(1,19,19,19),col="black",
	ylim=c(0,max(Ix)*1.5),xlab="time (day)",
	xlim=c(0,max(TD1)),
	ylab="intensity")

COR<-t(sapply(NDm[2:length(NDm)],function(x){
	NDmx=x
	In<-Ix[2:4]
	zXn<-subset(TDo, TDo <NDmx)
	zYn<-subset(zY, TDo <NDmx)

	
	zXn<-c(zXn, ND)
	zYn <-c(zYn, In)
	df2<-data.frame(x= zXn,y= zYn)
	
	
	#4) Fit a lognormal model on the combined data (fitted for <NDm days and real for > NDm days)
	opt2 <- optim(opt$par, function(p) 
		sum((dlnorm(df2$x, 
		p[1], p[2]) * p[3] - df2$y)^2))
	y=opt2$par[3] * dlnorm(TD1, 
		opt2$par[1], opt2$par[2])
	
	
	#plot the fitted values
	points(TD1, y,type="l",
		col=rgb(0.5,0.5,0.5,0.05),
		lwd=0.5,lty=1)
	#return mean of difference between fitted and real values, day util which predicted values must be considered, day of max fitted intensity, difference between max fitted intensity and max real intensity
	return(c(
		mean(abs(y[iTD]- In)),x,
		subset(TD1,y==max(y)),
		
		max(y)-max(In))) 

}))	

#ASSUME THE MAXIMUM INTENSITY CANNOT BE REACHED WITHIN Nmin days
COR<-subset(COR,COR[,3]>= Nmin)

NDmax<-subset(COR[,2],COR[,1]==min(COR[,1]))[1]

	Ix<-Ix[2:4]
	zX<-subset(TDo, TDo <NDmax)
	zY<-subset(zY, TDo <NDmax)

	
	zX<-c(zX, ND)
	zY <-c(zY, Ix)
	df2<-data.frame(x= zX,y= zY)
	
	
	#4) Fit a lognormal model on the combined data (fitted for <NDm days and real for > NDm days)
	opt2 <- optim(opt$par, function(p) 
		sum((dlnorm(df2$x, p[1], p[2]) * p[3] - df2$y)^2))
	y=opt2$par[3] * dlnorm(TD1, opt2$par[1], opt2$par[2])
	
	
	
points(TD1, y,type="l",col="red",lwd=2)

legend(x = "topright",cex=0.8, horiz=F,
	legend= c("test","best"),
	text.col= "black",box.col=NA,lwd=c(0.5,1.5),
	border=NA,col= c("grey","red"),title="final log normal")
	
	dev.off()

return(c(max(y),
		mean(subset(TD1,y==max(y))),
		max(Ix),NDmax,max(COR[,1])))

}))

setwd(pD)

colnames(RateYieldNew)<-c("YIELD","LAG","maxInt","corDay","CORest")

write.table(RateYieldNew,"NewMilag.txt")

RateYieldNew<-read.table("NewMilag.txt",header=T)

#Days where the maximum intensity is observed
OBSlag<-sapply(c(1:length(CHARx[,1])),
	function(x){
	#Real values for this replicate and this treatment (ad 0 at time 0)
	Ix= as.numeric(STTn[x,])
	return(subset(ND,Ix==max(Ix))[1])
})	

pdf("compariron-OBS-PRED-yield&lag.pdf",height=3,width=6)

par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))
plot(RateYieldNew[,1], RateYieldNew[,3],
	cex.axis=0.8,las=1,
	xlab="pred. yield (log scale)",
	ylab="obs. max. intensity (log scale)",
	pch=19,cex=0.5,log="xy",
	col=rgb(0,0,0,0.2))
	
plot(as.factor(OBSlag), RateYieldNew[,2],
	cex.axis=0.8,las=1,
	xlab="day with max. intensity",
	ylab="pred. log+lag",
	pch=19,cex=0.5,log="",
	col=rgb(0,0,0,0.2))	

dev.off()

ix<-c(1:length(RateYieldNew[,1]))

RateYieldNewF<-RateYieldNew

#For each treatment and each replicate, do the same and return predicted growth curves

pdf("determination-MiLag-new2.pdf",height=5,width=6)

#Graphical parameters for X-axis
sw=5	
Xtreat<-rbind(NDo,           #P20M20
	(NDo)+max(TD1)+ sw,		#P20M30
	(NDo)+2*(max(TD1)+ sw),	#P30M20
	(NDo)+3*(max(TD1)+ sw))  #P30M30
	
par(bty="n",mar=c(0,4,1,1))
plot(-10,-10,xlim=c(0,(max(TD1)+5)*4),
	ylim=c(0,1),xlab="",ylab="",
	cex.axis=0.8,las=1,xaxt="n",yaxt="n")
	
text(apply(Xtreat[,c(1,4)],1,mean),
	c(0,0,0,0),TT,cex=1)
	
par(mar=c(4,4,1,1),new=T)
plot(-10,-10,xlim=c(0,(max(TD1)+5)*4),
	ylim=c(-2,max(STTn)*1.1),
	xlab="",ylab="Spot intensity",
	cex.axis=0.8,las=1,xaxt="n",log="")
	
axis(1, Xtreat[1,],NDo,cex.axis=0.6)
axis(1,Xtreat[2,],NDo,cex.axis=0.6)
axis(1,Xtreat[3,],NDo,cex.axis=0.6)
axis(1,Xtreat[4,],NDo,cex.axis=0.6)	


TDo<-seq(min(NDo),max(NDo),length.out=10001)
TD1<-seq(min(NDo),max(NDo)* fdx,length.out=40001)

MiLagP<-t(sapply(ix,
	function(x){
	print(x)
	#Real values for this replicate and this treatment (ad 0 at time 0)
	Ix= as.numeric(STTn[x,])
	Ix<-ifelse(Ix<Io,Io,Ix)
	Ix= c(0, Ix)
	
	#1) Polynomial model fitting these values
	model<-lm(Ix ~I(NDo[1:4]^2)+I(NDo[1:4]))
	COEF<-model $coefficients
	#Predicted values from the polynomial model
	y=  COEF[2]*(TDo ^2) + COEF[3]* TDo + COEF[1]
	df<-data.frame(x= TDo,y=y)


	#2) Fit a log normal model on the polynomial residual (will be use only to adjust predicted values <NDm days)
	opt <- optim(c(1, 1, 1), function(p) 
		sum((dlnorm(df$x, p[1], p[2]) * p[3] - df$y)^2))
	zY<-opt$par[3] * dlnorm(TDo, opt$par[1], opt$par[2])

	#3) Keep only predicted values in the range 0-(<NDm) days and combine with real values for > NDm days

	NDmax<-RateYieldNewF[x,4]
	Ix<-Ix[2:4]
	zX<-subset(TDo, TDo <NDmax)
	zY<-subset(zY, TDo <NDmax)
	zX<-c(zX, ND)
	zY <-c(zY, Ix)
	df2<-data.frame(x= zX,y= zY)
	
	
	#4) Fit a lognormal model on the combined data (fitted for <NDm days and real for > NDm days)
	opt2 <- optim(opt$par, function(p) 
		sum((dlnorm(df2$x, p[1], p[2]) * p[3] - df2$y)^2))
	y=opt2$par[3] * dlnorm(TD1, opt2$par[1], opt2$par[2])
	
	#Clade color codes
	clx<-as.vector(CHARx[x,1])
	clx <-c(subset(COLcla[,2],COLcla[,1]==clx),NA)[1]
	clx2 <-rgb(col2rgb(clx)[1,]/255,
		col2rgb(clx)[2,]/255,
		col2rgb(clx)[3,]/255,0.1)
	#X-axis shift (for temperature treatment)	
	trx<-paste("P",CHARx[x,3],"M",CHARx[x,4],sep="")
	trx <-min(subset(Xtreat,TT== trx))
	
	#Plot the predicted curve (step 4) colored according to clade	
	points(TD1+ trx, y,type="l",lwd=0.25,col= clx2)
	#Plot the maximum intensity (yield) and lag time
	points(mean(subset(TD1,
		y==max(y)))+ trx,max(y), 
		pch=19,col=clx,cex=0.3)	
	return(y)

}))

dev.off()

#Remove replicates that did not reach maximal intensity after 36 days
ix <-subset(ix, 
	RateYieldNewF[,2]<max(RateYieldNewF[,2]))
RateYieldNewF<-RateYieldNewF[ix,]
MiLagP <-MiLagP[ix,]
CHARx <-CHARx[ix,]

#Remove negative controls
RateYieldNewF<-subset(RateYieldNewF, CHARx[,2]!="CN1")
MiLagP <-subset(MiLagP, CHARx[,2]!="CN1")
CHARx <-subset(CHARx, CHARx[,2]!="CN1")
RateYieldNewF<-subset(RateYieldNewF, CHARx[,2]!="CN2")
MiLagP <-subset(MiLagP, CHARx[,2]!="CN2")
CHARx <-subset(CHARx, CHARx[,2]!="CN2")

####Summary: average Mi and LAG and number of replicate left per strain and treatment

#Strains
ST<-sort(unique(CHARx[,2]))
#Temperature treatments
Tx<-c("20_20","20_30","30_20","30_30")

STt<-as.vector(unique(paste(CHARx[,2],CHARx[,3],CHARx[,4],sep="_")))

MiLagSum<-t(sapply(STt,function(x){
	X=x
	MiLagx<-subset(RateYieldNewF, 
		paste(CHARx[,2],CHARx[,3],CHARx[,4],
		sep="_")==X)
	#Number of replicates per temperature treatment that reached MI within 24 days, that didn't reach MI and average + sd LAG for those who did
	return(c(length(MiLagx[,1]),
		mean(MiLagx[,1],na.rm=T),
		mean(MiLagx[,2],na.rm=T),
		1/mean(MiLagx[,2],na.rm=T)))	
}))


CHARM<-as.data.frame(t(sapply(STt,function(x){
	CHARxx<-t(subset(CHARx,
		paste(CHARx[,2], CHARx[,3], CHARx[,4],sep="_")==x)[1,])[,1]
	return(as.vector(CHARxx[c(1:4)]))
})))

####Retrieve forest (SITE), plot (SS) host tree species (SPE), temperature of isolation (TI) date of sampling (D)

META<-t(sapply(CHARM[,2],function(x){
	Sx<-paste(unlist(strsplit(x,split="-"))[c(1:2)],collapse="-")
	METAx<-subset(GROUP,as.vector(GROUP[,1])== Sx)
	
	return(c(METAx[1,3], METAx[1,4],METAx[1,6],
		METAx[1,9],METAx[1,10]))
}))

#Remove strain that have less than nr replicates for at least one temperature treatment

nr=1

META<-subset(META,MiLagSum[,1]>= nr)
CHARM<-subset(CHARM,MiLagSum[,1]>= nr)
MiLagSum <-subset(MiLagSum ,MiLagSum[,1]>= nr)

##ANOVA

ST<-as.vector(CHARM[,2])
TP<-as.vector(CHARM[,3])
TM<-as.vector(CHARM[,4])
CL<-as.vector(CHARM[,1])
SITE<-as.vector(META[,1])
SS <-as.vector(META[,2])
SPE <-as.vector(META[,3])
TI <-as.numeric(META[,4])
D <-as.numeric(META[,5])

RATE<-as.numeric(MiLagSum[,4])
YIELD<-as.numeric(MiLagSum[,2])

write.table(cbind(CHARM,META,YIELD,RATE),paste("New-Yield-Rate-summary-nr=",nr,"-ST=",length(unique(ST)),".txt",sep=""))

#Linear Models on log values
MODr<-lm(log(RATE) ~ SITE*SPE*D*TM*CL*TP)
AIC(MODr)
ANOr<-anova(MODr)

MODy<-lm(log(YIELD) ~ SITE*SPE*D*TM*CL*TP)
AIC(MODy)
ANOy<-anova(MODy)

SUMano<-cbind(round(ANOr[,2]/sum(ANOr[,2]),4),round(ANOy[,2]/sum(ANOy[,2]),4),round(ANOr[,5],6),round(ANOy[,5],6))
rownames(SUMano)<-rownames(ANOr)
colnames(SUMano)<-c("Rate","Yield","pRate","pYield")

#Only keep significant factor and residuals + sum of non-significant factors
SUMano <-rbind(subset(SUMano,apply(SUMano[,c(3:4)],1,min)<=0.05),apply(subset(SUMano,apply(SUMano[,c(3:4)],1,min)>0.05),2,sum),subset(SUMano,rownames(SUMano)=="Residuals"))
SUMano<-ifelse(SUMano>1,NA, SUMano)

write.table(SUMano,paste("New-Yield-Rate-anova-nr=",nr,"-ST=",length(unique(ST)),".txt",sep=""))

#######COLOR CODES
#Color code for clades
colc<-as.vector(sapply(CL,function(x){
	return(unique(subset(COLcla[,2], COLcla[,1]==x)))
}))

colcu<-sapply(sort(unique(CL)),function(x){return(unique(subset(colc, CL==x)))})

#Coloc codes for dates
Du<-sort(unique(D))
coldu<-(Du-min(Du))
coldu <-coldu/max(coldu)
coldu <-rgb(1-coldu,1-coldu,0)
cold<-as.vector(sapply(D,function(x){
	return(subset(coldu,Du==x))
}))

#Month categories
DM<-c("June","June","July","Aug.","Aug.","Sept.","Sept.","Oct.")
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


#Coloc codes for temperature
Tu<-c("20","30")
coltu<-c("cyan2","red")

coltm<-as.vector(sapply(TM,function(x){
	return(subset(coltu,Tu==x))
}))

coltp<-as.vector(sapply(TP,function(x){
	return(subset(coltu,Tu==x))
}))

colti<-as.vector(sapply(TI,function(x){
	return(subset(coltu,as.numeric(Tu)==x))
}))

#Color codes for sites
cols<-ifelse(SITE=="MSH","green2","orange")

#Color codes for host tree species
colsp<-ifelse(SPE=="ACSA","purple2","violet")


pdf(paste("New-Figure-Main-nr=",nr,"-ST=",length(unique(ST)),".pdf",sep=""),height=6,width=6)

zones<-matrix(c(1,3,2,4),ncol=2)
layout(zones)
#layout.show(max((zones)))
par(bty="n")


#Average curve per clade


#Clade
par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"a",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

plot(-10,-10,xlim=c(0,max(NDo)*1.2),ylim=c(0,0.75*max(YIELD)),ylab="Average growth",xlab="Time (days)",las=1,cex.axis=0.8,log="",xaxt="n")

axis(1, NDo,NDo,cex.axis=0.8)

CLu<-sort(unique(CL))

sapply(CLu,function(x){
	MiLagPx <-subset(MiLagP, CHARx[,1]==x)
	MLx<-apply(MiLagPx,2,mean)
	SDx<-apply(MiLagPx,2,sd)	/3
	clx <-subset(colcu,CLu==x)
	clx2 =rgb(col2rgb(clx)[1]/255,
				col2rgb(clx)[2]/255,
				col2rgb(clx)[3]/255,0.2)
	polygon(c(TD1, TD1[length(TD1):1]),
		c(MLx+SDx,(MLx-SDx)[length(TD1):1]),
		border=NA,col=clx2,lwd=0.5)
})

sapply(CLu,function(x){
	MiLagPx <-subset(MiLagP, CHARx[,1]==x)
	MLx<-apply(MiLagPx,2,mean)
	clx <-subset(colcu,CLu==x)
	points(TD1,MLx,col=clx,type="l",lwd=1)
	points(subset(TD1,MLx==max(MLx)),
		subset(MLx,MLx==max(MLx)),
		col=clx,pch=21,cex=1,bg="white",lwd=1)	
})

legend(x = "topright",cex=0.8, horiz=F,
	legend= CLu,
	text.col= colcu,box.col=NA,lwd=1,pch=21,pt.bg="white",
	border=NA,col= colcu,title="Clade",title.col="black")

#Clade
par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"b",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

plot(RATE, YIELD,col=colc,pch=19,cex=0.3,lwd=0.5,log="",las=1,xlab="Growth rate (r)",ylab="Yield (Y)",cex.axis=0.8,xlim=c(min(RATE),max(RATE)*1),ylim=c(min(YIELD),max(YIELD)*1))

ordiellipse(cbind(RATE, YIELD),groups= CL,
	border= colcu,col= colcu,label=F,alpha=50,
	draw="polygon",lwd=1,kind="sd",conf=0.3)	
	
legend(x = "topright",cex=0.8, horiz=F,
	legend= sort(unique(CL)),
	box.col=NA,pch=19,
	border="white",col= colcu,text.font=1,
	text.col="black",title="Clade")

#Time of sampling 

par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"c",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

plot(D+sample(seq(-1,1,length.out=length(D))), RATE,col=cols,pch=19,cex=0.3,lwd=1,log="y",las=1,ylab="Growth rate (log scale)",xlab="Sampling date",cex.axis=0.8,xaxt="n",xlim=c(min(D)-5,max(D)+5))

segM<-sapply(Mu,function(x){
	return(mean(subset(D,M==x)))
})

par(cex.axis=0.8)
axis(1,at= segM,labels=Mu)

#dnsity (violin plot) parameters
DENS=100
SVP=5

sapply(sort(unique(SITE)),function(x){
	SITEx=x
	Dx<-subset(D, SITE== SITEx)
	Rx<-subset(RATE, SITE== SITEx)
	colx=subset(cols, SITE ==SITEx)
	
	sapply(sort(unique(Dx)),function(x){
		X=x
		Rxx<-subset(Rx,Dx==X)
		XX<-density(Rxx,n= DENS,na.rm=T)$y/SVP
		YY<-density(Rxx,n= DENS,na.rm=T)$x
		xx=subset(XX,abs(median(Rxx)-YY)==min(abs(median(Rxx)-YY)))
		yy=subset(YY,abs(median(Rxx)-YY)==min(abs(median(Rxx)-YY)))
		polygon(
			c(XX+X,(X-XX)[DENS:1]),
			c(YY,YY[DENS:1]),
			col=rgb(col2rgb(colx)[1]/255,
				col2rgb(colx)[2]/255,
				col2rgb(colx)[3]/255,0.3),border=colx,lwd=0.5)
			
		segments(xx+X,yy,X-xx,yy,col=rgb(col2rgb(colx)[1]/255,
				col2rgb(colx)[2]/255,
				col2rgb(colx)[3]/255,1),lwd=2.5,lend=2)	
				
	})
})

legend(x = "bottomright",cex=0.8, horiz=F,
	legend= unique(SITE),
	box.col=NA,pch=19,
	border="white",col= unique(cols),text.font=1,
	text.col="black",title="Forest")


#Monitoring temperature (normalized by clade to show negative trade-off between rate and yield)
par(mar=c(0,0,0,0))
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")
text(0.05,0.95,"d",cex=1.5,font=2)

par(mar=c(4,5,2,1),new=T)

RATEc<-lm(RATE ~ CL)$residuals
YIELDc<-lm(YIELD ~ CL)$residuals

plot(RATEc, YIELDc,col=coltm,pch=19,cex=0.3,lwd=0.5,log="",las=1,xlab="Residual lm(r~clade)",ylab="Residual lm(Y~clade)",cex.axis=0.8,xlim=c(min(RATEc),max(RATEc)*1),ylim=c(min(YIELDc),max(YIELDc)*1))

ordiellipse(cbind(RATEc, YIELDc),groups= TM,
	border= coltu,col= coltu,label=F,alpha=50,
	draw="polygon",lwd=1,kind="sd",conf=0.3)	
	
legend(x = "topright",cex=0.8, horiz=F,
	legend= paste(sort(unique(TM)),"°C"),
	box.col=NA,pch=19,
	border="white",col= coltu,text.font=1,
	text.col="black",title="Monitoning temperature")
	

dev.off()		
SUMmean<-rbind(
t(sapply(sort(unique(CL)),function(x){
	rx=subset(RATE,CL==x)
	yx=subset(YIELD,CL==x)
	return(c(length(rx),mean(rx),sd(rx),
		mean(yx),sd(yx)))
})),
t(sapply(Du,function(x){
	rx=subset(RATE,D==x)
	yx=subset(YIELD,D==x)
	return(c(length(rx),mean(rx),sd(rx),
		mean(yx),sd(yx)))
})),
t(sapply(Tu,function(x){
	rx=subset(RATE,TM==x)
	yx=subset(YIELD, TM ==x)
	return(c(length(rx),mean(rx),sd(rx),
		mean(yx),sd(yx)))
})),
t(sapply(sort(unique(SITE)),function(x){
	rx=subset(RATE,SITE==x)
	yx=subset(YIELD, SITE ==x)
	return(c(length(rx),mean(rx),sd(rx),
		mean(yx),sd(yx)))
})))

colnames(SUMmean)<-c("n","rate","rate.sd","yield","yield.sd")
rownames(SUMmean)<-c(sort(unique(CL)), Du, Tu,sort(unique(SITE)))


write.table(SUMmean ,paste("New-Yield-Rate-values-nr=",nr,"-ST=",length(unique(ST)),".txt",sep=""))