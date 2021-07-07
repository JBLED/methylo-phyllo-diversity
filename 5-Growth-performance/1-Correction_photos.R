library("caTools")

pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/5-Growth-performance/data"

#Petri name
NAME<-c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24") 

#Temperature
TP<-c("T20","T30")

#Date of pictures
TM<-c("22-02-2019","28-02-2019","11-03-2019")

######

SI=50 #minimum intensity allowed
MI=200 #maximum intensity allowed
v=30 #border to remove
t=0.95
WS=500 #number of pixel to reajust the picture
ps=15000 #number of pixel to sample for correction
WM=50 #window size for mean calculation

	
sapply(TM,function(x){
	TMx<-x
		
	sapply(TP,function(x){
		TX<-x
			
		sapply(NAME,function(x){	
			NAMEx<-x
			print(paste(TMx,TX,NAMEx))	
			P1<-paste(pD,"/Pic_",TMx,"/",TX,"/",NAMEx,sep="")
			setwd(P1)	
			data<-read.csv("BW.csv")
			data<-as.matrix(data[c(2:length(data[,1])),c(2:length(data[1,]))])  ##A
			data<-subset(data,apply(data,1,max)!=0)
			data<-t(subset(t(data),apply(data,2,max)!=0))
			data<-ifelse(data<=SI,NA,data)
			data<-ifelse(data>=MI,NA,data) ##B
			data.N<-read.csv("back.csv")
			data.N <-as.matrix(data.N[c(2:length(data.N[,1])),c(2:length(data.N[1,]))]) # C
			data.N <-subset(data.N,apply(data.N,1,max)!=0)
			data.N <-t(subset(t(data.N),apply(data.N,2,max)!=0))
			data.N<-ifelse(data.N<=SI,NA,data.N)
			data.N<-ifelse(data.N>=MI,NA,data.N) ##D
			

			X<-round(seq(0,length(data[,1]),length.out=(WS+1)),0)
			
			WSy<-round(WS*length(data[1,])/length(data[,1]),0)
			
			Y<-round(seq(0,length(data[1,]),length.out=(WSy +1)),0)

			#COR<-read.table("COR.txt")
			POS<-read.csv("POS.csv")
						
			
		#Average pixel intensity in WS slidding windows
			
			datam<-t(sapply(c(1:WS),function(x){
				#print(x)
				X1<-X[x]+1
				X2<-X[x+1]
				return(sapply(c(1: WSy),function(x){
					Y1<-Y[x]+1
					Y2<-Y[x+1]
					return(round(mean(data[c(X1:X2),c(Y1:Y2)],na.rm=T),2))
				}))
			})) ###E
			
			data.Nm<-t(sapply(c(1:WS),function(x){
				#print(x)
				X1<-X[x]+1
				X2<-X[x+1]
				return(sapply(c(1: WSy),function(x){
					Y1<-Y[x]+1
					Y2<-Y[x+1]
					return(round(mean(data.N[c(X1:X2),c(Y1:Y2)],na.rm=T),2))
				}))
			}))	###F
			
			#heatmap(datam,Rowv=NA,Colv=NA,
			#	labCol="",labRow="",scale="none")
			
			POSx<-round(POS[,6]*WSy/length(data[1,]),0)
			POSy<-round(POS[,7]*WS/length(data[,1]),0)			
			
			write.table(cbind(POSx,POSy),("POSm.txt"),
				row.names=F,col.names=F)
		
			data.N<-data.Nm
			data<-datam
		
			x1<-c((v+1):(length(data.N[1,])-v))
			x2<-c((v+1):(length(data.N[,1])-v))
			x0<-sapply(x1,function(x){
				z<-as.vector(data.N[x2,x])
				return(length(subset(x2,is.na(z)==T)))
			})
			x1<-subset(x1,x0!=0)
			x0<-sapply(x2,function(x){
				z<-as.vector(data.N[x,x1])
				return(length(subset(x1,is.na(z)==T)))
			})
			x2<-subset(x2,x0!=0)
			COORD<-unlist(sapply(x1,function(x){
				z<-as.vector(data.N[x2,x])
				return(paste(x,subset(x2,
					is.na(z)==T),sep="_"))
			}))
			COORD<-COORD[sample(c(1:length(COORD)),ps)]
			data.m<-sapply(COORD,function(x){
				z<-as.numeric(strsplit(x,split="_")[[1]])
				return(round(mean(data.N[c((z[2]-v):(z[2]+v)), c((z[1]-v):(z[1]+v))],na.rm=T),0))	
			})
			COORD <-matrix(unlist(strsplit(COORD,split="_")),ncol=2,byrow=T)
			COORD<-cbind(as.numeric(COORD[,1]),as.numeric(COORD[,2]),data.m)
			COORD<-subset(COORD,COORD[,3]!="NaN")

			m1=min(data.m,na.rm=T)-5
			m2=max(data.m,na.rm=T)+5
			M1=min(COORD,na.rm=T)
			M2=max(COORD,na.rm=T)

			data.X<-sapply(c(1:length(data.N[1,])),
			function(x){
				COORDx<-subset(COORD,
					as.numeric(COORD[,1])==x)	
				data.x<-data.N[,x]	
				data.x <-c(data.x[
					setdiff(c(1:length(data.x)), 
					COORDx[,2])],COORDx[,3])[
					order(c(setdiff(c(1:length(data.x)), 
					COORDx[,2]),COORDx[,2]))]	
				return(runmean(data.x,WM,endrule="NA"))	
			}) ###G

			data.Y<-sapply(c(1:length(data.N[,1])),
			function(x){
				COORDx<-subset(COORD,
					as.numeric(COORD[,2])==x)	
				data.x<-data.N[x,]	
				data.x <-c(data.x[
					setdiff(c(1:length(data.x)), 
					COORDx[,1])],COORDx[,3])[
					order(c(setdiff(c(1:length(data.x)), 
					COORDx[,1]),COORDx[,1]))]	
				return(runmean(data.x,WM,endrule="NA"))	
			}) ###H

			data.XY<-(data.X+t(data.Y))/2 ###I
			data.XY<-data-data.XY ###J
			data.XY<-ifelse(data.XY<(-50),NA, data.XY) ###K

			write.table(round(data.XY,2),
				"CORm.txt",col.names=F,row.names=F)		
		
		
		
		})

	})	
})




######

sapply(TM,function(x){
	TMx<-x
	
	sapply(TP,function(x){
		TX<-x
		
		sapply(NAME,function(x){
			NAMEx<-x
			print(paste(TMx,TX,NAMEx))	
			P1<-paste(pD,"/Pic_",TMx,"/",TX,"/",NAMEx,sep="")
			setwd(P1)	
			data<-read.csv("BW.csv")
			data<-as.matrix(data[c(2:length(data[,1])),c(2:length(data[1,]))])
			data<-subset(data,apply(data,1,max)!=0)
			data<-t(subset(t(data),apply(data,2,max)!=0))
			data<-ifelse(data<=SI,NA,data)
			data<-ifelse(data>=MI,NA,data)
			data.N<-read.csv("back.csv")
			data.N <-as.matrix(data.N[c(2:length(data.N[,1])),c(2:length(data.N[1,]))])
			data.N <-subset(data.N,apply(data.N,1,max)!=0)
			data.N <-t(subset(t(data.N),apply(data.N,2,max)!=0))
			data.N<-ifelse(data.N<=SI,NA,data.N)
			data.N<-ifelse(data.N>=MI,NA,data.N)
			x1<-c((v+1):(length(data.N[1,])-v))
			x2<-c((v+1):(length(data.N[,1])-v))
			x0<-sapply(x1,function(x){
				z<-as.vector(data.N[x2,x])
				return(length(subset(x2,is.na(z)==T)))
			})
			x1<-subset(x1,x0!=0)
			x0<-sapply(x2,function(x){
				z<-as.vector(data.N[x,x1])
				return(length(subset(x1,is.na(z)==T)))
			})
			x2<-subset(x2,x0!=0)
			COORD<-unlist(sapply(x1,function(x){
				z<-as.vector(data.N[x2,x])
				return(paste(x,subset(x2,is.na(z)==T),sep="_"))
			}))
			COORD<-COORD[sample(c(1:length(COORD)),50000)]
			data.m<-sapply(COORD,function(x){
				z<-as.numeric(strsplit(x,split="_")[[1]])
				return(round(mean(data.N[c((z[2]-v):(z[2]+v)), c((z[1]-v):(z[1]+v))],na.rm=T),0))	
			})
			COORD <-matrix(unlist(strsplit(COORD,split="_")),ncol=2,byrow=T)
			COORD<-cbind(as.numeric(COORD[,1]),as.numeric(COORD[,2]),data.m)
			COORD<-subset(COORD,COORD[,3]!="NaN")

			m1=min(data.m,na.rm=T)-5
			m2=max(data.m,na.rm=T)+5
			M1=min(COORD,na.rm=T)
			M2=max(COORD,na.rm=T)

			data.X<-sapply(c(1:length(data.N[1,])),function(x){
				COORDx<-subset(COORD,as.numeric(COORD[,1])==x)	
				data.x<-data.N[,x]	
				data.x <-c(data.x[setdiff(c(1:length(data.x)), COORDx[,2])],COORDx[,3])[order(c(setdiff(c(1:length(data.x)), COORDx[,2]),COORDx[,2]))]	
				return(runmean(data.x,200,endrule="NA"))	
			}) ##L

			data.Y<-sapply(c(1:length(data.N[,1])),function(x){
				COORDx<-subset(COORD,as.numeric(COORD[,2])==x)	
				data.x<-data.N[x,]	
				data.x <-c(data.x[setdiff(c(1:length(data.x)), COORDx[,1])],COORDx[,3])[order(c(setdiff(c(1:length(data.x)), COORDx[,1]),COORDx[,1]))]	
				return(runmean(data.x,200,endrule="NA"))	
			}) ##M

			data.XY<-(data.X+t(data.Y))/2 #N
			data.XY<-data-data.XY #O
			data.XY<-ifelse(data.XY<(-50),NA, data.XY)

			setwd(P1)

			write.table(round(data.XY,2),"COR.txt",col.names=F,row.names=F)

		})
	})
})

