pD<-"/Users/jean-baptisteleducq/Desktop/Github-R-Scripts-Methylo-Adaptation/5-Growth-performance/data"

NAME<-c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24") 

TP<-c("T20","T30")

TM<-c("22-02-2019","28-02-2019","11-03-2019")

######

Rx=15 #size of the window
mp=200 #minimum number of pixels required to allow a t-test 
Lp=1e-50 #minimum pvalue for t-test

######

sapply(TM,function(x){
	TMx<-x
	Pc<-paste(pD,"/Pic_",TMx,"/COR_files",sep="")
	
	sapply(TP,function(x){
		TX<-x
		
		sapply(NAME,function(x){
			NAMEx<-x
			print(paste(TMx,TX,NAMEx))		
			P1<-paste(pD,"/Pic_",TMx,"/",TX,"/",NAMEx,sep="")
			setwd(P1)
			COR<-read.table("CORm.txt")
			COR <-as.matrix(COR[,c(1:length(COR[1,]))])
	
			POS<-read.table("POSm.txt")
		
			pdf("Spot_profile.pdf")
			par(mfrow=c(6,6),mar=c(2,2,1,1))


			IM<-t(sapply(c(1:length(POS[,1])),function(x){
				X=POS[x,2]
				Y=POS[x,1]
				X1<-X-2*Rx
				X2<-X+2*Rx
				Y1<-Y-2*Rx
				Y2<-Y+2*Rx
				#heatmap(COR[c(X1:X2),c(Y1:Y2)],Rowv=NA,Colv=NA)
				distx<-matrix(sapply(c(X1:X2),function(x){
					Dx<-x
					return(sapply(c(Y1:Y2),function(x){
						Dy<-x
						return(c(COR[Dx,Dy],sqrt(((Dx-X)*(Dx-X))+((Dy-Y)*(Dy-Y)))))
					}))
				}),ncol=2,byrow=T)
				distx<-subset(distx,is.na(distx[,1])==F)
	
				IMx<-sapply(c(1:(Rx*2)),function(x){	
					return(mean(subset(distx[,1],round(distx[,2],0)==x),na.rm=T))	
				})
				zz<-subset(c(1:(Rx*2)), IMx!="NaN")
				IMx<-subset(IMx, IMx!="NaN")
	
				sumx<-sapply(zz,function(x){	
					#print(x)
					return(c(length(subset(distx[,1],round(distx[,2],0)<=x)),length(subset(distx[,1],round(distx[,2],0)>x))))
				})	

				zz<-subset(zz, apply(sumx,2,min)>=mp)
				IMx<-subset(IMx,apply(sumx,2,min)>=mp)	
	
				TTx<-sapply(zz,function(x){	
					#print(x)
					return(c(t.test(subset(distx[,1],round(distx[,2],0)<=x),subset(distx[,1],round(distx[,2],0)>=x))$statistic,t.test(subset(distx[,1],round(distx[,2],0)<=x),subset(distx[,1],round(distx[,2],0)>=x))$p.value)	)
				})	
	
				plot(zz, IMx,pch=19,col=ifelse(TTx[1,]==max(TTx[1,]),"red","black"),cex=ifelse(TTx[1,]==max(TTx[1,]),1,0.2),main=paste(x, ",",X,",",Y))
	
				return(subset(cbind(zz,t(TTx)),TTx[1,]==max(TTx[1,])))
	
			}))
			dev.off()
		
		
			Tx<-ifelse(IM[,3]>Lp,mean(subset(IM[,1],IM[,3]<=Lp)),IM[,1])

		
			pdf("spot_distribution.pdf")
			par(mfrow=c(6,6),mar=c(2,2,1,1))
			SPOT<-t(sapply(c(1:length(POS[,1])),function(x){
				X=POS[x,2]
				Y=POS[x,1]
				X1<-X-2*Rx
				X2<-X+2*Rx
				Y1<-Y-2*Rx
				Y2<-Y+2*Rx
	
				distx<-matrix(sapply(c(X1:X2),function(x){
					Dx<-x
					return(sapply(c(Y1:Y2),function(x){
						Dy<-x
						return(c(COR[Dx,Dy],sqrt(((Dx-X)*(Dx-X))+((Dy-Y)*(Dy-Y)))))
					}))
				}),ncol=2,byrow=T)
	
				MM<-mean(subset(distx[,1], distx[,2]<=Tx[x]),na.rm=T)
				SM<-sd(subset(distx[,1], distx[,2]<=Tx[x]),na.rm=T)

				distc<-subset(distx, distx[,2]>Tx[x])

				MMc<-mean(distc[,1],na.rm=T)
				SMc<-sd(distc[,1],na.rm=T)

				bx=round((max(distx[,1],na.rm=T)-min(distx[,1],na.rm=T))/2,0)

				hist(distc[,1],breaks= bx,col=rgb(1,0,0,1),border=NA,xlim=c(min(distx[,1],na.rm=T),max(distx[,1],na.rm=T)),main=x)
				par(new=T)

				hist(subset(distx[,1], distx[,2]<=Tx[x]),breaks= bx,col=rgb(0,0,1,0.5),border=NA,xlim=c(min(distx[,1],na.rm=T),max(distx[,1],na.rm=T)),main="",xaxt="n",yaxt="n",xlab="",ylab="")
	
				return(c(x,X,Y,Tx[x],MM,SM,MMc,SMc))
	
			}))

			dev.off()

			INTc<-SPOT[,5]-SPOT[,7]

			z=abs(INTc)/max(INTc)
			colz<-ifelse(z>0.5,"white","black")

			pdf("spot_position.pdf")

			plot(SPOT[,2], SPOT[,3],pch=21,cex=Tx/3, bg=rgb(1-z,1-z,1-z,1),xlim=c(0,length(COR[,1])),ylim=c(0,length(COR[1,])),xlab="X",ylab="Y",las=1,cex.axis=0.6,main=paste(paste(TX, NAMEx)))

			text(SPOT[,2], SPOT[,3],SPOT[,1],cex=1, col= colz)

			dev.off()


			colnames(SPOT)<-c("spot","X_position","Y_position","radius","av_intensity","sd_intensity","av_intensity_back","sd_intensity_back")

		write.table(SPOT,"spot_intensities.txt",row.names=F)
		
		
		})
	})
})


######
