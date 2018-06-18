source("fANCOVA_all.r")
#to get the overall QTL threshold 
#Create a sine curve and add some noise:
tdf<-read.csv("allBinres.csv")
x<-1:nrow(tdf)
y<-tdf[,3]
#Plot the points on this noisy sine curve:
#plot(x,y, main="Sine Curve + 'Uniform' Noise")
#mtext("showing loess smoothing (local regression smoothing)")

y.loess2 <- loess.as( x,y, plot=FALSE)
#Compute loess smoothed values for all points along the curve:
y.predict2 <- predict(y.loess2, data.frame(x=x))
#Plot the loess smoothed curve along with the points that were already plotted:
#lines(x,y.predict2,col="red")
#to obtain the threshold of QTL regions
Qthred<-quantile(y.predict2,c(1-0.05/nrow(tdf)))
Qthred
tmax=max(y)
tmin=min(y)
#lines(x,rep(Qthred,nrow(tdf)),col="red")

#to get the detla value for the confience intervals
peak <- optimize(function(x, model)
            predict(model, data.frame(x=x)),
            c(min(x),max(x)),
            maximum=TRUE,
            model=y.loess2) 
#points(peak$maximum,peak$objective,
#            pch=FILLED.CIRCLE<-19, col="green")

tright<-as.integer(peak$maximum)
delta<-peak$objective-quantile(y.predict2,c(.95))

chrarr<-unique(tdf$chr)
i=1
pdf("QTG_Seq_plots.pdf")
while(i<=length(chrarr)){
	ttdf<-subset(tdf,chr==chrarr[i])
	if(nrow(ttdf)>20){
		tx<-1:nrow(ttdf)
		ty<-ttdf[,3]
		plot(ttdf$position/1000000,ty, main="QTG-Seq For QTG Fine-mapping",xlab="Physical Coordinate (Mb)",ylab=colnames(tdf)[3],ylim=c(0,tmax))
		mtext(chrarr[i])
		ty.loess2 <- loess.as( tx,ty, plot=FALSE)
		#Compute loess smoothed values for all points along the curve:
		ttdf$predict2 <- predict(ty.loess2, data.frame(x=tx))
		#Plot the loess smoothed curve along with the points that were already plotted:
		lines(ttdf$position/1000000,ttdf$predict2,col="red")
		lines(ttdf$position/1000000,rep(Qthred,nrow(ttdf)),col="red")
		ttdf$Qthred=Qthred
		ttdf$delta=delta
		ttdf[ttdf$predict2>=Qthred,4]=1
		tname<-paste("QTG_Seq_R_summary_",chrarr[i],".csv",sep="")
		write.csv(ttdf,tname,col.names = FALSE, row.names = FALSE)
		write.table(tname,"QTG_Seq_R_summaryfile.txt",append=T,quote = F,row.names = F,col.names =F)
	}
	i=i+1
}
dev.off()

