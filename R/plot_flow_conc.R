



plot_flow_conc <- function( date, discharge, conc ) {

  par(mar=c(4, 4, 4, 4) + 0.1)
  plot(date, discharge, ylim=c(0,max(discharge, na.rm=TRUE)), xlab="", ylab="",type="l",
       col="blue", main="")
  points(date,discharge,pch=20,col="blue",cex=0.3)
  #axis(2, ylim=c(0,max(discharge, na.rm=TRUE)),col="black",lwd=1.5)
  mtext(2,text="Discharge, in cfs",line=2)

  par( new=TRUE )
  plot(date, conc, axes=F, ylim=c(0,max(conc, na.rm=TRUE)), xlab="", ylab="",type="l",
       col="green", main="")
  points(date,conc,pch=21,col="green", cex=0.3)
  axis(4, ylim=c(0,max(conc, na.rm=TRUE)),col="black")
  mtext(4,text="Specific conductance, umhos/cm",line=2)

}
