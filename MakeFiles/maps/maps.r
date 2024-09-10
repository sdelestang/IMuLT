library(gplots)
library(PBSmapping)

WAcoast<-read.table("maps/WAcoastPointsNew.txt", header=T)
WAislands<-read.table(paste(getwd(),"/maps/WAislandsPointsNew.txt",sep=''), header=T)
Shelf<-read.table(paste(getwd(),"/maps/shelf.txt",sep=''), header=T)
Towns<-read.table(paste(getwd(),"/maps/towns.txt",sep=''), header=T)
D40 <- dget(paste(getwd(),"/maps/Depth40",sep=''))
D100 <- dget(paste(getwd(),"/maps/Depth100",sep=''))
D200 <- dget("maps/Depth200")
EEZ <- data.frame(x=c(114.9211, 110.1427, 110.0848, 110.0848, 110.0558, 109.9400, 109.7952, 109.6504, 109.5056, 109.3608,
                      109.3029, 109.3029, 109.3897, 109.6214, 109.8242, 110.1138, 110.6061, 110.9246, 111.2722, 111.3590,
                      111.4749, 111.5038, 111.2142, 110.9826, 110.7509, 110.4613, 117.0642),
                  y=c(-21.79, -21.79,       -22.19474, -22.39979, -22.65609, -22.96365, -23.34811, -23.73256, -24.32206,
                      -24.83466, -25.52668, -26.29559, -27.14139, -27.73088, -28.29475, -28.85861, -29.49937, -30.16576,
                      -30.72962, -31.29349, -31.83172, -32.39559, -33.01072, -33.36954, -33.98467, -34.38, -34.38))

plottowns <- function(x,CEX=0.8){
  points(Towns$lon, Towns$lat, pch=16, col='darkgreen')
  text(Towns$lon, Towns$lat, Towns$loc, pos=4, cex=CEX)}

plotBBaAbr <- function(x){
  arrows(112.6242+0.5, -27.18997, 113.3444+0.5,-27.04339 , length=0.1, code=1)
  text(113.3444+0.5,-27.04339, 'Big Bank', pos=4)
  arrows(113.2415+0.7, -28.80240, 113.8515+0.7,-28.36265 , length=0.1, code=1)
  text( 113.8515+0.7,-28.36265, 'Abrolhos Islands', pos=4)
  }

plotislands <- function(x=grey(0.7)){
  for (j in 1:length(unique(WAislands$ID))){
    polygon(WAislands$Longitude[WAislands$ID==unique(WAislands$ID)[j]],WAislands$Latitude[WAislands$ID==unique(WAislands$ID)[j]], col=x)}}

plotbath <- function(x){
  lines(Shelf$lon,Shelf$lat, lty=3)
  unique(WAislands$ID)
}

plotzones <- function(x){
  lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
  #lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
  lines(c(112.0,114.9), c(-30,-30), lty=2)
  text(113.7,-29.2,labels="A Zone",cex=0.8)
  text(113.2,-27,labels="B Zone",cex=0.8)
  text(113.2,-32.0,labels="C Zone",cex=0.8)
  #text(112.8,-27.3,labels="Big Bank")
  }

plotzones2 <- function(x,sz=0.8){
  lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
  #lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
  lines(c(113.0,115.9), c(-30.5,-30.5), lty=2)
  text(113.5,-29.3,labels="Abrolhos Zone", cex=sz)
  text(113.5,-30,labels="North Zone", cex=sz)
  text(113.4,-31.0,labels="South Zone", cex=sz)
}

plotmodelmap <- function(long, lat, nam="", all.zones=FALSE, zones.num=FALSE){
        par(las=1)
        plot(long, lat, type='n', main='', axes=F, xlab=" ", ylab=" ") 
        polygon(c(100, 100, 130, 130), c(-20, -40, -40, -20), col='#F0F8FF') 
        polygon(WAcoast$Longitude,WAcoast$Latitude, col='#FFE7BA')
        for (j in 1: length(unique(WAislands$ID))){
            polygon(WAislands$Longitude[WAislands$ID==unique(WAislands$ID)[j]],WAislands$Latitude[WAislands$ID==unique(WAislands$ID)[j]], col='#FFE7BA')}
        par(new=T) 
        plot(long, lat, type='n', main=nam, xlab="Longitude", ylab="Latitude")
        lines(D40$x,D40$y, lty=3, lwd=0.1, col='grey50')
        lines(D200$x,D200$y, lty=3, lwd=0.1, col='grey50')
        
        lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
        lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
        lines(c(113.0,114.9), c(-30,-30), lty=2)        
        if(all.zones==TRUE) {text(113.6,-29.2,labels="A Zone")
                             text(112.4,-29,labels="B Zone")
                             text(113.4,-31.0,labels="C Zone")
                             text(112.6,-27.2,labels="Big Bank")  }
        points(( 115.03022), -30.30578, pch=19); text(( 115.03022), -30.30578, labels="Jurien", pos=4)
        points(( 114.92433), -29.25258, pch=19); text(( 114.92433), -29.25258, labels="Dongara", pos=4)
        points(( 114.16472), -27.71286, pch=19); text(( 114.16472), -27.71286, labels="Kalbarri", pos=4)
        points(( 115.34065), -31.02539, pch=19); text(( 115.34065), -31.02539, labels="Lancelin", pos=4)
        points(( 115.75642),-32.03915, pch=19); text(( 115.75642), -32.03915, labels="Fremantle", pos=4)
        if(zones.num==TRUE) {text(c(115.5, 115.25, 115.0, 114.8, 114.9, 114.59, 114.05, 113.3, 113.7, 114.0, 113.4),
            c(-31.6, -32.3, -30.5, -30.5, -29.5, -29.4, -27.4, -27.0, -28.5, -29, -27.8),
            labels=1:11,font=2)}
       lines(c(114.85,115.3),c(-31, -31), lty=1)
       lines(c(114.5,114.97),c(-30, -30), lty=1)
       lines(c(113.82, 113.01),c(-28-1/6, -28-1/6), lty=1)
       lines(c(113.75, 114.15),c(-28, -28), lty=1)       
                      }

addlines <- function(x){
        lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
        lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
        lines(c(113.0,114.9), c(-30,-30), lty=2)
        text(113.4,-29.2,labels="A Zone")
        text(112.5,-29,labels="B Zone")
        text(113.4,-31.0,labels="C Zone")
        text(112.8,-27.3,labels="Big Bank")
                      }

plotmodelmapNew <- function(long=c(112,117), lat=c(-27,-33), nam="", all.zones=FALSE, zones.num=FALSE){
  par(las=1)
  plot(long, lat, type='n', main='', axes=F, xlab=" ", ylab=" ") 
  polygon(c(100, 100, 130, 130), c(-20, -40, -40, -20), col='#F0F8FF') 
  polygon(WAcoast$Longitude,WAcoast$Latitude, col='#FFE7BA')
  for (j in 1: length(unique(WAislands$ID))){
    polygon(WAislands$Longitude[WAislands$ID==unique(WAislands$ID)[j]],WAislands$Latitude[WAislands$ID==unique(WAislands$ID)[j]], col='#FFE7BA')}
  par(new=T) 
  plot(long, lat, type='n', main=nam, xlab="Longitude", ylab="Latitude")
  lines(D40$x,D40$y, lty=3, lwd=0.1, col='grey50')
  lines(D200$x,D200$y, lty=3, lwd=0.1, col='grey50')
  
  lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
  lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
  lines(c(113.0,114.9), c(-30,-30), lty=2)        
  if(all.zones==TRUE) {text(113.6,-29.2,labels="A Zone")
    text(112.4,-29,labels="B Zone")
    text(113.4,-31.0,labels="C Zone")
    text(112.6,-27.2,labels="Big Bank")  }
  points(( 115.03022), -30.30578, pch=19); text(( 115.03022), -30.30578, labels="Jurien", pos=4)
  points(( 114.92433), -29.25258, pch=19); text(( 114.92433), -29.25258, labels="Dongara", pos=4)
  points(( 114.16472), -27.71286, pch=19); text(( 114.16472), -27.71286, labels="Kalbarri", pos=4)
  points(( 115.34065), -31.02539, pch=19); text(( 115.34065), -31.02539, labels="Lancelin", pos=4)
  points(( 115.75642),-32.03915, pch=19); text(( 115.75642), -32.03915, labels="Fremantle", pos=4)
#  if(zones.num==TRUE) {
    
text(c(114.6022,114.3350,114.2407,114.0364,114.0992,113.8242,113.2742,112.8498,112.9363,113.3999,112.9048)+0.9,c(-30.35128 ,-30.47526 ,-29.62806 ,-29.60740 ,-28.68788 ,-28.50191 ,-28.11964 ,-28.28495 ,-27.45842 ,-26.98316 ,-26.80753)-1.5,  labels=(areas$new+1), font=2)
    }
    
  locator(11)
    # text(c(115.0475, 114.8037, 114.5377, 114.1941, 113.7286, 113.1522, 113.5623, 113.0414)+.25,
    #      c(-30.97891, -31.23014, -29.18439, -28.93316, -27.24631, -27.00405, -28.49350, -28.85241),
    #      labels=tareas$new8, font=1, cex=0.8, pos=4, col=3)

    text(c(115.0475, 114.8037, 114.6264, 114.4047, 114.5377, 114.1941, 113.7286, 113.1522, 113.5623, 113.0414, 113.0746)+.25,
         c(-30.97891, -31.23014, -30.04576, -29.89323, -29.18439, -28.93316, -27.24631, -27.00405, -28.49350, -28.85241, -27.59624),
         labels=areas$old11, font=1, cex=0.8, pos=1, col=2)
    
    legend(115.5,-27.5, text.col=1:3, legend=c('Areas in current model', 'Old 11 areas', 'New 8 areas'))
    
    }
  lines(c(113.75, 114.15),c(-28, -28), lty=1)       
}

addlines <- function(x){
  lines(c(113,113.5,114.5,113,113), c(-27.5,-27.5,-29.5,-29.5,-27.5), lty=2)
  lines(c(113.5,112.5), c(-27.5,-25.5), lty=2)
  lines(c(113.0,114.9), c(-30,-30), lty=2)
  text(113.4,-29.2,labels="A Zone")
  text(112.5,-29,labels="B Zone")
  text(113.4,-31.0,labels="C Zone")
  text(112.8,-27.3,labels="Big Bank")
}




