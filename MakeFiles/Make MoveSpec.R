library(dplyr)
library(magrittr)

migrate999 <- migrate %>% filter(Season==999)
dat <- expand.grid(age=(1:ages)-1, area=sort(unique(areas$newarea))-1, step=sort(unique(times$tstep))-1)
dat2 <- matrix(0, nrow=nrow(dat), ncol=length(startseason:endseason))
for(i in 1:nrow(migrate999)){
  dat2[dat$age==(migrate999$Age[i]-1) & dat$area==(migrate999$Source[i]-1) & dat$step==(migrate999$tstep[i]-1),] <- i
}

## See if there were any specific years when migration was to be different and implement this.
yrs2change <- unique(migrate$Season[migrate$Season!=999])
if(length(yrs2change)>0){
  for(y in 1:length(yrs2change)){
    tmpdf <- migrate[migrate$Season==yrs2change[y],]
    dat2[dat$age==(tmpdf$Age-1) & dat$area==(tmpdf$Source-1) & dat$step==(tmpdf$tstep-1),(startseason:endseason)==tmpdf$Season] <- tmpdf$pointer
  }}

dat <- cbind(dat,dat2)
dat %<>% arrange(age, area, step)

migratesum <- migrate %>% group_by() 

tmp <- list()
tmp <- c(tmp, "# Movement section\n# Number of movement patterns\n",nrow(migrate)+1,"\n")
tmp <- c(tmp, "# Pattern\tType\tDest\tExtra\t(Type 0: none; 1 constant [prespecified or estimated; 1 parameter]; 2 knife-eded-specific [pre-specified or estimated; 2 parameters])\n")
tmp <- c(tmp,"\t",paste(c(0,0,0,0),collapse = "\t"),"\n")
for(i in 1:nrow(migrate)){ tmp <- c(tmp,"\t",paste(c(i,1,migrate$Dest[i]-1,0),collapse = "\t"),"\n")}

tmp <- c(tmp, "# Movement specifications\n")
tmp <- c(tmp, "#Age\tArea\tTStep\t", paste(startseason:endseason,collapse = "\t"),"\n")
for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

tmp <- c(tmp, "# Movement parameters\n# Lower\tUpper\tEstimate\tPhase\t\tSource to Dest & Age\n")
Mpar <- migrate %>% dplyr::select(lower, upper, est, Phase)
for(i in 1:nrow(Mpar)){ tmp <- c(tmp,"\t",paste(Mpar[i,],collapse = "\t"),paste("\t\t #", migrate$Source[i],"to", migrate$Dest[i],"&",migrate$Age[i],"\n"))}

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/MOVESPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
