## Compund growthinside or outside (faster) of model.  0 is outside 
cpinside <- 0
# ## Open up fleetcode what a fleet represents is a mixture of area and vessel type (e.g. commercial or recreational)
growth <- readWorkbook(wb,sheet='Growth', startRow = 2)
growth999 <- growth %>% filter(season==999)
growthdf <- expand.grid(sex=sexs, season = startseason:(endseason+projectseason), area=sort(unique(areas$newarea)),tstep=sort(unique(times$tstep)),grow=0)
growthdf %<>% mutate(grow = growth999$grow[match(paste(sex,area,tstep),paste(growth999$sex-1, growth999$area,growth999$tstep))],grow=ifelse(is.na(grow),0,grow)) %>% mutate(gwtharea=areas$growthmat[match(area, areas$newarea)])

## See if there were any specific years when growth was to be different and implement this.
yrs2change <- unique(growth$season[growth$season!=999])
if(length(yrs2change)>0){
  for(y in 1:length(yrs2change)){
    tmpgrowthdf <- growthdf[growthdf$season==yrs2change[y],]
    tmpgrowthdf$grow <- 0
    dum <-  growth[growth$season==yrs2change[y],]
    tmpgrowthdf %<>% mutate(grow = dum$grow[match(paste(season,area,tstep),paste(dum$season, dum$area, dum$tstep))],grow=ifelse(is.na(grow),0,grow))
    growthdf[growthdf$season==yrs2change[y],] <- tmpgrowthdf
      }}

allmoult <- growthdf %>% filter(grow!=0) %>% group_by(sex,gwtharea) %>% summarise(grow=mean(grow))

nstm <- (length(unique(areas$growthmat)))*length(unique(growthdf$grow[growthdf$grow!=0]))*length(sexs)
gspec <- data.frame(Pattern=0:(nstm-1), Type=1, Sex=rep(0:1,each=nstm/2),Extra=0,Pointer=0:(nstm-1),Mpower=1,tsteps=paste(unique(growthdf$tstep[growthdf$grow>0]),collapse = " "), growthareas=sort(unique(areas$growthmat)), Compound=unique(growthdf$grow[growthdf$grow!=0]))

## to swap around compounding growth inside or outside of model 
if(cpinside==1) {gspec$Mpower <- allmoult$grow[match(paste(gspec$Sex, gspec$growthareas),paste(allmoult$sex,allmoult$gwtharea))]
  gspec$Compound <- 1} else { gspec$Mpower <- 1
  gspec$Compound <- allmoult$grow[match(paste(gspec$Sex, gspec$growthareas),paste(allmoult$sex,allmoult$gwtharea))]}

dat <- expand.grid(sex=sexs, age=(1:ages)-1, area=sort(unique(areas$newarea))-1, step=sort(unique(times$tstep))-1)
dat2 <- matrix(-1, nrow=nrow(dat), ncol=length(startseason:endseason))
for(r in 1:nrow(dat)){
  tmparea  <- dat$area[r]+1
  tmpsex   <-  dat$sex[r]
  tmptstep   <-  dat$step[r]+1
  for(c in 1:ncol(dat2)){    
    sea <- (startseason:endseason)[c]
    tmpdat  <- growthdf[growthdf$sex==tmpsex & growthdf$area==tmparea & growthdf$season==sea & growthdf$tstep==tmptstep,]
    if(tmpdat$grow>0){
      point <- gspec$Pointer[gspec$Sex==tmpsex & gspec$growthareas==tmpdat$gwtharea]
         dat2[r,c]  <- point  }}}
    
dat <- cbind(dat,dat2)
dat %<>% arrange(sex,age,area,step)
tmp <- list()
tmp <- c(tmp, "# Growth specification\n\n# Number of growth Patterns (Mpower is legacy and needed for power function on growth)\n",nrow(gspec),"\n")
tmp <- c(tmp, "# ", paste(colnames(gspec),collapse='\t'),"\n")
for(i in 1:nrow(gspec)){ tmp <- c(tmp,'\t',paste(gspec[i,1:6],collapse = "\t"),'\t# ',paste(gspec[i,7:9],collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Specifications for growth\t\n# Sex\tAge\tArea\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")
for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Number prespecified size-transition\n",length(unique(gspec$Pointer)),"\n")
tmp <- c(tmp, "# Sex of matrices\n",paste(gspec$Sex,collapse = "\t"),"\n")

tmp <- c(tmp, "# Prespecified size-transition\n")
STM <- growth <- readWorkbook(wb,sheet='SizeTransMatrices', startRow = 2, colNames = F)
pos <- c(grep('#',STM[,1],fixed=T),1+nrow(STM))

for(g in 1:nrow(gspec)){
  nm <- gspec$Compound[g]
  mat <- as.matrix(STM[(pos[g]+1):(pos[g]+length(lens)),1:length(lens)])
  mat <- apply(mat, c(1,2), as.numeric)
  toXL(mat)
  tmat <- mat
  ## Make sum to 1
  tmat <- apply(tmat, 2, function(x) x/sum(x))
  if(nm>1){
    for(r in 1:nrow(mat)){
      sz <- rep(0, ncol(mat)); sz[r] <- 1
      for(m2 in 2:nm){ sz <- (mat)%*%sz }  # currently represents one moult, compound it up to the correct number of moults
        tmat[,r] <- sz    }}
  tmp <- c(tmp, "\n# Matrix ",STM[pos[g],1] ,"\n")
  #Nname <- paste(sx,'stm',nu,"_months_",nm,".csv", sep = "")
  #write.table(tmat, Nname, sep=',',row.names = F, col.names = F)
  for(i in 1:nrow(mat)){ tmp <- c(tmp,paste(round(tmat[i,],7),collapse = "\t"),"\n")}}

tmp <- c(tmp, "\n# Growth parameters\n#LB\tUP\tEstimate\tPhase\n")

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/GROWTHSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
