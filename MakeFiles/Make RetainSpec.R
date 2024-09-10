library(dplyr)
library(magrittr)

hgrad <- readWorkbook(wb,sheet='HighGrading', startRow = 2) %>% mutate(tstep=times$tstep[match(paste(month,hm),paste(times$month,times$hm))]) %>% group_by(season, oldarea, tstep) %>% summarise(prop=mean(prop)) %>% mutate(propfl=1-trunc(prop/0.02)*0.02)
hglist <- sort(unique(c(hgrad$propfl,1)))

hgrad99 <- hgrad %>% filter(oldarea==99)

dat <- expand.grid(sex=sexs, age=(1:ages)-1, fleet=sort(unique(fleets$fleet))-1, step=sort(unique(times$tstep))-1)
dat2 <- matrix(length(hglist)-1, nrow=nrow(dat), ncol=length(startseason:endseason))
for(r in 1:nrow(hgrad)){
  if(hgrad$oldarea[r]==99) fl <- fleets$fleet[fleets$group=='comm']-1
  if(hgrad$oldarea[r]!=99) fl <- fleets$fleet[fleets$group=='comm' & fleets$newarea==hgrad$oldarea[r]]-1
  dat2[dat$fleet%in%fl & dat$step==(hgrad$tstep[r]-1),(startseason:endseason)==hgrad$season[r]]  <- which(hglist==hgrad$propfl[r])-1}
  
dat <- cbind(dat,dat2)
dat %<>% arrange(sex, age, fleet, step)
tmp <- list()
tmp <- c(tmp, "# Retain specification (This represents the proportion of LEGAL lobster retained - (1-high-graded due to low value))\n\n# Number Retain Patterns\n",length(hglist),"\n")
tmp <- c(tmp, "# Pattern\tType\tSex\tExtra\tPointer","\n")
for(i in 1:length(hglist)){ tmp <- c(tmp,"\t",paste(c((i-1),1,-1,0,i-1),collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Specifications for retention. Fleet = ", paste(fleets$group,collapse=' '),"\n# Sex\tAge\tFleet\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")

for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# retention - Proportion of LEGAL lobster retained - NOT high-graded\n",length(hglist),"\n")
for(i in 1:length(hglist)){ tmp <- c(tmp,paste(rep(hglist[i],length(lens)),collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Retention parameters\n# Lower\tUpper\tEstimate\tPhase\n")

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/RETAINSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
