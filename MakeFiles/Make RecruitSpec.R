library(dplyr)
library(magrittr)

tmp <- list()
tmp <- c(tmp, "# Recruitment specifications\n\n# Number of sex_area_allocation options\n",1,"\n\n")
tmp <- c(tmp, "# Allocate_yearxarea (form 0 is estimate sex split and then area split; form 1 is estimate sex*area split)\n")
tmp <- c(tmp, "# Type\tForm\tAlloc_sex_area\n",paste(c(0,0,0),collapse = '\t'),'\n\n')

dum <- data.frame(type=(1:length(unique(areas$recruitarea)))-1, sex1=0, sex2=0, Extra1=-1, Extra2=-1)
dum$sex1 <- -seq(1,(2*nrow(dum))-1,2)
dum$sex2 <- dum$sex1-1

tmp <- c(tmp, "# Number of length_allocation options (size compoisition of recruits by area Not sex)\n",length(unique(areas$recruitarea)),"\n\n")
tmp <- c(tmp, "# Allocate_length\n# Type\tAllocate_sex1\tAllocate_sex2\tExtra1\tExtra2\t: 2nd column sex 1 and 3rd col for sex 2 (Offset for PreSpecRecFrac = -1*col.value-1)\n")
for(r in 1:nrow(dum)) {
  tmp <- c(tmp, "\t", paste(dum[r,],collapse = "\t"),paste("\t# Areas ", paste(unique(areas$aname[areas$recruitarea==r]), collapse=','), "\n"))
           }
tmp <- c(tmp, "\n# Recruit by year and time-step\n# TimeStep\t", paste(startseason:(endseason+projectseason+10),collapse = "\t"),"\n")
dat <- matrix(-1, nrow=length(unique(times$tstep)), ncol=length(startseason:(endseason+projectseason+10)))
dat[unique(times$tstep[times$recruit==1]),] <- 0
dat <- cbind(data.frame(ts=sort(unique(times$tstep))-1),dat)
for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

rec <- areas %>% group_by(newarea) %>% summarise(area=mean(recruitarea)-1)
tmp <- c(tmp, "\n# Recruit by area\n",paste(rec$area,collapse = "\t"),"\n")

rec <- readWorkbook(wb,sheet='Recruitment', startRow = 2)
usepar <- tail(rec$Use.Parameters,1)

tmp <- c(tmp, "\n# Recruitment fractions use prespecified (0) or calculate from parameters (1). If pre-specified need to code out parameters #\n")  
tmp <- c(tmp, paste(usepar, "\t# Use parameters 0 = No, 1 = Yes.\n"))

tmp <- c(tmp, "\n# Number of pre-specified recruitment functions\n",max(areas$recruitarea)*length(sexs),"\n")
tmp <- c(tmp, "# Prespecified recruitment fractions (based on mean CL + SD of lobster at the start of the year age 3)\n")

recdist <- rec %>% dplyr::select(starts_with('Prespecified'))
recdist <- recdist[!is.na(recdist[,1]),]
for(i in 1:nrow(recdist)){ tmp <- c(tmp,paste(recdist[i,],collapse = "\t"),"\n")  }

tmp <- c(tmp, "\n#  Recuitment1 parameters\n#Lower\tUpper\tEstimate\tPhase::\t Number of recruitment fraction parameters must match Number of pre-specified recruitment functions above.\n")

rec %<>% dplyr::select(Use.Parameters,lower,upper,est,Phase,description) %>% mutate(Use.Parameters=ifelse(Use.Parameters==1,'','#'), description =paste('#', description ))
for(a in 1:nrow(rec)){ tmp <- c(tmp, paste(rec[a,],collapse='\t'), "\n")  }

# Bias ramp
tmp <- c(tmp, "\n# Bias ramp - insert description here\n",paste(c(startseason, startseason+2, endseason-5, endseason-1),collapse = "\t"),"\t#description\tdescription\tdescription\tdescription\n")

puerpar <- readWorkbook(wb,sheet='PuerulusPar', startRow = 2) %>% mutate(description=paste('#', description))

tmp <- c(tmp, "\n# Puerulus Power for puerulus to recruit relationship\n",nrow(puerpar),"\n#LB\tUP\tEstimate\tPhase\n")
for(a in 1:nrow(puerpar)){ tmp <- c(tmp, paste(puerpar[a,],collapse = "\t"), "\n") }

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/RECRUITSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)









