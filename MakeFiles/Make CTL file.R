setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(magrittr)

tmp <- list()
tmp <- c(tmp, "# WRL fishery\n\n# weight-at-length (W=aL^b) (kg) \n")
biol <- readWorkbook(wb,sheet='lengthweight', startRow = 2)
female.wat <- round(biol$a[biol$sex=='female']*lens^biol$b[biol$sex=='female'] ,3)
male.wat <- round(biol$a[biol$sex=='male']*lens^biol$b[biol$sex=='male'] ,3)
tmp <- c(tmp, paste(male.wat,collapse = "\t"),"\n", paste(female.wat,collapse = "\t"),"\n")

tmp <- c(tmp, "\n# Maturity (by area) \n")
matdat <- readWorkbook(wb,sheet='maturity', startRow = 2)

code <- matrix(0, nrow=max(area$newarea), ncol=length(lens))
for(i in 1:nrow(code)){ code[i,] <- round(1/(1+exp((lens-matdat$a[matdat$cat=='mature'][i])/matdat$b[matdat$cat=='mature'][i])),3)}
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Egg production (by area)\n")
fecundity <- matdat$a[matdat$cat=='fecundity']*lens^matdat$b[matdat$cat=='fecundity']
code2 <- matrix(0, nrow=max(area$newarea), ncol=length(lens))
for(i in 1:nrow(code)){ code2[i,] <- round(1/(1+exp((lens-matdat$a[matdat$cat=='dspawn'][i])/matdat$b[matdat$cat=='dspawn'][i])),3)}
for(i in 1:nrow(code)){ code[i,] <- round(  (code[i,]+code2[i,])*fecundity,1)}
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Egg time step (Assume this is when to determine egg production?? time step 6 [6-1])\n",1,'\n')

tmp <- c(tmp, "\n# Fleet specification\n# Fleet Area\n\t")
fleetcode2 <- fleets %>% group_by(fleet) %>% summarise(area=unique(newarea), group=unique(group))
for(i in 1:nrow(fleetcode2)){ tmp <- c(tmp, paste(fleetcode2[i,c('fleet','area')]-1, collapse = "\t"),"\n\t")}

tmp <- c(tmp, "\n# Number of Zones\n", zones, "\n# Areas in each Zone\n", paste(as.numeric(sapply(zoneareas, length)), collapse = " "),'\n')
tmp <- c(tmp, "\n# The Zones\n")
for(z in 1:length(zoneareas)){
  tmp <- c(tmp, paste(zoneareas[[z]]-1, collapse = " "),'\n')}

tmp <- c(tmp, "\n# Discard mortality\n#Age\tFleet\ttstep\t",paste(startseason:endseason,collapse = "\t"),'\n')
dat <- expand.grid(age=(1:ages)-1, fleet=sort(unique(fleets$fleet))-1, step=sort(unique(times$tstep))-1)
dat %<>% arrange(age,fleet,step)
dat2 <- matrix(0.05, nrow=nrow(dat), ncol=length(startseason:endseason))
dat <- cbind(dat,dat2)
for(i in 1:nrow(dat)){ tmp <- c(tmp, paste(dat[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Recruitment_deviations\n",startseason, "\t\t\t# First year with estimated recruitment deviations\n")
tmp <- c(tmp, "#", endseason, "\t\t\t# last year with estimated recruitment deviations\n")
tmp <- c(tmp, 1, "\t\t\t# Phase for recruitment deviations\n")

tmp <- c(tmp, "\n# Spatial_deviations_in_recruitment\n",startseason, "\t\t\t# First year with estimates spatial recruitment deviations\n")
tmp <- c(tmp, "#", endseason, "\t\t\t# Last year with estimates spatial recruitment deviations\n")
tmp <- c(tmp, 3, "\t\t\t# Phase for spatial recruitment deviations\n\n")

tmp <- c(tmp, "# Prespecify_rec_devs :  dev # Year\n",1,"\t\t\t# 1 = rec_devs are to be pre-specified\n")
dat <- data.frame(rec_dev=0, year=startseason:(endseason+projectseason+burnin))
for(i in 1:nrow(dat)){ tmp <- c(tmp, dat[i,1],"\t#\t", dat[i,2],"\n")}

tmp <- c(tmp, "\n# Prespecify_spatial_rec_devs : dev # Year Area\n",1,"\t\t\t# 1 = Spatial_rec_devs are to be pre-specified\n")
dat <- data.frame(Spatial_rec_dev=0, expand.grid(year=startseason:(endseason+projectseason+burnin),area=(2:length(unique(area$newarea)))))
for(i in 1:nrow(dat)){ tmp <- c(tmp, dat[i,1],"\t#\t", dat[i,2],"\t", dat[i,3],"\n")}

tmp <- c(tmp, "\n# Weights on the data (simple)\n")
tmp <- c(tmp, 1.0, "\t\t# Weight on CPUE data\n")
tmp <- c(tmp, 1.0, "\t\t# Weight on catch-numbers data\n")
tmp <- c(tmp, 0.05, "\t\t# Weight on Length-frequency data\n")
tmp <- c(tmp, 1.0, "\t\t# Weight on larval data\n")
tmp <- c(tmp, 0.0, "\t\t# Weight on Tag1 data\n")
tmp <- c(tmp, 0.0, "\t\t# Weight on Tag2 data\n")
tmp <- c(tmp, 0.0, "\t\t# Weight on initial numbers\n")
tmp <- c(tmp, 0.0, "\t\t# Weight on initial penalty (InitOpt=3 or 5)\n")


tmp <- c(tmp, "\n# Weights by fleet (Type: 1=Cpue,2=Numbers,3=Length;4=Larva;) - Can tweak by Fleet and index.  See values to -1 for them to encompass all options.  e.g. time-step set to -1 covers all timesteps. our fleets\n# Type\tFleet\tTime-step\tsex\n")
tmp <- c(tmp, 0, "\t\t\t# set to >0 for yes\n")
dat <- data.frame(type=c(1,1,2,3,3,3,3,3,3),fleet=c(0,0,0,0,2,3,0,2,3),Tstep=c(1,2,1,-1,-1,-1,-1,-1,-1),sex=c(-1,-1,-1,0,0,0,1,1,1),weight=c(1,1,1,0.04,0.05,0.05,0.02,0.01,0.01),id=rep(c('cpue','numbers','lengths'),c(2,1,6)))
for(i in 1:nrow(dat)){ tmp <- c(tmp,paste("#",paste(dat[i,1:5],collapse = "\t")),"\t#\t", dat[i,6],"\n")}

tmp <- c(tmp, "\n# Basic parameters (lower, upper, estimate, phase)\n")

Mpar1 <- 12
tmp <- c(tmp, paste(c(10,25,Mpar1,1),collapse = "\t"), "\t\t\t# Rbar\n")
Mpar2 <- 0.12 
for(i in 1:length(unique(area$newarea))) {tmp <- c(tmp, paste(c(0.1,0.2,Mpar2,-1),collapse="\t"), paste("\t\t\t# M area",i,"\n"))}
Mpar3 <- rep(1,ages) 
for(i in 1:ages){ 
  tmp <- c(tmp, paste(c(0,10,Mpar3[i],-1),collapse="\t"), paste("\t\t\t# Multiplier for age ",i-1,"\n")) }
tmp <- c(tmp, paste(c(1,5,3,-1),collapse="\t"), "\t\t\t# WhiteScaleM\n")
tmp <- c(tmp, paste(c(0,1,0.025,-1),collapse="\t"), "\t\t\t# MVirgin (additive).  Only one value not area specific\n")
tmp <- c(tmp, paste(c(-10,-1,-2.5,-1),collapse="\t"), "\t\t\t# Minflection % of virgin biomass in that area. In log space. Only one value not area specific.\n")
tmp <- c(tmp, paste(c(0,1,0.3,2),collapse="\t"), "\t\t\t# RedScaleQ\n")
tmp <- c(tmp, paste(c(0.001,1,0.5,-1),collapse="\t"), "\t\t\t# SigmaR\n")

for(a in 1:length(unique(area$newarea))){ 
  tmp <- c(tmp, paste(c(10,25,15,-1),collapse="\t"), paste("\t\t\t# Rintro - initial estimated recruitment. Area",a,"\n"))
}
tmp <- c(tmp, paste(c(-6,1,-0.2,-1),collapse="\t"), "\t\t\t# Initial F - what F was occuring back at the start of the model\n")

tmp <- c(tmp, "\n# Initial_dev_option\n")
tmp <- c(tmp, 0, "\t\t\t\t# 0=convetional; 1=alternative; 2=Something; 3=Something else; 4=Yet another option; 5 as for 3 but with initial values for Rinitial by area\n")
tmp <- c(tmp, 1, "\t\t\t\t# Initial value options (0=default; 1=same for all)\n")
tmp <- c(tmp, "# Initial size parameters\n")
tmp <- c(tmp, paste(-100,100,0,1,collapse='\t'), "\t\t\t# 0=convetional; 1=alternative\n")

tmp <- c(tmp, "\n# Q parameters\n")
Qpar1 <- 1
nqs <- length(unique(fleets$fleet[fleets$group=='comm']))*length(unique(times$tstep))
for(i in 1:nqs){ tmp <- c(tmp, paste(c(-100,100,Qpar1,-1),collapse="\t"), "\t\t\t# Multiplier for Environmental index - Keep to 1\n")}

tmp <- c(tmp, "\n# Efficiency parameters (lower, upper, estimate, phase) \n")
ecpar1 <- effic[,c('lower','upper','est','Phase')]
nECvec <- max(fleets$effic.creep)
nECpar <- floor((endseason-startseason+1)/effic$temporal.cover)
tmp <- c(tmp, sum(nECpar),"\t# Number of Efficiency parameters \n")
for(nv in 1:nECvec){ 
  for(np in 1:nECpar[nv]){
    tmp <- c(tmp, paste(ecpar1[nv,],collapse="\t"), paste("\t\t\t# Efficiency creep par - Pointer", effic$pointer[nv]," One par every",effic$temporal.cover[nv],"years.\n")   )}}

tmp <- c(tmp, "\n# variance specification parameters (1=SSB;2=Recruitment;3legal biomass by area;4:exploitation rate by zone)\n")
tmp <- c(tmp, "# Number of variance specifications\n",4)
tmp <- c(tmp, "\n# variance components\n",paste(1:4,collapse = '\t'),"\n")

tmp <- c(tmp, "\n1 # use the pin file for specifying parameters (ADMB)")
tmp <- c(tmp, "\n1 # last function call (ADMB)\n")

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/CONTROL.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)



## MAke Starter file

tmp <- list()
tmp <- c(tmp, "DATA.DAT                              # General specifications file \n")
tmp <- c(tmp, "CONTROL.DAT                              # Control file\n")
tmp <- c(tmp, "SELEXSPEC.DAT                        # Specifications for selectivity \n")
tmp <- c(tmp, "RETAINSPEC.DAT                       # Specifications for retention \n")
tmp <- c(tmp, "RECRUITSPEC.DAT                      # Specifications for recruitment \n")
tmp <- c(tmp, "GROWTHSPEC.DAT                       # Specifications for growth \n")
tmp <- c(tmp, "MOVESPEC.DAT                         # Specifications in movement \n")
tmp <- c(tmp, "TAGFILE.TXT                          # Specifications for tags \n")
tmp <- c(tmp, "PROPORTIONS.TXT                      # Proportions \n")
tmp <- c(tmp, "PROJECTIONS.DAT                      # Projections file \n\n\n")
tmp <- c(tmp, "1                                    # Stop after this phase \n")

write.table(tmp, paste(floc,'/STARTER.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)



