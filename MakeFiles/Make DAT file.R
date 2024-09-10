setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(magrittr)
library(openxlsx)
library(emmeans)

lens <- seq(dynamics$value[dynamics$object=='lblwr'],dynamics$value[dynamics$object=='lbupr'],dynamics$value[dynamics$object=='lbgap'])
lensPlus1 <- c(lens, (lens[2]-lens[1])+tail(lens,1))

## Make data input of new RL model
tmp <- list()
tmp <- c(tmp, "# Lobster Data Set", "\n", "# First year of the assessment")

#get latest catch and effort data
dat <- readWorkbook(wb,sheet='CatchnEffort', startRow = 2)
dat %<>% mutate(newarea=areas$newarea[match(oldarea,areas$oldarea)], tstep=times$tstep[match(paste(month,hm),paste(times$month,times$hm))]) %>% mutate(finyr=ifelse(month>=7,year,year-1))
dat %>% filter(is.na(tstep))
cpuedat <- dat  ## keep all info before summing by area and tstep.
dat %<>% filter(finyr%in%startseason:endseason) 
tmp <- c(tmp,"\n", startseason, "\n# Last year of the assessment\n", endseason,"\n# Maximum projection years\n", projectseason+5)

bin <- areas %>% group_by(newarea) %>% summarise(av=floor(mean(burn_in)))

tmp <- c(tmp, "\n# Burn-in\n", paste(bin$av,collapse=' '), 
         "\n# Time steps per year\n", max(times$tstep), 
         "\n# Number of areas of data included in the file\n", max(areas$newarea),
         "\n# Number of sexes\n",length(unique(sexs)),
         "\n# Number of ages\n", ages,
         "\n# Number of fleets\n",max(fleets$fleet))

## Lower Length bins
femaleLB <- maleLB <- lensPlus1
tmp <- c(tmp, "\n# Number of size-classes (males then females)\n", length(maleLB)-1, " ", length(femaleLB)-1)
# Timesteps
tstmp <- times %>% group_by(tstep) %>% mutate(prop=length(month)/nrow(times)) %>% summarise(prop=mean(prop))
tmp <- c(tmp, "\n# The Time steps\n# Starts in July\n", paste(round(tstmp$prop,5), collapse = "\t"))

tmp <- c(tmp, "\n# Loop counter for initial conditions\n", 10,"\n# Years over which to tune (one per area)\n", dynamics$value[dynamics$object=='tune_years'])
# Length Bins cont.
tmp <- c(tmp, "\n# Lower Length Bins (one more than number of size-classes)\n", paste(maleLB, collapse = "\t"), "\n", paste(femaleLB, collapse = "\t"))																										

## Catch data
dat %<>% mutate(icpue=catch/effort) %>% filter(catch>0) %>% group_by(finyr, Newtstep=tstep, Newarea=newarea) %>% summarise(Catch=round(sum(as.numeric(catch)),6), Effort=round(sum(as.numeric(effort)),6), cpue=round(mean(icpue),3), sdcpue=round(sd(icpue),3)) %>% mutate(group='comm') %>% mutate(fleet = fleets$fleet[match(paste(group,Newarea), paste(fleets$group,fleets$newarea))]) %>% mutate(Modtstep=Newtstep, Modfleet=fleet) 

## Recreational Catch All in timestep 7 and split 75:25 south : north
rec <- readWorkbook(wb,sheet='RecCatch', startRow = 2)
rec %<>% mutate(year=as.numeric(substr(Season,1,4)), month=12, area=min(areas$newarea[areas$oldzone=='C']), Adjusted.catch=Adjusted.catch) %>% dplyr::select(finyr=year, month, area, Adjusted.catch) %>% filter(finyr%in%startseason:endseason)
recBZone <- rec %>% mutate(Adjusted.catch=Adjusted.catch*0.25) %>% mutate(area=min(areas$newarea[areas$oldzone=='B']))
recCZone <- rec %>% mutate(Adjusted.catch=Adjusted.catch*0.75)
rec <- rbind(recBZone,recCZone) %>% filter(!is.na(Adjusted.catch)) %>% mutate(Step=times$tstep[match(paste(month,0),paste(times$month,times$hm))])
rfc <- fleets[fleets$group=='rec',]
rec %<>% mutate(fleet=rfc$fleet[match(area, rfc$newarea)]) %>% arrange(finyr,fleet,Step) 
tmp <- c(tmp, "\n# Catch data (kg) - Number of observations\n", nrow(dat)+nrow(rec), "\n#Year\tstep\tfleet\tcatch\n")

for(i in 1:nrow(dat)){
  tmp <- c(tmp, paste(dat[i,c('finyr', 'Modtstep', 'Modfleet', 'Catch')], collapse = "\t"),"\n")}
for(i in 1:nrow(rec)){
  tmp <- c(tmp, paste(rec[i,c('finyr', 'Step', 'fleet', 'Adjusted.catch')], collapse = "\t"),"\n")}

##  Catch Rate Indices / CPUE
lbdat <- readWorkbook(wb,sheet='LBCPUE', startRow = 2) %>% mutate(source='lb')
cdrdat <- readWorkbook(wb,sheet='CDRCPUE', startRow = 2) %>% mutate(source='cdr')
cpue <- rbind(lbdat,cdrdat)
cpue %<>% mutate(newarea=areas$newarea[match(oldarea,areas$oldarea)], tstep=times$tstep[match(paste(month,hmonth),paste(times$month,times$hm))]) %>% mutate(finyr=ifelse(month>=7,year,year-1), U=catch/effort, yrts=paste(finyr,tstep), farea=as.factor(newarea)) %>% filter(catch>0) %>% mutate(lfb1=substr(LFB,1,1), lfb2=substr(LFB,2,nchar(LFB)), lfb2=numchar(lfb2), lfb=paste(lfb1,lfb2,sep=''))
#head(cpue)
#summary(cpue$U)

cpueout <- data.frame(newarea=NA, finyear=NA, tstep=NA, est=NA, se=NA)
for(a in sort(unique(cpue$newarea))){
  tcpue <- cpue %>% filter(newarea==a) %>% group_by(lfb) %>% mutate(num=length(U)) %>% filter(num>50)
  #length(table(tcpue$lfb))
  lm1 <- lm(log(U)~yrts+LFB, data=tcpue)
  out <- summary(emmeans(lm1, ~yrts, type='response',rg.limit = 200000))
  out %<>% mutate(finyear=as.numeric(substr(yrts,1,4)), tstep=as.numeric(do.call('rbind',strsplit(as.character(yrts), ' '))[,2]), newarea=a)  %>% dplyr::select(newarea, finyear, tstep, est=response, se=SE) 
  cpueout <- rbind(cpueout,out)}
  
cpueout %<>% mutate(group='comm') %>% mutate(fleet = fleets$fleet[match(paste(group,newarea), paste(fleets$group,fleets$newarea))],Qid = fleets$commonq[match(paste(group,newarea), paste(fleets$group,fleets$newarea))]) %>% mutate(Modtstep=tstep, Modfleet=fleet, cpue=est, finyr=finyear, Newtstep=tstep, Newarea=newarea) 

comdat <- cpueout %>% filter(!is.na(cpue) & cpue>0) %>% mutate(sex=0, cpuecv=se/cpue) %>% group_by(finyr) %>% filter(finyr%in%startseason:endseason) %>% dplyr::select(Qid,Modfleet,Newarea,sex,finyr,Modtstep,cpue,cpuecv) %>%  as.data.frame() %>% arrange(Modfleet,sex,Modtstep)  

if(dynamics$value[dynamics$object=='change_q_bymonth']==1) {
  comdat %<>% mutate(Qid=Qid+(Modtstep/10)) %>% as.data.frame() }
comdat %<>% mutate(Qid=as.numeric(as.factor(Qid))-1) %>% arrange(Modfleet,sex,Modtstep)

qidarea <- comdat %>% dplyr::select(Qid,Modfleet,Modtstep) %>%  group_by(Qid,Modfleet,Modtstep) %>% summarise(n=length(Qid)) %>% as.data.frame()

comdat %<>% dplyr::select(Qid,Modfleet,sex,finyr,Modtstep,cpue,cpuecv) %>% mutate(Uid=as.numeric(factor(as.character(paste(Modfleet,Modtstep)),levels=unique(as.character(paste(Modfleet,Modtstep)))))-1)

Qs <- comdat %>% group_by(Uid,Modfleet,Modtstep) %>% summarise(Qid=mean(Qid)) %>% mutate(area=fleets$newarea[match(Modfleet,fleets$fleet)]) %>% as.data.frame()

comdat %<>% dplyr::select(Qid=Uid,Modfleet,sex,finyr,Modtstep,cpue,cpuecv) %>% arrange(Qid,Modfleet,sex,finyr,Modtstep)  

## Get IBSS & ISS indices
ibss <- readWorkbook(wb,sheet='IBSS')
ibss %<>% mutate(group='ibss',month=11, cv=round(SE/response,3), mn=round(response,3),fleet=fleets$fleet[fleets$group=='ibss'][match(Ibssarea, fleets$newarea[fleets$group=='ibss'])],Qid = fleets$commonq[match(paste(group,Ibssarea), paste(fleets$group,fleets$newarea))]) %>% mutate(tstep=times$tstep[match(paste(month,0),paste(times$month,times$hm))], Qid=as.numeric(as.factor(Qid))) %>% filter(year%in%startseason:endseason) %>% group_by(Qid, fleet,sex,year,tstep) %>% summarise(mn=mean(mn), cv=mean(cv)) %>% dplyr::select(Qid, fleet, sex, year, tstep, mn, cv) %>% ungroup() %>%  mutate(Qid=(Qid-min(Qid))+max(comdat$Qid)+1)

#ibss$Qid

ibssqs <- ibss %>% group_by(Qid, fleet, tstep) %>% summarise(num=length(mn)) %>% mutate(area=fleets$newarea[match(fleet,fleets$fleet)],qs=Qs$Uid[match(paste(area,tstep),paste(Qs$area,Qs$Modtstep))])

iss <- readWorkbook(wb,sheet='ISS')
iss %<>% mutate(group='iss',month=3, cv=round(SE/response,3), mn=round(response,3), fleet=fleets$fleet[fleets$group=='iss'][match(Issarea, fleets$newarea[fleets$group=='iss'])],Qid = fleets$commonq[match(paste(group,Issarea), paste(fleets$group,fleets$newarea))]) %>% mutate(tstep=times$tstep[match(paste(month,0),paste(times$month,times$hm))], Qid=as.numeric(as.factor(Qid))) %>% group_by(Qid, fleet,sex,finyr,tstep) %>% summarise(mn=mean(mn), cv=mean(cv)) %>% filter(finyr%in%startseason:endseason) %>% dplyr::select(Qid, fleet, sex, year=finyr, tstep, mn, cv) %>% ungroup() %>% mutate(Qid=(Qid-min(Qid))+max(ibss$Qid+1))

issqs <- iss %>% group_by(Qid, fleet, tstep) %>% summarise(num=length(mn)) %>% mutate(area=fleets$newarea[match(fleet,fleets$fleet)],qs=Qs$Uid[match(paste(area,tstep),paste(Qs$area,Qs$Modtstep))])

cpuenumbers <- (c(Qs$Qid,unique(ibssqs$qs),unique(issqs$qs)))

## Comm=1-16, IBSSa2=17, IBSSa4=18, IBSSa5=19, IBSSa6=20, IBSSa8=21, ISSa1=5, ISSa3=6, ISSa5=7, ISSa7=8
tmp <- c(tmp, "\n# Index data (# Comm=0-15 - 0-7 = reds, 8-15 = whites, IBSS=16+ is Ibss and Iss: unknown sex = 0)\n#Number of cpue datasets\n", length(cpuenumbers))
tmp <- c(tmp, "\n# Type of index (1=weight;2=numbers)\n", paste(c(rep(1,length(unique(comdat$Qid))),rep(2,length(c(unique(ibss$Qid),unique(iss$Qid))))), collapse=" "))
tmp <- c(tmp, "\n# Treatment of sigma (unique value represents a unique SS for the series)\n", paste((1:length(cpuenumbers))-1, collapse=" "))  ## Fixsigma
tmp <- c(tmp, "\n# Treatment of q (a value represents a unique q for that series)\n", paste(cpuenumbers, collapse=" "))
#tmp <- c(tmp, "\n# Environmental Index / Fishing Efficiency (value points to index, 0 = no index)\n", paste(c(unname(tapply(qidarea$Newarea, list(qidarea$Qid), mean)),rep(0,100))[1:length(cpuenumbers)], collapse=" "))    
tmp <- c(tmp, "\n# Environmental Index (value points to index, 0 = no index)\n", paste(c(1:length(unique(comdat$Qid)),rep(0,100))[1:length(cpuenumbers)], collapse=" "))  
## Efficiency creep
tmp1 <- comdat %>% group_by(Qid) %>% summarise(mn=mean(Modfleet))
tmp2 <- ibss %>% group_by(Qid) %>% summarise(mn=mean(fleet))
tmp3 <- iss %>% group_by(Qid) %>% summarise(mn=mean(fleet))
tmp4 <- rbind(tmp1,tmp2,tmp3)
effcreep <- fleets$effic.creep[match(tmp4$mn,fleets$fleet)]
tmp <- c(tmp, "\n# Efficiency creep (value points to index, 0 = no index)\n", paste(effcreep, collapse=" "))
tmp <- c(tmp, "\n# Efficiency creep year lag (each par compounds for this many years until next par starts\n", paste(effic$temporal.cover, collapse=" "))
tmp <- c(tmp, "\n# Minimum sigma\n", 0.05)
#Size of cpue data																																
tmp <- c(tmp,"\n# The cpue data\n", (nrow(comdat)+nrow(ibss)+nrow(iss)))	

#CpueInd	#Fleet	#Sex	#Year	#Step	#Index	#CV
tmp <- c(tmp, "\n#CpueInd Fleet\tSex\tYear\tStep\tIndex\tCV\n")
for(i in 1:nrow(comdat)){
  tmp <- c(tmp, paste(comdat[i,], collapse = "\t"),"\n")}
for(i in 1:nrow(ibss)){
  tmp <- c(tmp, paste(ibss[i,], collapse = "\t"),"\n")}
if(nrow(iss)>0) {for(i in 1:nrow(iss)){
  tmp <- c(tmp, paste(iss[i,], collapse = "\t"),"\n")}}

tmp <- c(tmp,"\n# Numbers data (recreational fishery only)\n# Number of catch data sets\n",0)	
tmp <- c(tmp,"\n# Treatment of catch in numbers series\n#", 0, "\n# Minimum sigma\n",0.05)

## Build the numbers data
tmp <- c(tmp,'\n# The numbers data\n',0)	
tmp <- c(tmp,'\n#Group  Fleet  Year  Step  Catch  CV\n')	
# mutate(fleet=fleet-1, sex=sex-1, tstep=tstep-1) %>% 

## Get the length data
#Commercial  file:  M:\Fisheries Research\Invertebrates Fisheries\Rock Lobster\Length Monitoring\Analysis/mk.lf.input(51.121).r
#len <- read.csv(direct('Data/Com monitoring/NewModel_Lenfreq.csv'))
len <- readWorkbook(wb,sheet='LengthFreqComm', startRow = 2)
#toXL(len[1,])
len %<>% mutate(fleet=fleets$fleet[fleets$group=='comm_monitor'][match(area, fleets$newarea[fleets$group=='comm_monitor'])], fleet=ifelse(is.na(fleet),min(fleets$fleet[fleets$group=='comm_monitor']),fleet)) %>% filter(season%in%startseason:endseason & !is.na(season)) %>% mutate(sex=sex, tstep=tstep) %>% dplyr::select(fleet,sex,season,tstep,numb,starts_with('X')) %>% group_by(fleet, sex, season, tstep) %>% summarise_at(c('numb', colnames(len)[grepl('X',colnames(len))]),sum)

#IBSS and ISS   file:///C:/Users/snd/Rock Lobster/Data/IBSS/mk.lf.input(51.121).r

ibss <- readWorkbook(wb,sheet='LengthFreqIBSS', startRow = 2)%>% mutate(area=NA, season=ifelse(month>=7,year,year-1)) %>%  filter(!grepl('30deg', loc, ignore.case = T)) %>% filter(!grepl('Marshalls', loc, ignore.case = T)) %>% mutate(loc=ifelse(grepl('deep',loc,ignore.case=T), 'deep',loc))
for(i in 1:nrow(areas)){  ## Add newarea based on location name
  Name <- areas$ibssname[i]  
  if(!is.na(Name)){
    Names <- do.call('rbind', strsplit(Name, ','))
    for(n in 1:length(Names)){
      tNames <- Names[n]
      tNames <- gsub(" ","",tNames)
      ibss$area[grepl(tNames, ibss$loc, ignore.case = T)] <- areas$newarea[i] 
    } }  }

#table(ibss$area, ibss$loc)

fleetibss <- fleets[fleets$group=='ibss',]
ibss %<>% mutate(fleet=fleetibss$fleet[match(area, fleetibss$newarea)]) %>% filter(season%in%startseason:endseason & !is.na(season)) %>% mutate(sex=ifelse(sex=='F',1,2), tstep=times$tstep[match(paste(month,hm),paste(times$month,times$hm))]) %>% dplyr::select(fleet,sex,season,tstep,numb,starts_with('X')) %>% group_by(fleet, sex, season, tstep) %>% summarise_at(c('numb', colnames(len)[grepl('X',colnames(len))]),sum)

iss <- readWorkbook(wb,sheet='LengthFreqISS', startRow = 2) %>% mutate(area=NA, season=ifelse(month>=7,year,year-1))
for(i in 1:nrow(areas)){
  Name <- areas$issname[i]  
  if(!is.na(Name)){
    Names <- do.call('rbind', strsplit(Name, ','))
    for(n in 1:length(Names)){
      tNames <- Names[n]
      tNames <- gsub(" ","",tNames)
      iss$area[grepl(tNames, iss$loc, ignore.case = T)] <- areas$newarea[i] 
    } }  }

#tapply(iss$loc, list(iss$loc,iss$area), length)
fleetiss <- fleets[fleets$group=='iss',]
iss %<>% mutate(fleet=fleetiss$fleet[match(area, fleetiss$newarea)]) %>% filter(season%in%startseason:endseason & !is.na(season)) %>% mutate(sex=ifelse(sex=='F',1,2), tstep=times$tstep[match(paste(month,hm),paste(times$month,times$hm))]) %>% dplyr::select(fleet,sex,season,tstep,numb,starts_with('X')) %>% group_by(fleet, sex, season, tstep) %>% summarise_at(c('numb', colnames(len)[grepl('X',colnames(len))]),sum)

tmp <- c(tmp,"\n# Length compostion\n", (nrow(len)+nrow(ibss)), '\n#Fleet\tSex\tSEASON\ttstep\tInd\t',paste(maleLB, collapse = "\t"),'\n')	
for(i in 1:nrow(len)){
  tmp <- c(tmp, paste(len[i,], collapse = "\t"),"\n")}
for(i in 1:nrow(ibss)){
  tmp <- c(tmp, paste(ibss[i,], collapse = "\t"),"\n")}
if(nrow(iss)>0) {for(i in 1:nrow(iss)){
  tmp <- c(tmp, paste(iss[i,], collapse = "\t"),"\n")}
}

#Puerulus index																																	
tmp <- c(tmp,"\n# Larval index (puerulus)\n# Likelihood for larval (puerulus) data (0=lognormal, else normal\n",0)	
tmp <- c(tmp,"\n# Delay from puerulus to entering the model (years)\n",3)	
puer <- readWorkbook(wb,sheet='Puerulus')  %>%  mutate(area=area,mn=round(mn,1), sd=round(cv*mn,2))  %>% dplyr::select(area,season, mn, sd) %>% filter(season %in% startseason:endseason, area %in% areas$newarea)
tmp <- c(tmp,"\n# Number of data points (puerulus samples)\n",nrow(puer),'\n# Area\tYear\tIndex\tSD\n')	
for(i in 1:nrow(puer)){ tmp <- c(tmp, paste(puer[i,], collapse = "\t"),"\n")}

#	Environmental	Data   / Fishing efficiency
#fe <- read.csv(direct('Lobster/Minor stuff/Efficiency/2020/Fish_eff_estimates.csv'))
fe <- readWorkbook(wb,sheet='CommEfficiency')  
fe2 <- expand.grid(year=startseason:endseason, area=fe$newarea, creep=1) %>% mutate(pwr=year-startseason)
for(a in unique(fe$newarea)){
  fe2$creep[fe2$area==a] <- fe2$creep[fe2$area==a] * fe$annualE[fe$newarea==a]^fe2$pwr[fe2$area==a]}

dat <- expand.grid(year=startseason:endseason, Qid=sort(unique(comdat$Qid)), step=sort(unique(times$tstep)))
dat %<>% arrange(year,Qid,step)
dat %<>% mutate(area=qidarea$Newarea[match(Qid,qidarea$Qid)])
dat$effcreep <- round(log(fe2$creep[match(paste(dat$year, dat$area), paste(fe2$year,fe2$area))]),5)
dat$effcreep[is.na(dat$effcreep)] <- 0

#nseries <- dat %>% group_by(fleet) %>% summarise(num=length(effcreep)) 
tmp <- c(tmp,"\n# Environmental Data - ln(Efficiency creep)\n# Number of environmental series (commercial efficiency creep for each area)\n",length(unique(comdat$Qid)))	
tmp <- c(tmp,"\n# Years of data per series\n",paste(unname(table(dat$Qid)), collapse = "\t"),'\n')	
tmp <- c(tmp,"# Year\tTstep\tln(creep)\tSeries\n")	
for(s in 1:length(unique(dat$Qid))){
  tdat <- dat[dat$Qid==sort(unique(dat$Qid))[s],]
  for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,c('year', 'step', 'effcreep')], collapse = "\t"), "\t#\t",s,"\n")}}

#write.table(tmp,  paste(getwd(),'/Simon/Northa.dat',sep=''), sep="", row.names = F, col.names = F, quote=F)

# # Movement Data
# move <- read.csv(direct('Lobster/Minor stuff/Growth/All_move.growth.data.csv'))
# move %<>% filter(!is.na(Lsex), Lsex!='U' , LCl>=51, LCl<=200) %>% filter(libm>2)%>% filter(Cloc!=0)%>% filter(Ccl >= 51 , Ccl <=200)
# move %<>% mutate(sex=ifelse(Lsex=='M',1,2), RelLB=as.numeric(cut(LCl, c(lens[1:(length(lens)-1)]-0.1,180))), RecLB=as.numeric(cut(Ccl, c(lens[1:(length(lens)-1)]-0.1,180)))) %>%
#   mutate(group=ifelse(source=='comm', 'comm', 'ibss'), fleet=fleetcode$fleet[match(paste(Cloc,group),paste(fleetcode$area,fleetcode$group))]) %>%
#   select(Lloc,sex,Lseason,Ltstep,RelLB,fleet,Cseason,Ctstep,RecLB) %>% filter(!is.na(fleet)) %>%
#   group_by(Lloc,sex,Lseason,Ltstep,RelLB,fleet,Cseason,Ctstep,RecLB) %>% summarise(Num=length(sex))
# unique(move$fleet)
# tmp <- c(tmp,"\n# Movement and growth data from tagging\n")
# tmp <- c(tmp,"# Length of Movement data\n", nrow(move), '\n#RelArea   Sex   RelSeason   RelTstep   RelLB   RecFleet   RecSeason   RecTstep   RecLB   NumObs\n')
# for(i in 1:nrow(move)){ tmp <- c(tmp, paste(move[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/DATA.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
getwd()
