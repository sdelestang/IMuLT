
egap <- readWorkbook(wb,sheet='Escapegaps', startRow = 2)
egappar <- egap %>% filter(!is.na(gap))
#iswhite <- readWorkbook(wb,sheet='migrate', startRow = 2)
iswhite <- readWorkbook(wb,sheet='IsWhite', startRow = 2)

tmp <- list()
tmp <- c(tmp, "# Selex specification\n# Number Selex Patterns\n",nrow(egappar)*2)
tmp <- c(tmp, "\n# Pattern Type Sex Extra Pointer\n")
for(i in 1:(nrow(egappar)*2)){ tmp <- c(tmp,paste(c((i-1),1,-1,0,i-1),collapse = "\t"),"\n")}

tmp <- c(tmp, "# Specifications for selectivity (Selectivity of the pots [escape gaps, females then males], 0 = None, 1 = 54 mm, 2 = 55 mm) \n", "# Sex\tAge\tFleet\tStep: ",paste(startseason:endseason,collapse = "\t"),"\n")
leg <- expand.grid(sex=sexs, age=0:(ages-1), fleet=unique(fleets$fleet-1), step=sort(unique(times$tstep))-1)
leg %<>% arrange(sex, age, fleet, step)
id <- paste(leg[,1],leg[,2],leg[,3],leg[,4], sep="-")
code <- matrix(0, nrow=nrow(leg), ncol=length(startseason:endseason), dimnames = list(pat=id,year=paste('Y',startseason:endseason,sep='')))
code[leg$fleet %in% (fleets$fleet[fleets$group %in% c('comm_old','comm','comm_monitor','rec')]-1) , startseason:endseason%in%egap$season[egap$gaplink==1]] <- 1
code[leg$fleet %in% (fleets$fleet[fleets$group %in% c('comm_old','comm','comm_monitor')]-1) , startseason:endseason%in%egap$season[egap$gaplink==2]] <- 2
code[leg$sex==1,] <- code[leg$sex==1,] + 3  ## adds offset for males to use a different selectivity
code <- cbind(leg,code)
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "# selectivity\n", nrow(egappar)*2,"\n")
code <- matrix(0, nrow=(nrow(egappar)*2), ncol=length(lens), dimnames = list(pat=paste(rep(c('f','m'),each=nrow(egappar)),egappar$gap,sep=''),year=paste('Y',lens,sep='')))
# This adjust for pre-determined q of lobster based on size
for(i in 1:nrow(egappar)){                                     # 1.097297e+02  1.364472e+02  5.571002e+00  5.224848e+00
  qselect <- 1/(1+exp(((lens)-136)/5.22))  ## Females
  code[i,] <- round(1/(1+exp(((lens+1)-egappar$inflect[i])/egappar$slope[i])) * qselect,4)
  qselect <- 1/(1+exp(((lens)-110)/5.57)) ## Males   
  code[(i+nrow(egappar)),] <- round(1/(1+exp(((lens)-egappar$inflect[i])/egappar$slope[i])) * qselect,4)
   }
code[,lens<75] <- 0                          ## Remove all below size at legal as selectivity is dodgy down here
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

## Move onto Legal patterns
gauge <- readWorkbook(wb,sheet='LegalID', startRow = 2) %>% mutate(hash='#', type=1, Extra=0, pointer=pos-1, pos=pointer) %>% dplyr::select(pos, type, Extra, pointer, hash, id)
tmp <- c(tmp, "\n# Number Legal patterns (What is legal and can be retained [e.g. setose, berried, maxsize, minsize] or in a survey all can be caught)\n", nrow(gauge),"\n# Pattern\tType\tExtra\tPointer\t#  Decade Sex Zone Depth Tstep min guage, maxguage biocontrol\n")
#ids2 <- apply(as.matrix(ids),2,function(x) as.character(x))
for(i in 1:nrow(gauge)){ tmp <- c(tmp, paste(gauge[i,], collapse = "\t"), "\n")    }  

gauge2 <- readWorkbook(wb,sheet='LegalPattern', startRow = 2) %>% dplyr::select(starts_with('lb')) 

tmp <- c(tmp, "\n# legal patterns\n", nrow(gauge2),"\n")
for(i in 1:nrow(gauge2)){ tmp <- c(tmp, paste(round(as.numeric(gauge2[i,]),4), collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Specifications for Fleet legal assignment\n#Sex\tAge\tFleet\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")

gauge3 <- readWorkbook(wb,sheet='LegalList', startRow = 2) 
agemat <- dynamics$value[dynamics$object=='agemat']
gauge3$point <- gauge$pointer[match(gauge3$id, gauge$id)]
gauge3 %<>% mutate(Sex=ifelse(sex=='F',0,1))
fleets %<>% mutate(zone=areas$newzone[match(newarea, areas$newarea)])

leg <- expand.grid(Sex=sexs, Age=0:(ages-1),  Fleet=sort(unique(fleets$fleet))-1,Step=sort(unique(times$tstep))-1)
id <- paste(leg[,1],leg[,2],leg[,3],leg[,4], sep="-")
code <- matrix(0, nrow=nrow(leg), ncol=length(startseason:endseason), dimnames = list(pat=id,year=paste('Y',startseason:endseason,sep='')))
for(r in 1:nrow(code)){
  if(fleets$Is.Survey[fleets$fleet==(leg$Fleet[r]+1)]==1) {
  code[r,] <- gauge$pos[gauge$id=='Survey']
  } else {
    for(c in 1:ncol(code)){
      Sex <- leg$Sex[r] ; Sex <- ifelse(leg$Age[r]<(agemat-1), 1, Sex)  ## This is for WRL to set legal the same as males when age < maturity
      if(Sex==0) {
        Season <- (startseason:endseason)[c]
        Tstep <- leg$Step[r]+1
        Zone <- which(LETTERS==fleets$zone[fleets$fleet==(leg$Fleet[r]+1)])
        Depth <- areas$depth[areas$newarea==fleets$newarea[fleets$fleet==(leg$Fleet[r]+1)]][1]
        point <- gauge3 %>% filter(sex=='F', sea==Season, ts==Tstep, Depth==depth, zone==as.character(Zone)) %>% dplyr::select(point) %>% as.numeric()
      } else {
        Season <- (startseason:endseason)[c]
        Tstep <- leg$Step[r]+1
        point <- gauge3 %>% filter(sex=='M', sea==Season, ts==Tstep, depth=='X', zone=='X') %>% dplyr::select(point) %>% as.numeric()
        }
   code[r,c] <- point   
    }}}

code <- cbind(leg, code)
code %<>% arrange(Sex, Age, Fleet, Step)
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}


tmp <- c(tmp, "\n# Specifications for legal biomass\n#Sex\tAge\tArea\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")
leg <- expand.grid(Sex=sexs, Age=(1:ages)-1,  Area=sort(unique(areas$newarea))-1,Step=sort(unique(times$tstep))-1)
id <- paste(leg[,1],leg[,2],leg[,3],leg[,4], sep="-")
code <- matrix(0, nrow=nrow(leg), ncol=length(startseason:endseason), dimnames = list(pat=id,year=paste('Y',startseason:endseason,sep='')))
for(r in 1:nrow(code)){
  id1 <- leg[r,]
  (s1 <- ifelse(id1$Sex==1,'M','F'))  # pull out sex 
  s1 <- ifelse(id1$Age<(agemat-1), 'M', s1)  ## This is for WRL to set legal the same as males when age < maturity
  if(s1=='F'){
  (z1 <- areas$newzone[areas$newarea==(id1$Area+1)][1]); z1 <- which(LETTERS==z1)
  (a1 <- id1$Area)#trunc((id1$Area+1)/2)*2)
  (t1 <- id1$Step)
  tgauge <- gauge3 %>% filter(sex==s1, zone==z1, ts==(t1+1)) 
  point <- tgauge$point[match(startseason:endseason, tgauge$sea)]
  code[r,] <- point}
  if(s1=='M'){
    (z1 <- 'X')
    (a1 <- id1$Area)#trunc((id1$Area+1)/2)*2)
    (t1 <- id1$Step)
    tgauge <- gauge3 %>% filter(sex==s1, zone==z1, ts==(t1+1)) 
    point <- tgauge$point[match(startseason:endseason, tgauge$sea)]
    code[r,] <- point
  }
  }

code <- cbind(leg, code)
code %<>% arrange(Sex, Age, Area, Step)
for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Reference selectivity pattern (This is to set a constant Legal definition)\n")
## Set Base LegalBiomass to Legal definition of a male in 1992 which is a min CL of 76 mm
tmp <- c(tmp, paste(gauge2[gauge3$point[gauge3$sea==1992 & gauge3$sex=='M' & gauge3$min==76][1]+1,], collapse = "\t"))

tmp <- c(tmp, "\n\n# IsRed specifications - assignment of unique life stage quality\n")
dat <- expand.grid(sex=sexs,age=(1:ages)-1, area=sort(unique(areas$newarea))-1, step=sort(unique(times$tstep))-1, state=1)
if(nrow(iswhite)>0) {
  for(i in 1:nrow(iswhite)){
   dat$state[dat$age==(iswhite$Age[i]-1) & dat$area==(iswhite$Area[i]-1) & dat$step==(iswhite$Tstep[i]-1)] <- 0
 }}
dat %<>% arrange(sex, age, area, step)
dat %>%  filter(state==0)
tmp <- c(tmp, "#Sex\tAge\tArea\tTStep\tState","\n")
for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Final check\n123456")

write.table(tmp, paste(floc,'/SELEXSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
