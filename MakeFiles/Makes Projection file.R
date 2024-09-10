find <- function(KeyWord, DataFile, Offset){
  KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
  if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,1]))+Offset
  if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])))+Offset
  if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,3])))+Offset
  if(length(KeyWord)==4) pos1 <- which(4==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,3])+grepl(KeyWord[4],DataFile[,4])))+Offset
  return(pos1)}


tmp <- list()
tmp <- c(tmp, "# Notes\n")
tmp <- c(tmp, "# Recruitment, movement, growth propotions are as for the last year\n")
tmp <- c(tmp, "\n# Number of projection years (must be less than the maximum number of projection years)\n")
tmp <- c(tmp,projectseason,'\n')
tmp <- c(tmp, "# Selectivity\n")
tmp <- c(tmp, "# Specifications for selectivity (Selectivity of the pots [escape gaps, males then females], 0 = None, 1 = 54 mm, 2 = 55 mm)\n")

tdat  <- read.table(paste(floc,'/SELEXSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
pos1 <- find(c("#",'Specifications','for','selectivity'), tdat, 2)
pos2 <- find(c("#",'selectivity'), tdat, -1)
tdat <- tdat[pos1:pos2,c(1:4, sum(!is.na(tdat[pos1,])))  ]
for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  } 
colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n# Retention\n")
tmp <- c(tmp, "# Specifications for retention. Fleet = comm comm comm comm comm comm comm comm comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor rec rec rec rec ibss ibss ibss ibss ibss iss iss iss iss\n")
tdat  <- read.table(paste(floc,'/RETAINSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
pos1 <- find(c("#",'Specifications','for','retention'), tdat, 2)
pos2 <- find(c("#",'retention',"-"), tdat, -1)
tdat <- tdat[pos1:pos2,c(1:4, sum(!is.na(tdat[pos1,])))  ]
for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  } 
colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "#\n# Specifications for Fleet legal assignment\n")
tdat  <- read.table(paste(floc,'/SELEXSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
pos1 <- find(c("#",'Specifications','for','Fleet'), tdat, 2)
pos2 <- find(c("#",'Specifications','for','legal'), tdat, -1)
tdat <- tdat[pos1:pos2,c(1:4, sum(!is.na(tdat[pos1,])))  ]
for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  } 
colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

tmp <- c(tmp, "\n#	Discard	mortality\n")
tdat  <- read.table(paste(floc,'/CONTROL.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
pos1 <- find(c("#",'Discard','mortality'), tdat, 2)
pos2 <- find(c("#",'Recruitment_deviations'), tdat, -1)
tdat <- tdat[pos1:pos2,c(1:3, sum(!is.na(tdat[pos1,])))  ]
for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  } 
colnames(tdat) <- c("Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

# tmp <- c(tmp, "\n#\n# Recruitment index\n# Number of data points (puerulus samples)")
# tdat <- expand.grid(Area=sort(unique(areas$narea)), Year=(endseason+1):(endseason+projectseason))
# puer <- read.csv(direct('Data/Puerulus/Data analysis/Puerulus for New model.csv'))
# puer %<>% mutate(area=area,mn=round(mn,1), sd=round(cv*mn,2))  %>% dplyr::select(area,season, mn, sd) %>% 
#   filter(season>endseason, area %in% areas$narea)
# tdat$Index[tdat$Year<=(min(tdat$Year)+1)] <- puer$mn[match(paste(tdat$Area[tdat$Year<=(min(tdat$Year)+1)], tdat$Year[tdat$Year<=(min(tdat$Year)+1)]), paste(puer$area,puer$season))]
# tdat$SD[tdat$Year<=(min(tdat$Year)+1)] <- puer$sd[match(paste(tdat$Area[tdat$Year<=(min(tdat$Year)+1)], tdat$Year[tdat$Year<=(min(tdat$Year)+1)]), paste(puer$area,puer$season))]
# tdat$Index[is.na(tdat$Index)] <- -1.0
# tdat$SD[is.na(tdat$SD)] <- -1.0
# tmp <- c(tmp, "\n# ",paste(colnames(tdat), collapse = "\t"), "\n")
# for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

tdat  <- read.table(paste(floc,'/DATA.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
pos1 <- find(c("#",'Catch','data'), tdat, 3)
pos2 <- find(c("#",'Index','data'), tdat, -1)
Names <- tdat[pos1-1,c(1:4)]
tdat <- tdat[pos1:pos2,c(1:4)]
names(tdat) <- Names 
tdat %<>% filter(`#Year`==max(`#Year`)) %>% mutate(prop=as.numeric(catch)/sum(as.numeric(catch))) %>% 
  mutate(catch=round(projectcatch*prop,1)) %>% dplyr::select(-prop)

tmp <- c(tmp, "\n# Specifications for projections (1=catch;2=effort;3=?)\n",1,"\n#\n# Catch data (kg) - Number of observations\n", projectseason*nrow(tdat), "\n")
tmp <- c(tmp,paste(colnames(tdat), collapse = "\t"), "\n")

for(i in 1:projectseason){ 
  tdat$`#Year` <- endseason+i 
  for(j in 1:nrow(tdat)){
    tmp <- c(tmp, paste(tdat[j,], collapse = "\t"),"\n")
  }}

write.table(tmp, paste(floc,'/PROJECTIONS.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

 