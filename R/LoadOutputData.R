LoadOutputData <- function(is95){
  library(dplyr)
  library(magrittr)
  library(reshape2)
  library(ggplot2)
  library(tidyr)
  library(openxlsx)

  dat  <- read.table("Output.RL",comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  echo  <- read.table("Echo.out",comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  selx <- read.table(paste("../SELEXSPEC.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  lbin1  <- read.table(paste("../DATA.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  ctl1  <- read.table(paste("../CONTROL.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  mov1  <- read.table(paste("../MOVESPEC.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
  wb <- loadWorkbook(file="../../MakeFiles/ModelStructure.xlsx")

  #This is the location of the data input files and their associated parameters
  dynamics <- readWorkbook(wb,sheet='dynamics', startRow = 2)
  startseason <- as.numeric(dynamics$value[dynamics$object=='startseason'])
  endseason <- as.numeric(dynamics$value[dynamics$object=='endseason'])
  yr <- startseason:endseason
  projectseason <- as.numeric(dynamics$value[dynamics$object=='projectedseason'])
  projectcatch <- as.numeric(dynamics$value[dynamics$object=='projectedcatch'])
  burnin <- as.numeric(dynamics$value[dynamics$object=='burnin'])
  ages <- as.numeric(dynamics$value[dynamics$object=='ages'])
  sexs <- 0:as.numeric(dynamics$value[dynamics$object=='sexs'])
  areas <- readWorkbook(wb,sheet='area', startRow = 2)
  nareas <- length(unique(areas$newarea))
  times <- readWorkbook(wb,sheet='times', startRow = 2)
  fleets <- readWorkbook(wb,sheet='fleetcode', startRow = 2)

  #mov2 <- findNclean(c('#','Movement'), mov1, 1,0, char=T)
  # KeyWord <- c('Cpue','data'); DataFile <- echo; Offset<- 1; char=T
  #KeyWord <- c('#','Parameter','Par') ; char = T; DataFile <- dat; Offset <- 1; convert=0 lb <- findNclean(, dat, 2)
  findNclean <- function(KeyWord, DataFile, Offset=1, convert=0, char=F){
    hash <- c(which(grepl('#',DataFile[,1])),nrow(DataFile))
    if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,1]))
    if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])))
    if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,3])))
    if(length(KeyWord)==4) pos1 <- which(4==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,3])+grepl(KeyWord[4],DataFile[,4])))
    if(!(pos1+1)%in%hash) {
      pos2 <- hash[hash>(pos1+Offset)][1]-1
      adj <- 0 } else {
        pos2 <- hash[hash>(pos1+1+Offset)][1]-1
        adj <- 1}
    if(pos2>pos1){
      tmp <- DataFile[(pos1+Offset):pos2,]
      tmp <- tmp[!grepl('#', tmp[,1],fixed = T),]
      tmp <- tmp[tmp[,1]!='',]
      if(nrow(tmp)>0){
        rname <- DataFile[(pos1+adj),]
        rname <- gsub('#','',rname)
        rname <- rname[!is.na(rname) & rname!='' & rname!='NA']
        maxcol <- max(which(!is.na(tmp) & tmp!='', arr.ind = T)[,2])
        tmp <- tmp[!is.na(tmp[,1]),1:maxcol]
        tmp[tmp=='NaN'] <- 0
        #Add tweak here
        if(!is.null(dim(tmp))) if(length(rname)!=ncol(tmp)) { rname <- c(rname, paste('a',1:200,sep=''))[1:ncol(tmp)]}
        if(is.null(ncol(tmp))){ return(as.numeric(tmp)) } else {
          convert1 <- 1:ncol(tmp)
          convert <- convert1[!convert1%in%convert]
          tmp <- data.frame(tmp)
          chartmp <- tmp
          if(nrow(tmp)==1)  suppressWarnings(tmp[convert] <- (apply(as.matrix(tmp[,convert]),2,function(q) as.numeric(as.character(q)))))
          if(nrow(tmp)>1)   suppressWarnings(tmp[,convert] <- data.frame(apply(as.matrix(tmp[,convert]),2,function(q) as.numeric(as.character(q)))))
          tmp <- tmp[!is.na(tmp[,1]),!is.na(tmp[1,])]
          if(char==T) { return(chartmp) }
          if(char==F) {if(!is.null(dim(tmp))) colnames(tmp) <- rname[1:length(colnames(tmp))]
          return(tmp) }}} else { return(NA) } } else { return(NA) }
  }
  # find <- function(KeyWord, DataFile, Offset){
  #   KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
  #   if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,1]))+Offset
  #   if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])))+Offset
  #   if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,3])))+Offset
  #   return(pos1)}

  find <- function(KeyWord, DataFile, Offset){
    KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
    pos1 <- list(NA)
    for(K in 1:length(KeyWord)){
      pos1[[K]] <- which(grepl(KeyWord[K],DataFile[,K]))
    }
    if(length(KeyWord)>1) {ids <- NULL
    for(K in 1:(length(KeyWord)-1)){
      ids <- c(ids,pos1[[K]][which(pos1[[K]]%in%pos1[[K+1]])] )
    }
    ids <- table(ids)
    return(as.numeric(names(which.max(ids)))+Offset)} else {
      return(as.numeric(pos1[[1]])) }}



  ##Get parameters
  p1 <- find("# parameter table",dat,0)
  p2 <- find("#Total estimated parameters",dat,0)
  pout <- dat[(p1+2):(p2-1),c(1,3,4)]
  names(pout) <- c('name', 'estimated', 'value')
  pout$estimated <- ifelse(is.na(pout$estimated), 0, 1)

  ## GetSDReport
  sdr <- read.delim(paste("SDReport.RL",sep=''), sep=' ')
  names(sdr) <- c('name','Estimate','SE')
  sdr %<>% mutate(SE=ifelse(is.na(SE),0,SE)) %>% filter(!is.na(name), nchar(name)>0) %>% mutate(up95=Estimate+SE*1.96,low95=Estimate-SE*1.96,upr68=Estimate+SE,low68=Estimate-SE, cv=SE/Estimate, lwr=low95, upr=up95)
  if(!is95) sdr %<>% mutate(lwr=low68,upr=upr68)

  Fdims <- function(x){
    dims <- c(1,1)
    if(x==2)   dims <- c(1,2)
    if(x%in%3:4)   dims <- c(2,2)
    if(x%in%5:6)   dims <- c(2,3)
    if(x%in%7:9)   dims <- c(3,3)
    if(x%in%9:12)   dims <- c(3,4)
    if(x>12)   {     nr <- ceiling(x/4)
    dims <- c(nr,4)
      }
    return(dims)}

  ##define lenbin
  lbin <- findNclean(c('#','Lower'), lbin1, 1, convert=0)
  lbin <- as.numeric(unname(lbin[1,]))
  lbin <- lbin + diff(lbin)[1]/2
  lbinl <- lbin
  lbin <- lbin[1:(length(lbin)-1)]

}
