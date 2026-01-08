library(dplyr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(tidyr)

if(!exists("fls")) fls <- list.files(pattern = 'Run')[1]
dat  <- read.table(paste("../",fls,"/Output/Output.RL",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)

selx <- read.table(paste("../",fls,"/SELEXSPEC.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)

lbin1  <- read.table(paste("../",fls,"/DATA.dat",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)

ctl1  <- read.table(paste("../",fls,"/CONTROL.dat",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)

mov1  <- read.table(paste("../",fls,"/MOVESPEC.dat",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)

#mov2 <- findNclean(c('#','Movement'), mov1, 1,0, char=T)
# KeyWord <- c('#','Years','over'); DataFile <- lbin1; Offset<- 1; char=T
#KeyWord <- c('#','Movement','parameters'); DataFile <- mov1; Offset <- 1; convert=0
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
      tmp <- tmp[!is.na(tmp[,1]),!is.na(tmp[1,]) & tmp[1,]!='']
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
find <- function(KeyWord, DataFile, Offset){
  KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
  if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,1]))+Offset
  if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])))+Offset
  if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,1])+grepl(KeyWord[2],DataFile[,2])+grepl(KeyWord[3],DataFile[,4])))+Offset
  return(pos1)}
##Get parameters
p1 <- which(dat[,1]=="#" & dat[,2]=='parameter' & dat[,3]=='table')
p2 <- which(dat[,1]=="#Total" & dat[,2]=='estimated' & dat[,3]=='parameters:')
pout <- dat[(p1+2):(p2-1),c(1,3,4)]
names(pout) <- c('name', 'estimated', 'value')
pout$estimated <- ifelse(is.na(pout$estimated), 0, 1)

## GetSDReport
sdr <- read.delim(paste("../",fls,"/Output/SDReport.RL",sep=''), sep=' ')
#sdr <- read.delim(paste("SDReport.RL",sep=''), sep=' ')
names(sdr) <- c('name','Estimate','SE')
sdr %<>% mutate(SE=ifelse(is.na(SE),0,SE)) %>% filter(!is.na(name), nchar(name)>0) %>% mutate(up95=Estimate+SE*1.96,low95=Estimate-SE*1.96, cv=SE/Estimate)

Fdims <- function(x){
  dims <- c(1,1)
  if(x==2)   dims <- c(1,2)
  if(x%in%3:4)   dims <- c(2,2)
  if(x%in%5:6)   dims <- c(2,3)
  if(x%in%7:9)   dims <- c(3,3)
  if(x%in%9:12)   dims <- c(3,4)
  if(x>12)   dims <- c(4,4)
  return(dims)}

##define lenbin
lbin <- findNclean(c('#','Lower'), lbin1, 1, convert=0)
lbin <- as.numeric(unname(lbin[1,]))
lbin <- lbin + diff(lbin)[1]/2
lbinl <- lbin
lbin <- lbin[1:(length(lbin)-1)]
