################################################################################

SetInitialAndPhases <- function(ParOld,parameters,InitialVars,CurrPhase)
{

  # Pointer to old parameters
  Ipnt <- 0
  SinglePhase<-function(Est,Bnd,Phase,CurrPhase)
  {
    map <- NULL
    estVec <- NULL; lowBnd <- NULL; uppBnd <- NULL
    Npar <- length(Est)

    for (Ipar in 1:Npar)
    {
      if (Phase[Ipar] > 0 & Phase[Ipar] <= (CurrPhase-1))
      {
        Ipnt <<- Ipnt + 1
        Est[Ipar] <- ParOld[Ipnt]
      }
      if (Phase[Ipar] > 0 & Phase[Ipar] <= CurrPhase)
      {
        map <- c(map,Ipar)
        estVec <- c(estVec,Est[Ipar])
        lowBnd <- c(lowBnd,Bnd[Ipar,1])
        uppBnd <- c(uppBnd,Bnd[Ipar,2])
      }
      else
        map <- c(map,as.factor(NA))
    }
    ReturnObj <- NULL
    ReturnObj$map <- as.factor(map)
    ReturnObj$estVec <- estVec
    ReturnObj$lowBnd <- lowBnd
    ReturnObj$uppBnd <- uppBnd
    return(ReturnObj)
  }  # SinglePhase

  map <- list()
  estvec <- NULL; lowBnd <- NULL; uppBnd <- NULL
  for (ParName in names(parameters))
  {
    if (ParName != "dummy")
    {
      if (length(InitialVars[[ParName]]$Initial) >0)
      {
        ThePar <- InitialVars[[ParName]]
        PhaseOut <- SinglePhase(ThePar$Initial,ThePar$Bnd,ThePar$Phase,CurrPhase)
        map <- append(map,list(ParName=PhaseOut$map))
        estvec <- c(estvec,PhaseOut$estVec)
        lowBnd <- c(lowBnd,PhaseOut$lowBnd)
        uppBnd <- c(uppBnd,PhaseOut$uppBnd)
      }
      else
      {
        estvec <- c(estvec,0)
        map <- append(map,list(ParName=factor(NA)))
      }
    }
    else
    {
      if(max(as.vector(sapply(map, function(x) max(as.numeric(!is.na(x))))))==0) {
        map <- append(map,list(dummy=factor(1)))
        estvec <- c(estvec,0)
        lowBnd <- c(lowBnd,-1)
        uppBnd <- c(uppBnd,1)
      } else { map <- append(map,list(dummy=factor(NA)))}
    }
  }
  names(map) <- names(parameters)
  # dummy
  ReturnObj <- NULL
  ReturnObj$map = map
  ReturnObj$parameters <- parameters
  ReturnObj$EstVec <- estvec
  ReturnObj$lowBnd <- lowBnd
  ReturnObj$uppBnd <- uppBnd
  return(ReturnObj)
}


getPhase <- function(InitialVars){
  mX <- 1
  num <- length(InitialVars)
  for(n in 1:num){
    tmp <- InitialVars[[n]]
    mX <- max(c(mX, tmp$Phase))
  }
  return(mX)
}

csemod <- function(x){
  mod <-  as.numeric(dlg_input(c('Choose a model:',paste(1:length(x), x, sep=(" : ") )), 1)$res)
  if (!length(mod)) {# The user clicked the 'cancel' button
    cat(paste("OK, the default model is",x[1],"\n"))
  } else {
    cat(paste("Model", x[mod], "has been chosen"), "\n")
  }
  return(x[mod])}

choose_model <- function(pattern = 'Run') {
  x <- list.files(pattern = pattern)

  mod <- as.numeric(dlg_input(c('Choose a model:',
                                paste(1:length(x), x, sep = " : ")),
                              1)$res)

  if (!length(mod)) {
    cat(paste("OK, the default model is", x[1], "\n"))
    selected <- x[1]
  } else {
    cat(paste("Model", x[mod], "has been chosen"), "\n")
    selected <- x[mod]
  }

  setwd(file.path(getwd(), selected))
  invisible(selected)  # Returns but doesn't auto-print
}


UpdatePars <- function(todo=' '){
  if(todo==' ') {todo <- dlg_list(c('Yes','No','                '),  title=c('Update Parameters?                    '))$res  }
  ## Get estimated parameters  KeyWord <- locs$id[i]
  if(todo=='Yes'){
find <- function(KeyWord, DataFile, Offset){
    KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
    if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,2]))+Offset
    if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,2])+grepl(KeyWord[2],DataFile[,3])))+Offset
    if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,2])+grepl(KeyWord[2],DataFile[,3])+grepl(KeyWord[3],DataFile[,4])))+Offset
    return(pos1)}

  pout <- read.delim('Output/model final.par',sep='\t',stringsAsFactors =F)
  unique(pout$name)
  locs <-data.frame(par=c('MainPars', 'RecruitPars','PuerPowPars','RecDevs','Qpars','efpars','RecSpatDevs','MovePars','SelPars'), file=c('CONTROL.DAT', 'RECRUITSPEC.DAT', 'RECRUITSPEC.DAT','CONTROL.DAT','CONTROL.DAT','CONTROL.DAT','CONTROL.DAT','MOVESPEC.DAT','SELEXSPEC.DAT'), id=c('Basic parameters','Recuitment1 parameters','Puerulus Power for','Prespecify_rec_devs','Q parameters','Efficiency parameters','Prespecify_spatial_rec_devs','Movement parameters','Selectivity Parameters'), off=c(1,2,3,2,1,2,2,2,2), col=c(3,3,3,1,3,3,1,3,3))

  for (i in 1:nrow(locs)){
    DataFile <- read.table(locs$file[i],comment.char = "?",fill=T, blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
    ptmp <- pout$est[grepl(locs$par[i], pout$name)]
    roff <- locs$off[i]
    coff <- locs$col[i]
    pos <- find(c(unlist(strsplit(locs$id[i],' '))), DataFile, roff)
    DataFile[pos:(pos+length(ptmp)-1),coff] <- ptmp
    write.table(DataFile, locs$file[i], na=" ", sep=" ", row.names = F, col.names = F, quote=F)
    print(paste("Parameters upated: ", locs$par[i]))
  }
}}



SolveModelNew <- function(phit=500,lphit=1000, mxph=MaxPhase, PrintLag = 50, report=F){
  MaxPhase=ifelse(mxph==0,1,mxph)
  for (CurrPhase in 1:MaxPhase) {
    MaXeVaL <- ifelse(CurrPhase<MaxPhase, phit, lphit)
    parameters <- list(MainPars=InitialVars$MainPars$Initial,RecruitPars=InitialVars$RecruitPars$Initial,PuerPowPars=InitialVars$PuerPowPars$Initial,SelPars=InitialVars$SelPars$Initial,RetPars=InitialVars$RetPars$Initial,RecDevs=InitialVars$RecDevs$Initial,Qpars=InitialVars$Qpars$Initial,efpars=InitialVars$efpars$Initial,InitPars=InitialVars$InitPars$Initial,RecSpatDevs=InitialVars$RecSpatDevs$Initial,MovePars=InitialVars$MovePars$Initial,GrowthPars=InitialVars$GrowthPars$Initial,dummy=0)
    RunSpecs <- SetInitialAndPhases(ParOld,parameters,InitialVars,CurrPhase=CurrPhase)  # Set parameters and mapping
    ## Make model
    cat("Making model object that will solve for",sum(!is.na(unlist(RunSpecs$map))) ,"parameters.","Phase =",CurrPhase,"\n")
    pnames <- names(unlist(RunSpecs$map)[!is.na(unlist(RunSpecs$map))]);
    nam <- stringr::str_extract(pnames, "[\\p{Letter}]+")
    num <- stringr::str_extract(pnames, "\\d+$")
    unnam <- nam[!duplicated(nam)]
    for(iii in 1:length(unnam))  { print(paste(unnam[iii], length(num[nam==unnam[iii]]),'parameters'))  }
    model <- MakeADFun(Data, parameters, map=RunSpecs$map, DLL="IMuLT",silent=T)
    BestFn <- model$fn()
    initBestFn <- BestFn
    FnCallNo <<- 0;
    model$fn_Orig <- model$fn
    yy <- 1e+10
    model$fn <- function(x)  {
      tyy <- model$fn_Orig(x)
      yy <<- ifelse(is.na(tyy),yy,tyy)
      FnCallNo <<- FnCallNo + 1
      if(BestFn>yy){
        if ((FnCallNo %% PrintLag)==0) {
          cat("Phase ", CurrPhase," ",FnCallNo," -LogLike / Delta: ",yy," / ",round(100*(1-(yy/BestFn)),6),"%; npar = ",length(x),"\n",sep="")
          BestFn <<- yy }
      }
      return(tyy);  }
    model$par <- RunSpecs$EstVec  # Assign new parameters associated with the correct phase
    # Run model
    BestFn <- model$fn(model$par)
    initBestFn <- BestFn
    mout<-nlminb(model$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL))
    pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par;
    Grad <- model$gr(mout$par)
    grad <- round(max(abs(Grad)),6)
    cat("Likelihood: ",round(initBestFn,6),' to ' ,round(mout$objective,6),"Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","Max Gradient:",round(max(abs(model$gr(mout$par))),6)," Interations:",mout$iterations,"Evalutions:",mout$evaluations,"\n")
    # Store and save parameters
    pout <- unlist(parameters); pout[names(pout)%in%names(pars)] <- pars; suffix <- ifelse(CurrPhase==MaxPhase," final", CurrPhase); write.table(pout, paste("Output/model",suffix,".par",sep=""), sep='\t', col.names = c('name\test'), quote=F)
  if(report==T & CurrPhase==MaxPhase){       cat("Making report object.\n")
      print("Loading report")
      Report <- model$report()
      best <- mout$par
      print("Loading SD report (can take quite a long time)")
      SDrep <- sdreport(model)
      fullrep <- summary(SDrep)
      BigSave <-NULL
      BigSave$Report <- Report
      BigSave$SDrep <- SDrep
      BigSave$map <- RunSpecs$map
      BigSave$Data <- Data
      BigSave$fullrep <- fullrep
      BigSave$parameters <- parameters
      BigSave$pin  <- pout
      BigSave$best <- best
      BigSave$Gradient <- Grad
      #BigSave$lowlike <- model$fn(best)
      if(max(list.files()=='Output')==1) { setwd(paste(getwd(), "/Output",sep=""))  }
      save(BigSave,file="BigSave.lda")
      print("making Output.RL")
      WriteOutput(Report,SDrep,fullrep,parameters,pout,GeneralSpecs,ControlSpecs,TheData,CurrPhase=0,best=best,grad=Grad)
      }
    }
  }



LoadPars <- function(aask=''){
  for(i in 1:length(Data)){ isnafunc(Data[[i]],i)}
  InitialVars <<- ReadInitialValues(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
  if(aask=='test')dlg_message("Check Console for summary of parameter inputs")
  if(Data$Narea==1){## need to trick SetInitialAndPhases because only one area
    InitialVars$RecSpatDevs$Initial <<- 0
    InitialVars$RecSpatDevs$Bnd <<- c(-15,15)
    InitialVars$RecSpatDevs$Phase <<- -1
  }
  Parssolved(InitialVars) ## Records which parameters were used to solve the model for output file
  parameters <- list(MainPars=NULL,RecruitPars=NULL,PuerPowPars=NULL,SelPars=NULL,RetPars=NULL,RecDevs=NULL,Qpars=NULL,efpars=NULL,InitPars=NULL,RecSpatDevs=NULL,MovePars=NULL,GrowthPars=NULL,dummy=0)
  (MaxPhase <<- getPhase(InitialVars));ParOld <<- NULL;CurrPhase <<- 1
}

MakeDiagReport <- function(is95=T) {
  ## Run and output diagnostics file
  print("Making Diagnostics report")
  if(max(list.files()%in%'Output')==1) {  setwd(makehtml::filenametopath(getwd(),'Output'))}
  #source('../../R files/MakeOutPut.R')
  MakeOutPut(is95)
}

AdjustPhase <- function(dum=' '){
  MaxPhase <- 0
  for(i in 1:length(names(InitialVars))) {
    if(dum!='dummy') {
      todo <- dlg_list(c('All','Individual','Skip','                '),  title=paste('Phase for',names(InitialVars)[i],'                             ') )$res

      ## Get estimated parameters  KeyWord <- locs$id[i]
      if(todo=='All'){
        nphase <- dlg_list(c('-1','1','2','3','4','                '),  title=paste('Phase for all',names(InitialVars)[i],'                             ') )$res
        InitialVars[[i]]$Phase <- rep(as.numeric(nphase), length(InitialVars[[i]]$Phase))
      }

      if(todo=='Individual'){
        for(ip in 1:length(InitialVars[[i]]$Initial)){
          nphase <- dlg_list(c('-1','1','2','3','4','                '),  title=paste('Phase',names(InitialVars)[i],'par #',ip,'                             ') )$res
          InitialVars[[i]]$Phase[ip] <- as.numeric(nphase)
        }
      }
      if(max(InitialVars[[i]]$Phase)>MaxPhase) MaxPhase <- max(InitialVars[[i]]$Phase)
    }
    if(dum=='dummy') {
      InitialVars[[i]]$Phase <- rep(-1, length(InitialVars[[i]]$Phase))
      tmp <- InitialVars$Growth; tmp$Phase <- 1
      InitialVars$dummy <- tmp
    }
  }
  InitialVars <<- InitialVars
  MaxPhase <<- MaxPhase
}

UpdateLFWeights <- function(todo='No'){
  library(dplyr); library(magrittr)
  if(todo=='Yes'){
    tmp  <- read.csv(paste0(getwd(),"/Output/Summary/result/Tuning.csv"))
    # open Control file and find weightings
    DataFile <- read.table('CONTROL.DAT',comment.char = "?",fill=T, blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
    ## Find the length freqs weights
    pos1 <- which(1==(grepl('length',DataFile[,8])))
    if(length(pos1)!=nrow(tmp)) {
      print("Predetermined weights do not match length compoitions")
      break  }
    for(i in 1:length(pos1)){
      ttmp <- DataFile[pos1[i],1:10]
      Scale <- tmp$Multiscale[(tmp$Fleet-1)==as.numeric(ttmp[,2]) & (tmp$Sex-1)==as.numeric(ttmp[,4])]
      DataFile[pos1[i],5] <-  as.numeric(DataFile[pos1[i],5]) * as.numeric(Scale)
    }
    write.table(DataFile, 'CONTROL.DAT',na=" ", sep=" ", row.names = F, col.names = F, quote=F)
  }
  }

