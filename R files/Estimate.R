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
      map <- append(map,list(dummy=factor(1)))
      estvec <- c(estvec,0)
      lowBnd <- c(lowBnd,-1)
      uppBnd <- c(uppBnd,1)
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

UpdatePars <- function(){
  todo <- dlg_list(c('Yes','No','                '),  title=c('Update Parameters?                    '))$res  
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
  locs <-data.frame(par=c('MainPars', 'RecruitPars','PuerPowPars','RecDevs','Qpars','efpars','RecSpatDevs','MovePars'), file=c('Control.DAT', 'RECRUITSPEC.dat', 'RECRUITSPEC.dat','Control.DAT','Control.DAT','Control.DAT','Control.DAT','MOVESPEC.dat' ), id=c('Basic parameters','Recuitment1 parameters','Puerulus Power for','Prespecify_rec_devs','Q parameters','Efficiency parameters','Prespecify_spatial_rec_devs','Movement parameters'), off=c(1,2,3,2,1,2,2,2), col=c(3,3,3,1,3,3,1,3))
  
  for (i in 1:nrow(locs)){
    DataFile <- read.table(locs$file[i],comment.char = "?",fill=T, blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
    ptmp <- round(pout$est[grepl(locs$par[i], pout$name)],5)
    roff <- locs$off[i]
    coff <- locs$col[i]
    pos <- find(c(unlist(strsplit(locs$id[i],' '))), DataFile, roff)
    DataFile[pos:(pos+length(ptmp)-1),coff] <- ptmp
    write.table(DataFile, locs$file[i], na=" ", sep=" ", row.names = F, col.names = F, quote=F)
    print(paste("Parameters upated: ", locs$par[i]))
  }
}}


SolveModel <- function(phit=500,lphit=5000, mxph=MaxPhase){
  MaxPhase=mxph
  for (CurrPhase in 1:MaxPhase) {
    MaXeVaL <- ifelse(CurrPhase<MaxPhase, phit, lphit)
    RunSpecs <- SetInitialAndPhases(ParOld,parameters,InitialVars,CurrPhase=CurrPhase)  # Set parameters and mapping
    parameters <- list(MainPars=InitialVars$MainPars$Initial,RecruitPars=InitialVars$RecruitPars$Initial,PuerPowPars=InitialVars$PuerPowPars$Initial,SelPars=InitialVars$SelPars$Initial,RetPars=InitialVars$RetPars$Initial,RecDevs=InitialVars$RecDevs$Initial,Qpars=InitialVars$Qpars$Initial,efpars=InitialVars$efpars$Initial,InitPars=InitialVars$InitPars$Initial,RecSpatDevs=InitialVars$RecSpatDevs$Initial,MovePars=InitialVars$MovePars$Initial,GrowthPars=InitialVars$GrowthPars$Initial,dummy=0)
    ## Make model
    cat("Making model object that will solve for",sum(!is.na(unlist(RunSpecs$map))) ,"parameters.","Phase =",CurrPhase,"\n")
    pnames <- names(unlist(RunSpecs$map)[!is.na(unlist(RunSpecs$map))]); print(pnames); map <- RunSpecs$map
    model <- MakeADFun(Data, parameters, map=map, DLL="Model1",silent=T)
    model$par <- RunSpecs$EstVec  # Assign new parameters associated with the correct phase
    # Run initial model
    system.time(mout<-nlminb(model$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL, trace=50)))
    # re-run model
    mout<<-nlminb(mout$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL, trace=50))
    pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par; 
    cat("Convergence:",round(mout$objective,4),"Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","Interations:",mout$iterations,"Evalutions:",mout$evaluations,"\n")
    # Store and save parameters
    pout <- unlist(parameters); pout[names(pout)%in%names(pars)] <- pars; suffix <- ifelse(CurrPhase==MaxPhase," final", CurrPhase); write.table(pout, paste("Output/model",suffix,".par",sep=""), sep='\t', col.names = c('name\test'), quote=F)
  }
}

LoadPars <- function(){
  for(i in 1:length(Data)){ isnafunc(Data[[i]])}
  InitialVars <- ReadInitialValues(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
  dlg_message("Check Console for summary of parameter inputs")
  if(Data$Narea==1){## need to trick SetInitialAndPhases because only one area
    InitialVars$RecSpatDevs$Initial <- 0
    InitialVars$RecSpatDevs$Bnd <- c(-15,15)
    InitialVars$RecSpatDevs$Phase <- -1
  }
  Parssolved(InitialVars) ## Records which parameters were used to solve the model for output file
  parameters <- list(MainPars=NULL,RecruitPars=NULL,PuerPowPars=NULL,SelPars=NULL,RetPars=NULL,RecDevs=NULL,Qpars=NULL,efpars=NULL,InitPars=NULL,RecSpatDevs=NULL,MovePars=NULL,GrowthPars=NULL,dummy=0)
  
  (MaxPhase <- getPhase(InitialVars));ParOld <- NULL;CurrPhase <- 1
}


LoadData <- function() {
  # read in files
  StarterFile <- read.table("Starter.dat",comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100) 
  Starter <- ReadStarterFile(StarterFile)
  DataFile <- read.table(Starter$DataFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  ControlFile <- read.table(Starter$ControlFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  SelexFile <- read.table(Starter$SelexFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  RetenFile <- read.table(Starter$RetainFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  RecruitFile <- read.table(Starter$RecruitFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  GrowthFile <- read.table(Starter$GrowthFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  MoveFile <- read.table(Starter$MoveFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  #TagFile <- read.table(Starter$TagFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  #PropnFile <- read.table(Starter$PropFFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  ProjFile <- read.table(Starter$ProjectionsFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  
  EchoFile <- "Output/Echo.out"
  write("The is an echo file",EchoFile)
  
  # Read in the data file
  GeneralSpecs <- ReadGeneralFile(DataFile)
  Data <- GeneralSpecs
  
  # Read in the data file
  TheData <- ReadDataFile(DataFile,GeneralSpecs)
  Data <- append(Data,TheData)
  GeneralSpecs$NQpars <- Data$NQparPass
  
  # # Read in the tag file
  #TagSpecs <- ReadTagFile(TagFile,PropnFile,GeneralSpecs,TheData)
  #Data <- append(Data,TagSpecs)
  
  # Read in the Control file
  ControlSpecs <- ReadControlFile(ControlFile,GeneralSpecs,TheData)
  Data <- append(Data,ControlSpecs)
  # Read in the projections file
  ProjectSpecs <- ReadProjFile(ProjFile,GeneralSpecs,Data$Phi1) 
  Data <- append(Data,ProjectSpecs)
  # Read in the Selectivity file
  SelexSpecs <- ReadSelexFile(SelexFile,GeneralSpecs)
  Data <- append(Data,SelexSpecs)
  # Read in the Retention file
  RetenSpecs <-ReadRetenFile(RetenFile,GeneralSpecs)
  Data <- append(Data,RetenSpecs)
  # Read in the Recruiment file
  RecruitSpecs <-ReadRecruitFile(RecruitFile,GeneralSpecs)
  Data <- append(Data,RecruitSpecs)
  # Read in the Growth file
  GrowthSpecs <-ReadGrowthFile(GrowthFile,GeneralSpecs)
  Data <- append(Data,GrowthSpecs)
  # Read in the Movement file
  MoveSpecs <- ReadMoveFile(MoveFile,GeneralSpecs)
  Data <- append(Data,MoveSpecs)
}

RunReport <- function(){
  print("Loading report")
  Report <- model$report()
  best <- mout$par
  print("Loading SD report (can take quite a long time)")
  rep <- sdreport(model)  
  fullrep <- summary(rep)
  BigSave <-NULL
  BigSave$Report <- Report
  BigSave$rep <- rep
  BigSave$map <- RunSpecs$map
  BigSave$Data <- Data
  BigSave$fullrep <- fullrep
  BigSave$parameters <- parameters
  BigSave$best <- best
  #BigSave$lowlike <- model$fn(best)
  save(BigSave,file="Output/BigSave.lda")
  print("making Output.RL")
  WriteOutput(Report,rep,fullrep,parameters,GeneralSpecs,ControlSpecs,TheData, CurrPhase=CurrPhase,best=best)
  ## Run and output diagnostics file
  print("making Diagnostics report")
  source("../R files/MakeDiagnostics.R")
}
