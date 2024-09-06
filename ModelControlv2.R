setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("R files\\Utilities.r")  
source("R files\\ReadMaterial.r")
source("R files\\WriteMaterial.r")
source("R files\\Estimate.r")
source("R files\\WriteDat.r")
library(readxl)
library(reshape2)
library('svDialogs')
library(stats4)
library(TMB)
#################################################################################

compile("Model1.cpp")#, flags="-Wno-ignored-attributes")
dyn.load(dynlib("Model1"))

################################################################################
# Choose Model run
fls <- csemod(list.files(pattern = 'Run')) 
setwd(paste(getwd(), "\\",fls,sep=""))

## Update starting parameters to last run 
UpdatePars()

# New boolean to write more stuff out
FullOutput <- F

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


################################################################################
# TMB
################################################################################

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

# Set the projection flag (0=n0;1=Yes)
Data$DoProject = 0;

## Check the input data for NA and/or missing data and Read in the initial values
for(i in 1:length(Data)){ isnafunc(Data[[i]])}
InitialVars <- ReadInitialValues(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,
                                 GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
dlg_message("Check Console for summary of parameter inputs")

if(Data$Narea==1){## need to trick SetInitialAndPhases because only one area
  InitialVars$RecSpatDevs$Initial <- 0
  InitialVars$RecSpatDevs$Bnd <- c(-15,15)
  InitialVars$RecSpatDevs$Phase <- -1
  }
Parssolved(InitialVars) ## Records which parameters were used to solve the model for output file
parameters <- list(MainPars=NULL,RecruitPars=NULL,PuerPowPars=NULL,SelPars=NULL,RetPars=NULL,RecDevs=NULL,Qpars=NULL,efpars=NULL,InitPars=NULL,RecSpatDevs=NULL,MovePars=NULL,GrowthPars=NULL,dummy=0)

(MaxPhase <- getPhase(InitialVars));ParOld <- NULL;CurrPhase <- 1
#EvalMax <- rep(100,MaxPhase)
#StoreReport <-NULL

for (CurrPhase in 1:MaxPhase) {
  MaXeVaL <- ifelse(CurrPhase<MaxPhase, 500, 5000)
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
  mout<-nlminb(mout$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL, trace=50))
  pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par; 
  cat("Convergence:",round(mout$objective,4),"Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","Interations:",mout$iterations,"Evalutions:",mout$evaluations,"\n")
  # Store and save parameters
  pout <- unlist(parameters); pout[names(pout)%in%names(pars)] <- pars; suffix <- ifelse(CurrPhase==MaxPhase," final", CurrPhase); write.table(pout, paste("Output/model",suffix,".par",sep=""), sep='\t', col.names = c('name\test'), quote=F)
}

model$fn()
Report <- model$report()
lb <- Report$LegalBio76
yrs <- Data$Year1:Data$Year2
par(mfrow=c(3,3), mar=c(4,5,0.5,2),las=1)
for(i in 1:8){
  tlb <- lb[,i]/1000000
  mxY <- max(tlb)
  plot(yrs, tlb, type='o', axes=F, pch=16, cex=0.9,ylab='Total Legal Biomass (1000s t)', xlab='Season', main=paste('Area ',i), xlim=c(1975,2025), ylim=c(0,mxY))
  axis(1)
  axis(2)}

best <- mout$par
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
BigSave$lowlike <- model$fn(best)
save(BigSave,file="Output/BigSave.lda")
WriteOutput(Report,rep,fullrep,parameters,GeneralSpecs,ControlSpecs,TheData, CurrPhase=CurrPhase,best=best)
# plot((Data$Year1-Data$BurnIn):(Data$Year2+Data$MaxProjYr), rowSums(Report$LegalBioAll)/1000, type='o', ylim=c(0,1000))
# abline(v=1978, col=2)

## Run and output diagnostics file
source("../R files/MakeDiagnostics.R")
