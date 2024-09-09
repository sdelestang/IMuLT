WriteOutput <- function(Report,stdrep,fullrep,parameters,GeneralSpecs,ControlSpecs,TheData,CurrPhase=2,best=rep(-1,1000)){

  files <- list.files()
  if(sum(files=='Output')>0) { setwd(paste(getwd(), "/Output",sep=""))  } 
  
  write.table(fullrep,'SDReport.RL',sep=' ', quote=F)
  stdrep <- summary(stdrep)
  OutputFile <- "Output.RL"
  ParFileName <- "Model.ppp"

  write("# Likelihood summary",OutputFile)
  write("Catch penalty",OutputFile,append=T)
  write(paste("Total objective function",Report$neglogL),OutputFile,append=T)
  write(paste("Catch likelihood",Report$CatchLike),OutputFile,append=T)
  write(paste("Cpue likelihood",Report$CpueLike, " (Weighted ", Report$Weighted_CpueLike, ")"),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLike, " (Weighted ", Report$Weighted_NumbersLike, ")"),OutputFile,append=T)
  write(paste("Length likelihood",Report$LengthLike, " (Weighted ", Report$Weighted_LengthLike, ")"),OutputFile,append=T)
  write(paste("Larval likelihood",Report$LarvalLike, " (Weighted ", Report$Weighted_LarvalLike, ")"),OutputFile,append=T)
 # write(paste("Tag1 Likelihood",sum(Report$TagLike1)," (Weighted ", Report$Weighted_TagLike1,")"),OutputFile,append=T)
 # write(paste("Tag2 Likelihood",sum(Report$TagLike2)," (Weighted ", Report$Weighted_TagLike2,")"),OutputFile,append=T)
  write(paste("Initial N penalty",Report$Initial_pen),OutputFile,append=T)
  write(paste("Recruitment penalty",Report$Rec_Penal),OutputFile,append=T)
  write(paste("Recruitment smooth penalty",Report$Rec_Penal_Smooth),OutputFile,append=T)

  write("\n# Likelihood by fleet",OutputFile,append=T)
  write(paste("Cpue likelihood",paste(Report$CpueLikeComp[])),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLikeComp),OutputFile,append=T)
  write(paste("Length likeliood",Report$LengthLikeComps),OutputFile,append=T)
  write(paste("Larval likeliood",Report$LarvalLikeComps),OutputFile,append=T)

  write("\n# Likelihood by group",OutputFile,append=T)
 # write(paste("Tag Likelihood 1",Report$TagLike1),OutputFile,append=T)
 # write(paste("Tag Likelihood 2",Report$TagLike2),OutputFile,append=T)


  #print("AEP IS STILL WORKING ON THIS")
  write("\n# Labels",OutputFile,append=T)
  write("#Simplified N-matrix",OutputFile,append=T)
  write("#F-matrix",OutputFile,append=T)
  write("#Recruitment patterns by timestep, area and sex",OutputFile,append=T)
  write("#Recruitment patterns by sex and size",OutputFile,append=T)

  
  write("\n#Initiation Option",OutputFile,append=T)
  write(Data$InitOpt,OutputFile,append=T)
  
  
  write("\n# parameter table",OutputFile,append=T)
  write("# Parameter Par_cnt Estpar_cnt Estimate SD lwrBound uprBound",OutputFile,append=T)
  ParName <- names(parameters)
  Ipnt <- 0; Iqnt <- 0
  write("# parameters",ParFileName)
  write("# dummy",ParFileName,append=T)
  write("0",ParFileName,append=T)
  for (ParName in names(parameters))
   {
    if (ParName != "dummy") 
     {
      print(ParName)
      ThePar <- InitialVars[[ParName]]
      for (Ipar in 1:length(ThePar$Initial))
       {  
        Iqnt <- Iqnt + 1
        write(paste("#",ParName,"_",Ipar," ",Iqnt," ",sep=""),ParFileName,append=T)
        if (ThePar$Phase[Ipar] > 0 & ThePar$Phase[Ipar] <= CurrPhase)  
         {
          Ipnt <- Ipnt + 1;
          xx <- paste(ParName,"_",Ipar," ",Iqnt," " ,Ipnt," ", stdrep[Ipnt,1]," ",stdrep[Ipnt,2]," ",ThePar$Bnd[Ipar,1]," ",ThePar$Bnd[Ipar,2],sep="")
          write(best[Ipnt],ParFileName,append=T)
         }  
        else
         {  
          xx <- paste(ParName,"_",Ipar," ",Iqnt," NA ", ThePar$Initial[Ipar],sep="")
          write(ThePar$Initial[Ipar],ParFileName,append=T)
         }  
        write(xx,OutputFile,append=T)
       }   
     } # Parameters within ParName
    } # ParNames
   write(paste("#Total estimated parameters:",Ipnt),OutputFile,append=T)  
  
  IvarPnt <- Ipnt + 1
    
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1
  
  #### ===================================================================================== 
  write("\n#Egg Production",OutputFile,append=T) 
  write("#Year Total By_area (SD?)",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+GeneralSpecs$MaxProjYr
  Ncol <- 1+1+GeneralSpecs$Narea
  SSBOut <- matrix(0,nrow=Nyears,ncol=Ncol)
  SSBOut[,1] <-GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1+1:Nyears
  SSBOut[,2] <- Report$MatBio[1:Nyears]
  for (Iarea in 1:GeneralSpecs$Narea)
   SSBOut[,2+Iarea] <- Report$MatBioArea[Iarea,1:Nyears]  
  for (Ivar in 1:ControlSpecs$NvarTypes)
   if (ControlSpecs$VarTypes[Ivar]==1)
    {
     SSBVar <- rep(0,Nyears)
     Ncol <- Ncol + 1
     for (IvarYr in 1:Nyears)
      {
       IvarPnt <- IvarPnt + 1
       SSBVar[IvarYr] <- fullrep[IvarPnt,2]
      }
     SSBOut <- cbind(SSBOut,SSBVar)
     IvarPnt <- IvarPnt+1
    }
   
  write(t(SSBOut),ncol=Ncol,OutputFile,append=T)
     
  write("\n#Total Recruitment",OutputFile,append=T) 
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  write(t(cbind(Years,Report$Recruits)),ncol=2,OutputFile,append=T)

  write("\n#Recruitment by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
   {
    Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1+1 
    Ncol <- 3
    RecOut <- cbind(rep(Iarea,Nyears+1),Years,Report$RecruitmentByArea[Iarea,])
    for (Ivar in 1:ControlSpecs$NvarTypes)
     if (ControlSpecs$VarTypes[Ivar]==2)
      {
       RecVar <- rep(0,Nyears+1)
       Ncol <- Ncol + 1
       for (IvarYr in 1:(Nyears+1))
        {
         IvarPnt <- IvarPnt + 1
         RecVar[IvarYr] <- fullrep[IvarPnt,2]
        }
       RecOut <- cbind(RecOut,RecVar)
      }
    write(t(RecOut),ncol=Ncol,OutputFile,append=T)
   } 

  write("\n#Recruitment Fractions",OutputFile,append=T) 
  lens <- Data$LowLenBin[1,1:dim(Report$RecruitFrac)[2]]
  write(t(cbind(lens,t(Report$RecruitFrac))),ncol=(dim(Report$RecruitFrac)[1]+1),OutputFile,append=T)

  write("\n#Legal Biomass by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
   {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$sLegalBio[,Iarea])
    for (Ivar in 1:ControlSpecs$NvarTypes)
     if (ControlSpecs$VarTypes[Ivar]==3)
      {
       LBVar <- rep(0,GeneralSpecs$Nyear)
       Ncol <- Ncol + 1
       for (IvarYr in 1:GeneralSpecs$Nyear)
        {
         IvarPnt <- IvarPnt + 1
         LBVar[IvarYr] <- fullrep[IvarPnt,2]
        }
       LBOut <- cbind(LBOut,LBVar)
     }
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
   } 

  write("\n#Legal 76 Biomass by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$LegalBio76[,Iarea])
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } 

  write("\n#Extended Legal Biomass by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    #Years <- 1:(GeneralSpecs$BurnIn+GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Years <- (GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)):(GeneralSpecs$Year2+GeneralSpecs$MaxProjYr)
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,length(Years)),Years,Report$LegalBioAll[,Iarea])
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } 
  ## Mtempts(Narea, Nage, BurnIn+Nyear+MaxProjYr+1)
  write("\n#Natural Mortality by Area, Age and Year",OutputFile,append=T) 
  for (Iage in 1:GeneralSpecs$Nage){
    for (Iarea in 1:GeneralSpecs$Narea){
    Years <- (GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)):(GeneralSpecs$Year2+GeneralSpecs$MaxProjYr)
    ## BurnIn+Nyear+MaxProjYr+1
    Ncol <- 4
    LBOut <- cbind(rep(Iarea,length(Years)),rep(Iage,length(Years)),Years,Report$Mtempts[Iarea,Iage,1:length(Years)])
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } }
  
  write("\n#Virgin Biomass",OutputFile,append=T) 
      write(Report$VirginBio,ncol=1,OutputFile,append=T)
  
  
  write("\n#Harvest rate by zone (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:ControlSpecs$Nzone)
   {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate[,Iarea])
    for (Ivar in 1:ControlSpecs$NvarTypes)
     if (ControlSpecs$VarTypes[Ivar]==3)
      {
       HRVar <- rep(0,GeneralSpecs$Nyear)
       Ncol <- Ncol + 1
       for (IvarYr in 1:GeneralSpecs$Nyear)
        {
         IvarPnt <- IvarPnt + 1
         HRVar[IvarYr] <- fullrep[IvarPnt,2]
        }
       HROut <- cbind(HROut,HRVar)
      }
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
   } 

  write("\n#Harvest 76 rate by zone (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:ControlSpecs$Nzone)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate76[,Iarea])
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
  } 
  
  
  #### ===================================================================================== 
  write("\n#Catches",OutputFile,append=T)
  write("#Fleet Area Year Step Observed Predicted Fishing_mortality",OutputFile,append=T)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
  {
    TotalCatch <- sum(TheData$Catch[,,Ifleet]) 
    if (TotalCatch > 0)
    {  
      Iarea = ControlSpecs$Fleet_area[Ifleet]
      for (Iyear in 1:GeneralSpecs$Nyear)
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          Summ <- c(Ifleet,Iarea+1,Iyear+GeneralSpecs$Year1-1,Istep,TheData$Catch[Iyear,Istep,Ifleet],Report$CatchCheck[Iyear,Istep,Ifleet],Report$Hrate[Iyear+max(GeneralSpecs$BurnIn),Istep,Ifleet])     
          write(Summ,ncol=7,OutputFile,append=T)
        } 
    }
  }  
  
  #### =====================================================================================   
  
  write("\n#Fishing Efficiency Commercial",OutputFile,append=T)
  write("#id Year Predicted",OutputFile,append=T)
  fcreep <- Report$CpueEcreep
  for (Iid in 1:ncol(fcreep)) {
    tmp <- cbind(rep(Iid, nrow(fcreep)), Data$Year1:(Data$Year1+nrow(fcreep)-1), fcreep[,Iid])
    write(t(tmp),ncol=3,OutputFile,append=T) 
    #for(wte in 1:nrow(fcreep)) write(tmp[wte,],ncol=3,OutputFile,append=T)  
      }  
  
  #### ===================================================================================== 
  write("\n#Index data",OutputFile,append=T)
  write("#Data_set Fleet Sex Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResCpue <- matrix(0,nrow=TheData$Ncpue,ncol=9)
  for (Ipnt in 1:TheData$Ncpue)  
   {
    ResCpue[Ipnt,] <- c(TheData$IndexI[Ipnt,]+c(1,1,0,GeneralSpecs$Year1,1),TheData$IndexR[Ipnt,],Report$PredCpue[Ipnt,])
   } 
  write(t(ResCpue),OutputFile,ncol=9,append=T)
  write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
   {
    IndexPoints <- which(TheData$FixSigmaCpue==IdataSet)-1
    if (length(IndexPoints) >0)
     {  
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaCpue[IdataSet+1],Report$CpueLikeComp[IdataSet+1],Report$LambdaCpue2[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }  
   }  
  write("\n#Data_set N Q",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
  {
    IndexPoints <- which(TheData$TreatQcpue==IdataSet)-1
    if (length(IndexPoints) >0)
    {  
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$CpueQ[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }  
  }  
    
  #### --------------------------------------------------------------------------------------
  write("\n#Numbers data",OutputFile,append=T)
  write("#Index Fleet Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResNumbers <- matrix(0,nrow=TheData$Nnumbers,ncol=8)
  if (TheData$Nnumbers > 0)
   {  
    for (Ipnt in 1:TheData$Nnumbers)  
     {
      ResNumbers[Ipnt,] <- c(TheData$NumbersI[Ipnt,]+c(1,1,GeneralSpecs$Year1-1,1),TheData$NumbersR[Ipnt,],Report$PredNumbers[Ipnt,])
     } 
    write(t(ResNumbers),OutputFile,ncol=8,append=T)
    write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
    for (IdataSet in (0:(Data$NcatchDataSeries-1)))
     {  
      Use <- TheData$NumbersI[,1] == IdataSet
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaNumbers[IdataSet+1],Report$NumbersLikeComp[IdataSet+1],Report$LambdaNumbers2[IdataSet+1])
      write(Summ,OutputFile,append=T)
     }
  }

  write("\n#Growth Curves",OutputFile,append=T)   ## (Nyear,Narea,Nsex,Nage,MaxLen)
  write(c("#year area sex age lengthbins"),OutputFile,append=T)
  Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
  tempout <- expand.grid(year=Years, area=1:GeneralSpecs$Narea, sex=1:GeneralSpecs$Nsex, age=1:GeneralSpecs$Nage)
  Mat <- matrix(NA, nrow=nrow(tempout), ncol=GeneralSpecs$Nlen[1])
  for (i in 1:nrow(tempout))   {
    GROut <- c(as.numeric(tempout[i,]),  Report$GrowthOut[which(Years==tempout$year[i]),tempout$area[i],tempout$sex[i],tempout$age[i],]  )
    write(GROut,ncol=length(GROut), OutputFile,append=T) }
  
  #### ===================================================================================== 
#   write("\n#Tagging data - Tag numbers by ",OutputFile,append=T)
#   write("#sex Relgrp area year step repSplit Obs Est",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)    
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (Iyear in 1:TagSpecs$NyearTags)
#       for (Istep in 1:GeneralSpecs$Nstep)
#        for (IrepSplit in 1:TagSpecs$NrepSplit)
#         if (TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep] > 0)
#          {
#           VAL1<- TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           VAL2 <- Report$RecapNum[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           Summ <- paste(Isex,Igrp,Iarea,Iyear+TagSpecs$Year1Tag[1]-1,Istep,IrepSplit,VAL1,VAL2)
#           write(Summ,OutputFile,append=T)
#          }
#     
#   write("\n#Tagging length data - Tag numbers by ",OutputFile,append=T)    
#   write("#sex Relgrp area repSplit size ObsTot ObsProp EstProp",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)    
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (IrepSplit in 1:TagSpecs$NrepSplit)
#       if (TagSpecs$FitTagSizes[IrepSplit] == 1)
#        {
#         ObsSS = 0
#         for (Iyear in 1:TagSpecs$NyearTags)
#          for (Istep in 1:GeneralSpecs$Nstep)
#           ObsSS = ObsSS + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1]
#         if (ObsSS > 0)
#          {
#           Vec <- rep(0,GeneralSpecs$Nlen[Isex])
#           for (Isize in 1:GeneralSpecs$Nlen[Isex])
#            {
#             ObsEE <- 0
#             for (Iyear in 1:TagSpecs$NyearTags)
# 	     for (Istep in 1:GeneralSpecs$Nstep)
# 	      ObsEE = ObsEE + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1+Isize]
# 	     Vec[Isize] <- ObsEE 
# 	     PredEE <- Report$PredTagSize[Isex,Igrp,Iarea,Isize]
#              Summ <- paste(Isex,Igrp,Iarea,IrepSplit,Isize,ObsSS,Vec[Isize],PredEE)
#              write(Summ,OutputFile,append=T)
#             }
#         
#          }
#        }   

  
    
  #### ===================================================================================== 
  write("\n#Length data",OutputFile,append=T)
  #TagSpecs 
  
  # Calculate Francis weights and McAllister-Ianelli mode
  write("#Fleet Sex Npnts Francis_Multiplier",OutputFile,append=T)
  EffN <- rep(0,TheData$NlenComp)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
   for (Isex in 1:GeneralSpecs$Nsex)  
     {
      Use <-  TheData$LenCompI[,1]+1 == Ifleet & TheData$LenCompI[,2]+1 == Isex
      Indexes <- c(1:TheData$NlenComp)[Use]
      Residuals <- NULL
      for (Index2 in 1:length(Indexes))
       {
        Index <- Indexes[Index2]
        Top <- 0; Bot <- 0
        Nlens <- GeneralSpecs$Nlen[Isex]
        Length <- GeneralSpecs$MidLenBin[Isex,1:Nlens]
        ObsProp <- TheData$LenCompR[Index,1:Nlens]/sum(TheData$LenCompR[Index,1:Nlens])
        PredProp <- Report$PredLengthComp[Index,1:Nlens]
        for (Ilen in 1:Nlens) 
         {
          Bot <- Bot + (ObsProp[Ilen] - PredProp[Ilen])^2
          Top <- Top + PredProp[Ilen]*(1.0-PredProp[Ilen])
         }  
        MeanOL <- sum(ObsProp*Length)
        MeanPL <- sum(PredProp*Length)
        SD <- sqrt((sum(PredProp*Length^2)-MeanPL^2)/TheData$Stage1W[Index])
        Residual <- (MeanOL-MeanPL)/SD
        Residuals <- c(Residuals,Residual)
        EffN[Index] <- Top/Bot
       }   
      if (!is.na(Residuals[1]))
       {  
        LenWghtMultipliers <- 1.0/var(Residuals)
        Summ <- c(Ifleet,Isex,length(Residuals),LenWghtMultipliers)
        write(Summ,OutputFile,append=T,ncol=5)
       }  
     }  
  write("\n#Obs/Pred Fleet Sex Year Step Nsamp EffN proportions",OutputFile,append=T)
  for (Ipnt in 1:TheData$NlenComp)
   {
    Ifleet <-  TheData$LenCompI[Ipnt,1]+1
    Isex <- TheData$LenCompI[Ipnt,2]+1
    Iyear <-  TheData$LenCompI[Ipnt,3]+GeneralSpecs$Year1-1
    Istep <- TheData$LenCompI[Ipnt,4]+1
    Obs <- c(round(as.vector(TheData$LenCompR[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("O",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Obs)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
    Pred <- c(round(as.vector(Report$PredLengthComp[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("P",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Pred)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
   }  

  #### --------------------------------------------------------------------------------------
  write("\n#Larval data",OutputFile,append=T)
  write("#Area Year Observed SD Predicted Residual",OutputFile,append=T)
  ResLarval <- matrix(0,nrow=TheData$NLarvalData,ncol=6)
  if (TheData$NLarvalData > 0)
   {
    for (Ipnt in 1:TheData$NLarvalData)  
     {
      ResLarval[Ipnt,] <- c(TheData$Lar_dataI[Ipnt,]+c(1,GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)),TheData$Lar_dataR[Ipnt,],Report$PredLarval[Ipnt,])
     } 
   # print(ResLarval)
    write(t(ResLarval),OutputFile,ncol=6,append=T)
   }


  #### =====================================================================================
  write("\n#Full Selectivity",OutputFile,append=T) 
  NselexPatterns <-length(Report$ActSelex[,1])
  for (Ipattern in 1:NselexPatterns)
   {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActSelex[Ipattern,],5)))) 
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
   }  
  write("\n#Retention",OutputFile,append=T) 
  NretenPatterns <-length(Report$ActReten[,1])
  for (Ipattern in 1:NretenPatterns)
   {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActReten[Ipattern,],5)))) 
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
   }  

  #### =====================================================================================
  write("\n#Recruitment",OutputFile,append=T) 

  write("#Recruitment patterns by timestep, area and sex ",OutputFile,append=T) 
  write("#year, tstep, sex area1, area2, area3, area4, area5, area6, area7, area8",OutputFile,append=T) 
  for (Iyear in 1:dim(Report$ActRecruitAreaSexDist)[1]){
    for (Istep in 1:dim(Report$ActRecruitAreaSexDist)[2]){
        for (Isex in 1:dim(Report$ActRecruitAreaSexDist)[4]){
     write(paste(Iyear,Istep,Isex,Report$ActRecruitAreaSexDist[Iyear,Istep, ,Isex]),OutputFile,append=T)
    }}}
  
  write("#Recruitment patterns by sex and size ",OutputFile,append=T) 
  write("#Pattern, sex, size ",OutputFile,append=T) 
  for (Ipattern in 1:dim(Report$ActRecruitLenDist)[1]){
   for (Isex in 1:dim(Report$ActRecruitLenDist)[2]){
     write(paste(Ipattern,Iarea,Isex,Report$ActRecruitLenDist[Ipattern,Isex,]),OutputFile,append=T)
    }}
  
  #### =====================================================================================
  write("\n#Movement",OutputFile,append=T) 
  
  write("#Movement patterns by size ",OutputFile,append=T) 
  write("#Pattern, size",OutputFile,append=T) 
  NmovePatterns <- length(Report$ActMove[,1])
  for (Ipattern in 1:NmovePatterns)
   write(paste(Ipattern,Report$ActMove[Ipattern,]),OutputFile,append=T)
  
  #### ===================================================================================== 
  write("\n#Initial N-matrix\n#Area\tSex\tAge",OutputFile,append=T) 
  Nout <-rep(0,3+GeneralSpecs$MaxLen)
  
  for (Iarea in 1:GeneralSpecs$Narea)
   for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))
     {       
      Nout[1:3]  <- c(Iarea,Isex,Iage)
      Nout[4:(3+GeneralSpecs$MaxLen)] <- round(Report$Ninit[Iarea,Isex,Iage,],2)
      write(t(Nout),OutputFile,append=T,ncol=3+GeneralSpecs$MaxLen)
     }  
  
  #### ===================================================================================== 
  write("\n#Simplified N-matrix",OutputFile,append=T) 
  ncol <- GeneralSpecs$MaxLen+4
  for (Iarea in 1:GeneralSpecs$Narea)
   for (Isex in 1:GeneralSpecs$Nsex)
    { 
     write("#Area Sex Year(raw) Year(act) Lengths",OutputFile,append=T,ncol=100)
     Nrow = (Nyears)
     Nout <- matrix(0,Nrow,4+GeneralSpecs$MaxLen)
     Ipnt <- 0
     for (Iyear in 1:Nyears)
      {
       Ipnt <- Ipnt + 1
       Nout[Ipnt,1:4]  <- c(Iarea,Isex,Iyear-max(GeneralSpecs$BurnIn)+1,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1)
       for (Iage in 1:(GeneralSpecs$Nage))
        Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))] <- Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))]+ round(Report$N[Iarea,Iyear,1,Isex,Iage,],2)
      }
     write(t(Nout),OutputFile,append=T,ncol=ncol)
   } 
  

  #### ===================================================================================== 
  write("\n#F-matrix",OutputFile,append=T) 
  write("#Year Time-step F_by_fleet",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)
  Nrow = GeneralSpecs$Nstep*Nyears
  Fout <- matrix(0,Nrow,2+GeneralSpecs$Nfleet)
  Ipnt <- 0
  for (Iyear in 1:Nyears)
    for (Istep in 1:GeneralSpecs$Nstep)
    {
     Ipnt <- Ipnt + 1
     Fout[Ipnt,1:2]  <- c(Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep-1)
     Fout[Ipnt,(3:(2+GeneralSpecs$Nfleet))] <- round(Report$Hrate[Iyear,Istep,],5)
    }
  write(t(Fout),OutputFile,append=T,ncol=GeneralSpecs$Nfleet+2)
  
  write("\n#Full N-matrix",OutputFile,append=T) 
  write("#Area Sex Age Year Time-step Lengths",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1
  Nrow = (GeneralSpecs$Narea*GeneralSpecs$Nsex*GeneralSpecs$Nage*Nyears*GeneralSpecs$Nstep)
  Nout <- matrix(NA,Nrow,5+GeneralSpecs$MaxLen)
  ncol <- GeneralSpecs$MaxLen+5
  Ipnt <- 0
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
      for (Iage in 1:GeneralSpecs$Nage)
       for (Iyear in 1:Nyears)
        for (Istep in 1:GeneralSpecs$Nstep)
         {
          Ipnt <- Ipnt + 1
          Nout[Ipnt,1:5]  <- c(Iarea,Isex,Iage,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep)
          Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- 0
          Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- Report$N[Iarea,Iyear,Istep,Isex,Iage,]
         }
  write(t(round(Nout,1)),OutputFile,append=T,ncol=ncol)
  setwd(dirname(getwd()))
  
}

WriteOutput.noSD <- function(Report,parameters,GeneralSpecs,ControlSpecs,TheData,CurrPhase=2,best=rep(-1,1000))
{
  OutputFile <- "Output.RL"
  ParFileName <- "Model.ppp"
  write("# Likelihood summary",OutputFile)
  write("Catch penalty",OutputFile,append=T)
  write(paste("Total objective function",Report$neglogL),OutputFile,append=T)
  write(paste("Catch likelihood",Report$CatchLike),OutputFile,append=T)
  write(paste("Cpue likelihood",Report$CpueLike, " (Weighted ", Report$Weighted_CpueLike, ")"),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLike, " (Weighted ", Report$Weighted_NumbersLike, ")"),OutputFile,append=T)
  write(paste("Length likelihood",Report$LengthLike, " (Weighted ", Report$Weighted_LengthLike, ")"),OutputFile,append=T)
  write(paste("Larval likelihood",Report$LarvalLike, " (Weighted ", Report$Weighted_LarvalLike, ")"),OutputFile,append=T)
 # write(paste("Tag1 Likelihood",sum(Report$TagLike1)," (Weighted ", Report$Weighted_TagLike1,")"),OutputFile,append=T)
 #write(paste("Tag2 Likelihood",sum(Report$TagLike2)," (Weighted ", Report$Weighted_TagLike2,")"),OutputFile,append=T)
  write(paste("Initial N penalty",Report$Initial_pen),OutputFile,append=T)
  write(paste("Recruitment penalty",Report$Rec_Penal),OutputFile,append=T)
  write(paste("Recruitment smooth penalty",Report$Rec_Penal_Smooth),OutputFile,append=T)
  
  write("\n# Likelihood by fleet",OutputFile,append=T)
  write(paste("Cpue likelihood",paste(Report$CpueLikeComp[])),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLikeComp),OutputFile,append=T)
  write(paste("Length likeliood",Report$LengthLikeComps),OutputFile,append=T)
  write(paste("Larval likeliood",Report$LarvalLikeComps),OutputFile,append=T)
  
  write("\n# Likelihood by group",OutputFile,append=T)
#  write(paste("Tag Likelihood 1",Report$TagLike1),OutputFile,append=T)
#  write(paste("Tag Likelihood 2",Report$TagLike2),OutputFile,append=T)

  write("\n# parameter table",OutputFile,append=T)
  write("# Parameter Par_cnt Estpar_cnt Estimate SD lwrBound uprBound",OutputFile,append=T)
  write("# parameters",ParFileName)
  write("# dummy",ParFileName,append=T)
  write("0",ParFileName,append=T)
  ParName <- names(parameters)
  Ipnt <- 0; Iqnt <- 0
  for (ParName in names(parameters))
  {
    if (ParName != "dummy") 
    {
       print(ParName)
       ThePar <- InitialVars[[ParName]]
      for (Ipar in 1:length(ThePar$Initial))
      {  
        Iqnt <- Iqnt + 1
        write(paste("#",ParName,"_",Ipar," ",Iqnt," ",sep=""),ParFileName,append=T)
        if (ThePar$Phase[Ipar] > 0 & ThePar$Phase[Ipar] <= CurrPhase)  
        {
          Ipnt <- Ipnt + 1;
          xx <- paste(ParName,"_",Ipar," ",Iqnt," " ,Ipnt," ", best[Ipnt]," ","stdrep[Ipnt,2]"," ",ThePar$Bnd[Ipar,1]," ",ThePar$Bnd[Ipar,2],sep="")
          write(best[Ipnt],ParFileName,append=T)
        }  
        else
        {  
          xx <- paste(ParName,"_",Ipar," ",Iqnt," NA ", ThePar$Initial[Ipar],sep="")
          write(ThePar$Initial[Ipar],ParFileName,append=T)
        }  
        write(xx,OutputFile,append=T)
      }   
    } # Parameters within ParName
  } # ParNames
  write(paste("#Total estimated parameters:",Ipnt),OutputFile,append=T)  
  
  IvarPnt <- Ipnt + 1
  
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1
  
  #### ===================================================================================== 
  write("\n#Egg Production",OutputFile,append=T) 
  write("#Year Total By_area (SD?)",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)
  Ncol <- 1+1+GeneralSpecs$Narea
  SSBOut <- matrix(0,nrow=Nyears,ncol=Ncol)
  SSBOut[,1] <-GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1+1:Nyears
  SSBOut[,2] <- Report$MatBio[1:Nyears]
  for (Iarea in 1:GeneralSpecs$Narea)
    SSBOut[,2+Iarea] <- Report$MatBioArea[Iarea,1:Nyears]  
  for (Ivar in 1:ControlSpecs$NvarTypes)
    if (ControlSpecs$VarTypes[Ivar]==1)
    {
      SSBVar <- rep(0,Nyears)
      Ncol <- Ncol + 1
      for (IvarYr in 1:Nyears)
      {
        IvarPnt <- IvarPnt + 1
        SSBVar[IvarYr] <- "?"#fullrep[IvarPnt,2]
      }
      SSBOut <- cbind(SSBOut,SSBVar)
      IvarPnt <- IvarPnt+1
    }
  
  write(t(SSBOut),ncol=Ncol,OutputFile,append=T)
  
  write("\n#Total Recruitment",OutputFile,append=T) 
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  write(t(cbind(Years,Report$Recruits)),ncol=2,OutputFile,append=T)
  
  write("\n#Recruitment by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1+1 
    Ncol <- 3
    RecOut <- cbind(rep(Iarea,Nyears+1),Years,Report$RecruitmentByArea[Iarea,])
    for (Ivar in 1:ControlSpecs$NvarTypes)
      if (ControlSpecs$VarTypes[Ivar]==2)
      {
        RecVar <- rep(0,Nyears+1)
        Ncol <- Ncol + 1
        for (IvarYr in 1:(Nyears+1))
        {
          IvarPnt <- IvarPnt + 1
          RecVar[IvarYr] <- "?"#fullrep[IvarPnt,2]
        }
        RecOut <- cbind(RecOut,RecVar)
      }
    write(t(RecOut),ncol=Ncol,OutputFile,append=T)
  } 
  
  write("\n#Legal Biomass by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$LegalBio[,Iarea])
    for (Ivar in 1:ControlSpecs$NvarTypes)
      if (ControlSpecs$VarTypes[Ivar]==3)
      {
        LBVar <- rep(0,GeneralSpecs$Nyear)
        Ncol <- Ncol + 1
        for (IvarYr in 1:GeneralSpecs$Nyear)
        {
          IvarPnt <- IvarPnt + 1
          LBVar[IvarYr] <- "?"#fullrep[IvarPnt,2]
        }
        LBOut <- cbind(LBOut,LBVar)
      }
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } 
  
  write("\n#Legal 76 Biomass by area (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$LegalBio76[,Iarea])
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } 
  
  write("\n#Harvest rate by zone (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:ControlSpecs$Nzone)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate[,Iarea])
    for (Ivar in 1:ControlSpecs$NvarTypes)
      if (ControlSpecs$VarTypes[Ivar]==3)
      {
        HRVar <- rep(0,GeneralSpecs$Nyear)
        Ncol <- Ncol + 1
        for (IvarYr in 1:GeneralSpecs$Nyear)
        {
          IvarPnt <- IvarPnt + 1
          HRVar[IvarYr] <- "?"#fullrep[IvarPnt,2]
        }
        HROut <- cbind(HROut,HRVar)
      }
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
  } 
  
  write("\n#Harvest 76 rate by zone (1+SD?)",OutputFile,append=T) 
  for (Iarea in 1:ControlSpecs$Nzone)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate76[,Iarea])
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
  } 
  
  
  #### ===================================================================================== 
  write("\n#Catches",OutputFile,append=T)
  write("#Fleet Area Year Step Observed Predicted Fishing_mortality",OutputFile,append=T)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
  {
    TotalCatch <- sum(TheData$Catch[,,Ifleet]) 
    if (TotalCatch > 0)
    {  
      Iarea = ControlSpecs$Fleet_area[Ifleet]
      for (Iyear in 1:GeneralSpecs$Nyear)
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          Summ <- c(Ifleet,Iarea+1,Iyear+GeneralSpecs$Year1-1,Istep,TheData$Catch[Iyear,Istep,Ifleet],Report$CatchCheck[Iyear,Istep,Ifleet],Report$Hrate[Iyear+max(GeneralSpecs$BurnIn),Istep,Ifleet])     
          write(Summ,ncol=7,OutputFile,append=T)
        } 
    }
  }  
  
  #### ===================================================================================== 
  write("\n#Index data",OutputFile,append=T)
  write("#Data_set Fleet Sex Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResCpue <- matrix(0,nrow=TheData$Ncpue,ncol=9)
  for (Ipnt in 1:TheData$Ncpue)  
  {
    ResCpue[Ipnt,] <- c(TheData$IndexI[Ipnt,]+c(1,1,0,GeneralSpecs$Year1,1),TheData$IndexR[Ipnt,],Report$PredCpue[Ipnt,])
  } 
  write(t(ResCpue),OutputFile,ncol=9,append=T)
  write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
  {
    IndexPoints <- which(TheData$FixSigmaCpue==IdataSet)-1
    if (length(IndexPoints) >0)
    {  
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaCpue[IdataSet+1],Report$CpueLikeComp[IdataSet+1],Report$LambdaCpue2[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }  
  }  
  write("\n#Data_set N Q",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
  {
    IndexPoints <- which(TheData$TreatQcpue==IdataSet)-1
    if (length(IndexPoints) >0)
    {  
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$CpueQ[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }  
  }  
  
  #### --------------------------------------------------------------------------------------
  write("\n#Numbers data",OutputFile,append=T)
  write("#Index Fleet Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResNumbers <- matrix(0,nrow=TheData$Nnumbers,ncol=8)
  if (TheData$Nnumbers > 0)
  {  
    for (Ipnt in 1:TheData$Nnumbers)  
    {
      ResNumbers[Ipnt,] <- c(TheData$NumbersI[Ipnt,]+c(1,1,GeneralSpecs$Year1-1,1),TheData$NumbersR[Ipnt,],Report$PredNumbers[Ipnt,])
    } 
    write(t(ResNumbers),OutputFile,ncol=8,append=T)
    write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
    for (IdataSet in (0:(Data$NcatchDataSeries-1)))
    {  
      Use <- TheData$NumbersI[,1] == IdataSet
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaNumbers[IdataSet+1],Report$NumbersLikeComp[IdataSet+1],Report$LambdaNumbers2[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }
  }
   #### ===================================================================================== 
#   write("\n#Tagging data - Tag numbers by ",OutputFile,append=T)
#   write("#sex Relgrp area year step repSplit Obs Est",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)    
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (Iyear in 1:TagSpecs$NyearTags)
#       for (Istep in 1:GeneralSpecs$Nstep)
#        for (IrepSplit in 1:TagSpecs$NrepSplit)
#         if (TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep] > 0)
#          {
#           VAL1<- TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           VAL2 <- Report$RecapNum[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           Summ <- paste(Isex,Igrp,Iarea,Iyear+TagSpecs$Year1Tag[1]-1,Istep,IrepSplit,VAL1,VAL2)
#           write(Summ,OutputFile,append=T)
#          }
#     
#   write("\n#Tagging length data - Tag numbers by ",OutputFile,append=T)    
#   write("#sex Relgrp area repSplit size ObsTot ObsProp EstProp",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)    
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (IrepSplit in 1:TagSpecs$NrepSplit)
#       if (TagSpecs$FitTagSizes[IrepSplit] == 1)
#        {
#         ObsSS = 0
#         for (Iyear in 1:TagSpecs$NyearTags)
#          for (Istep in 1:GeneralSpecs$Nstep)
#           ObsSS = ObsSS + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1]
#         if (ObsSS > 0)
#          {
#           Vec <- rep(0,GeneralSpecs$Nlen[Isex])
#           for (Isize in 1:GeneralSpecs$Nlen[Isex])
#            {
#             ObsEE <- 0
#             for (Iyear in 1:TagSpecs$NyearTags)
# 	     for (Istep in 1:GeneralSpecs$Nstep)
# 	      ObsEE = ObsEE + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1+Isize]
# 	     Vec[Isize] <- ObsEE 
# 	     PredEE <- Report$PredTagSize[Isex,Igrp,Iarea,Isize]
#              Summ <- paste(Isex,Igrp,Iarea,IrepSplit,Isize,ObsSS,Vec[Isize],PredEE)
#              write(Summ,OutputFile,append=T)
#             }
#         
#          }
#        }   

  
 
  #### ===================================================================================== 
  write("\n#Length data",OutputFile,append=T)
  
  # Calculate Francis weights and McAllister-Ianelli mode
  write("#Fleet Sex Npnts Francis_Multiplier",OutputFile,append=T)
  EffN <- rep(0,TheData$NlenComp)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Isex in 1:GeneralSpecs$Nsex)  
    {
      Use <-  TheData$LenCompI[,1]+1 == Ifleet & TheData$LenCompI[,2]+1 == Isex
      Indexes <- c(1:TheData$NlenComp)[Use]
      Residuals <- NULL
      for (Index2 in 1:length(Indexes))
      {
        Index <- Indexes[Index2]
        Top <- 0; Bot <- 0
        Nlens <- GeneralSpecs$Nlen[Isex]
        Length <- GeneralSpecs$MidLenBin[Isex,1:Nlens]
        ObsProp <- TheData$LenCompR[Index,1:Nlens]/sum(TheData$LenCompR[Index,1:Nlens])
        PredProp <- Report$PredLengthComp[Index,1:Nlens]
        for (Ilen in 1:Nlens) 
        {
          Bot <- Bot + (ObsProp[Ilen] - PredProp[Ilen])^2
          Top <- Top + PredProp[Ilen]*(1.0-PredProp[Ilen])
        }  
        MeanOL <- sum(ObsProp*Length)
        MeanPL <- sum(PredProp*Length)
        SD <- sqrt((sum(PredProp*Length^2)-MeanPL^2)/TheData$Stage1W[Index])
        Residual <- (MeanOL-MeanPL)/SD
        Residuals <- c(Residuals,Residual)
        EffN[Index] <- Top/Bot
      }   
      if (!is.na(Residuals[1]))
      {  
        LenWghtMultipliers <- 1.0/var(Residuals)
        Summ <- c(Ifleet,Isex,length(Residuals),LenWghtMultipliers)
        write(Summ,OutputFile,append=T,ncol=5)
      }  
    }  
  write("\n#Obs/Pred Fleet Sex Year Step Nsamp EffN proportions",OutputFile,append=T)
  for (Ipnt in 1:TheData$NlenComp)
  {
    Ifleet <-  TheData$LenCompI[Ipnt,1]+1
    Isex <- TheData$LenCompI[Ipnt,2]+1
    Iyear <-  TheData$LenCompI[Ipnt,3]+GeneralSpecs$Year1-1
    Istep <- TheData$LenCompI[Ipnt,4]+1
    Obs <- c(round(as.vector(TheData$LenCompR[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("O",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Obs)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
    Pred <- c(round(as.vector(Report$PredLengthComp[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("P",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Pred)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
  }  
  
  #### --------------------------------------------------------------------------------------
  write("\n#Larval data",OutputFile,append=T)
  write("#Area Year Observed SD Predicted Residual",OutputFile,append=T)
  ResLarval <- matrix(0,nrow=TheData$NLarvalData,ncol=6)
  if (TheData$NLarvalData > 0)
  {
    for (Ipnt in 1:TheData$NLarvalData)  
    {
      ResLarval[Ipnt,] <- c(TheData$Lar_dataI[Ipnt,]+c(1,GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)),TheData$Lar_dataR[Ipnt,],Report$PredLarval[Ipnt,])
    } 
    # print(ResLarval)
    write(t(ResLarval),OutputFile,ncol=6,append=T)
  }
  
  
  #### =====================================================================================
  write("\n#Full Selectivity",OutputFile,append=T) 
  NselexPatterns <-length(Report$ActSelex[,1])
  for (Ipattern in 1:NselexPatterns)
  {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActSelex[Ipattern,],5)))) 
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
  }  
  write("\n#Retention",OutputFile,append=T) 
  NretenPatterns <-length(Report$ActReten[,1])
  for (Ipattern in 1:NretenPatterns)
  {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActReten[Ipattern,],5)))) 
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
  }  
  
  #### =====================================================================================
  write("\n#Recruitment",OutputFile,append=T) 
  
  write("#Recruitment patterns by timestep, area and sex ",OutputFile,append=T) 
  write("#year, tstep, sex area1, area2, area3, area4, area5, area6, area7, area8",OutputFile,append=T) 
  for (Iyear in 1:dim(Report$ActRecruitAreaSexDist)[1]){
    for (Istep in 1:dim(Report$ActRecruitAreaSexDist)[2]){
      for (Isex in 1:dim(Report$ActRecruitAreaSexDist)[4]){
        write(paste(Iyear,Istep,Isex,Report$ActRecruitAreaSexDist[Iyear,Istep, ,Isex]),OutputFile,append=T)
      }}}
  
  write("#Recruitment patterns by sex and size ",OutputFile,append=T) 
  write("#Pattern, sex, size ",OutputFile,append=T) 
  for (Ipattern in 1:dim(Report$ActRecruitLenDist)[1]){
    for (Isex in 1:dim(Report$ActRecruitLenDist)[2]){
      write(paste(Ipattern,Iarea,Isex,Report$ActRecruitLenDist[Ipattern,Isex,]),OutputFile,append=T)
    }}
  
  #### =====================================================================================
  write("\n#Movement",OutputFile,append=T) 
  
  write("#Movement patterns by size ",OutputFile,append=T) 
  write("#Pattern, size",OutputFile,append=T) 
  NmovePatterns <- length(Report$ActMove[,1])
  for (Ipattern in 1:NmovePatterns)
    write(paste(Ipattern,Report$ActMove[Ipattern,]),OutputFile,append=T)
  
  #### ===================================================================================== 
  write("\n#Initial N-matrix\n#Area\tSex\tAge",OutputFile,append=T) 
  Nout <-rep(0,3+GeneralSpecs$MaxLen)
  
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
      for (Iage in 1:(GeneralSpecs$Nage))
      {       
        Nout[1:3]  <- c(Iarea,Isex,Iage)
        Nout[4:(3+GeneralSpecs$MaxLen)] <- round(Report$Ninit[Iarea,Isex,Iage,],2)
        write(t(Nout),OutputFile,append=T,ncol=3+GeneralSpecs$MaxLen)
      }  
  
  #### ===================================================================================== 
  write("\n#Simplified N-matrix",OutputFile,append=T) 
  ncol <- GeneralSpecs$MaxLen+4
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
    { 
      write("#Area Sex Year(raw) Year(act) Lengths",OutputFile,append=T,ncol=100)
      Nrow = (Nyears)
      Nout <- matrix(0,Nrow,4+GeneralSpecs$MaxLen)
      Ipnt <- 0
      for (Iyear in 1:Nyears)
      {
        Ipnt <- Ipnt + 1
        Nout[Ipnt,1:4]  <- c(Iarea,Isex,Iyear-max(GeneralSpecs$BurnIn)+1,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1)
        for (Iage in 1:(GeneralSpecs$Nage))
          Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))] <- Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))]+ round(Report$N[Iarea,Iyear,1,Isex,Iage,],2)
      }
      write(t(Nout),OutputFile,append=T,ncol=ncol)
    } 
  
  
  #### ===================================================================================== 
  write("\n#F-matrix",OutputFile,append=T) 
  write("#Year Time-step F_by_fleet",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)
  Nrow = GeneralSpecs$Nstep*Nyears
  Fout <- matrix(0,Nrow,2+GeneralSpecs$Nfleet)
  Ipnt <- 0
  for (Iyear in 1:Nyears)
    for (Istep in 1:GeneralSpecs$Nstep)
    {
      Ipnt <- Ipnt + 1
      Fout[Ipnt,1:2]  <- c(Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep-1)
      Fout[Ipnt,(3:(2+GeneralSpecs$Nfleet))] <- round(Report$Hrate[Iyear,Istep,],5)
    }
  write(t(Fout),OutputFile,append=T,ncol=GeneralSpecs$Nfleet+2)
  
  write("\n#Full N-matrix",OutputFile,append=T) 
  write("#Area Sex Age Year Time-step Lengths",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1
  Nrow = (GeneralSpecs$Narea*GeneralSpecs$Nsex*GeneralSpecs$Nage*Nyears*GeneralSpecs$Nstep)
  Nout <- matrix(NA,Nrow,5+GeneralSpecs$MaxLen)
  ncol <- GeneralSpecs$MaxLen+5
  Ipnt <- 0
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
      for (Iage in 1:GeneralSpecs$Nage)
        for (Iyear in 1:Nyears)
          for (Istep in 1:GeneralSpecs$Nstep)
          {
            Ipnt <- Ipnt + 1
            Nout[Ipnt,1:5]  <- c(Iarea,Isex,Iage,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep)
            Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- 0
            Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- Report$N[Iarea,Iyear,Istep,Isex,Iage,]
          }
  write(t(round(Nout,1)),OutputFile,append=T,ncol=ncol)
  setwd(dirname(getwd()))
  
  }

  