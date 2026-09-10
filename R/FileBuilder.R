#' Build IMuLT Model Input Files from Excel Template
#'
#' Creates all required input files (.DAT files) for running the IMuLT stock assessment
#' model. Reads model structure and parameters from a ModelStructure.xlsx workbook and
#' generates a complete set of properly formatted input files in a new run directory.
#'
#' @details
#' The function requires a ModelStructure.xlsx file in the current working directory.
#' This Excel workbook must contain the following sheets:
#' \itemize{
#'   \item Dynamics - Model dynamics and temporal settings
#'   \item Area - Spatial structure and management zones
#'   \item Times - Time step specifications
#'   \item Fleetcode - Fleet definitions
#'   \item Catch - Catch data by year, time step, and fleet
#'   \item CPUE - Catch-per-unit-effort indices
#'   \item LengthFreq - Length frequency data
#'   \item Puerulus - Larval settlement index
#'   \item Growth, lengthweight, maturity, etc. - Biological parameters
#' }
#'
#' The function creates the following output files in the run directory:
#' \itemize{
#'   \item STARTER.DAT - File paths and run settings
#'   \item DATA.DAT - Catch, CPUE, and biological data
#'   \item CONTROL.DAT - Model control parameters and biological specifications
#'   \item GROWTHSPEC.DAT - Growth transition matrices and parameters
#'   \item MOVESPEC.DAT - Movement/migration specifications
#'   \item RECRUITSPEC.DAT - Recruitment specifications and parameters
#'   \item RETAINSPEC.DAT - Retention (high-grading) specifications
#'   \item SELEXSPEC.DAT - Gear selectivity and legal size specifications
#'   \item PROJECTIONS.DAT - Future projection settings
#' }
#'
#' @return NULL (invisibly). Creates a new directory containing all model input files.
#'   The directory name follows the pattern: `[nAreas]`Area`[nAges]`AgeRun`[startYY]`_`[endYY]`
#'
#' @examples
#' \dontrun{
#' # Ensure ModelStructure.xlsx is in current directory
#' setwd("path/to/model/directory")
#'
#' # Build all input files
#' BuildInputFiles()
#'
#' # Files will be created in a new directory like: "8Area8AgeRun92_23"
#' }
#'
#' @seealso
#' See the IMuLT User Guide for ModelStructure.xlsx template format and specifications.
#'
#' @export
BuildInputFiles <- function(end_override = NULL){

  # Ensure the wd is set to the same location as ModelStructure.xls
  MSdir <- gsub('/ModelStructure.xlsx','',find_model_file())
  if(getwd()!=MSdir) { message("Changing working directory back to .../ModelStructure.xls")
    setwd(MSdir)
  }

  suppressPackageStartupMessages({
    library(dplyr, quietly = T)
    library(magrittr, quietly = T)
    library(openxlsx, quietly = T)})

  ## Make function that adjusts sex definations loaded through the excel file.
  adjsex <- function(x,nsex,section='CPUE'){
    if(nsex==1) { xout <- rep(1, length(x)) }
    if(nsex==2) { xout <- ifelse(toupper(x)%in%c('B','C'),0,ifelse(toupper(x)=='F',1,ifelse(toupper(x)=='M',2,x)))  }
    invalid <- unique(setdiff(x, c('B','C','F','M')))
    if(any(!x[!is.na(x)] %in% c('B','C','F','M'))) {
      warning("Unusual sex definition in ",section,":", paste(invalid, collapse=', '), '.\n',
              call. = FALSE)   }
    if(is.numeric(x)){ xout <- x }
    return(as.numeric(xout))
  }

  # safe_loadWorkbook <- function(file) {
  #   tryCatch(
  #     suppressWarnings(loadWorkbook(file = file)),
  #     error = function(e) {
  #       stop("'", file, "' appears to be open in Excel. Close it first and retry.", call. = FALSE)
  #     }
  #   )
  # }
  #
  # ## Open up file with all info
  # wb <- safe_loadWorkbook("ModelStructure.xlsx")

  wb <- load_model_structure()

  #This is the location of the data input files and their associated parameters
  dynamics <- readWorkbook(wb,sheet='Dynamics', startRow = 2)
  if(length(dynamics$value[dynamics$object=='estimateVariance'])>0){
    Varspos <- as.numeric(strsplit(dynamics$value[dynamics$object=='estimateVariance'], ",")[[1]])
    dynamics <- dynamics[dynamics$object!='estimateVariance', ]
    dynamics$value <- as.numeric(dynamics$value)
  }
  startseason <- as.numeric(dynamics$value[dynamics$object=='startseason'])
  endseason <- as.numeric(dynamics$value[dynamics$object=='endseason'])
  if (!is.null(end_override)) endseason <- end_override
  projectseason <- as.numeric(dynamics$value[dynamics$object=='projectedseason'])
  projectcatch <- as.numeric(dynamics$value[dynamics$object=='projectedcatch'])
  burnin <- as.numeric(dynamics$value[dynamics$object=='burnin'])
  ages <- as.numeric(dynamics$value[dynamics$object=='ages'])
  sexs <- 0:(as.numeric(dynamics$value[dynamics$object=='sexs'])-1)
  nsex <- length(sexs)
  areas <- readWorkbook(wb,sheet='Area', startRow = 2)
  times <- readWorkbook(wb,sheet='Times', startRow = 2)
  fleets <- readWorkbook(wb,sheet='Fleetcode', startRow = 2)
  effic <- readWorkbook(wb,sheet='EfficiencyCreep', startRow = 2)
  migrate <- readWorkbook(wb,sheet='Migrate', startRow = 2)
  zones <- length(unique(areas$ManageZone))
  area <- areas %>% group_by(ManageZone) %>% reframe(newarea=unique(AreaCode))
  (zoneareas <- split(area$newarea, area$ManageZone))
  lens <- seq(dynamics$value[dynamics$object=='lblwr'],dynamics$value[dynamics$object=='lbupr'],dynamics$value[dynamics$object=='lbgap'])+1

  ## Create a new folder for the model if one does not exist
  (files <- list.files(pattern = 'AgeRun'))
  (nfile <- paste(length(unique(area$newarea)),'Area',ages,'AgeRun',substr(startseason,3,4),"_",substr(endseason,3,4),sep=''))
  if(!(nfile%in%files)) {
    dir.create(paste(nfile,sep='')) ;  dir.create(paste(nfile,"/Output",sep=''))
  }
  (fls <- nfile)
  floc <- paste(getwd(),fls,sep='/')  # location of data files
  print(paste("Making new directory:: ", floc))
  print("#################################################################################################")

  #### Starter Filer ####
  ## Make Starter file
  print("Building Starter File")

  tmp <- list()
  tmp <- c(tmp, "DATA.DAT                              # General specifications file \n")
  tmp <- c(tmp, "CONTROL.DAT                              # Control file\n")
  tmp <- c(tmp, "SELEXSPEC.DAT                        # Specifications for selectivity \n")
  tmp <- c(tmp, "RETAINSPEC.DAT                       # Specifications for retention \n")
  tmp <- c(tmp, "RECRUITSPEC.DAT                      # Specifications for recruitment \n")
  tmp <- c(tmp, "REPOSPEC.DAT                         # Specifications for maturity \n")
  tmp <- c(tmp, "GROWTHSPEC.DAT                       # Specifications for growth \n")
  tmp <- c(tmp, "MOVESPEC.DAT                         # Specifications in movement \n")
  tmp <- c(tmp, "TAGSPEC.DAT                          # Specifications for tags \n")
  tmp <- c(tmp, "TAGPROP.DAT                          # Proportions \n")
  tmp <- c(tmp, "PROJECTIONS.DAT                      # Projections file \n\n\n")
  tmp <- c(tmp, "1                                    # Stop after this phase \n")

  write.table(tmp, paste(floc,'/STARTER.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)


  #### Data File ####
  print("Building Data File")

  lens <- seq(dynamics$value[dynamics$object=='lblwr'],dynamics$value[dynamics$object=='lbupr'],dynamics$value[dynamics$object=='lbgap'])
  lensPlus1 <- c(lens, (lens[2]-lens[1])+tail(lens,1))

  ## Make data input of new RL model
  tmp <- list()
  tmp <- c(tmp, "# Data Set", "\n", "# First year of the assessment")

  #get latest catch and effort data
  dat <- readWorkbook(wb,sheet='Catch', startRow = 2)
  dat %<>% filter(year%in%startseason:endseason) %>% arrange(year,step,fleet) ## ensure matches the model structure
  tmp <- c(tmp,"\n", startseason, "\n# Last year of the assessment\n", endseason,"\n# Maximum projection years\n", projectseason)

  bin <- areas %>% group_by(AreaCode) %>% summarise(av=floor(mean(burn_in)))

  tmp <- c(tmp, "\n# Burn-in whole model\n", burnin,
           "\n# Burn-in for F by area\n", paste(bin$av,collapse=' '),
           "\n# Time steps per year\n", max(times$tstep),
           "\n# Number of areas of data included in the file\n", max(areas$AreaCode ),
           "\n# Number of sexes\n",length(unique(sexs)),
           "\n# Number of ages\n", ages,
           "\n# Number of fleets\n",max(fleets$fleet))

  ## Lower Length bins
  femaleLB <- maleLB <- lensPlus1
  if(max(sexs)==1) tmp <- c(tmp, "\n# Number of size-classes (females then males)\n", length(lensPlus1)-1, " ", length(lensPlus1)-1)
  if(max(sexs)==0) tmp <- c(tmp, "\n# Number of size-classes (one sex)\n", length(lensPlus1)-1)

  # Timesteps
  tstmp <- times %>% group_by(tstep) %>% mutate(prop=length(month)/nrow(times)) %>% summarise(prop=mean(prop))
  tmp <- c(tmp, "\n# The Time steps\n# \n", paste(round(tstmp$prop,5), collapse = "\t"))

  # tmp <- c(tmp, "\n# Loop counter for initial conditions\n", 10,"\n# Years over which to tune (one per area)\n", dynamics$value[dynamics$object=='tune_years'])

  tmp <- c(tmp, "\n# Loop counter for initial conditions\n", 10,"\n# Years over which to tune (one per area)\n", paste(areas$tune_years,collapse=" "))

  # Length Bins cont.
  if(max(sexs)==0) tmp <- c(tmp, "\n# Lower Length Bins (one more than number of size-classes)\n", paste(lensPlus1, collapse = "\t"))
  if(max(sexs)==1) tmp <- c(tmp, "\n# Lower Length Bins (one more than number of size-classes)\n", paste(lensPlus1, collapse = "\t"), "\n", paste(lensPlus1, collapse = "\t"))

  ## Catch data
  tmp <- c(tmp, "\n# Catch data (kg) - Number of observations\n", nrow(dat), "\n#Year\tstep\tfleet\tcatch\n")

  for(i in 1:nrow(dat)){tmp <- c(tmp, paste(dat[i,], collapse = "\t"),"\n")}

  ##  Catch Rate Indices / CPUE - ensure that a cutfof does not leave just one obs!
  Udat <- readWorkbook(wb,sheet='CPUE', startRow = 2) %>% filter(Year%in%startseason:endseason) %>% group_by(Fleet) %>% mutate(nobs=length(unique(Year))) %>% filter(nobs>1) %>% select(-nobs) %>% mutate(CpueInd=as.factor(as.character(CpueInd)), CpueInd=as.numeric(CpueInd)) %>% ungroup() %>% mutate( CpueInd=CpueInd-min(CpueInd)) %>% arrange(CpueInd)
  Udat %<>% mutate(Sex=adjsex(Sex, nsex,section='CPUE'))
  SigmaCpueCeiling <- EstimateCpueSigmaCeiling(Udat)
  cpuenumbers <- unique(Udat$CpueInd)

  ## Comm=1-16, IBSSa2=17, IBSSa4=18, IBSSa5=19, IBSSa6=20, IBSSa8=21, ISSa1=5, ISSa3=6, ISSa5=7, ISSa7=8
  tmp <- c(tmp, "\n# Index data \n#Number of cpue datasets\n", length(cpuenumbers))
  tmp <- c(tmp, "\n# Type of index (1=weight;2=numbers)\n", paste(c(rep(1,length(cpuenumbers))), collapse=" "))
  tmp <- c(tmp, "\n# Treatment of sigma (unique value represents a unique SS for the series)\n", paste((1:length(cpuenumbers))-1, collapse=" "))  ## Fixsigma
  tmp <- c(tmp, "\n# Treatment of q (a value represents a unique q for that series)\n", paste(cpuenumbers, collapse=" "))
  tmp <- c(tmp, "\n# Environmental Index (value points to index, 0 = no index)\n", paste(rep(0,100)[1:length(cpuenumbers)], collapse=" "))
  ## Efficiency creep
  tmp1 <- Udat %>% group_by(CpueInd) %>% summarise(mn=mean(Fleet))
  effcreep <- fleets$effic.creep[match(tmp1$mn,fleets$fleet)]
  tmp <- c(tmp, "\n# Efficiency creep (value points to index, 0 = no index)\n", paste(effcreep, collapse=" "))
  tmp <- c(tmp, "\n# Efficiency creep year lag (each par compounds for this many years until next par starts\n", paste(effic$temporal.cover, collapse=" "))
  tmp <- c(tmp, "\n# Minimum sigma\n", 0.05)
  tmp <- c(tmp, "\n# Maximum sigma\n", round(SigmaCpueCeiling, 3))
  #Size of cpue data
  tmp <- c(tmp,"\n# The cpue data\n", nrow(Udat))

  #CpueInd	#Fleet	#Sex	#Year	#Step	#Index	#CV
  tmp <- c(tmp, "\n#CpueInd Fleet\tSex\tYear\tStep\tIndex\tCV\n")
  for(i in 1:nrow(Udat)){ tmp <- c(tmp, paste(Udat[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp,"\n# Numbers data \n# Number of catch data sets\n",0)
  tmp <- c(tmp,"\n# Treatment of catch in numbers series\n#", 0, "\n# Minimum sigma\n",0.05)

  ## Build the numbers data
  tmp <- c(tmp,'\n# The numbers data\n',0)
  tmp <- c(tmp,'\n#Group  Fleet  Year  Step  Catch  CV\n')

  ## Get the length data
  len <- readWorkbook(wb,sheet='LengthFreq', startRow = 2) %>% filter(Season%in%startseason:endseason)
  len %<>% mutate(Sex=adjsex(Sex, nsex,section='Length Freq.'))
  if(min(len$Sex) == 0) { len %<>% mutate(Sex = Sex + 1) }
  tmp <- c(tmp,"\n# Length compostion\n", nrow(len), '\n#Fleet\tSex\tSEASON\ttstep\tInd\t',paste(lensPlus1, collapse = "\t"),'\n')
  for(i in 1:nrow(len)){ tmp <- c(tmp, paste(len[i,], collapse = "\t"),"\n")}

  #Puerulus index
  tmp <- c(tmp,"\n# Larval index (puerulus)\n# Likelihood for larval (puerulus) data (0=lognormal, else normal\n",0)
  tmp <- c(tmp,"\n# Delay from puerulus to entering the model (years)\n",3)
  puer <- readWorkbook(wb,sheet='Puerulus',startRow = 2)  %>%  mutate(area=area,mn=round(mn,1), sd=round(cv*mn,2))  %>% dplyr::select(area,season, mn, sd) %>% filter(season %in% startseason:endseason, area %in% areas$AreaCode)
  if(dim(puer)[1]==0) { tmp <- c(tmp,"\n# Number of data points (puerulus samples)\n",'0\n# Area\tYear\tIndex\tSD\n') } else {
    tmp <- c(tmp,"\n# Number of data points (puerulus samples)\n",nrow(puer),'\n# Area\tYear\tIndex\tSD\n')
    for(i in 1:nrow(puer)){ tmp <- c(tmp, paste(puer[i,], collapse = "\t"),"\n")} }

  tmp <- c(tmp,"\n# Environmental Data - ln(Efficiency creep)\n# Number of environmental series (commercial efficiency creep for each area)\n",0)

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/DATA.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

  #### Control File ####
  print("Building Control File")

  tmp <- list()
  tmp <- c(tmp, "# Fishery\n\n# weight-at-length (W=aL^b) (kg) \n")
  biol <- readWorkbook(wb,sheet='Biologicals', startRow = 2) %>% filter(Type==4)
  female.wat <- round(biol$par_a[biol$Sex=='F']*lens^biol$par_b[biol$Sex=='F'] ,3)
  male.wat <- round(biol$par_a[biol$Sex=='M']*lens^biol$par_b[biol$Sex=='M'] ,3)
  if(length(female.wat)>0 & length(male.wat)>0){ tmp <- c(tmp, paste(male.wat,collapse = "\t"),"\n", paste(female.wat,collapse = "\t"),"\n") }
  if(length(female.wat)==0 & length(male.wat)>0){ tmp <- c(tmp, paste(male.wat,collapse = "\t"),"\n") }
  if(length(female.wat)>0 & length(male.wat)==0){ tmp <- c(tmp, paste(female.wat,collapse = "\t"),"\n") }

  tmp <- c(tmp, "\n# Egg time step (This is when to determine egg production)\n",0,'\n')

  tmp <- c(tmp, "\n# Biomass time step (This is when to determine Biomass)\n",0,'\n')

  tmp <- c(tmp, "\n# Biomass target, threshold, limit (one per area)\n")

  bio <- readWorkbook(wb,sheet='Area', startRow = 2) %>% select(starts_with('biomass'))
  for(i in 1:nrow(bio)){ tmp <- c(tmp, paste(bio[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Fleet specification\n# Fleet Area Name\n\t")
  fleetcode2 <- fleets %>% group_by(fleet) %>% summarise(area=unique(newarea), group=description) %>% mutate(fleet=fleet-min(fleet), area=area-min(area))
  for(i in 1:nrow(fleetcode2)){ tmp <- c(tmp, paste(fleetcode2[i,c('fleet','area','group')], collapse = "\t"),"\n\t")}
  #as.data.frame(fleetcode2)
  tmp <- c(tmp, "\n# Number of Zones\n", zones, "\n# Areas in each Zone\n", paste(as.numeric(sapply(zoneareas, length)), collapse = " "),'\n')
  tmp <- c(tmp, "\n# The Zones\n")
  for(z in 1:length(zoneareas)){
    tmp <- c(tmp, paste(zoneareas[[z]]-1, collapse = " "),'\n')}

  tmp <- c(tmp, "\n# Discard mortality\n#Age\tFleet\ttstep\t",paste(startseason:endseason,collapse = "\t"),'\n')
  dat <- expand.grid(age=(1:ages)-1, fleet=sort(unique(fleets$fleet))-1, step=sort(unique(times$tstep))-1)
  dat %<>% arrange(age,fleet,step)

  # Get discard mortality rate
  gauge <- readWorkbook(wb,sheet='Retention', startRow = 2) %>%  dplyr::select(StartSeason,EndSeason,Fleet,Age,DiscardMortality) %>% mutate(StartSeason=ifelse(StartSeason=='X', startseason, StartSeason),EndSeason=ifelse(EndSeason=='X', endseason, EndSeason), Age=ifelse(Age=='X',paste(ages,collapse=','),Age), primary=1)
  gauge2 <- gauge %>% filter(Fleet=='X')
  if(nrow(gauge2)>0){
    gauge %<>% filter(Fleet!='X')
    gauge2 %<>% mutate(Fleet=ifelse(Fleet=='X', paste(fleets$fleet,collapse=','), Fleet), primary=2)
    gauge <- rbind(gauge, gauge2)
  }
  gauge %<>% tidyr::separate_rows(Fleet, sep = ",", convert = TRUE)
  gauge %<>% group_by(across(all_of(c("StartSeason", "EndSeason", "Fleet", "Age")))) %>%
    filter(n() == 1 | primary == 1) %>% ungroup() %>% dplyr::select(-primary)

  dat2 <- matrix(gauge$DiscardMortality[1], nrow=nrow(dat), ncol=length(startseason:endseason))
  for(r in seq_len(nrow(gauge))){
    dat2[dat$age==(as.numeric(gauge$Age[r])-1) & dat$fleet==(as.numeric(gauge$Fleet[r])-1),  (startseason:endseason)%in%(gauge$StartSeason[r]:gauge$EndSeason[r])] <- gauge$DiscardMortality[r]
  }
  dat <- cbind(dat,dat2)
  for(i in 1:nrow(dat)){ tmp <- c(tmp, paste(dat[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Recruitment_deviations\n",startseason, "\t\t\t# First year with estimated recruitment deviations\n")
  tmp <- c(tmp, "#", endseason+projectseason-1, "\t\t\t# last year with estimated recruitment deviations\n")
  rec <- readWorkbook(wb,sheet='Recruitment', startRow = 2)
  tmp <- c(tmp, rec$Phase[1], "\t\t\t# Phase for recruitment deviations\n")

  tmp <- c(tmp, "\n# Spatial_deviations_in_recruitment\n",startseason, "\t\t\t# First year with estimates spatial recruitment deviations\n")
  tmp <- c(tmp, "#", endseason+projectseason-1, "\t\t\t# Last year with estimates spatial recruitment deviations\n")
  tmp <- c(tmp, rec$Phase[2], "\t\t\t# Phase for spatial recruitment deviations\n\n")

  tmp <- c(tmp, "# Prespecify_rec_devs :  dev # Year\n",1,"\t\t\t# 1 = rec_devs are to be pre-specified\n")
  dat <- data.frame(rec_dev=0, year=(startseason-max(areas$burn_in)):(endseason+projectseason+5))
  for(i in 1:nrow(dat)){ tmp <- c(tmp, dat[i,1],"\t#\t", dat[i,2],"\n")}

  tmp <- c(tmp, "\n# Prespecify_spatial_rec_devs : dev # Year Area\n",1,"\t\t\t# 1 = Spatial_rec_devs are to be pre-specified\n")
  dat <- data.frame(Spatial_rec_dev=0, expand.grid(year=(startseason-burnin):(endseason+projectseason+5),area=(2:length(unique(area$newarea)))))
  for(i in 1:nrow(dat)){ tmp <- c(tmp, dat[i,1],"\t#\t", dat[i,2],"\t", dat[i,3],"\n")}

  wei <- readWorkbook(wb,sheet='Weights', startRow = 2)
  wei %<>% mutate(sex=adjsex(sex, nsex,section='Weights')-1,)
  tmp <- c(tmp, "\n# Weights on the data (simple)\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='cpue'], "\t\t# Weight on CPUE data\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='numbers'], "\t\t# Weight on catch-numbers data\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='length'], "\t\t# Weight on Length-frequency data\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='larvae'], "\t\t# Weight on larval data\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='tag1'], "\t\t# Weight on Tag1 data\n")
  tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='tag2'], "\t\t# Weight on Tag2 data\n")
  # tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='initialN'], "\t\t# Weight on initial numbers\n")
  # tmp <- c(tmp, wei$value[wei$form=='global' & wei$type=='initialPen'], "\t\t# Weight on initial penalty (InitOpt=3 or 5)\n")

  wei %<>% filter(form=='individual')  ## These are pre specified in the doc.  also need to get printout if this exists
  wei %<>% mutate(type2=case_when(type=='cpue'~1,type=='numbers'~2,type=='length'~3,type=='larvae'~4), codeout = '#', id2 = paste(id,type)) %>% select(type2, fleet, tstep, sex, value, codeout, id2) %>% mutate(fleet=fleet-1)
  ## Add a line for every LF by sex and fleet for francis weightings later
  lenw <- len %>% group_by(Fleet,Sex) %>% summarise(num=sum(Samp),.groups = "drop_last") %>% mutate(type2=3, fleet=Fleet-1, tstep=-1, sex=Sex-1, value=1, codeout='#')
  lenw$fleettype <- fleets$group[match(lenw$Fleet, fleets$fleet)]
  lenw %<>% mutate(id2=paste(fleettype , Fleet, 'ifreq')) %>% ungroup() %>% dplyr::select(type2,fleet,tstep,sex,value,codeout,id2)
  ## Look for Length tuning file
  if(max("Length_freq_tunings.csv" %in%list.files(path = paste0(floc,'/Output')))==1){
    lenw <- read.csv(paste0(floc,'/Output/Length_freq_tunings.csv')) %>% mutate(type2=3,tstep= -1,codeout='#',id2='Length-Freq from tuning file',value=scale, fleet=fleet-1,sex=sex-1) %>% dplyr::select(type2,fleet,tstep,sex, value,codeout,id2)
  }
  if(exists('lenw')) {
    wei %<>% filter(type2!=3)
    wei <- rbind(wei,lenw)
  }

  tmp <- c(tmp, "\n# Weights by fleet (Type: 1=cpue,2=numbers,3=length;4=larvae;) - Can tweak by Fleet and index.  See values to -1 for them to encompass all options.  e.g. time-step set to -1 covers all timesteps. our fleets\n# Type\tFleet\tTime-step\tsex\n")
  tmp <- c(tmp, nrow(wei), "\t# set to number of individual weights defined below - 0 would define no individual weights\n")
  if(nrow(wei)>0) for(i in 1:nrow(wei)){ tmp <- c(tmp,paste(paste(wei[i,],collapse = "\t")),"\n")}

  tmp <- c(tmp, "\n# Basic parameters (lower, upper, estimate, phase, link, prior(0=no, 1=normal, 2=gamma, 3=lognormal), prior.mean, prior.sd) - (link will use same par for multiple areas)\n")

  mainpar <- readWorkbook(wb,sheet='MainParameters', startRow = 2)
  NumMPars <- 1+nrow(areas)+ages+3

  n_need <- NumMPars
  vals   <- nrow(mainpar)
  if (vals < n_need){ stop("Not enough parameters have been provided in the MainParameters tab ","(", vals, " supplied, ", n_need, " required).", call. = FALSE)}

  Mpar1 <- mainpar %>% filter(name=='MeanRecruitment') %>% dplyr::select(-name)
  tmp <- c(tmp, paste(Mpar1,collapse = "\t"), "\n")
  Mpar2 <- mainpar[2:(1+length(unique(area$newarea))),] %>% dplyr::select(-name)
  for(i in 1:nrow(Mpar2)) {tmp <- c(tmp, paste(Mpar2[i,],collapse="\t"),"\n")}
  off <- 2+length(unique(area$newarea))
  Mpar3 <- mainpar[off:(off+ages-1),] %>% dplyr::select(-name)
  for(i in 1:nrow(Mpar3)){ tmp <- c(tmp, paste(Mpar3[i,],collapse="\t"),"\n") }
  Mpar4 <- mainpar %>% filter(grepl('Morph1M', name, ignore.case=T)) %>% dplyr::select(-name)
  tmp <- c(tmp, paste(Mpar4,collapse="\t"), "\n")
  Mpar5 <- mainpar %>% filter(grepl('Morph0Q', name, ignore.case=T)) %>% dplyr::select(-name)
  tmp <- c(tmp, paste(Mpar5,collapse="\t"), "\n")
  Mpar6 <- mainpar %>% filter(grepl('Sigma', name, ignore.case=T)) %>% dplyr::select(-name)
  tmp <- c(tmp, paste(Mpar6,collapse="\t"), "\n")

  tmp <- c(tmp, "\n# Q parameters\n")
  Qpar1 <- 1
  nqs <- length(unique(fleets$fleet[fleets$group=='comm']))*length(unique(times$tstep))
  for(i in 1:nqs){ tmp <- c(tmp, paste(c(-100,100,Qpar1,-1),collapse="\t"), "\t\t\t# Multiplier for Environmental index - Keep to 1\n")}

  tmp <- c(tmp, "\n# Efficiency parameters (lower, upper, estimate, phase, link, useprior, prior, priorsd) \n")
  ecpar1 <- effic[,c('lower','upper','est','Phase','Link','useprior','prior','priorsd')] %>%  dplyr::select(lower,upper,est,Phase,Link,useprior,prior,priorsd)
  nECvec <- max(fleets$effic.creep)

  ## Test number of parameters supplied is correct
  n_need <- nECvec
  vals   <- nrow(ecpar1)
  if (vals < n_need){ stop("Not enough parameters have been provided in the Efficiency tab ","(", vals, " supplied, ", n_need, " required).", call. = FALSE)}

  nECpar <- floor((endseason-startseason+1)/effic$temporal.cover)
  tmp <- c(tmp, sum(nECpar),"\t# Number of Efficiency parameters \n")
  for(nv in 1:nECvec){
    for(np in 1:nECpar[nv]){
      tmp <- c(tmp, paste(ecpar1[nv,],collapse="\t"), paste("\t\t\t# Efficiency creep par - Pointer", effic$pointer[nv]," One par every",effic$temporal.cover[nv],"years.\n")   )}}

  tmp <- c(tmp, "\n# variance specification parameters (1=Egg Production; 2=Egg Production x area;3=Recruitment x area; 4=Legal Biomass x area;5=Harvest Rate;6=Catch rates;7=Fishing efficiency;8=Unspecified;9=Unspecified;10=Unspecified)\n")
  tmp <- c(tmp, "# Number of variance specifications\n",10)

  Vars <- rep(0, 10)
  if(exists('Varspos')){
    Vars[Varspos] <- 1     }
  tmp <- c(tmp, "\n# variance components\n",paste(Vars,collapse = '\t'),"\n")

  tmp <- c(tmp, "\n1 # use the pin file for specifying parameters (ADMB)")
  tmp <- c(tmp, "\n1 # last function call (ADMB)\n")

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/CONTROL.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

  #### Growth file ####
  print("Building Growth File")
  growth <- readWorkbook(wb,sheet='Growth', startRow = 2)
  ## Check if pars provided
  IsGPars <- ifelse(length(growth$est[!is.na(growth$est)])>0, 1, 0)
  if(IsGPars) {
    gpars <- growth %>% dplyr::select(lower,upper,est,Phase,Link,UsePrior,Prior,Priorsd,description)
  }
  growth <- growth %>% dplyr::select(startseason,endseason,sex,area,tstep,grow,matrix) %>% filter(!is.na(startseason))
  growth %<>% mutate(sex=adjsex(sex,nsex,section='Growth')) %>% rowwise() %>% mutate(Years=paste(startseason, endseason, sep='-')) %>% tidyr::separate_rows(area, sep = ",", convert = TRUE) %>% arrange(startseason,  sex, area, tstep) %>% as.data.frame()
  umat <- unique(growth$matrix)
  nstm <- length(umat)
  growth2 <- growth %>% filter(!is.na(startseason)) %>% mutate(area = as.numeric(area) - min(as.numeric(area))) %>% group_by(sex, tstep, grow, matrix) %>% summarise(areas = paste(area,collapse = ","),Years = paste(Years,collapse = ","), .groups = "drop")
  Sex <- growth2$sex-1
  Tstep <- as.numeric(growth2$tstep)
  Pointer <- 0:(nrow(growth2)-1)
  compound <- as.numeric(growth2$grow)
  Area <- growth2$areas
  Years <- growth2$Years

  gspec <- data.frame(Pattern=0:(nstm-1), Type=1, Sex=Sex, Extra=0,Pointer=Pointer,Mpower=1,hash='#',tsteps=Tstep, growthareas=Area, Years=Years, Compound=compound)
  if(IsGPars) { gspec$Type=2 } ## Change to estimatable

  dat <- expand.grid(sex=sexs, age=(1:ages)-1, area=sort(unique(areas$AreaCode))-1, step=sort(unique(times$tstep))-1)
  dat2 <- matrix(-1, nrow=nrow(dat), ncol=length(startseason:endseason))
  seasons <- startseason:endseason
  for(r in 1:nrow(growth)){
    dat2[dat$sex==growth$sex[r]-1  & dat$area==as.numeric(growth$area[r])-1 & dat$step==as.numeric(growth$tstep[r])-1, seasons%in%growth$startseason[r]:growth$endseason[r]] <-  which(growth$matrix[r]==umat)-1     }

  dat <- cbind(dat,dat2)
  dat %<>% arrange(sex,age,area,step)
  ## Check for complete data - does every sex age area moult at least once every year?
  tdat <- dat %>% tidyr::pivot_longer(cols = -c(sex, age, area, step), names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year)) %>% group_by(sex, age, area, year) %>%
    summarise(has_valid = any(value != -1), .groups = "drop") %>% filter(!has_valid) %>% as.data.frame()
  if (nrow(tdat) > 0) {
    missing_str <- paste(
      apply(tdat[, c("sex", "age", "area", "year")], 1, function(x)
        paste0("sex=", x["sex"], " age=", x["age"], " area=", x["area"], " year=", x["year"])
      ),
      collapse = "\n"
    )
    warning("The following sex/age/area/year combinations have no growth assigned:\n", missing_str)
  }

  tmp <- list()
  tmp <- c(tmp, "# Growth specification\n\n# Number of growth Patterns (Mpower is legacy and needed for power function on growth)\n",nrow(gspec),"\n")
  tmp <- c(tmp, "# ", paste(colnames(gspec),collapse='\t'),"\n")
  for(i in 1:nrow(gspec)){ tmp <- c(tmp,'\t',paste(gspec[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for growth\t\n# Sex\tAge\tArea\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")
  for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Number prespecified size-transition\n",length(unique(gspec$Pointer)),"\n")
  tmp <- c(tmp, "# Sex of matrices\n",paste(gspec$Sex,collapse = "\t"),"\n")

  tmp <- c(tmp, "# Prespecified size-transition\n")

  ## If IsFPars then make STMs otherwise load them from the file pre-specified
  if(IsGPars) {
    lbinL <- lens
    lbinM <- lens+(dynamics$value[dynamics$object=='lbgap']/2)
    lbinU <- lens+(dynamics$value[dynamics$object=='lbgap'])

    ## Make STM and then load it to the file
    for(st in 1:nstm) {
      srt <- (st-1)*8+1
      Pins <- gpars[srt:(srt+7),]

      ## Make Growth Vector
      Amax <- exp(Pins$est[1])
      P2   <- Pins$est[2]
      P1   <- exp(Pins$est[3])
      P3   <- exp(Pins$est[4])
      P5   <- exp(Pins$est[5])
      scale =  -1 / Pins$est[8]
      loc   = -Pins$est[7] / Pins$est[8]

      xdev  <- lbinM - P2
      grow1 <- 1 / (1 + exp(xdev / P1))
      grow2 <- 1 / (1 + exp(xdev / P3))
      swap1 <- 1 / (1 + exp(xdev / P5))
      swap2 <- 1 - swap1
      growthvec <- Amax * (grow1 * swap1 + grow2 * swap2)

      ## Make STM
      nlbin <- length(lens)
      STM <- matrix(0, ncol=nlbin, nrow=nlbin)
      for (fm in 1:nlbin) {
        mn_growth <- growthvec[fm]
        sd_growth <- exp(Pins$est[6]) * mn_growth
        Pmoult  <- 1/(1+exp((lbinM[fm]-loc)/scale))

        probs <- rep(0, nlbin)
        for (k in fm:(nlbin - 1)) {
          probs[k] <- pnorm(lbinU[k], lbinM[fm] + mn_growth, sd_growth) - pnorm(lbinL[k], lbinM[fm] + mn_growth, sd_growth)
        }
        probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbinM[fm] + mn_growth, sd_growth)
        probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])
        STM[fm:nlbin, fm] <- Pmoult * probs_norm
        STM[fm, fm] <- STM[fm, fm] + (1 - Pmoult)
      }

      ## Load each STM to the file
      tmp <- c(tmp,paste("\n# Matrix #",umat[st],"\n"))
      for(rr in 1:nlbin){  tmp <- c(tmp,paste(round(as.numeric(STM[rr,]),10),collapse = " "),"\n")      }
    } }


  if(!IsGPars){
    STM <- growth <- readWorkbook(wb,sheet='SizeTransMatricesNew', startRow = 2, colNames = F)
    STM[is.na(STM)] <- ''

    lines <- unlist(STM$X1)  # labels are in X1
    label_idx <- which(grepl("^# ", lines))
    stm_list <- vector("list", length(umat))
    names(stm_list) <- umat
    nlbin <- length(lens)
    for (r in 1:nstm) {
      label_idx <- which(grepl(paste0("^# *", umat[r]), lines))
      if (length(label_idx) == 0)
      {
        stop(paste("Mismatch between STM labels in Growth and SizeTransMatricesNew sheets. Label: ",umat[r]))
      }

      start <- label_idx + 1
      mat_rows <- STM[start:(start + nlbin - 1), 1:nlbin]
      tmp <- c(tmp,paste("\n# Matrix #",umat[r],"\n"))
      for(rr in 1:nrow(mat_rows)){  tmp <- c(tmp,paste(round(as.numeric(mat_rows[rr,]),10),collapse = " "),"\n")
      }}

  }
  tmp <- c(tmp, "\n# Growth parameters\n","# Lower, Upper, Estimate, Phase, Link, Prior(0=no, 1=normal, 2=gamma, 3=lognormal), prior.mean, prior.sd, ID\n")
  if(IsGPars) {
    for(rr in 1:nrow(gpars)){  tmp <- c(tmp,paste(gpars[rr,],collapse = " "),"\n") }
  }

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/GROWTHSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

  #### Reproduction file ####
  print("Building Reproduction File")
  repo <- readWorkbook(wb,sheet='Biologicals', startRow = 2) %>% mutate(Sex=adjsex(Sex,nsex,section='Reproduction')) %>% rowwise() %>% mutate(Years=paste(Startseason, Endseason, sep='-'))
  nrepo <- repo %>% group_by(Type) %>% summarise(num=length(Sex))
  nstm <- nrow(repo)

  ## Maturity
  mature <- repo %>% filter(Type==1) %>% ungroup() %>% mutate(pos = row_number()) %>% tidyr::separate_rows(area, sep = ",", convert = TRUE) %>% arrange(Startseason, Sex, area) %>% as.data.frame()
  matspec <- data.frame(Pointer=0:(nrepo$num[nrepo$Type==1]-1),Par_a=repo$par_a[repo$Type==1],Par_b=repo$par_b[repo$Type==1],hash='#',Area=repo$area[repo$Type==1], Years=repo$Years[repo$Type==1])
  dat <- expand.grid(age=(1:ages)-1, area=sort(unique(areas$AreaCode))-1)
  dat2 <- matrix(-1, nrow=nrow(dat), ncol=length(startseason:endseason))
  seasons <- startseason:endseason
  for(r in 1:nrow(mature)){
    dat2[dat$age>=mature$Minage[r]-1 & dat$area==as.numeric(mature$area[r])-1, seasons%in%mature$Startseason[r]:mature$Endseason[r]] <- mature$pos[r]-1     }
  dat <- cbind(dat,dat2)
  dat %<>% arrange(age,area)
  ## Check for complete data - does every sex age area moult at least once every year?
  tdat <- dat %>%
    tidyr::pivot_longer(cols = -c(age, area), names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year)) %>%
    group_by(area, year) %>%
    summarise(has_valid = any(value != -1), .groups = "drop") %>%
    filter(!has_valid) %>%
    as.data.frame()

  if (nrow(tdat) > 0) {
    missing_str <- paste(
      apply(tdat[, c("area", "year")], 1, function(x)
        paste0("area=", x["area"], " year=", x["year"])
      ),
      collapse = "\n"
    )
    warning("The following area/year combinations have no maturity assigned:\n", missing_str)
  }
  matdat <- dat
  minmaturity <- mature %>% group_by(area) %>% summarise(minage=min(Minage))
  if(min(minmaturity$minage)>(ages-1)) minmaturity$minage <- ages-1

  ## Multiple Spawn
  spawn <- repo %>% filter(Type==3) %>% ungroup() %>% mutate(pos = row_number()) %>% tidyr::separate_rows(area, sep = ",", convert = TRUE) %>% arrange(Startseason, Sex, area) %>% as.data.frame()
  if(nrow(spawn)>0) {
    spawnspec <- data.frame(Pointer=0:(nrepo$num[nrepo$Type==3]-1),Par_a=repo$par_a[repo$Type==3],Par_b=repo$par_b[repo$Type==3],Par_c=repo$par_c[repo$Type==3],hash='#',Area=repo$area[repo$Type==3], Years=repo$Years[repo$Type==3])
  }
  if(nrow(spawn)==0) {
    spawn <- data.frame(Type=3, Sex=0, Minage=0, Startseason=startseason, Endseason=endseason, area=sort(unique(areas$AreaCode)), par_a=1, par_b= -0.1, par_c=1, Years=paste0(startseason,'-',endseason), pos=1)
    spawnspec <- data.frame(Pointer=0:(nrepo$num[nrepo$Type==1]-1),Par_a=1,Par_b=-1,Par_c=1,hash='#',Area=repo$area[repo$Type==1], Years=repo$Years[repo$Type==1])
  }
  dat <- expand.grid(age=(1:ages)-1, area=sort(unique(areas$AreaCode))-1)
  dat2 <- matrix(-1, nrow=nrow(dat), ncol=length(startseason:endseason))
  for(r in 1:nrow(spawn)){
    dat2[dat$age>=spawn$Minage[r]-1 & dat$area==as.numeric(spawn$area[r])-1, seasons%in%spawn$Startseason[r]:spawn$Endseason[r]] <- spawn$pos[r]-1  }
  dat <- cbind(dat,dat2)
  dat %<>% arrange(age,area)
  ## Check for complete data - does every age area double spawn at least once every year?
  tdat <- dat %>%
    tidyr::pivot_longer(cols = -c(age, area), names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year)) %>%
    group_by(area, year) %>%
    summarise(has_valid = any(value != -1), .groups = "drop") %>%
    filter(!has_valid) %>%
    as.data.frame()

  if (nrow(tdat) > 0) {
    missing_str <- paste(
      apply(tdat[, c("area", "year")], 1, function(x)
        paste0("area=", x["area"], " year=", x["year"])
      ),
      collapse = "\n"
    )
    warning("The following area/year combinations have no multiple spawning assigned:\n", missing_str)
  }
  spawndat <- dat

  ## Fecundity
  fec <- repo %>% filter(Type==2) %>% ungroup() %>% mutate(pos = row_number()) %>% tidyr::separate_rows(area, sep = ",", convert = TRUE) %>% arrange(Startseason, Sex, area) %>% as.data.frame()
  fecspec <- data.frame(Pointer=0:(nrepo$num[nrepo$Type==2]-1),Par_a=repo$par_a[repo$Type==2],Par_b=repo$par_b[repo$Type==2],hash='#',Area=repo$area[repo$Type==2], Years=repo$Years[repo$Type==2])
  dat <- expand.grid(age=(1:ages)-1, area=sort(unique(areas$AreaCode))-1)
  dat2 <- matrix(-1, nrow=nrow(dat), ncol=length(startseason:endseason))
  for(r in 1:nrow(fec)){
    dat2[dat$age>=fec$Minage[r]-1 & dat$area==as.numeric(fec$area[r])-1, seasons%in%fec$Startseason[r]:fec$Endseason[r]] <- fec$pos[r]-1  }
  dat <- cbind(dat,dat2)
  dat %<>% arrange(age,area)
  ## Check for complete data - does every age area double spawn at least once every year?
  tdat <- dat %>%
    tidyr::pivot_longer(cols = -c(age, area), names_to = "year", values_to = "value") %>%
    mutate(year = as.integer(year)) %>%
    group_by(area, year) %>%
    summarise(has_valid = any(value != -1), .groups = "drop") %>%
    filter(!has_valid) %>%
    as.data.frame()

  if (nrow(tdat) > 0) {
    missing_str <- paste(
      apply(tdat[, c("area", "year")], 1, function(x)
        paste0("area=", x["area"], " year=", x["year"])
      ),
      collapse = "\n"
    )
    warning("The following area/year combinations have no fecundity assigned:\n", missing_str)
  }
  fecdat <- dat

  tmp <- list()
  tmp <- c(tmp, "# Reproductive specification\n# Age at maturity ", 'Area ',paste(minmaturity$area, collapse = ' '),'\n')
  tmp <- c(tmp, paste(minmaturity$minage,collapse = " "),"\n")
  tmp <- c(tmp, "\n# Number of maturity patterns\n",nrow(matspec),"\n")
  tmp <- c(tmp, "# ", paste(colnames(matspec),collapse='\t'),"\n")
  for(i in 1:nrow(matspec)){ tmp <- c(tmp,'\t',paste(matspec[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for maturity\t\n# Age\tArea\t",paste(startseason:endseason,collapse = "\t"),"\n")
  for(i in 1:nrow(matdat)){ tmp <- c(tmp,paste(matdat[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Number of multiple spawning patterns\n",nrow(spawnspec),"\n")
  tmp <- c(tmp, "# ", paste(colnames(spawnspec),collapse='\t'),"\n")
  for(i in 1:nrow(spawnspec)){ tmp <- c(tmp,'\t',paste(spawnspec[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for multiple spawning\t\n# Age\tArea\t",paste(startseason:endseason,collapse = "\t"),"\n")
  for(i in 1:nrow(spawndat)){ tmp <- c(tmp,paste(spawndat[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Number of fecundity patterns\n",nrow(fecspec),"\n")
  tmp <- c(tmp, "# ", paste(colnames(fecspec),collapse='\t'),"\n")
  for(i in 1:nrow(fecspec)){ tmp <- c(tmp,'\t',paste(fecspec[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for fecundity\t\n# Age\tArea\t",paste(startseason:endseason,collapse = "\t"),"\n")
  for(i in 1:nrow(fecdat)){ tmp <- c(tmp,paste(fecdat[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/REPOSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

  #### Movement file ####
  print("Building Migration File")

  migrate999 <- migrate %>% filter(Season==999)
  dat <- expand.grid(age=(1:ages)-1, area=sort(unique(areas$AreaCode))-1, step=sort(unique(times$tstep))-1)
  dat2 <- matrix(0, nrow=nrow(dat), ncol=length(startseason:endseason))
  for(i in 1:nrow(migrate999)){
    dat2[dat$age==(migrate999$Age[i]-1) & dat$area==(migrate999$Source[i]-1) & dat$step==(migrate999$tstep[i]-1),] <- i
  }
  migrateNOT999 <- migrate %>% filter(Season<=100)
  if(nrow(migrateNOT999)>0){
    nyears <- ncol(dat2)
    parpoint <- max(dat2)+1
    for(i in 1:nrow(migrateNOT999)){
      for(y in 1:ceiling(nyears/migrateNOT999$Season[i])){
        pos <- ((y-1)*migrateNOT999$Season[i])+(1:migrateNOT999$Season[i])
        pos <- pos[pos<=nyears]
        dat2[dat$age==(migrateNOT999$Age[i]-1) & dat$area==(migrateNOT999$Source[i]-1) & dat$step==(migrateNOT999$tstep[i]-1),pos] <- parpoint
        tmpmigrateNOT999 <- migrateNOT999[i,]
        tmpmigrateNOT999$pointer <- parpoint
        migrate999 <- rbind(migrate999,tmpmigrateNOT999)
        parpoint <- parpoint + 1
      }}
  }
  # ## See if there were any specific years when migration was to be different and implement this.
  # yrs2change <- unique(migrate$Season[migrate$Season!=999 & migrate$Season>1000])
  # if(length(yrs2change)>0){
  #   for(y in 1:length(yrs2change)){
  #     tmpdf <- migrate[migrate$Season==yrs2change[y],]
  #     dat2[dat$age==(tmpdf$Age-1) & dat$area==(tmpdf$Source-1) & dat$step==(tmpdf$tstep-1),(startseason:endseason)==tmpdf$Season] <- tmpdf$pointer
  #   }}

  dat <- cbind(dat,dat2)
  dat %<>% arrange(age, area, step)

  #migratesum <- migrate %>% group_by()

  tmp <- list()
  tmp <- c(tmp, "# Movement section\n# Number of movement patterns\n",nrow(migrate999)+1,"\n")
  tmp <- c(tmp, "# Pattern\tType\tDest\tExtra\t(Type 0: none; 1 constant [prespecified or estimated; 1 parameter]; 2 knife-eded-specific [pre-specified or estimated; 2 parameters])\n")
  tmp <- c(tmp,"\t",paste(c(0,0,0,0),collapse = "\t"),"\n")
  Nmigratepars <- nrow(migrate999)
  if(Nmigratepars>0) {
    for(i in 1:Nmigratepars) {
      tmp <- c(tmp,"\t",paste(c(i,1,migrate999$Dest[i]-1,0),collapse = "\t"),"\n")}}

  tmp <- c(tmp, "# Movement specifications\n")
  tmp <- c(tmp, "#Age\tArea\tTStep\t", paste(startseason:endseason,collapse = "\t"),"\n")
  for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

  tmp <- c(tmp, "# Movement parameters\n# Lower\tUpper\tEstimate\tPhase\tLink\tUsePrior\tPrior_mn\tPrior_sd\tSource to Dest & Age\n")
  if(Nmigratepars>0) {
    Mpar <- migrate999 %>% dplyr::select(lower, upper, est, Phase, Link, UsePrior, Prior, Priorsd)
    for(i in 1:nrow(Mpar)){
      tmp <- c(tmp,"\t",paste(Mpar[i,],collapse = "\t"),paste("\t\t #", migrate999$Source[i],"to", migrate999$Dest[i],"&",migrate999$Age[i],"\n"))}
  }

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/MOVESPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)


  #### TagRecapture file ####
  tag <- readWorkbook(wb,sheet='TagRecapture', startRow = 2,colNames = FALSE)
  LoadTdata <- tag[tag[,1]=='Load tag data',2]

  if(LoadTdata==0){
    tmp <- list()
    tmp <- c(tmp, "# IsTagData\n",0,"\n")
    write.table(tmp, paste(floc,'/TAGSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
    tmp <- c(tmp,0,"\t# Number proportion observations\n", "# Year Tstep area type1 type2\n")
    write.table(tmp, paste(floc,'/TAGPROP.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
  }

  if(LoadTdata==1){
    print("Building Tag Recapture File")

    tmp <- list()
    tmp <- c(tmp, "# IsTagData\n",tag[tag[,1]=='Load tag data',2],"\n")
    tmp <- c(tmp, "# Initial tagloss\n",tag[tag[,1]=='Initial tag loss',2],"\n")
    tmp <- c(tmp, "# Longterm tagloss\n",tag[tag[,1]=='Long-term tag loss',2],"\n")
    Nreptype <- as.numeric(tag[tag[,1]=='Num of reporting types',2])
    tmp <- c(tmp, "# Number of reporting types\n",tag[tag[,1]=='Num of reporting types',2],"\n")
    tmp <- c(tmp, "# Reporting rates of recapture types\n",paste(tag[tag[,1]=='Reporting rates (x type)',2:(Nreptype+1)], collapse=" "),"\n")
    tmp <- c(tmp, "# Use size (0=N, 1=Y)\n",paste(tag[tag[,1]=='Use recapture size (x type)',2:(Nreptype+1)],collapse=" "),"\n")
    tmp <- c(tmp, "# Number of tsteps to ignore\n",tag[tag[,1]=='Timesteps to ignore',2],"\n")
    tmp <- c(tmp, "# Release areas\n",tag[tag[,1]=='Num release areas',2],"\n")

    pos1 <- which(tag[,1]=='Release by lbin')+1; pos2 <- which(tag[,1]=='Recaptures by lbin')
    release <- data.frame(tag[(pos1+1):(pos2-1),])
    colnames(release) <- tag[pos1,]
    release <- release[,!is.na(release[1,])] %>% filter(!is.na(Total))
    release %<>% mutate(Sex=adjsex(Sex, nsex,section='Tag_release'))
    relyr <- release %>% group_by(RelArea) %>% summarise(Min=min(Year))
    tmp <- c(tmp, "# First release year\n",paste(relyr$Min, collapse=' '),"\n")

    tmp <- c(tmp, "# Release by lbin\n",nrow(release),"\t# number release observations\n")
    tmp <- c(tmp, "# Sex	Group	Area	Year	Tstep	Total ", paste0("lbin",lens, collapse=' '),"\n")
    for(i in 1:nrow(release)){ tmp <- c(tmp,paste(release[i,],collapse = " "),"\n")}

    pos1 <- which(tag[,1]=="Recaptures by lbin")+1; pos2 <- which(tag[,1]=="Proportions of effort between types")
    recap <- data.frame(tag[(pos1+1):(pos2-1),])
    #tail(recap)
    colnames(recap) <- tag[pos1,]
    recap <- recap[,!is.na(recap[1,])]
    recap %<>% mutate(Sex=adjsex(Sex, nsex,section='Tag_recapture'))
    tmp <- c(tmp, "# Recaptures by lbin\n",nrow(recap),"\t# Number recapture observations\n")
    tmp <- c(tmp, "# Sex RelArea RecArea Type Year Tstep Total ", paste0("lbin",lens, collapse=' '),"\n")
    for(i in 1:nrow(recap)){ tmp <- c(tmp,paste(recap[i,],collapse = " "),"\n")}

    recap2 <- recap %>% dplyr::select(Sex, RecArea, RelArea, RecType, Year, Tstep, Total) %>% mutate(Total=as.numeric(Total), Tstep=as.numeric(as.character(Tstep)), Tstep=paste0('ts', Tstep) ) %>% tidyr::pivot_wider(names_from = Tstep, values_from = Total,names_sort = TRUE)
    recap2[is.na(recap2)] <- 0

    tmp <- c(tmp, "# Recaptures by timestep\n",nrow(recap2),"\t# Number recapture observations\n")
    tmp <- c(tmp, "#", paste(colnames(recap2), collapse=" "),"\n")
    for(i in 1:nrow(recap2)){ tmp <- c(tmp,paste(recap2[i,],collapse = " "),"\n")}

    tmp <- c(tmp, "\n# Final check\n123456")
    write.table(tmp, paste(floc,'/TAGSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

    tmp <- list()
    pos1 <- which(tag[,1]=="Proportions of effort between types")+1
    prop <- data.frame(tag[(pos1+1):nrow(tag),])
    colnames(prop) <- tag[pos1,]
    prop <- prop[,!is.na(prop[1,])]
    tmp <- c(tmp,nrow(prop),"\t# Number proportion observations\n", "# Year Tstep area type1 type2\n")
    for(i in 1:nrow(prop)){ tmp <- c(tmp,paste(prop[i,],collapse = " "),"\n")}

    tmp <- c(tmp, "\n# Final check\n123456")
    write.table(tmp, paste(floc,'/TAGPROP.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
  }
  #### Recruitment file ####
  print("Building Recruitment File")

  tmp <- list()
  tmp <- c(tmp, "# Recruitment specifications\n\n# Number of sex_area_allocation options\n",1,"\n\n")
  tmp <- c(tmp, "# Allocate_yearxarea (form 0 is estimate sex split and then area split; form 1 is estimate sex*area split)\n")
  tmp <- c(tmp, "# Type\tForm\tAlloc_sex_area\n",paste(c(0,0,0),collapse = '\t'),'\n\n')

  dum <- data.frame(type=(1:length(unique(areas$recruitarea)))-1, sex1=0, sex2=0, Extra1=-1, Extra2=-1)
  dum$sex1 <- -seq(1,(2*nrow(dum))-1,2)
  dum$sex2 <- dum$sex1-1

  tmp <- c(tmp, "# Number of length_allocation options (size compoisition of recruits by area Not sex)\n",length(unique(areas$recruitarea)),"\n\n")
  tmp <- c(tmp, "# Allocate_length\n# Type\tAllocate_sex1\tAllocate_sex2\tExtra1\tExtra2\t: 2nd column sex 1 and 3rd col for sex 2 (Offset for PreSpecRecFrac = -1*col.value-1)\n")
  for(r in 1:nrow(dum)) {
    tmp <- c(tmp, "\t", paste(dum[r,],collapse = "\t"),paste("\t# Areas ", paste(unique(areas$aname[areas$recruitarea==r]), collapse=','), "\n"))
  }
  tmp <- c(tmp, "\n# Recruit by year and time-step\n# TimeStep\t", paste(startseason:(endseason+projectseason+10),collapse = "\t"),"\n")
  dat <- matrix(-1, nrow=length(unique(times$tstep)), ncol=length(startseason:(endseason+projectseason+10)))
  dat[unique(times$tstep[times$recruit==1]),] <- 0
  dat <- cbind(data.frame(ts=sort(unique(times$tstep))-1),dat)
  for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = "\t"),"\n")}

  rec <- areas %>% group_by(AreaCode) %>% summarise(area=mean(recruitarea)-1)
  nsizecomp <- length(unique(rec$area))*length(sexs)
  tmp <- c(tmp, "\n# Recruit by area\n",paste(rec$area,collapse = "\t"),"\n")

  rec <- readWorkbook(wb,sheet='Recruitment', startRow = 2)
  usepar <- tail(rec$Use.Parameters,1)

  tmp <- c(tmp, "\n# Recruitment fractions use prespecified (0) or calculate from parameters (1). If pre-specified need to code out parameters #\n")
  tmp <- c(tmp, paste(usepar, "\t# Use parameters 0 = No, 1 = Yes.\n"))

  tmp <- c(tmp, "\n# Number of pre-specified recruitment functions\n",max(areas$recruitarea)*length(sexs),"\n")
  tmp <- c(tmp, "# Prespecified recruitment fractions (based on mean CL + SD of lobster at the start of the year age 3)\n")

  recdist <- rec %>% dplyr::select(starts_with('Prespecified'))
  recdist <- recdist[!is.na(recdist[,1]),]
  ## Trim in case more info has been added but we want only have minimal size dist for recruitment, i.e. same size comp for all areas
  if(nrow(recdist)!=nsizecomp) warning("Number of predetermined size at recruitment does not match recruitment areas defined in Area tab. \nThey have been truncated.")
  recdist <- recdist[1:nsizecomp,]
  for(i in 1:nrow(recdist)){ tmp <- c(tmp,paste(recdist[i,],collapse = "\t"),"\n")  }

  tmp <- c(tmp, "\n#  Recuitment1 parameters\n# lower, upper, estimate, phase, link, prior(0=no, 1=normal, 2=gamma, 3=lognormal), prior.mean, prior.sd::\t Number of recruitment size parameter pairs (mean+sd) must match Number of pre-specified recruitment functions above.\n")

  rec %<>% dplyr::select(Use.Parameters,lower,upper,est,Phase,Link,useprior,mnprior,sdprior,description) %>% mutate(Use.Parameters=ifelse(Use.Parameters==1,'','#'), description =paste('#', description ))
  ## Remove the Rec deviations
  rec <- rec[3:nrow(rec),]
  npars <- (length(unique(areas$AreaCode)))+nsizecomp*2
  if(nrow(rec)!=npars) warning("Number of recruitment pars for size at recruitment does not match recruitment areas defined in Area tab. \nThey have been truncated.")
  rec <- rec[1:npars,]
  for(a in 1:nrow(rec)){ tmp <- c(tmp, paste(rec[a,],collapse='\t'), "\n")  }

  # Bias ramp
  tmp <- c(tmp, "\n# Bias ramp - insert description here\n",paste(c(startseason, startseason, endseason, endseason),collapse = "\t"),"\t#description\tdescription\tdescription\tdescription\n")

  if(suppressWarnings(!is.null(readWorkbook(wb,sheet='PuerulusPar', startRow = 3)))){
    puerpar <- readWorkbook(wb,sheet='PuerulusPar', startRow = 2) %>% mutate(description=paste('#', description))
    tmp <- c(tmp, "\n# Puerulus Power for puerulus to recruit relationship\n",nrow(puerpar),"\n#LB\tUP\tEstimate\tPhase\n")
    for(a in 1:nrow(puerpar)){ tmp <- c(tmp, paste(puerpar[a,],collapse = "\t"), "\n") }} else {
      tmp <- c(tmp, "\n# Puerulus Power for puerulus to recruit relationship\n",0,"\n#LB\tUP\tEstimate\tPhase\n")
    }

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/RECRUITSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)


  #### Retainment file ####
  print("Building Retain/Discard File")

  hgrad <- readWorkbook(wb,sheet='HighGrading', startRow = 2) %>% group_by(season, area, tstep) %>% summarise(prop=mean(prop), .groups = 'drop') %>% mutate(propfl=1-trunc(prop/0.01)*0.01)
  hglist <- sort(unique(c(hgrad$propfl,1)))
  hgrad99 <- hgrad %>% filter(area==99)

  dat <- expand.grid(sex=sexs, age=(1:ages)-1, fleet=sort(unique(fleets$fleet))-1, step=sort(unique(times$tstep))-1)
  dat2 <- matrix(length(hglist)-1, nrow=nrow(dat), ncol=length(startseason:endseason))
  for(r in 1:nrow(hgrad)){
    if(hgrad$area[r]==99) fl <- fleets$fleet[fleets$group=='comm']-1
    if(hgrad$area[r]!=99) fl <- fleets$fleet[fleets$group=='comm' & fleets$newarea==hgrad$area[r]]-1
    dat2[dat$fleet%in%fl & dat$step==(hgrad$tstep[r]-1),(startseason:endseason)==hgrad$season[r]]  <- which(hglist==hgrad$propfl[r])-1}

  dat <- cbind(dat,dat2)
  dat %<>% arrange(sex, age, fleet, step)
  tmp <- list()
  tmp <- c(tmp, "# Retain specification (This represents the proportion of LEGAL animals retained - (1-high-graded due to low value))\n\n# Number Retain Patterns\n",length(hglist),"\n")
  tmp <- c(tmp, "# Pattern\tType\tSex\tExtra\tPointer","\n")
  for(i in 1:length(hglist)){ tmp <- c(tmp," ",paste(c((i-1),1,-1,0,i-1),collapse = " "),"\n")}

  tmp <- c(tmp, "\n# Specifications for retention. Fleet = ", paste(fleets$group,collapse=' '),"\n# Sex Age Fleet Step ",paste(startseason:endseason,collapse = " "),"\n")

  for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = " "),"\n")}

  tmp <- c(tmp, "\n# retention - Proportion of LEGAL animals retained - NOT high-graded. (values have been rounded to 1%)\n",length(hglist),"\n")
  for(i in 1:length(hglist)){ tmp <- c(tmp,paste(rep(hglist[i],length(lens)),collapse = " "),"\n")}

  tmp <- c(tmp, "\n# Retention parameters\n# Lower\tUpper\tEstimate\tPhase\n")

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/RETAINSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)


  #### Selection file ####
  print("Building Gear selectivity File")
  ## Retention
  egap <- readWorkbook(wb,sheet='Selectivity', startRow = 2)
  egappar <- egap %>% filter(!is.na(yearlink)) %>% mutate(uniq=paste(sex,yearlink)) %>% dplyr::select(!starts_with('fleet'))%>% mutate(Sex=ifelse(sex=='F',0,1), Sex=Sex-min(Sex))

  ## Look at number of pars
  fleetyr <- egap %>% dplyr::select(starts_with('fleet'))
  if(ncol(fleetyr)==1) {tfleetyr <- fleetyr;
  colnames(tfleetyr) <- 'fleet999';
  fleetyr <- cbind(fleetyr, tfleetyr)  }
  fleetyr <- fleetyr[egap$season%in%startseason:endseason,]
  pars <- sort(unique(as.vector(as.matrix(fleetyr))))

  egappar_sum <- egappar %>% mutate(Sex=ifelse(sex=='F',0,1), Sex=Sex-min(Sex)) %>% group_by(yearlink,Sex,form,uniq) %>% summarise(num=length(Sex), .groups = 'drop') %>% ungroup() %>% mutate(type=case_when(form=='logistic'~3, form=='doublelogistic'~9, form=='knife'~4)) %>% as.data.frame() %>% arrange(Sex) %>% ungroup() %>% mutate(pattern=as.numeric(rownames(.))-1, Pointer=pattern) %>% dplyr::select(pattern, type, Sex, num, Pointer, uniq, yearlink)

  egappar_sumog <- egappar_sum
  for(p in pars){
    if(!p%in%egappar_sum$yearlink){
      reps <- which(egappar_sumog$yearlink==floor(p))
      tmpe <- egappar_sum[reps,]
      tmpe$uniq <- paste(substr(tmpe$uniq,1,1), p)
      tmpe$pattern <- max(egappar_sum$pattern) + (1:nrow(tmpe))
      tmpe$yearlink <- p
      egappar_sum <- rbind(egappar_sum, tmpe)
    }
  }
  egappar_sum %<>% mutate(Pointer=pattern)
  nes <- nrow(egappar_sum)
  negappar <- nrow(egappar_sum)


  Selid <- expand.grid(Sex=sort(unique(sexs)), Age=sort(unique(1:ages))-1, Fleet=sort(unique(fleets$fleet))-1, Step=sort(unique(times$tstep))-1) %>% arrange(Sex, Age, Fleet)
  Semat <- matrix(0, nrow=nrow(Selid), ncol=length(startseason:endseason))
  SFleets <- as.numeric(gsub('fleet','',colnames(fleetyr)))
  Sseasons <- egap$season[egap$season%in%startseason:endseason]

  if(!length(unique(Selid$Fleet))==length(SFleets)) warning(paste0("Number of fleets in Fleet tab (",length(unique(Selid$Fleet)),") do not match the columns of fleets in the Retention tab (",length(SFleets),")"))

  for(i in 1:nrow(Selid)){
    Yrlinks <- fleetyr[,(SFleets-1)==Selid$Fleet[i]]
    Sexegappar_sum <- egappar_sum %>% filter(Sex==Selid$Sex[i])
    Pointers <- Sexegappar_sum$Pointer[match(Yrlinks,Sexegappar_sum$yearlink)]
    if(length(Pointers)==length(startseason:endseason)) {  Semat[i,] <- Pointers } else { warning('Seasons needed for fleet assignment are not all present in the year link assignment on the selectivity tab. Ensure the entire timescale is represented'); break  }
  }

  iswhite <- readWorkbook(wb,sheet='IsMorph', startRow = 2)

  tmp <- list()
  tmp <- c(tmp, "# Selex specification\n# Number Selex Patterns\n",negappar)
  tmp <- c(tmp, "\n# Pattern Type Sex Npars Pointer # Type PRESPECIFIED 1, COEFFICIENTS 2, LOGISTIC 3, KNIFE 4, DOUBLELOG 9\n")
  for(i in 1:nrow(egappar_sum)){ tmp <- c(tmp,paste(egappar_sum[i,1:5],collapse = " "),"\n")}

  ids <- paste(0:7, unique(egappar$comment), sep='=', collapse = ", ")

  tmp <- c(tmp, paste0("# Specifications for selectivity (for example escape gaps),", ids,". \n"), "# Sex Age Fleet Step: ",paste(startseason:endseason,collapse = " "),"\n")

  code <- cbind(Selid,Semat)
  for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = " "),"\n")}

  tmp <- c(tmp, "# Selectivity Parameters\n", "# Lower, Upper, Estimate, Phase, Link, Prior(0=no, 1=normal, 2=gamma, 3=lognormal), prior.mean, prior.sd, ID\n")

  tegappar <- ExpandSelectPars(wb, startseason, endseason)

  for(i in 1:nrow(tegappar)){
    tmp <- c(tmp, paste(tegappar[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "# selectivity\n", negappar,"\n")
  for(i in 1:negappar){
    tmpegappar <- tegappar[tegappar$uniq==unique(tegappar$uniq)[i],]
    if(unique(tmpegappar$form)=='logistic')
      qselect <- 1.0/(1.0+exp(-tmpegappar$par[tolower(tmpegappar$id)=='p2']*(lens-tmpegappar$par[tolower(tmpegappar$id)=='p1'])))
    if(unique(tmpegappar$form)=='doublelogistic'){
      qselect <- (1.0/(1.0+exp(-tmpegappar$par[tolower(tmpegappar$id)=='p2']*(lens-tmpegappar$par[tolower(tmpegappar$id)=='p1']))))*(1.0/(1.0+exp(-tmpegappar$par[tolower(tmpegappar$id)=='p4']*(lens-tmpegappar$par[tolower(tmpegappar$id)=='p3']))))
      qselect <- qselect/max(qselect)}
    qselect <- round(qselect,4)
    tmp <- c(tmp, paste(qselect, collapse = "\t"),"\n")}

  gauge <- readWorkbook(wb,sheet='Retention', startRow = 2) %>% mutate(hash='#', type=1, Extra=0, pointer=pos-1, pos=pointer) %>% dplyr::select(pos, type, Extra, pointer, hash, id)

  tmp <- c(tmp, "\n# Number Legal patterns (What is legal and can be retained [e.g. above Min Legal length, not egg bearing] or in a survey all can be caught)\n", nrow(gauge),"\n# Pattern\tType\tExtra\tPointer\t#  Description\n")

  for(i in 1:nrow(gauge)){ tmp <- c(tmp, paste(gauge[i,], collapse = "\t"), "\n")    }
  gauge <- readWorkbook(wb,sheet='Retention', startRow = 2)
  Pos <- which(colnames(gauge)=='DiscardMortality')+1
  gauge2 <- gauge[,Pos:ncol(gauge)]
  tmp <- c(tmp, "\n# legal patterns\n", nrow(gauge2),"\n")
  for(i in 1:nrow(gauge2)){ tmp <- c(tmp, paste(round(as.numeric(gauge2[i,]),4), collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for Fleet legal assignment\n#Sex\tAge\tFleet\tStep\t",paste(startseason:endseason,collapse = "\t"),"\n")

  gauge3 <- gauge[,1:which(colnames(gauge)=='IsConstantLegal')] %>% mutate(Sex=ifelse(Sex=='F',1, ifelse(Sex=='M',2,Sex)))
  if(length(sexs)==1) gauge3$Sex <- 1 ## Handle if only a one sex model

  ## Expand gauge3 if there are common fleets/ages or timesteps
  gauge3 <- gauge3 %>%  tidyr::separate_rows(Fleet, sep = ",", convert = TRUE) %>%  tidyr::separate_rows(Age, sep = ",", convert = TRUE) %>%  tidyr::separate_rows(TimeStep, sep = ",", convert = TRUE) %>% arrange(pos)

  # Count how many X values and then sort by this.
  gauge3$nX <- apply(as.matrix(gauge3), 1, function(x)  length(x[x=='X']))
  gauge3 %<>% arrange(desc(nX))

  leg <- expand.grid(Sex=sexs, Age=0:(ages-1),  Fleet=sort(unique(fleets$fleet))-1,Step=sort(unique(times$tstep))-1)
  id <- paste(leg[,1],leg[,2],leg[,3],leg[,4], sep="-")
  code <- matrix(0, nrow=nrow(leg), ncol=length(startseason:endseason), dimnames = list(pat=id,year=paste('Y',startseason:endseason,sep='')))
  Yrs <- startseason:endseason
  for(r in 1:nrow(gauge3)){
    tgau <- gauge3[r,]
    if(tgau$StartSeason=='X') {SS <- startseason } else { SS <- tgau$StartSeason}
    if(tgau$EndSeason=='X') {ES <- endseason} else {ES <- tgau$EndSeason}
    if(tgau$Sex=='X') {Sx <- sexs}else{Sx <- as.numeric(tgau$Sex)-1}
    if(tgau$Fleet=='X') {Ft <- fleets$fleet-1}else{Ft <- as.numeric(tgau$Fleet)-1}
    if(tgau$TimeStep=='X') {Ts <- sort(unique(times$tstep))-1} else{Ts <- as.numeric(tgau$TimeStep)-1}
    if(tgau$Age=='X') {Ag <-  sort(unique(1:ages))-1}else{Ag <- as.numeric(tgau$Age)-1}
    code[leg$Sex%in%Sx & leg$Age%in%Ag & leg$Fleet%in%Ft & leg$Step%in%Ts, SS<=Yrs & ES>=Yrs] <- tgau$pos-1
  }

  code <- cbind(leg, code)
  code %<>% arrange(Sex, Age, Fleet, Step)
  for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Specifications for legal biomass\n#Sex Age Area Step ",paste(startseason:endseason,collapse = " "),"\n")
  leg <- expand.grid(Sex=sexs, Age=(1:ages)-1,  Area=sort(unique(areas$AreaCode ))-1,Step=sort(unique(times$tstep))-1)
  id <- paste(leg[,1],leg[,2],leg[,3],leg[,4], sep="-")
  code <- matrix(0, nrow=nrow(leg), ncol=length(startseason:endseason), dimnames = list(pat=id,year=paste('Y',startseason:endseason,sep='')))
  if(length(unique(gauge3$Fleet))==1) {
    if(unique(gauge3$Fleet)=='X') {
      gauge3$UseArea <- 0
      for(f in unique(fleets$newarea)){
        tgauge <- gauge3[1,]
        tgauge$Fleet <- fleets$fleet[fleets$newarea==f][1]
        tgauge$UseArea <- 1
        gauge3 <- rbind(gauge3, tgauge)
      }
    }}
  gauge4 <- gauge3 %>% filter(UseArea==1) %>% mutate(Area=fleets$newarea[match(Fleet,fleets$fleet)])
  for(r in 1:nrow(gauge4)){
    tgau <- gauge4[r,]
    if(tgau$StartSeason=='X') {SS <- startseason}else{SS <- tgau$StartSeason}
    if(tgau$EndSeason=='X') {ES <- endseason} else {ES <- tgau$EndSeason} # ES <-1950
    if(tgau$Sex=='X') {Sx <- sexs}else{Sx <- as.numeric(tgau$Sex)-1}
    if(tgau$Area=='X') {Ar <- sort(unique(fleets$newarea))-1} else {Ar <- as.numeric(tgau$Area)-1}
    if(tgau$TimeStep=='X') {Ts <- sort(unique(times$tstep))-1} else {Ts <- as.numeric(tgau$TimeStep)-1}
    if(tgau$Age=='X') {Ag <-  sort(unique(1:ages))-1} else {Ag <- as.numeric(tgau$Age)-1}
    code[leg$Sex%in%Sx & leg$Age%in%Ag & leg$Area%in%Ar & leg$Step%in%Ts, SS<=Yrs & ES>=Yrs] <- tgau$pos-1
  }

  code <- cbind(leg, code)
  code %<>% arrange(Sex, Age, Area, Step)
  for(i in 1:nrow(code)){ tmp <- c(tmp, paste(code[i,], collapse = " "),"\n")}

  tmp <- c(tmp, "\n# Reference selectivity pattern (This is to set a constant Legal definition)\n")
  ## Set Base LegalBiomass to Legal definition of a male in 1992 which is a min CL of 76 mm
  conslb <- gauge2[gauge3$pos[gauge3$IsConstantLegal==1],][1,]

  if(nsex==2) {
    conslb <- rbind(conslb,conslb)
    Fem <- gauge3$pos[gauge3$IsConstantLegal==1 & gauge3$Sex==1][1]
    Mal <- gauge3$pos[gauge3$IsConstantLegal==1 & gauge3$Sex==2][1]
    if(is.na(Fem)) Fem <- Mal
    if(is.na(Mal)) Mal <- Fem
    if(!is.na(Fem) & !is.na(Mal)) conslb <- gauge2[c(Fem,Mal),]
  }

  for(r in 1:nrow(conslb)){ tmp <- c(tmp, paste(conslb[r,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n\n# IsMorph specifications - assignment of unique life stage quality\n")
  dat <- expand.grid(sex=sexs,age=(1:ages)-1, area=sort(unique(areas$AreaCode ))-1, step=sort(unique(times$tstep))-1, state=1)
  if(nrow(iswhite)>0) {
    for(i in 1:nrow(iswhite)){
      dat$state[dat$age==(iswhite$Age[i]-1) & dat$area==(iswhite$Area[i]-1) & dat$step==(iswhite$Tstep[i]-1)] <- 0
    }}
  dat %<>% arrange(sex, age, area, step)
  tmp <- c(tmp, "#Sex Age Area TStep State","\n")
  for(i in 1:nrow(dat)){ tmp <- c(tmp,paste(dat[i,],collapse = " "),"\n")}

  tmp <- c(tmp, "\n# Final check\n123456")

  write.table(tmp, paste(floc,'/SELEXSPEC.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)

  #### Projection file ####
  print("Building Projection File")

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
  tmp <- c(tmp, "# Specifications for gear selectivity (for example impact of escape gaps)\n")

  tdat  <- read.table(paste(floc,'/SELEXSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  pos1 <- find(c("#",'Specifications','for','selectivity'), tdat, 2)
  pos2 <- find(c("#",'Selectivity'), tdat, -1)
  tdat <- tdat[pos1:pos2,c(1:4,sum(!is.na(tdat[pos1,]) & tdat[pos1,]!=''))  ]
  for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  }
  colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
  tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
  for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n# Retention\n")
  tmp <- c(tmp, "# Specifications for retention. Fleet = comm comm comm comm comm comm comm comm comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor comm_monitor rec rec rec rec ibss ibss ibss ibss ibss iss iss iss iss\n")
  tdat  <- read.table(paste(floc,'/RETAINSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  pos1 <- find(c("#",'Specifications','for','retention'), tdat, 2)
  pos2 <- find(c("#",'retention',"-"), tdat, -1)
  tdat <- tdat[pos1:pos2,c(1:4, sum(!is.na(tdat[pos1,]) & tdat[pos1,]!=''))  ]
  for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  }
  colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
  tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
  for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "#\n# Specifications for Fleet legal assignment\n")
  tdat  <- read.table(paste(floc,'/SELEXSPEC.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  pos1 <- find(c("#",'Specifications','for','Fleet'), tdat, 2)
  pos2 <- find(c("#",'Specifications','for','legal'), tdat, -1)
  tdat <- tdat[pos1:pos2,c(1:4, sum(!is.na(tdat[pos1,]) & tdat[pos1,]!=''))  ]
  for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  }
  colnames(tdat) <- c("Sex","Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
  tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
  for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

  tmp <- c(tmp, "\n#	Discard	mortality\n")
  tdat  <- read.table(paste(floc,'/CONTROL.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  tdat[tdat==''&!is.na(tdat)] <- NA
  pos1 <- find(c("#",'Discard','mortality'), tdat, 2)
  pos2 <- find(c("#",'Recruitment_deviations'), tdat, -1)
  tdat <- tdat[pos1:pos2,c(1:3, sum(!is.na(tdat[pos1,]) & tdat[pos1,]!=''))  ]
  for(proj in 2:projectseason){  tdat <- cbind(tdat, nm=tdat[,ncol(tdat)])  }
  colnames(tdat) <- c("Age","Fleet","Step:",(endseason+1):(endseason+projectseason))
  tmp <- c(tmp, "#",paste(colnames(tdat),collapse = "\t"),"\n")
  for(i in 1:nrow(tdat)){ tmp <- c(tmp, paste(tdat[i,], collapse = "\t"),"\n")}

  tdat  <- read.table(paste(floc,'/DATA.DAT',sep=''),comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  pos1 <- find(c("#",'Catch','data'), tdat, 3)
  pos2 <- find(c("#",'Index','data'), tdat, -1)
  Names <- tdat[pos1-1,c(1:4)]
  tdat <- tdat[pos1:pos2,c(1:4)]
  names(tdat) <- Names
  tdat %<>% filter(`#Year`==max(`#Year`)) %>% mutate(prop=as.numeric(catch)/sum(as.numeric(catch))) %>%
    mutate(catch=round(projectcatch*1000*prop,1)) %>% dplyr::select(-prop) %>% mutate(Hrate=dynamics$value[tolower(dynamics$object)=='projectedhr'])

  tmp <- c(tmp, "\n# Specifications for projections (1=Catch;2=HarvestRate)\n",dynamics$value[tolower(dynamics$object)=='whichproject'],"\n#\n# Catch data (kg) / Harvest Rate - Number of observations\n", projectseason*nrow(tdat), "\n")
  tmp <- c(tmp,paste(colnames(tdat), collapse = "\t"), "\n")

  for(i in 1:projectseason){
    tdat$`#Year` <- endseason+i
    for(j in 1:nrow(tdat)){
      tmp <- c(tmp, paste(tdat[j,], collapse = "\t"),"\n")
    }}

  write.table(tmp, paste(floc,'/PROJECTIONS.DAT',sep=''), sep="", row.names = F, col.names = F, quote=F)
  invisible(floc)
}


