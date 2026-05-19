#' Generate Comprehensive HTML Report from IMuLT Model Outputs
#'
#' Master function to create a complete HTML report summarizing IMuLT model results.
#' Reads model output files, creates diagnostic plots and tables for all model
#' components, and compiles them into an organized HTML document with navigation.
#'
#' @param is95 Logical indicating whether to use 95% confidence intervals (TRUE)
#'   or 68% confidence intervals (FALSE, approximately ±1 SE). Default is TRUE
#'
#' @param folder_name The name for the folder to contain the outputs (inside Summary).
#' If omitted or left blank it will revert to the default behaviour of storing outputs
#' in 'summary/result/'.
#'
#' @param openfile Whether to open the html file on completion. Default TRUE
#'
#' @return Invisibly returns NULL. Side effects include:
#' \itemize{
#'   \item Creates 'Summary/result/' directory with HTML report and all figures
#'   \item Optionally archives previous report to timestamped folder
#'   \item Opens HTML report in default browser
#'   \item Copies input DAT files and key output files to report directory
#' }
#'
#' @details
#' This function is the main post-processing workflow for IMuLT stock assessments.
#' It creates a comprehensive, organized HTML report containing:
#'
#' **Model Diagnostics**:
#' \itemize{
#'   \item Likelihood components (raw and weighted)
#'   \item Penalty terms (initial N, recruitment, recruitment smoothing)
#'   \item Estimated parameter table with gradients and bounds
#'   \item Data tuning statistics for size compositions (Francis multipliers)
#' }
#'
#' **Data and Model Specifications**:
#' \itemize{
#'   \item Data availability plot showing temporal coverage by fleet/area
#'   \item Fleet descriptions and area definitions
#'   \item Selectivity curves by sex and fleet
#'   \item Retention curves (legal size regulations) by sex, fleet, and time
#'   \item Growth curves derived from size-transition matrices
#' }
#'
#' **Model Fits to Observations**:
#' \itemize{
#'   \item Commercial catches (observed vs predicted by area and total)
#'   \item Abundance indices (CPUE) with confidence intervals
#'   \item Fishing efficiency trends (technological creep)
#'   \item Size composition fits (multiple visualizations: overlaid, residuals, time series)
#'   \item Puerulus settlement data fits
#' }
#'
#' **Population Dynamics Outputs**:
#' \itemize{
#'   \item Recruitment estimates by area and year
#'   \item Recruitment size distributions
#'   \item Movement/migration patterns between areas
#'   \item Legal biomass (B/B0) trajectories with reference points
#'   \item Harvest rates by management zone
#'   \item Fishing mortality by area and fleet
#'   \item Egg production (spawning stock) with breeding stock management area summaries
#'   \item Natural mortality (including density dependence if applicable)
#'   \item Virgin and initial population size structures
#' }
#'
#' **Report Organization**: Uses the makehtml/hplot package system to create
#' a navigable HTML interface with categories (sidebar menu) containing related
#' plots and tables. Each plot includes a descriptive caption.
#'
#' **Archiving**: Prompts user whether to archive previous report. If 'Y',
#' copies existing report to timestamped 'archive' folder before creating new report.
#'
#' **File Dependencies**: Requires the following files in working directory:
#' \itemize{
#'   \item Output.RL - Main model outputs
#'   \item SDReport.RL - Parameter estimates and standard errors
#'   \item Echo.out - Echo of specifications
#'   \item Parameters_solved.txt - Table of estimated parameters
#'   \item ../DATA.DAT, ../CONTROL.DAT, ../SELEXSPEC.DAT, ../MOVESPEC.DAT - Input files
#'   \item ../../ModelStructure.xlsx - Model configuration workbook
#' }
#'
#' **Special Processing**:
#' \itemize{
#'   \item Automatically detects 8-area model for Breeding Stock Management Area plots
#'   \item Handles variable sex structure (combined, female only, or both sexes)
#'   \item Adapts plot panels based on number of areas/fleets using Fdims() helper
#'   \item Converts between different data representations (numbers ↔ proportions)
#'   \item Handles missing standard errors gracefully
#' }
#'
#' The function loads required packages: makehtml, hplot, dplyr, magrittr, tidyr,
#' ggplot2, reshape2, and openxlsx.
#'
#' @note This function should be run from the model output directory (typically
#' 'Run/') after successful model estimation. It expects a specific directory
#' structure with input files one level up (../) and ModelStructure.xlsx two
#' levels up (../../).
#'
#' @examples
#' \dontrun{
#' # After running IMuLT model, generate report with 95% CIs
#' setwd("Run/")
#' MakeOutPut(is95 = TRUE)
#'
#' # Generate report with narrower 68% CIs
#' MakeOutPut(is95 = FALSE)
#' }
#'
#' @seealso \code{\link{LoadOutputData}} for loading outputs without report generation
#'
#' @export
MakeOutPut <- function(is95=TRUE,folder_name='',openfile=TRUE){

  library(makehtml, quietly = T)
  library(hplot, quietly = T) # for plotprep and parset; automates the use of png
  library(dplyr, quietly = T)
  library(magrittr, quietly = T)
  library(tidyr, quietly = T)
  library(ggplot2, quietly = T)
  library(dplyr, quietly = T)
  library(reshape2, quietly = T)
  library(openxlsx, quietly = T)

  options(dplyr.summarise.inform = FALSE) ## Removes useless dplyr warnings
  starttime <- as.character(Sys.time())

  SclErr <- ifelse(is95,1.96,0.842)

  ddir <- filenametopath(getwd(),"")
  indir <- filenametopath(ddir,"Summary")

  if (folder_name == '') {
    # Ask what to call this run
    run_name <- dlg_input("Name for this model run (blank = 'result'):")$res
    if (run_name == '') run_name <- 'result'

    rundir <- filenametopath(indir, run_name)

    # Only ask about archiving if that folder already exists with files
    if (dir.exists(rundir)) {
      f <- list.files(rundir, include.dirs = FALSE, full.names = TRUE, recursive = TRUE)
      if (length(f) > 0) {
        Archive <- toupper(dlg_input(
          paste0("'", run_name, "' already exists. Archive it first? (Y or N)")
        )$res)
        if (Archive == 'Y') {
          ctime <- gsub(' ', '', gsub(".", "", format(Sys.time(), '%Y.%m.%d %H.%M'), fixed = TRUE))
          rundirA <- paste0(indir, '/archive', ctime)
          invisible(dir.create(rundirA))
          invisible(file.copy(rundir, rundirA, recursive = TRUE))
          print("Archived old report")
        }
        suppressWarnings(invisible(file.remove(f)))
      }
    }
  } else {
    rundir <- filenametopath(indir, folder_name)
  }

  dirExists(rundir,verbose=TRUE)  ## This makes it

  # Copy correlation matrix if it exists
  corfile_src <- file.path(getwd(), "CorrelationMatrix.csv")
  #print(paste("Looking for correlation file at:", corfile_src))
  #print(paste("Exists:", file.exists(corfile_src)))
  if (file.exists(corfile_src)) {
    file.copy(corfile_src, filenametopath(rundir, "CorrelationMatrix.csv"), overwrite = TRUE)
  }

  analysis <- "IMuLT"
  resfile <- setuphtml(rundir=rundir) # creates resultTable.csv in rundir

  #### Open up all data ####
  dat  <- read.table("Output.RL",comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
  echo  <- read.table("Echo.out",comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
  selx <- read.table(paste("../SELEXSPEC.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
  lbin1  <- read.table(paste("../DATA.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
  ctl1  <- read.table(paste("../CONTROL.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
  mov1  <- read.table(paste("../MOVESPEC.DAT",sep=''),comment.char = "?",fill=T,blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)

  file_path <- "../../ModelStructure.xlsx"
  # Test workbook and inform user that it is open
  is_file_locked <- function(path) {
    tryCatch({
      con <- file(path, open = "a")  # "a" = append, requires exclusive access
      close(con)
      FALSE  # not locked
    }, error = function(e) {
      TRUE   # locked / open elsewhere
    })
  }

  if (is_file_locked(file_path)) {
    stop("ModelStructure.xlsx is currently open in another application. Please close it and try again.")
  } else {  wb <- loadWorkbook(file = file_path)  }

  #This is the location of the data input files and their associated parameters
  dynamics <- readWorkbook(wb,sheet='dynamics', startRow = 2)
  startseason <- as.numeric(dynamics$value[dynamics$object=='startseason'])
  endseason <- as.numeric(dynamics$value[dynamics$object=='endseason'])
  yr <- startseason:endseason
  projectseason <- as.numeric(dynamics$value[dynamics$object=='projectedseason'])
  projectcatch <- as.numeric(dynamics$value[dynamics$object=='projectedcatch'])
  burnin <- as.numeric(dynamics$value[dynamics$object=='burnin'])
  ages <- as.numeric(dynamics$value[dynamics$object=='ages'])
  sexs <- 1:as.numeric(dynamics$value[dynamics$object=='sexs'])
  areas <- readWorkbook(wb,sheet='area', startRow = 2)
  nareas <- length(unique(areas$AreaCode))
  times <- readWorkbook(wb,sheet='times', startRow = 2)
  fleets <- readWorkbook(wb,sheet='fleetcode', startRow = 2)

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
  Pout <- pout

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

  suppressWarnings(fleetarea <- fleets %>% mutate(areaname=do.call('rbind',strsplit(description,'_'))[,1],fleettype=group))

  ####
  ## Model run statistics
  ### Likelihoods ###
  print("Making Likelihoods")
  like <- dat[4:10,c(1,3,5)]
  colnames(like) <- c('Id','Raw LL', 'Weighted LL')
  like <- rbind(like,like[1,])
  like[1,3] <- like[1,2]
  like[nrow(like),] <- c('Total', NA,dat[3,4])
  suppressWarnings(like[,2] <- round(as.numeric(like[,2]),2))
  suppressWarnings(like[,3] <- round(as.numeric(like[,3]),2))
  row.names(like) <- 1:nrow(like)
  like[is.na(like)] <- ''
  filen <- "Likelihood.csv"  # csv files only
  addtable(intable=like,filen=filen,rundir=rundir,category="Like",
           caption="Raw and weighted likelihoods.")

  ### Penalties ####
  INP <- as.numeric(dat[11,4])
  RP <-  as.numeric(dat[12,3])
  RSP <-  as.numeric(dat[13,4])
  Prior_Main <-  as.numeric(dat[14,4])
  Prior_Rec <-  as.numeric(dat[15,4])
  Prior_Sel <-  as.numeric(dat[16,4])
  Prior_Eff <-  as.numeric(dat[17,4])
  dfram <- data.frame(id=c('Initial Nunbers','Recruitment Deviations', 'Recruitment Devs Smoother','Priors on Main Pars','Priors on Recruit Pars','Priors on Selectivity Pars','Priors on Efficiency Pars'),Value=c(round(INP,1),round(RP,1),round(RSP,1),round(Prior_Main,1),round(Prior_Rec,1),round(Prior_Sel,1),round(Prior_Eff,1)))
  filen <- "Penality.csv"  # csv files only
  addtable(intable=dfram,filen=filen,rundir=rundir,category="Like",
           caption="Penalities added to likelihoods.")

  ### Parameters: Main and those estimated ####
  mainP <- pout[grepl('MainPars', pout$name),]
  ctl2 <- findNclean(c('#','Basic','parameters'), ctl1, 1, char=T)
  mpnames <- ctl2$X10
  mainP$name <- mpnames[as.numeric(do.call('rbind', strsplit(mainP$name, '_'))[,2])]
  mainP$est <- round(mainP$est,2)
  #mainP <- mainP[c(1,nrow(mainP)),]
  #rownames(mainP) <- 1:2
  #pander(mainP)

  pout <- read.delim(paste("Parameters_solved.txt",sep=''),sep='\t',stringsAsFactors =F)
  pout$Parameter <- ifelse(pout$Parameter=='MainPars', mpnames[pout$Number],pout$Parameter)
  un_names <- unique(c(mainP$name,pout$Parameter))
  mxlen <- length(un_names)
  ntab <- data.frame(Parameter=rep(as.character('.'),mxlen), Value=as.character('.'), Estimated=as.character('N'),stringsAsFactors =F)
  for(i in 1:nrow(ntab)){
    ntab$Parameter[i] <- un_names[i]
    ntab$Value[i] <- ifelse(nrow(mainP)>=i, mainP$est[i],'.')
    ntab$Estimated[i] <- ifelse(un_names[i]%in%pout$Parameter, 'Y','N')
  }

  ### Model Run Comments ####
  #init <- findNclean('Initiation', dat, 0)
  #init <- init[!is.na(init)]
  #txt <- paste('Initiation option is ', init,'.\n', sep='')

  ttxt1 <- findNclean(c('#','Burn-in', 'whole'), lbin1, 1, T); txt2 <- paste('Initiate model years: ', paste(ttxt1, collapse=' '),'.\n', sep='')
  ttxt1.1 <- findNclean(c('#','Burn-in', 'for'), lbin1, 1, T); txt2.1 <- paste('Burn-in years: ', paste(ttxt1.1, collapse=' '),'.\n', sep='')

  ttxt2 <- findNclean(c('#','Loop', 'counter'), lbin1, 1); txt3 <- paste('Loops to refine initial F: ', ttxt2,'.\n', sep='')
  ttxt3 <- findNclean(c('#','Years','over'), lbin1, 1, T); txt4 <- paste('Number of years to base initial F on: ', paste(ttxt3, collapse=' '),'\n', sep='')

  #### Data in ####
  print("Making Data Summary")
  CAtch <- find('Catch data by', echo, 1); ECAtch <- find('Number of cpue', echo, -1)
  CAtch <- echo[CAtch:ECAtch,]
  if(dim(CAtch)[1]==1){
    ttmp <- rep(0, dim(CAtch)[2])
    CAtch <- rbind(CAtch, ttmp)  }
  CAtch <- CAtch[,!is.na(CAtch[1,])&CAtch[1,]!='']
  CAtch <- apply(CAtch,2,as.numeric)
  nyears <- 1+(as.numeric(echo[find('Year2', echo, 0),2]) - as.numeric(echo[find('Year1', echo, 0),2])) + as.numeric(echo[find('MaxProjYr', echo, 0),2])
  dimnames(CAtch)[[2]] <- (as.numeric(echo[find('Year1', echo, 0),2]):(as.numeric(echo[find('Year1', echo, 0),2])+(nyears-1)))
  ntstep <- as.numeric(echo[find('Number of time', echo, 0),5])
  nfleet <- as.numeric(echo[find('Number of fleets', echo, 0)[1],4])
  idmat <- expand.grid(fleet=1:nfleet, tstep=1:ntstep) %>% arrange(fleet)
  CAtch <- as.data.frame(cbind(idmat,CAtch)) %>% pivot_longer(!c(fleet,tstep),names_to = 'year') %>% mutate(group=fleets$group[match((fleet),fleets$fleet)], area=fleets$newarea[match((fleet),fleets$fleet)], time=round(as.numeric(year)+(tstep-1)/(max(tstep)),4), id=paste(group,area)) %>% group_by(time,id) %>% summarise(num=sum(value)) %>%  filter(num>0) %>% mutate(type=paste('Catch')) # %>% pivot_wider(names_from = 'time', , values_from = num)

  CPue <- find('Cpue data', echo, 1); ECPue <- find('Number of Q-related', echo, -1)
  CPue <- echo[CPue:ECPue,]
  CPue <- CPue[,!is.na(CPue[1,])&CPue[1,]!='']
  names(CPue) <- c('Ind','fleet','sex','year','tstep','obs','cv')
  CPue <- as.data.frame(apply(as.matrix(CPue),2,as.numeric))
  CPue %<>% mutate(year=as.numeric(year)+startseason, group=fleets$group[match((fleet+1),fleets$fleet)], area=fleets$newarea[match((fleet+1),fleets$fleet)], time=round(year+tstep/(max(tstep)+1),4), id=paste(group,area)) %>% group_by(time,id) %>% summarise(num=length(obs))  %>% mutate(type=paste('Index')) # %>% pivot_wider(names_from = 'time', , values_from = num)

  LEngth <- find('Size-composition data', echo, 1); ELEngth <- find('Puerulus data', echo, -1);
  LEngth <- echo[LEngth:ELEngth,]
  LEngth <- LEngth[,1:5]
  names(LEngth) <- c('fleet','sex','year','tstep','obs')
  LEngth <- as.data.frame(apply(as.matrix(LEngth),2,as.numeric))
  LEngth %<>% mutate(year=as.numeric(year)+startseason, group=fleets$group[match((fleet+1),fleets$fleet)], area=fleets$newarea[match((fleet+1),fleets$fleet)], time=round(year+tstep/(max(tstep)+1),4), id=paste(group,area)) %>% group_by(time,id) %>% summarise(num=length(obs))  %>% mutate(type=paste('Size-comp')) # %>% pivot_wider(names_from = 'time', , values_from = num)

  # NUmbers

  # PUerulus
  PUer <- find('Puerulus data', echo, 1); EPUer <- find('Environmental data', echo, -2)
  PUer <- echo[PUer:EPUer,]
  PUer <- PUer[,!is.na(PUer[1,])&PUer[1,]!='']
  if(nrow(PUer)>4){
    names(PUer) <- c('area','year','obs','cv')
    PUer <- as.data.frame(apply(as.matrix(PUer),2,as.numeric)) %>% mutate(area=area+1, year=year+(startseason-burnin-1))
    PUer %<>% mutate(time=year, id=paste('survey',area)) %>% group_by(time,id) %>% summarise(num=length(obs)) %>% mutate(type=paste('Larval')) #%>% pivot_wider(names_from = 'time', , values_from = num)
  }
  adat <- rbind(CAtch, CPue, LEngth, PUer) %>% mutate(id2 = paste(id,type)) %>% filter(!is.na(time)) %>% as.data.frame()
  unid <- unique(adat$id2)
  lab <- as.data.frame(do.call('rbind',strsplit(unid,' ')))
  names(lab) <- c('source','area','type')
  lab %<>% mutate(order=1:nrow(lab), unid=unid, id1 = paste(source,type)) %>% filter(!is.na(source))
  unid2 <- unique(lab$id1)
  lab %<>% mutate(Col=topo.colors(length(unid2), alpha = 1, rev = T)[match(id1,unid2)])
  adat$Col <- lab$Col[match(adat$id2, lab$unid)]
  adat$yax <- match(adat$id2,unid)
  adat %<>% mutate(id3 = paste(id,type))

  filename <- filenametopath(rundir,"DataIn.png")
  plotprep(width=8, height=8, filename=filename, cex=1.2, verbose=FALSE)
  parset(plots=c(1,1), margin=c(0.45, 1.6, 0.15, 1.6), cex=1.2)

  plot(adat$time, adat$yax, col=1, bg=adat$Col, pch=21, cex=1.8,
       axes=FALSE, xlab='', ylab='', xlim=range(adat$time, na.rm=TRUE))
  abline(h=1:max(adat$yax), col="grey85", lty=1)
  points(adat$time, adat$yax, col=1, bg=adat$Col, pch=21, cex=1.8)

  lab1 <- lab %>% group_by(source, type, id1) %>% summarise(mnpos1=mean(order))
  mtext(side=4, at=lab1$mnpos1, lab1$source, las=1, cex=0.9, font=2)
  mtext(side=2, at=lab$order, lab$area, las=1, cex=0.8)
  mtext(side=2, line=3, at=lab1$mnpos1, lab1$type, las=1, cex=0.9, font=2)
  axis(1, cex.axis=0.9)
  box()
  caption <- "Data loaded into the model as recorded in the echo file."
  addplot(filen=filename,rundir=rundir,category="Data",caption=caption)

  ## Add Fleet descriptions
  fleetareatmp <- fleetarea %>% dplyr::select(Area=areaname, Fleet=fleet, FleetType=fleettype)
  addtable(intable=fleetareatmp,filen="Fleets.csv",rundir=rundir,category="Data",caption="Fleet Descriptions")

  ### Selectivity
  print("Making Selectivity and Retention Plots")
  fleet2 <- findNclean(c('#', 'Sex','Age', 'Fleet'), selx, 1, char=F)
  sel <- findNclean(c('Full','Selectivity'), dat, 1)
  sel <- sel[,2:ncol(sel)]
  ids <- findNclean(c('#','Selectivity','Parameters'), selx, 2, char = T)
  infpos <- which(grepl('inflect',ids[1,]))
  ids$name <- paste(ids[,(infpos+1)],ids[,(infpos+2)])
  ids <- data.frame(link=NA, name=unique(ids$name))
  ids$link <- as.numeric(row.names(ids))-1
  fleet3 <- fleet2 %>% pivot_longer(!c(Sex, Age, Fleet, `Step:`), names_to = 'year', values_to = 'link') %>% left_join(ids, by='link')
  fleet4 <- fleet3 %>% group_by(Sex, Fleet, link,name) %>% summarise(minyr=min(year), tsteps=paste(unique(`Step:`),collapse='.')) %>% mutate(name2=paste(name, minyr))

  for(f in unique(fleet4$Fleet)){
    if(length(unique(fleet4$Sex))==1) { fleet4$sex <- 'Sex 1' } else { fleet4$sex <- c('F','M')[(fleet4$Sex+1)]}
    fleet5 <- fleet4 %>% filter(Fleet==f)
    filename <- filenametopath(rundir,paste0('Fleet.',(f+1),"_Selectivity.png"))
    plotprep(width=10,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=c(1,1))
    sel2 <- sel[(1+fleet5$link),]
    sel_long <- sel2 %>% rename(a0 = Selectivity) %>%                        # rename to bin 1 (or a0 as bin 1)
      mutate(name2 = fleet5$name2, sex=fleet5$sex) %>%
      pivot_longer(cols = starts_with("a"), names_to = "bin", values_to = "selectivity") %>%
      mutate(bin = as.numeric(gsub("a", "", bin)) + 1, bin2 = lbin[bin])
    p <- ggplot(sel_long, aes(x = bin2, y = selectivity, colour = name2)) +
      geom_line(linewidth = 1) +
      scale_colour_viridis_d(name = NULL) +
      #scale_x_continuous(breaks = seq(0, 40, 5)) +
      labs(x = "LengthBin (midpoint)", y = "Selectivity") +
      facet_wrap(~sex)+
      theme_bw() +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(nrow = 2))
    print(p)
    caption <- paste("Selectivity curves estimated by the model.")
    addplot(filen=filename,rundir=rundir,category="Selectivity_Retenion",caption=caption)
  }

  ####  Retention ####
  ret <- findNclean(c('#Legal','Selectivity','by','sex'), dat, 1)  #Legal Selectivity by
  names(ret) <- c('sex','age','fleet','year','tstep',paste0('lb',1:(ncol(ret)-5)))
  for(ft in sort(unique(ret$fleet))){
    for(sx in sort(unique(ret$sex))){
      for(ag in sort(unique(ret$age))){
        tfleet1 <- ret[ret$sex==sx & ret$fleet==ft & ret$age==ag,]
        tfleet1 <- tfleet1[!duplicated(apply(as.matrix(tfleet1[,6:ncol(tfleet1)]),1,paste0,collapse=' ')), ]
        if(nrow(tfleet1)>4){
          tfleet2 <- tfleet1 %>% tidyr::pivot_longer(cols = starts_with("lb"),names_to = "length_bin",names_prefix = "lb",names_transform = list(length_bin = as.integer),values_to = "proportion") %>% mutate(length_bin = lbin[length_bin], id=as.factor(year), Age=paste('Age',age), Tstep=paste('Tstep',tstep))
          filename <- filenametopath(rundir,paste("Sex",sx,"Fleet",ft,"Age",ag,"Retention",".png"))
          plotprep(width=10,height=14,filename=filename,cex=0.9,verbose=FALSE)
          parset(plots=Fdims(length(unique(tfleet2$tstep))), margin = c(0.5,0.5,0.25,0.05))
          suppressWarnings(print(ggplot(tfleet2, aes(x=length_bin,y=proportion,colour=id))+
                                   geom_line()+ylab('Proportion')+xlab("Length bin (mm)")+
                                   scale_color_discrete(name = "Year") + scale_y_continuous(limits = c(0, 1))+
                                   theme_bw() + ggtitle(paste('Fleet',ft,'Sex',sx,'Age',ag)) +
                                   facet_wrap(~Tstep)+
                                   theme(legend.position = "top", strip.background = element_rect(fill = "white"),
                                         panel.grid.minor = element_blank())))
          caption <- paste("Retention curves for fleet", ft,",sex", sx,"and age",ag,".")
          addplot(filen=filename,rundir=rundir,category="Selectivity_Retenion",caption=caption)

        } else {
          tfleet1 <- ret[ret$sex==sx & ret$fleet==ft,]
          tfleet1 <- tfleet1[!duplicated(apply(as.matrix(tfleet1[,6:ncol(tfleet1)]),1,paste0,collapse=' ')), ]
          tfleet2 <- tfleet1 %>% tidyr::pivot_longer(cols = starts_with("lb"),names_to = "length_bin",names_prefix = "lb",names_transform = list(length_bin = as.integer),values_to = "proportion") %>% mutate(length_bin = lbin[length_bin], id=as.factor(year), TstepAge=paste("Ts",tstep,"Age",age))
          filename <- filenametopath(rundir,paste("Sex",sx,"Fleet",ft,"Retention",".png"))
          plotprep(width=10,height=14,filename=filename,cex=0.9,verbose=FALSE)
          parset(plots=Fdims(length(unique(tfleet2$tstep))), margin = c(0.5,0.5,0.25,0.05))
          suppressWarnings(print(ggplot(tfleet2, aes(x=length_bin,y=proportion,colour=id))+
                                   geom_line()+ylab('Proportion')+xlab("Length bin (mm)")+
                                   scale_color_discrete(name = "Year") + scale_y_continuous(limits = c(0, 1))+
                                   theme_bw() + ggtitle(paste('Fleet',ft,'Sex',sx,'Age',ag)) +
                                   facet_wrap(~TstepAge)+
                                   theme(legend.position = "top", strip.background = element_rect(fill = "white"),
                                         panel.grid.minor = element_blank())))
          caption <- paste("Retention curves for fleet", ft,"and sex", sx,".")
          addplot(filen=filename,rundir=rundir,category="Selectivity_Retenion",caption=caption)

        }
      }}}

 #### Growth ####
  print("Making Growth Curves")
  #  grow <- findNclean(c('#Growth','Curves'), dat, 2)
  #  head(grow)
  #  num <- length(unique(grow$sex))*length(unique(grow$area))

  find_unique_growth_years <- function(GrowthPnt) {
    dims  <- dim(GrowthPnt)
    Narea <- dims[1]; Nsex <- dims[2]; Nyear <- dims[4]

    results <- list()
    for (a in 1:Narea) {
      for (s in 1:Nsex) {
        gmat         <- GrowthPnt[a, s, 1, , ]
        unique_years <- which(!duplicated(gmat))
        results[[paste0("a", a, "_s", s)]] <- unique_years
      }
    }
    results
  }

  compound_growth_trajectory <- function(Data, GrowthPnt, length_midpoints = Data$MidLenBin[1,1:Data$Nlen[1]], Nyear_traj = 30) {

    dims  <- dim(GrowthPnt)
    Narea <- dims[1]; Nsex <- dims[2]; Nstep <- dims[5]
    nlbin <- dim(Data$TransInp)[2]

    unique_years <- find_unique_growth_years(GrowthPnt)

    all_results <- list()

    for (a in 1:Narea) {
      for (s in 1:Nsex) {
        sexlabels <- c("Female", "Male")
        if(Nsex==1) sexlabels <- c("Male", "Male")
        key           <- paste0("a", a, "_s", s)
        pattern_years <- unique_years[[key]]

        for (py in pattern_years) {

          dist <- c(1, rep(0, nlbin - 1))

          mean_len  <- numeric(Nyear_traj)
          modal_len <- numeric(Nyear_traj)
          sd_len    <- numeric(Nyear_traj)

          for (y in 1:Nyear_traj) {

            for (t in 1:Nstep) {
              ptr <- GrowthPnt[a, s, 1, py, t]
              if (ptr >= 0) {
                stm  <- Data$TransInp[ptr + 1, , ]
                dist <- as.numeric(stm %*% dist)
              }
            }

            # Normalise just in case
            dist_norm <- dist / sum(dist)

            mean_len[y]  <- sum(dist_norm * length_midpoints)
            modal_len[y] <- length_midpoints[which.max(dist_norm)]
            sd_len[y]    <- sqrt(sum(dist_norm * (length_midpoints - mean_len[y])^2))
          }

          all_results[[paste0(key, "_py", py)]] <- data.frame(
            age          = 1:Nyear_traj,
            mean_len     = mean_len,
            modal_len    = modal_len,
            sd_len       = sd_len,
            lo_len       = mean_len - sd_len,
            hi_len       = mean_len + sd_len,
            area         = factor(a),
            sex          = factor(s, labels = sexlabels[s]),
            pattern_year = factor(py)
          )
        }
      }
    }

    do.call(rbind, all_results)
  }

  # Run
  growth_traj <- compound_growth_trajectory(Data, Data$GrowthPnt, Nyear_traj = 30)

  growth_traj %<>% mutate(year=factor((Data$Year1:Data$Year2)[as.numeric(as.character(pattern_year))]))

  for(s in unique(growth_traj$sex)){
    filename <- filenametopath(rundir,paste(s,"Growth_Curves1.png"))
    plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(1))
    # Plot with ribbon for +/- 1 SD
    print(ggplot(growth_traj[growth_traj$sex==s,], aes(x = age, colour = year, fill = year)) +
            geom_ribbon(aes(ymin = lo_len, ymax = hi_len), alpha = 0.15, colour = NA) +
            geom_line(aes(y = mean_len), linewidth = 0.8) +
            facet_wrap( ~ area) +
            labs(
              x      = "Age (years since recruitment)",
              y      = "Mean length (mm)",
              colour = "Year first seen",
              fill   = "Year first seen",
              title  = "Growth by area"
            ) +
            theme_bw())

    caption <- paste("Inputted growth trajectories by model areas for",s)
    addplot(filen=filename,rundir=rundir,category="Growth",caption=caption)}

  filename <- filenametopath(rundir,"Growth_Curves2.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(1))
  print(ggplot(growth_traj, aes(x = age, y = mean_len, colour = area, linetype = year)) +
          geom_line(linewidth = 0.8) +
          facet_wrap(~ sex) +
          labs(
            x            = "Age (years since 1st Length bin)",
            y            = "Mean length (mm)",
            colour       = "Model area",
            title        = "Growth by sex"
          ) + theme_bw())
  caption <- "Inputted growth trajectories between model areas."
  addplot(filen=filename,rundir=rundir,category="Growth",caption=caption)

  #### Fit to Data ####
  #### Commercial Catches ####
  print("Making Model fit to Catch")
  catch <- findNclean('#Catches', dat, 0)
  Zcatch <- catch %>% group_by(Year, Area) %>% summarise(obs=sum(Observed), est=sum(Predicted))
  graphrange <- range(Zcatch$Year)

  filename <- filenametopath(rundir,"Total Catches.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(1))
  tZcatch <- Zcatch %>% group_by(Year) %>% summarise(obs=sum(obs),est=sum(est))
  ymx <- pretty(trunc(range(c(tZcatch$obs,tZcatch$est))/10)*10)
  ymn<- ymx[1]; ymx <- max(ymx)
  suppressWarnings(with(tZcatch, plot(Year, obs/1000, type='o', pch=16,cex=0.7,ylim=c(ymn,ymx)/1000, axes=F, ylab='Total Catch (t)', xlab='Fishing Season', lty=1, main='Total')))
  with(tZcatch, lines(Year, est/1000, type='o',cex=0.7,col=2, lty=1))
  axis(1, seq(min(Zcatch$Year), max(Zcatch$Year),2))
  axis(2,seq(ymn/1000,ymx/1000,length=5), las=1)
  caption <- "Observed (black) and estimated (red 95% CI grey) total commercial catches."
  addplot(filen=filename,rundir=rundir,category="Catches",caption=caption)

  filename <- filenametopath(rundir,"Cumulative Catches.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(1))
  Zcatch %<>% mutate(Location=fleetarea$areaname[match(Area,fleetarea$newarea)], CatchT=obs/1000)
  print(ggplot(Zcatch, aes(fill=Location, y=CatchT, x=Year)) +
          viridis::scale_fill_viridis(discrete = T) +
          geom_bar(position="stack", stat="identity") +
          ylab('Catch (t)') +
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.0, angle = 45)))
  caption <- "Commercial catches by location"
  addplot(filen=filename,rundir=rundir,category="Catches",caption=caption)

  filename <- filenametopath(rundir,"Catches2.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  Zcatch %<>% mutate(arean=fleetarea$areaname[match(Area, fleetarea$newarea)])  %>% pivot_longer(col=c(obs,est), names_to = 'type') %>% mutate(`Catch (t)` = value/1000)
  print(ggplot(Zcatch, aes(x=Year, y=`Catch (t)`,colour=type))+
          geom_line()+geom_point()+
          facet_wrap(~arean)+
          scale_color_manual(values=c("red","black")) +
          scale_size_manual(values = c(1, 0.5)) +
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.0, angle = 45)))
  caption <- "Observed (black) and estimated (red 95% CI grey) commercial catches by Model area on same scale."
  addplot(filen=filename,rundir=rundir,category="Catches",caption=caption)

  #### Discards ####
  print("Making Model Discard plots")
  disc_raw <- findNclean(c('#','Discards'), dat, 0)
  colnames(disc_raw) <- c("Year","Step","Fleet","DiscardWt","DeadDiscardWt")
  disc_raw$Area <- fleetarea$newarea[disc_raw$Fleet]
  disc_raw$AreaName <- fleetarea$areaname[disc_raw$Fleet]

  # ── Global: Total Discards & Dead Discards ──
  disc_total <- disc_raw %>%
    group_by(Year) %>%
    summarise(Discard = sum(DiscardWt)/1000, DeadDiscard = sum(DeadDiscardWt)/1000)

  filename <- filenametopath(rundir,"Total Discards.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(1))
  ymx <- pretty(c(0, max(disc_total$Discard)))
  ymn <- ymx[1]; ymx <- max(ymx)
  with(disc_total, plot(Year, Discard, type='o', pch=16, cex=0.7,
                        ylim=c(ymn, ymx), axes=F, ylab='Discard Weight (t)',
                        xlab='Fishing Season', lty=1, main='Total Discards', col=1))
  with(disc_total, lines(Year, DeadDiscard, type='o', pch=16, cex=0.7, col=2, lty=1))
  axis(1, seq(min(disc_total$Year), max(disc_total$Year), 2))
  axis(2, seq(ymn, ymx, length=5), las=1)
  legend("topright", legend=c("Total Discards","Dead Discards"),
         col=c(1,2), lwd=2, pch=16, bty="n", cex=0.85)
  caption <- "Total discards (black) and dead discards (red) by fishing season."
  addplot(filen=filename,rundir=rundir,category="Discards",caption=caption)

  # ── Stacked bar: Discards by Area ──
  disc_area <- disc_raw %>%
    group_by(Year, AreaName) %>%
    summarise(DiscardT = sum(DiscardWt)/1000, .groups="drop")

  filename <- filenametopath(rundir,"Discards by Area stacked.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(disc_area, aes(fill=AreaName, y=DiscardT, x=Year)) +
          viridis::scale_fill_viridis(discrete = T) +
          geom_bar(position="stack", stat="identity") +
          ylab('Discard Weight (t)') +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.0, angle = 45)))
  caption <- "Total discards by location."
  addplot(filen=filename,rundir=rundir,category="Discards",caption=caption)

  # ── By Area: Discards & Dead Discards faceted ──
  disc_area2 <- disc_raw %>%
    group_by(Year, AreaName) %>%
    summarise(Discard = sum(DiscardWt)/1000,
              DeadDiscard = sum(DeadDiscardWt)/1000, .groups="drop") %>%
    pivot_longer(cols=c(Discard, DeadDiscard), names_to='Type', values_to='Weight_t')

  filename <- filenametopath(rundir,"Discards by Area.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(disc_area2, aes(x=Year, y=Weight_t, colour=Type)) +
          geom_line() + geom_point() +
          facet_wrap(~AreaName) +
          scale_color_manual(values=c("red","black"),
                             labels=c("Dead Discards","Total Discards")) +
          ylab('Discard Weight (t)') +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.0, angle = 45)))
  caption <- "Total discards (black) and dead discards (red) by area."
  addplot(filen=filename,rundir=rundir,category="Discards",caption=caption)

  # ── By Fleet: Discards & Dead Discards faceted ──
  disc_fleet <- disc_raw %>%
    group_by(Year, Fleet) %>%
    summarise(Discard = sum(DiscardWt)/1000,
              DeadDiscard = sum(DeadDiscardWt)/1000, .groups="drop") %>%
    mutate(FleetName = paste("Fleet", Fleet)) %>%
    pivot_longer(cols=c(Discard, DeadDiscard), names_to='Type', values_to='Weight_t')

  filename <- filenametopath(rundir,"Discards by Fleet.png")
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(disc_fleet, aes(x=Year, y=Weight_t, colour=Type)) +
          geom_line() + geom_point() +
          facet_wrap(~FleetName) +
          scale_color_manual(values=c("red","black"),
                             labels=c("Dead Discards","Total Discards")) +
          ylab('Discard Weight (t)') +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.0, angle = 45)))
  caption <- "Total discards (black) and dead discards (red) by fleet."
  addplot(filen=filename,rundir=rundir,category="Discards",caption=caption)

  #### Index data ####
  print("Making Model fit to Abundance Indices")
  tdat <- findNclean(c('Index','data'), dat, 2)
  cpuesd <- sdr[grepl('PredCpue', sdr$name),]

  if(nrow(cpuesd)>0) suppressWarnings(tdat <- cbind(tdat,cpuesd))
  if(nrow(cpuesd)==0) tdat1 <- tdat %>% mutate(yts=Year+(Time_step-1)/max(Time_step)) %>% group_by(Sex, Fleet,Year,Time_step) %>%  summarise(obs=mean(Observed), Lse=log(mean(Observed))*mean(Relative_CV), obsUp=exp(log(mean(Observed))+(Lse*SclErr)) ,obsLow=exp(log(mean(Observed))-(Lse*SclErr)), est=mean(Predicted), esd=NA,estlwr=est,estupr=est)
  if(nrow(cpuesd)>0)  tdat1 <- tdat %>% mutate(yts=Year+(Time_step-1)/max(Time_step)) %>% group_by(Sex, Fleet,Year,Time_step) %>%  summarise(obs=mean(Observed), Lse=log(mean(Observed))*mean(Relative_CV), obsUp=exp(log(mean(Observed))+(Lse*SclErr)) ,obsLow=exp(log(mean(Observed))-(Lse*SclErr)), est=mean(Predicted), estupr=exp(log(mean(Predicted))+(log(mean(Predicted))*mean(cv))*SclErr), estlwr=exp(log(mean(Predicted))-(log(mean(Predicted))*mean(cv))*SclErr))
  if(nrow(cpuesd)>0) tdat %<>% mutate(se=sum(SE))

  tdat2 <- tdat1 %>% pivot_longer(col=c(obs,est), names_to = 'type') %>% mutate(lwr=ifelse(type=='obs',obsLow,estlwr),upr=ifelse(type=='obs',obsUp,estupr))

  tdat2$Fleettype <- fleets$group[match(tdat2$Fleet, fleets$fleet)]
  tdat2$Area <- fleets$newarea[match(tdat2$Fleet, fleets$fleet)]
  tdat2$Descrip <- fleets$description[match(tdat2$Fleet, fleets$fleet)]
  tdat2 %<>% mutate(Areaname = unlist(strsplit(Descrip, '_'))[1])

  unfleet <- sort(unique(tdat2$Fleet))
  for(uf in unfleet)  {
    tdat3 <- tdat2 %>% filter(Fleet==uf) %>% mutate(aSex=recode_values(Sex, -1~'comb',0~'F',1~'M'))
    filename <- filenametopath(rundir,paste0("Fleet ",uf,".png"))
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=c(1,1))
    if(length(unique(tdat3$aSex))>1) {print(ggplot(tdat3, aes(x=Year, y=value, colour=type))+
                                              geom_line()+geom_point()+
                                              geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2)+
                                              facet_grid(Time_step~aSex)+
                                              scale_color_manual(values=c("red","black")) +
                                              scale_size_manual(values = c(0.5, 0.5)) +
                                              theme(panel.background = element_rect(fill = "white",colour = NA),
                                                    panel.border = element_rect(fill = NA, colour = "grey20"),
                                                    axis.text.x = element_text(vjust = 0.0, angle = 45),legend.position = 'bottom')+
                                              ylab('Catch rate (kg/pot)'))} else {
                                                {print(ggplot(tdat3, aes(x=Year, y=value, colour=type))+
                                                         geom_line()+geom_point()+
                                                         geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2)+
                                                         facet_wrap(~Time_step)+
                                                         scale_color_manual(values=c("red","black")) +
                                                         scale_size_manual(values = c(0.5, 0.5)) +
                                                         theme(panel.background = element_rect(fill = "white",colour = NA),
                                                               panel.border = element_rect(fill = NA, colour = "grey20"),
                                                               axis.text.x = element_text(vjust = 0.0, angle = 45),legend.position = 'bottom')+
                                                         ylab('Catch rate (kg/pot)'))}
                                              }
    caption <- paste(unique(tdat3$aSex), unique(tdat3$Areaname), "Observed (black) and estimated (red 95% CI grey) catch rates for each fleet and or timestep.")
    addplot(filen=filename,rundir=rundir,category="Index",caption=caption)
  }

  #### Fishing Efficiency ####
  print("Making Fishing Efficiency")
  filename <- filenametopath(rundir,paste0("Fishing_Efficiency.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  fcreep <- findNclean(c('Fishing','Efficiency'), dat, 2)
  fcreepsd <- sdr[grepl('CpueEcreep', sdr$name),]

  if(nrow(fcreepsd)>0) {fcreep <- cbind(fcreep,fcreepsd) %>% mutate(Match=paste(id,Year))} else{
    fcreep %<>% mutate(SE=0, upr=Predicted, lwr=Predicted, Match=paste(id,Year)) }
  tdat <- expand.grid(area=1:nareas, year=yr)
  fcre <- fleetarea %>% filter(group=='comm')
  tdat$fcreep <- fcre$effic.creep[match(tdat$area, fcre$newarea)]
  fcreep %<>% rename(est=Predicted, fcreep=id, year=Year) %>% dplyr::select(est,fcreep,year,lwr,upr)
  tdat %<>% left_join(fcreep, by=c('year', 'fcreep'))
  fleetarea1 <- fleetarea %>% group_by(newarea, areaname) %>% summarise(cnt=length(newarea))
  tdat %<>% mutate(locname=fleetarea1$areaname[match(area,fleetarea1$newarea)])
  suppressWarnings(print(ggplot(data=tdat, aes(year, est)) +
                           geom_errorbar(aes(ymin=lwr,ymax=upr),colour='grey70')+
                           geom_line() + geom_point(size=0.5)+
                           scale_color_manual(values = c("red")) +
                           #geom_errorbar(aes(ymin=mn-sd, ymax=mn+sd, color=name) , width=.2,position=position_dodge(0.05)) +
                           labs(x="Year",y="Fishing efficiency")+
                           facet_wrap(~locname) + theme_bw()))

  caption <- "Estimated (95% CI grey) compounding commercial fishing efficiency for each model area."
  addplot(filen=filename,rundir=rundir,category="Fishing_Efficiency",caption=caption)

  #### Recruitment ####
  ### Puerulus Data
  print("Making Recruitment")
  rec <- findNclean(c('Larval','data'), dat, 1)

  if(!is.null(dim(rec)))  {
    recsd <- sdr[grepl('Larval', sdr$name),]
    if(length(recsd$SE[!is.na(recsd$SE)])>0) { rec <- cbind(rec,recsd) } else {rec %<>% mutate(upr=Predicted, lwr=Predicted)}
    mxyr <- max(rec$Year)
    rec %<>% filter(Year<(mxyr-2)) %>% mutate(ObsUp=Observed+SD*SclErr, ObsLow=Observed-SD*SclErr, PredUp=upr, PredLow=lwr)
    filename <- filenametopath(rundir,paste0("Recruitmant_Index.png"))
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(length(unique(rec$Area))))
    for(i in sort(unique(rec$Area))){
      tdat2 <- rec[rec$Area==i,]
      Mx <- max(c(tdat2$Observed, tdat2$Predicted),na.rm=T)
      suppressWarnings(with(tdat2, plot(Year, Observed, type='o', axes=F, ylab='Predicted recruitment', xlab='Puerulus Settlement', lty=1, main=paste('Area',i),ylim=c(0,Mx), xlim=c(1970,graphrange[2])) ))
      suppressWarnings(with(tdat2, arrows(Year, ObsUp,y1=ObsLow,code=3,angle=90,length=0.025,col=grey(0.3,0.3))))
      with(tdat2, lines(Year, Predicted, col=2, pch=16, type='o'))
      if(sum(recsd$SE[!is.na(recsd$SE)])!=0) suppressWarnings(with(tdat2, arrows(Year, Estimate+SE,y1=Estimate-SE, code=3,angle=90,length=0.05,col=2)))
      axis(1); axis(2)
    }
    caption <- "Observed (black) and estimated (red 9% CI grey) puerulus levels in each area of the model."
    addplot(filen=filename,rundir=rundir,category="Recruitment",caption=caption)
  }
  # Mean recruitment by area
  rec <- findNclean(c('Recruitment','by'), dat, 1,0)
  if(!'se'%in%colnames(rec)) rec %<>% mutate(se=0)
  rec2 <- rec %>% filter(rec$Year==min(rec$Year)) %>% mutate(cv=se/est, Nme=fleetarea$areaname[match(area,fleetarea$newarea)]) %>% mutate(Prop=est/sum(est), se2=cv*Prop, lwr=Prop-se2*SclErr, upr=Prop+se2*SclErr, lwr=ifelse(lwr<0,0,lwr), upr=ifelse(upr>1,1,upr))
  # rec2 <- rec %>% filter(rec$Year==min(rec$Year)) %>% mutate(lwr=est-se*SclErr,upr=est+se*SclErr, Nme=fleetarea$areaname[match(area,fleetarea$newarea)], lwr=ifelse(lwr<0,0,lwr))
  nareas <- length(unique(rec2$area))
  Ylim <- c(0, max(rec2$upr))
  filename <- filenametopath(rundir,paste0("Recruitmant_By_Area1.png"))
  plotprep(width=5,height=5,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(
    ggplot(data=rec2, aes(x=area,y=Prop))+
      geom_errorbar(aes(ymin=lwr,ymax=upr), width=0.25, colour='grey60')+
      geom_point()+
      theme(panel.background = element_rect(fill = "white",colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey20"),
            axis.text.x = element_text(vjust = 0.0, angle = 0)) +
      scale_x_discrete('Area',breaks=1:nareas,labels=rec2$Nme, limits=as.character(c(1:6)))+
      ylab('Relative mean recruitment')
  )
  caption <- "Relative mean recruitment by area."
  addplot(filen=filename,rundir=rundir,category="Recruitment",caption=caption)

  rec2 <- rec %>% mutate(lwr=est-se*SclErr,upr=est+se*SclErr, Area=fleetarea$areaname[match(area,fleetarea$newarea)], lwr=ifelse(lwr<0,0,lwr)) %>% filter(Year%in%Data$Year1:Data$Year2)
  filename <- filenametopath(rundir,paste0("Recruitmant_By_Area2.png"))
  plotprep(width=5,height=5,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(
    ggplot(rec2, aes(x=Year, y=est, colour=Area))+
      geom_line(linewidth=1) +
      theme(panel.background = element_rect(fill = "white",colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey20"),
            axis.text.x = element_text(vjust = 0.0, angle = 0)) +
      ylab('Recruitment (numbers)')
  )
  caption <- "Annual mean recruitment by area."
  addplot(filen=filename,rundir=rundir,category="Recruitment",caption=caption)

  filename <- filenametopath(rundir,paste0("Recruitmant_By_Area3.png"))
  plotprep(width=5,height=5,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(
    ggplot(rec2, aes(x=Year, y=est))+
      geom_errorbar(aes(ymin=lwr, ymax=upr), colour='grey60')+
      geom_line(linewidth=1) +
      theme(panel.background = element_rect(fill = "white",colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey20"),
            axis.text.x = element_text(vjust = 0.0, angle = 0)) +
      ylab('Recruitment (numbers)')+
      facet_wrap(~Area)
  )
  caption <- "Annual mean recruitment by area with 95% CI."
  addplot(filen=filename,rundir=rundir,category="Recruitment",caption=caption)

  recfrac <- findNclean(c('Recruitment','Fractions'), dat, 1,0)

  filename <- filenametopath(rundir,paste0("Recruitmant_Size_Dis.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(dim(recfrac)[2]-1))
  for(cl in 2:dim(recfrac)[2]) {
    suppressWarnings(plot(recfrac[,1], recfrac[,cl], axes=F, pch=16,xlab='Lower length bin (mm)', ylab='Proportion', type='o',cex=0.7, ylim=c(0,max(recfrac[,2:ncol(recfrac)])), main=))
    lines(recfrac[,1], recfrac[,cl], type='o',cex=0.7, col=(cl-1))
    axis(1);axis(2)}
  caption <- "Recruiting size composition."
  addplot(filen=filename,rundir=rundir,category="Recruitment",caption=caption)


  #### Movement ####
  print("Making Movement")
  moveP <- Pout[grepl('MovePars', Pout$name),]
  mov2 <- findNclean(c('#','Movement','parameters'), mov1, 1,0, char=T)
  mov3 <- findNclean(c('#','Movement','specifications'), mov1, 1,0, char=T)
  mov4 <- findNclean(c('#','Pattern','Type'), mov1, 1,0, char=F)
  mov3 <- mov3[mov3[,5]>0,]
  if(nrow(mov3)>0){
    nms <- rep(NA,nrow(mov2)) ; stcol <- which(mov2[1,]=='#')+1
    for(r in 1:nrow(mov2)){ nms[r] <- paste(mov2[r,stcol:ncol(mov2)],collapse=' ') }
    filename <- filenametopath(rundir,paste0("Movement.png"))
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(nrow(mov3)))
    nm <- unique(nms)
    for(i in 1:nrow(mov3)){
      mtmp <- mov3[i,]
      Age <- as.numeric(mtmp[1]) + 1
      fromArea <- as.numeric(mtmp[2]) + 1
      TStep <- as.numeric(mtmp[3]) + 1
      vec <- as.numeric(mtmp[4:ncol(mtmp)])
      toArea <- unique(mov4[as.numeric(mov4[,1])%in%vec,3]+1)
      vec2 <- as.numeric(moveP$value[vec])
      yrs <- startseason:endseason
      suppressWarnings(plot(yrs, vec2[1:length(yrs)], ylim=c(0,1), type='o', axes=F,pch=16,xlab='Years', ylab='Proportion moving', main=paste(fromArea,'to',toArea)))
      axis(1,cex.axis=1)
      axis(2)
    }
    caption <- "Estimated (red 95% CI grey) migration between areas of the model."
    addplot(filen=filename,rundir=rundir,category="Movement",caption=caption)
  }

  #### Size compositions ####
  ### Virgin Size Composition
  print("Making Size Compositions")
  ### Population size distribution - pooled
  tdat <- findNclean(c('Obs/Pred','Fleet'), dat, 1, convert=1)
  cnames <- colnames(tdat)
  cnames[substr(cnames,1,1)=='a'] <- paste('prop', 1:sum(substr(cnames,1,1)=='a'))
  colnames(tdat) <- cnames
  end <- 3+length(which(grepl('prop',colnames(tdat))))
  nms <-  colnames(tdat)[which(grepl('prop',colnames(tdat)))]
  colnames(tdat)[colnames(tdat)=='Obs/Pred'] <- 'O.P'
  ## Convert proportions to numbers
  for (i in 1:length(nms)){   tdat[, nms[i]==names(tdat)] <- tdat$Nsamp * tdat[,nms[i]==names(tdat)] }
  tdat1 <- tdat  %>% group_by(O.P,Fleet,Sex) %>% summarise_at(c('Nsamp',nms), sum) %>% filter(!is.na(Nsamp)) %>% as.data.frame()
  ## Convert back to proportions
  for (i in 1:length(nms)){   tdat1[, nms[i]==names(tdat1)] <- tdat1[,nms[i]==names(tdat1)]/tdat1$Nsamp }
  pos <- which(grepl('prop', colnames(tdat1)))
  tdat1$tot <- apply(as.matrix(tdat1[,pos]),1,sum)
  tdat1[,pos] <- tdat1[,pos]/tdat1$tot
  colnames(tdat1)[pos] <- paste('lb',lbinl[1:(length(lbinl)-1)])
  for(isex in unique(tdat1$Sex)){
    filename <- filenametopath(rundir,paste0(isex, "Fitted_Size_Comp.png"))
    plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
    #parset(plots=Fdims(length(unique(tdat1$Fleet))),margin = c(0.45, 0.45, 0.1, 0.05))
    parset(plots=c(1,1))
    #par(mfrow=c(Fdims(length(unique(tdat1$Fleet)))))

    tdat2 <- tdat1 %>% filter(Sex==isex) %>% pivot_longer(starts_with("lb"),names_to='albin',values_to = 'prop') %>% mutate(lbin=as.numeric(gsub('lb ','',albin))) %>%  mutate(fname = fleetarea$description[match(Fleet,fleetarea$fleet)]) %>% mutate(source=ifelse(O.P=='O','Observed','Predicted'))
    suppressMessages(print(ggplot(tdat2, aes(x=lbin, y=prop,colour=source)) +
            geom_line()+geom_point(size = 0.9) +
            ylab('Proportion') +
            facet_wrap(~fname) +
            scale_color_manual(values=c("black","red")) +
            theme(panel.background = element_rect(fill = "white",colour = NA),
                  panel.border = element_rect(fill = NA, colour = "grey20"),
                  axis.text.x = element_text(vjust = 0.0, angle = 45))+
            xlab('Length bin (mm)')
    ))
    caption <- paste('Sex = ',isex, "Catches by fleet, summed over years and time-steps (weighted by observations - Obs, black vs Exp, red).")
    addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)
  }

  tdat <- findNclean(c('Obs/Pred','Fleet'), dat, 1, convert=1)
  cnames <- colnames(tdat)
  cnames[substr(cnames,1,1)=='a'] <- paste('prop', 1:sum(substr(cnames,1,1)=='a'))
  colnames(tdat) <- cnames
  end <- 3+length(which(grepl('prop',colnames(tdat))))
  nms <-  colnames(tdat)[which(grepl('prop',colnames(tdat)))]
  colnames(tdat)[colnames(tdat)=='Obs/Pred'] <- 'O.P'
  ## Convert proportions to numbers
  for (i in 1:length(nms)){   tdat[, nms[i]==names(tdat)] <- tdat$Nsamp * tdat[,nms[i]==names(tdat)] }
  tdat1 <- tdat  %>% group_by(O.P,Fleet,Sex,Year) %>% summarise_at(c('Nsamp',nms), sum) %>% filter(!is.na(Nsamp)) %>% as.data.frame()
  ## Convert back to proportions
  for (i in 1:length(nms)){   tdat1[, nms[i]==names(tdat1)] <- tdat1[,nms[i]==names(tdat1)]/tdat1$Nsamp }
  pos <- which(grepl('prop', colnames(tdat1)))
  tdat1$tot <- apply(as.matrix(tdat1[,pos]),1,sum)
  tdat1[,pos] <- tdat1[,pos]/tdat1$tot
  colnames(tdat1)[pos] <- paste('lb',lbinl[1:(length(lbinl)-1)])
  yrs <- c(min(tdat1$Year),median(tdat1$Year), max(tdat1$Year))
  for(isex in unique(tdat1$Sex)){
    for(ifleet in unique(tdat1$Fleet)){
      tdat2 <- tdat1 %>% filter(Fleet==ifleet & Sex==isex) %>% pivot_longer(starts_with("lb"),names_to='albin',values_to = 'prop') %>% mutate(lbin=as.numeric(gsub('lb ','',albin))) %>%  mutate(fname = fleetarea$description[match(Fleet,fleetarea$fleet)]) %>% mutate(source=ifelse(O.P=='O','Observed','Predicted'))
      fname <- unique(tdat2$fname)
      filename <- filenametopath(rundir,paste0(isex," ",ifleet,"Fitted_Size_Comp2.png"))
      plotprep(width=9,height=9,filename=filename,cex=0.9,verbose=FALSE)
      parset(plots=c(1,1), margin = c(.5,.5,.5,.2))
      suppressMessages(print(ggplot(tdat2, aes(x=lbin, y=prop,colour=source)) +
              geom_line()+geom_point(size = 0.9) +
              ylab('Proportion') +
              facet_wrap(~Year) +
              scale_color_manual(values=c("black","red")) +
              theme(panel.background = element_rect(fill = "white",colour = NA),
                    panel.border = element_rect(fill = NA, colour = "grey20"),
                    axis.text.x = element_text(vjust = 0.0, angle = 45))+
              xlab('Length bin (mm)') + labs(color = fname)
      ))
      caption <- paste('Sex =',isex, "Fleet =",fname, "Size compositions by fleet, area and year, summed over time-steps (weighted by observations - Obs, black vs Exp, red).")
      addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)

      ## Line graphs
      filename <- filenametopath(rundir,paste0(isex," ",ifleet,"Median_Size_Comp2.png"))
      plotprep(width=9,height=9,filename=filename,cex=0.9,verbose=FALSE)
      parset(plots=c(1,1))
      ttmp <- tdat1 %>% filter(Fleet==ifleet & Sex==isex)
      tid <- ttmp %>% select(!starts_with('lb')) %>% mutate(O.P=ifelse(O.P=='O','Observed','Predicted')) %>% as.data.frame()
      tmat <- ttmp %>% select(starts_with('lb')) %>% as.data.frame()
      lbinM <- lbin + (diff(lbin[1:2])/2)
      funcoutMn <- function(x)  mean(rep(lbinM,x*1000))
      funcoutSd <- function(x)  sd(rep(lbinM,x*1000))
      tid$`Mean width (mm)` <- apply(tmat, 1, funcoutMn)
      tid$sd <- apply(tmat, 1, funcoutSd)
      nms <- unique(fleetarea$description[match(tid$Fleet, fleetarea$fleet)])
      suppressMessages(print(ggplot(tid,aes(x=Year,y=`Mean width (mm)`, colour=O.P), )+
              ggtitle(paste(isex,nms))+
              scale_color_manual(values=c(1,2))+
              geom_errorbar(data=tid[tid$O.P=='Observed',], aes(ymin=`Mean width (mm)`-sd, ymax=`Mean width (mm)`+sd), width=c(0.2), linewidth=0.9)+
              geom_errorbar(data=tid[tid$O.P=='Predicted',], aes(ymin=`Mean width (mm)`-sd, ymax=`Mean width (mm)`+sd), width=c(0.2), colour=2, linewidth=0.6)+
              geom_line(linewidth = 1)+geom_point(size=2)+
              theme(panel.background = element_rect(fill = "white",colour = NA),
                    panel.border = element_rect(fill = NA, colour = "grey20"),
                    axis.text.x = element_text(vjust = 0.0, angle = 0))))

      caption <- paste('Sex =',isex, "Fleet =",ifleet, "Area =", nms, "Median size compositions by fleet/area and year, summed over time-steps (weighted by observations - Obs, black vs Exp, red).")
      addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)
    }
  }

  ### Commercial Catch Size Composition
  tdat <- findNclean(c('Obs/Pred','Fleet'), dat, 1, convert=1)
  cnames <- colnames(tdat)
  cnames[substr(cnames,1,1)=='a'] <- paste('prop', 1:sum(substr(cnames,1,1)=='a'))
  colnames(tdat) <- cnames
  end <- 3+length(which(grepl('prop',colnames(tdat))))
  nms <-  colnames(tdat)[which(grepl('prop',colnames(tdat)))]
  tryCatch(tdat1 <- tdat %>% group_by(Fleet,Sex,Year,Step) %>% summarise_at(nms, diff) %>%
    mutate(yrstep = Year+(Step-1)/(max(Step))) %>% as.data.frame(),
    error = function(e) {
      stop(paste("Error: likely caused by multiple LF lines for the same fleet, step and sex. Original error:\n ", conditionMessage(e)))
    })
  pos <- which(grepl('prop', colnames(tdat1)))
  colnames(tdat1)[pos] <- paste('lb',lbinl[1:(length(lbinl)-1)])
  scale <- 0.05
  for(isex in unique(tdat1$Sex)){
    Sex <- case_when( isex %in% -1 ~  'U', isex %in% 0 ~ 'F',isex %in% 1 ~ 'M')
    filename <- filenametopath(rundir,paste0(isex, "Fitted_Size_Comp3.png"))
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(1+length(unique(tdat1$Fleet))))
    for(ifleet in unique(tdat1$Fleet)){
      tdat2 <- data.frame(tdat1[tdat1$Fleet==ifleet &tdat1$Sex==isex,])
      fname <- fleetarea$description[fleetarea$fleet==ifleet]
      if(nrow(tdat2)>0){
        tdat2 <- reshape2::melt(tdat2, variable.name = "Lbin",value.name = "diff",measure.vars=pos)
        if(!length(tdat2$Lbin)>0) {
          tdat2$Lbin <- tdat2$variable
          tdat2$diff <- tdat2$value  }
        tdat2$lb <- as.numeric(substr(tdat2$Lbin,4,7))
        if(nrow(tdat2)>0){

          tmp <- tdat2 %>% filter(abs(diff)>0) %>% mutate(mn=min(lb), mx=max(lb)) %>% group_by(mn, mx) %>% summarise(n=length(mn))
          suppressMessages(suppressWarnings(with(tdat2[tdat2$diff>0,], symbols(yrstep, lb, circles = diff, bg=rgb(1,0,0,0.2), fg=rgb(1,0,0,0.2), inches=scale, ylab='Length Bin', xlab='Fishing Season', bty='l', xlim=c(1970,graphrange[2]),ylim=c(tmp$mn,tmp$mx), main=paste('Sex',isex,', ',fname),las=1))))
          with(tdat2[tdat2$diff<0,], symbols(yrstep, lb, circles = -diff,bg=rgb(0,0,1,0.2), fg=rgb(0,0,1,0.2), inches=scale, add = T))
        }}
    }
    suppressMessages(suppressWarnings(symbols(rep(1,4),seq(1.8,0.6,-0.4),circles=c(1,0.75,0.5,0.25), bg=rgb(1,0,0,0.2), fg=rgb(1,0,0,0.2), inches=scale, ylim=c(0,3), ylab='', xlab='', bty='n', xlim=c(0.5,1.5), axes=F)))
    text(1,2,paste('Size Composition\nSex',isex), cex=0.8, pos=3)
    text(rep(1.1,4),seq(1.8,0.6,-0.4),c(1,0.75,0.5,0.25), pos=4, cex=0.8)
    text(1,0.1,'Proportional difference', cex=0.8)

    caption <- paste('Sex = ',isex, "Positive (red) and negative (blue) residuals between the observed and predicted proportions across lengths bins from each model area.")
    addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)

  }

  Ilen <- findNclean(c('Initial','N-matrix'), dat, 2,0)
  colnames(Ilen)[is.na(colnames(Ilen))] <- paste('L',1:length(colnames(Ilen)[is.na(colnames(Ilen))]),sep='')
  filename <- filenametopath(rundir,paste0("Virgin_Size_Comp.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(Ilen$Area))*length(unique(Ilen$Sex))))
  for(i in sort(unique(Ilen$Area))){
    for(s in sort(unique(Ilen$Sex))){
      tdat2 <- Ilen[Ilen$Area==i & Ilen$Sex==s,]
      lens <- tdat2[, grepl('Len',colnames(tdat2)) | grepl('a',substr(colnames(tdat2),1,1))]
      mx <- sqrt(max(lens/1e+6))
      Col <- ifelse(s==1,rgb(1,0,0,0.3),rgb(0,0,1,0.3,0.3))
      suppressMessages(suppressWarnings(plot(lbin, sqrt(lens[1,]/1e+6), type='l', cex=0.8, pch=16, axes=F, ylab='Size composition (sqrt-millions)', xlab='Length Bin', ylim=c(0, mx), lty=1, col=Col, main=paste('Area',i))))
      polygon(c(lbin,lbin[c(length(lbin),1)]),  c(as.numeric(sqrt(lens[1,]/1e+6)),0,0), col=Col,border=NA)
      abline(v=76, lty=1)
      for(sa in 2:length(unique(tdat2$Age))){
        polygon(c(lbin,lbin[c(length(lbin),1)]),  c(as.numeric(sqrt(lens[sa,]/1e+6)),0,0), col=Col,border=NA)
      }
      axis(1,lbin); axis(2)
    }}
  caption <- "The size composition at the end of the burn-In.  Each plot represents one area in the model and the various lines are the different sex and age groups. Sex 1 is red and Sex 2 blue."
  addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)

  ## Initial size composition
  Ilen <- findNclean(c('Full','N-matrix'), dat, 2,0)
  colnames(Ilen)[is.na(colnames(Ilen))] <- paste('L',1:length(colnames(Ilen)[is.na(colnames(Ilen))]),sep='')
  Ilen %<>% filter(Year==startseason, `Time-step`==max(`Time-step`))

  filename <- filenametopath(rundir,paste0(startseason, "_Size_Comp.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(Ilen$Area))*length(unique(Ilen$Sex))))
  for(i in sort(unique(Ilen$Area))){
    for(s in sort(unique(Ilen$Sex))){
      tdat2 <- Ilen[Ilen$Area==i & Ilen$Sex==s,]
      lens <- tdat2[, grepl('Len',colnames(tdat2)) | grepl('a',substr(colnames(tdat2),1,1))]
      mx <- sqrt(max(lens/1e+6))
      Col <- ifelse(s==1,rgb(1,0,0,0.3),rgb(0,0,1,0.3,0.3))
      suppressMessages(suppressWarnings(plot(lbin, sqrt(lens[1,]/1e+6), type='l', cex=0.8, pch=16, axes=F, ylab='Size composition (sqrt-millions)', xlab='Length Bin', ylim=c(0, mx), lty=1, col=Col, main=paste('Area',i))))
      polygon(c(lbin,lbin[c(length(lbin),1)]),  c(as.numeric(sqrt(lens[1,]/1e+6)),0,0), col=Col,border=NA)
      abline(v=76, lty=1)
      for(sa in 2:length(unique(tdat2$Age))){
        polygon(c(lbin,lbin[c(length(lbin),1)]),  c(as.numeric(sqrt(lens[sa,]/1e+6)),0,0), col=Col,border=NA)
      }
      axis(1,lbin); axis(2)
    }}
  caption <- "The size composition at the start of the model time-series (Year 1, final time-step).  Each plot represents one area in the model and the various modes are the different sex and age groups. Sex 1 is red and Sex 2 blue."
  addplot(filen=filename,rundir=rundir,category="FittedSizeComp",caption=caption)

  ### Tuning Length Composition Sample Size
  tdat <- findNclean(c('Obs/Pred','Fleet'), dat, 1, convert=1)
  cnames <- colnames(tdat)
  cnames[is.na(cnames)] <- paste('prop', 1:sum(is.na(cnames)))
  colnames(tdat) <- cnames
  end <- 3+length(which(grepl('prop',colnames(tdat))))
  nms <-  colnames(tdat)[which(grepl('prop',colnames(tdat)))]
  colnames(tdat)[colnames(tdat)=='Obs/Pred'] <- 'O.P'
  tdat <- tdat[tdat$O.P=="P",]
  Sexs <- sort(unique(tdat$Sex))
  Fleets <- sort(unique(tdat$Fleet))
  tune <- findNclean(c('Length','data','tuning'), dat, 1, convert=1) %>% mutate(Multiscale=NA)
  didtune <- findNclean(c('#','Weights','by','fleet'), ctl1, 3)
  tdat$ScaleNsamp <- didtune[match(tdat$Fleet,(didtune$Fleet+1)),5]
  tdat %<>% mutate(ScaleNsamp=ifelse(is.na(ScaleNsamp),1,ScaleNsamp))

  filename <- filenametopath(rundir,paste0("Tuning_Size_Comp.png"))
  plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(Fleets))*length(unique(Sexs))))
  Tuneout <- expand.grid(sex=Sexs,fleet=Fleets,scale=NA)
  for (isex in 1:length(Sexs)){
    for (ifleet in 1:length(Fleets)){
      tmp <- tdat[tdat$Sex==Sexs[isex] & tdat$Fleet==Fleets[ifleet],c("ScaleNsamp","EffN")]
      rng <- c(0,max(tmp,na.rm=T))
      plot(tmp$ScaleNsamp, tmp$EffN,pch=16,col=grey(0.4,0.3), ylim=c(0,max(tmp$EffN)), xlim=c(0,max(tmp$ScaleNsamp)), xlab='Scaled Sample Size', ylab='EffN',bty='l')
      lm1 <- lm(EffN~ScaleNsamp-1, data=tmp)
      abline(a=0,b=1,lty=1,col="red",lwd=2)
      abline(a=0,b=coef(lm1),lty=1,col="green",lwd=2)
      Tuneout$scale[Tuneout$sex==Sexs[isex] & Tuneout$fleet==Fleets[ifleet]] <- round(coef(lm1),5)
      tune$Multiscale[tune$Sex==isex & tune$Fleet==Fleets[ifleet]] <- round(coef(lm1),5)
      mtext(paste("Sx= ",Sexs[isex]," Ft= ",Fleets[ifleet], " Multi= ",round(coef(lm1),5)," Francis ",round(tune$Francis_Multiplier[tune$Sex==isex & tune$Fleet==Fleets[ifleet]],5),sep=""),3,line = -0.5, cex = 0.7)
    } }

  caption <- "By sex and fleet, summed over years and time-steps (red is 1:1 and green is data regression)."
  addplot(filen=filename,rundir=rundir,category="TuningSizeComp",caption=caption)

  filen <- "Tuning.csv"  # csv files only
  addtable(intable=tune,filen=filen,rundir=rundir,category="TuningSizeComp",caption="Tuning values.")

  #### Fit to Tagging data ####
  print("Making fit to Tagging data")
  tag <- findNclean(c('#Tagging','data'), dat, 2, convert = 2)
  if(exists("tag") && is.data.frame(tag) && nrow(tag) > 0){
    tag %<>% dplyr::mutate(RLArea = paste0('Area', RelArea),
                           RCArea = as.character(RecArea),
                           pearson = (Obs - Est) / sqrt(Est + 1e-5),
                           yt = year + tstep / (max(tstep) + 1))

    # Plot 1: Aggregated obs vs est by release area x recapture area
    tag_agg <- tag %>%
      group_by(RLArea, RCArea) %>%
      summarise(TotalObs = sum(Obs, na.rm = TRUE),
                TotalEst = sum(Est, na.rm = TRUE),
                .groups = "drop")

    mx <- ceiling(max(log(c(tag_agg$TotalObs, tag_agg$TotalEst) + 1e-5)))
    mn <- floor(min(log(c(tag_agg$TotalObs, tag_agg$TotalEst) + 1e-5)))

    filename <- filenametopath(rundir, paste0("Tag_recapt_by_release.png"))
 plotprep(width = 7, height = 7, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1, 1))
  suppressWarnings(print(
    ggplot(tag_agg, aes(x = log(TotalObs + 1e-5), y = log(TotalEst + 1e-5), color = RCArea)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      geom_point(size = 3) +
      facet_wrap(~RLArea) +
      coord_equal(xlim = c(mn, mx), ylim = c(mn, mx)) +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey20")) +
      xlab("log(Total Observed)") + ylab("log(Total Estimated)")
  ))
  caption <- "Aggregated observed vs estimated tag recaptures by release and recapture area."
  addplot(filen = filename, rundir = rundir, category = "Tag-Recapture", caption = caption)

  # Plot 2: Pearson residuals by release area, faceted by recapture area
  for (a in as.numeric(sort(unique(tag$RelArea)))) {
    tmp <- tag[tag$RelArea == a, ]
    filename <- filenametopath(rundir, paste0("Tag_recapt_by_release", a, ".png"))
    plotprep(width = 7, height = 7, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1, 1))
    suppressWarnings(print(
      ggplot(tag_agg, aes(x = log(TotalObs + 1e-5), y = log(TotalEst + 1e-5), color = RCArea)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
        geom_point(size = 3) +
        facet_wrap(~RLArea) +
        coord_equal(xlim = c(mn, mx), ylim = c(mn, mx)) +
        theme(panel.background = element_rect(fill = "white", colour = NA),
              panel.border = element_rect(fill = NA, colour = "grey20")) +
        xlab("log(Total Observed)") + ylab("log(Total Estimated)")
    ))
    caption <- "Aggregated observed vs estimated tag recaptures by release and recapture area."
    addplot(filen = filename, rundir = rundir, category = "Tag-Recapture", caption = caption)
    }
    # Plot 2: Pearson residuals by release area, faceted by recapture area
    for (a in as.numeric(sort(unique(tag$RelArea)))) {
      tmp <- tag[tag$RelArea == a, ]
      filename <- filenametopath(rundir, paste0("Tag_recapt_by_release", a, ".png"))
      plotprep(width = 7, height = 7, filename = filename, cex = 0.9, verbose = FALSE)
      parset(plots = c(1, 1))
      suppressWarnings(print(
        ggplot(tmp, aes(x = yt, y = pearson)) +
          geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
          geom_point(alpha = 0.5) +
          facet_wrap(~RCArea, scales = "free_y") +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20")) +
          xlab("Year") + ylab("Pearson Residual  (O-E)/sqrt(E)") +
          ggtitle(paste('Release Area', a))
      ))
      caption <- paste("Pearson residuals of tag recaptures for release area", a, ".")
      addplot(filen = filename, rundir = rundir, category = "Tag-Recapture", caption = caption)
    }
    tag <- findNclean(c('#Tagging','length'), dat, 2, convert = 2)
    tag %<>% pivot_longer(
      cols = c(ObsProp, EstProp),
      names_to = "Type",
      values_to = "Proportion"    )
    for(s in as.numeric(sort(unique(tag$Sex)))){
      for(a in as.numeric(sort(unique(tag$RelArea)))){
        filename <- filenametopath(rundir,paste0("Tag_recapt_by_release",s,a,".png"))
        plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
        parset(plots=c(1,1))
        tmp <- tag %>% filter(RelArea==a & Sex==s) %>% group_by(Type,RecArea) %>% mutate(Proportion=Proportion/(sum(Proportion)+1e-7)) %>% mutate(lbin=lbin[Size])
        suppressWarnings(print(ggplot(tmp, aes(x=lbin, y=Proportion, colour=Type))+
                                 geom_point()+geom_line()+
                                 facet_wrap(~RecArea, scale='free_y')+
                                 scale_color_discrete(
                                   name = "Type", labels = c("Observed", "Estimated")) +
                                 theme(panel.background = element_rect(fill = "white",colour = NA),
                                       panel.border = element_rect(fill = NA, colour = "grey20"))+
                                 xlab("Length Bin (mm)")+ylab("Proportion")+ggtitle(paste('Sex',s,'Release Area',a))))
        caption <- paste("Size composition of observed and estimated tag recaptures by sex",s,"and release area",a,".")
        addplot(filen=filename,rundir=rundir,category="Tag-Recapture",caption=caption)
      }
    }
  }
  #### Model Outputs ####
  ### Relative Legal Biomass by area
  print("Making Legal Biomass")
  lb <- findNclean(c('#Legal','Biomass', 'by'), dat, 2, convert = 2)
  if(ncol(lb)==3) lb$se <- 0
  lb %<>% dplyr::select(area,year=Year,est,se)   ## Uses the reference legal statement >76
  #write.csv(lb,'egg.csv')
  virgin <- data.frame(est=findNclean(c('#Virgin','Legal'), dat, 1)) %>% mutate(area=1:nrow(.))  ## Uses the reference legal statement >76
  lb %<>% mutate(virgin=virgin$est[match(area,virgin$area)], `Legal Biomass (t)`=est/1e+3, LBlwr=(est-se*SclErr)/1e+3, LBupr=(est+se*SclErr)/1e+3, rel=est/virgin, rellwr=(est-se*SclErr)/virgin, relupr=(est+se*SclErr)/virgin) %>% filter(est>0)
  lb %<>% filter(year>=GeneralSpecs$Year1) %>% group_by(area) %>% mutate(`B/B0`=rel, lwr=rellwr, upr=relupr, lwr=ifelse(lwr<0,0,lwr), upr=ifelse(upr>1,1,upr))
  areas <- readWorkbook(wb,sheet='area', startRow = 2)
  lb$areaname <- areas$Name[match(lb$area,areas$AreaCode)]
  reflev <- findNclean(c('#', 'Biomass', 'target'), ctl1, 1)
  colnames(reflev) <- c('target','threshold','limit')
  styr <- findNclean(c('#','First', 'year'), lbin1, 1)

  ## Relative biomass ####
  alllb <- lb %>% group_by(year) %>% summarise(est=sum(est), vir=sum(virgin), se=sqrt(sum(se^2))) %>% mutate(rel=est/vir, lwr=(est-se*SclErr)/vir, upr=(est+se*SclErr)/vir) %>% mutate(`B/B0`=rel, lwr=lwr, upr=upr)
  filename <- filenametopath(rundir,paste0("Relative_Legal_Biom.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(alllb, aes(x=year,y=`B/B0`))+
          geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2,colour='grey70')+
          geom_line( linewidth=1)+
          geom_point(size=1)+
          geom_hline(yintercept = reflev$target, colour='green', linewidth=1)+
          geom_hline(yintercept = reflev$threshold, colour='orange', linewidth=1)+
          geom_hline(yintercept = reflev$limit, colour='red', linewidth=1)+
          geom_hline(yintercept = 0, colour='white')+
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.5, angle = 45))+
          ylab('B/B0')+xlab('Year'))
  caption <- "Annual estimates of all Biomass relative to Virgin."
  addplot(filen=filename,rundir=rundir,category="Biomass",caption=caption)

  ## Relative Biomass by area ####
  filename <- filenametopath(rundir,paste0("Relative_Legal_Biom_Area.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(lb, aes(x=year,y=`B/B0`))+
          geom_errorbar(aes(ymin=lb$lwr, ymax=lb$upr), width=.2,position=position_dodge(0.05), colour='grey70')+
          geom_line()+
          geom_point(size=0.75)+
          geom_hline(yintercept = reflev$target, colour='green')+
          geom_hline(yintercept = reflev$threshold, colour='orange')+
          geom_hline(yintercept = reflev$limit, colour='red')+
          facet_wrap(~areaname)+
          ylim(0,1.05)+
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.5, angle = 45)))

  caption <- "Annual estimates of all Biomass in each model area relative to Virgin."
  addplot(filen=filename,rundir=rundir,category="Biomass",caption=caption)

  ### Legal Biomass (legal size) by Area
  filename <- filenametopath(rundir,paste0("Legal_Biom.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  print(ggplot(lb, aes(x=year,y=`Legal Biomass (t)`))+
          geom_errorbar(aes(ymin=LBlwr, ymax=LBupr), width=.2,position=position_dodge(0.05), colour='grey70')+
          geom_line()+
          geom_point(size=0.75)+
          facet_wrap(~areaname)+
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.5, angle = 45)))

  caption <- "Annual estimates (95% CI) of Legal Biomass (assignment of legal based on reference selectivity) in each model area."
  addplot(filen=filename,rundir=rundir,category="Biomass",caption=caption)

  print("Making Harvest Rates")
  lb <- findNclean(c('#Harvest','rate'), dat, 2)
  colnames(lb) <- c('Zone','Year','est', 'sd', 'est76')[1:ncol(lb)]
  if(exists('lb$sd')) lb %<>% mutate(up95=est+sd*SclErr,lw95=est-sd*SclErr)
  filename <- filenametopath(rundir,paste0(Sex, "HarvestRate_Biom.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(lb$Zone))))
  for(i in 1:length(unique(lb$Zone))){
    tlb <- lb[lb$Zone==i,]
    mxY <- max(tlb$est)
    ymn <- 0; ymx <- ifelse(exists('lb$sd'),max(lb$up95,na.rm=T)+0.1, max(lb$est,na.rm=T)+0.1)
    suppressWarnings(plot(tlb$Year, tlb$est, type='o', axes=F, pch=16, cex=0.9,ylab='Harvest Rate', xlab='Season', main=paste('Zone ',i), xlim=c(graphrange), ylim=c(0,ymx)))
    if(exists('lb$sd')) polygon(c(tlb$Year, rev(tlb$Year)), c(tlb$up95, rev(tlb$lw95)), col='grey90', border = F)
    points(tlb$Year, tlb$est, type='o',pch=16)
    axis(1)
    axis(2)
  }

  caption <- "Annual estimates (95% CI) of Harvest Rate (assignment of legal based on the management arrangements during each season) in each management Zone."
  addplot(filen=filename,rundir=rundir,category="Biomass",caption=caption)

  ### Fishing Mortality ####
  ## Compute F
  print("Making Fishing Mortality")
  catch <- findNclean('Catches', dat, 2)
  catch%<>% mutate(fdate = Year + Step/7)
  areas <- unique(catch$Area)
  filename <- filenametopath(rundir,paste0(Sex, "F_Mort.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(areas))))
  for(i in 1:length(areas)){
    tdat <- catch[catch$Area==i,]
    fleets <- unique(tdat$Fleet)
    ymx <- max(c(tdat$Fishing_mortality))
    ymn<- 0; ymx <- 1
    with(tdat[tdat$Fleet==fleets[1],], plot(fdate, Fishing_mortality, type='o', pch=16,cex=0.7,ylim=c(ymn,ymx), axes=F, ylab='Fishing Mort', xlab='Fishing Season', lty=1, main=paste("Area",i)))
    if(length(fleets)>1){
      for (f in 2:length(fleets)){
        with(tdat[tdat$Fleet==fleets[f],], lines(Year, Fishing_mortality, type='o',cex=0.7,col=f, lty=1))
      }
    }
    axis(1, seq(min(Zcatch$Year), max(Zcatch$Year),2))
    axis(2,las=1)
  }

  caption <- "Estimate Fishing mortality by area (the various fleets are shown in different colours)."
  addplot(filen=filename,rundir=rundir,category="Biomass",caption=caption)

  ### Total Mature Biomass ####
  print("Making Egg Production")
  ###Mature Biomass sex by area
  lb <- findNclean(c('#Mature','Biomass', 'sex'), dat, 2, convert = 2) %>% filter(Year>=GeneralSpecs$Year1)
  if(length(unique(lb$sex))>1){
    lb$estM <- lb$est/1e+3
    filename <- filenametopath(rundir,"Mature.Biomass.Sex_by_Area.png")
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(length(unique(lb$area))))
    for(i in 1:length(unique(lb$area))){
      tlb <- lb[lb$area==i,]
      mxY <- ceiling(max(tlb$estM))
      tlb <- lb[lb$sex==c('F','M')[1] & lb$area==i,]
      suppressWarnings(plot(tlb$Year, tlb$estM, type='o', axes=F, pch=16, cex=0.9,ylab='Biomass >76 (t)', xlab='Season', main=paste('Area ',i), ylim=c(0,mxY), col='red'))
      tlb <- lb[lb$sex==c('F','M')[2] & lb$area==i,]
      suppressWarnings(lines(tlb$Year, tlb$estM, type='o', col='blue', pch=16, cex=0.9))
      axis(1)
      axis(2)}

    caption <- "Annual estimates of Mature Biomass (1000s t) (assignment of mature  is age-based) in each model area by sex."
    addplot(filen=filename,rundir=rundir,category="Egg_Production",caption=caption)
  }

  egg <- findNclean(c('#Egg','Production','Total'), dat, 2) %>% filter(Year>=GeneralSpecs$Year1)
  if(ncol(egg)==2) egg$se <- 0
  egg %<>% mutate(lwr=est-se*SclErr, upr=est+se*SclErr)
  filename <- filenametopath(rundir,paste0(Sex, "Breeding_Biomass1.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(1))

  print(
    ggplot(egg,aes(x=Year,y=est))+
      geom_errorbar(aes(ymin=lwr,ymax=upr), colour='grey60')+
      geom_line(linewidth=1) +
      geom_hline(yintercept = max(egg$est)*reflev$target[1], colour='green')+
      geom_hline(yintercept = max(egg$est)*reflev$threshold[1], colour='orange')+
      geom_hline(yintercept = max(egg$est)*reflev$limit[1], colour='red')+
      theme(panel.background = element_rect(fill = "white",colour = NA),
            panel.border = element_rect(fill = NA, colour = "grey20"),
            axis.text.x = element_text(vjust = 0.5, angle = 45))+
      ylab('Egg Production')
  )

  caption <- "Estimated relative egg production (95% CI) of the whole fishery."
  addplot(filen=filename,rundir=rundir,category="Egg_Production",caption=caption)

  ## By BMSA
  egg <- findNclean(c('#Egg','Production','by'), dat, 2) %>% filter(Year>=GeneralSpecs$Year1)
  names(egg) <-c('year','area','est','se')[1:ncol(egg)]
  nareas <- length(unique(egg$area))
  if(nareas==8){
    egg2 <- egg %>% mutate(Loc=case_when(area == 2 ~ "South", area %in% c(4,6) ~ 'Central', area == 8 ~ "North", area == 5 ~'Abrolhos', .default = "other")) %>% filter(Loc!='other') %>% mutate(Loc=factor(Loc, levels=c("North",'Abrolhos','Central',"South")))
    if(length(egg2$se)==0) egg2$se <- 0
    egg2 %<>% group_by(year, Loc) %>% summarise(est=sum(est), se=sqrt(sum(se^2)))
    filename <- filenametopath(rundir,paste0(Sex, "BSMA.png"))
    plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=Fdims(length(unique(egg2$Loc))))
    for(i in 1:length(unique(egg2$Loc))){
      tdat <- egg2 %>% filter(Loc==sort(unique(egg2$Loc))[i])
      thres <- 1983:1985
      if(unique(tdat$Loc)=='North')  thres <- 1990:1993
      threslevel <- mean(tdat$est[tdat$year%in%thres])
      limlevel <- threslevel * 0.8
      ymn<- 0; ymx <- 1
      with(tdat, plot(year, est/max(est), type='o', pch=16,cex=0.7,ylim=c(ymn,ymx), axes=F, ylab='Relative Egg Production', xlab='Fishing Season', lty=1, main=unique(tdat$Loc)))
      if(sum(tdat$se)>0) arrows(tdat$year, tdat$est-tdat$se*SclErr, y1=tdat$est-tdat$se*SclErr, code=3, angle=90, length=0.05)
      lines(range(tdat$year), rep(threslevel/max(tdat$est),2), lwd=1.5, col='orange')
      lines(range(tdat$year), rep(limlevel/max(tdat$est),2), lwd=1.5, col='red')
      axis(1)
      axis(2,las=1)
    }

    caption <- "Estimated Relative Egg Production by Breeding Stock Management Area."
    addplot(filen=filename,rundir=rundir,category="Egg_Production",caption=caption)
  }

  ## By Area
  egg2 <- egg
  filename <- filenametopath(rundir,paste0(Sex, "Breeding_Biomass2.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=Fdims(length(unique(egg2$area))))
  for(i in 1:length(unique(egg2$area))){
    tdat <- egg2 %>% filter(area==sort(unique(egg2$area))[i])
    ymn<- 0; ymx <- 1
    with(tdat, plot(year, est/max(est), type='o', pch=16,cex=0.7,ylim=c(ymn,ymx), axes=F, ylab='Relative Egg Production', xlab='Fishing Season', lty=1, main=paste('Area: ',unique(tdat$area))))
    axis(1)
    axis(2,las=1)
  }

  caption <- "Estimated Relative Egg Production by model area."
  addplot(filen=filename,rundir=rundir,category="Egg_Production",caption=caption)

  ####  Egg Production ####
  egg <- findNclean(c('Egg','Production'), echo, 3)
  egg <- egg[1:(which(is.na(egg[,20]))[1]-1),]
  names(egg) <- c('age','area','season', paste0('lb',1:(ncol(egg)-3)))

  # Pivot to long format
  egg <- egg %>% tidyr::pivot_longer(cols = starts_with("lb"), names_to = "lbin", values_to = "value") %>% mutate(lbin = as.numeric(gsub("lb", "", lbin)))
  # Create a signature for each year's curve shape (per area/age)
  curve_sig <- egg %>% group_by(area, age, season) %>%
    summarise(sig = paste(round(value, 6), collapse = "_"), .groups = "drop")
  # Find distinct curves and their year ranges
  curve_groups <- curve_sig %>%
    group_by(area, age, sig) %>%
    summarise(yr1 = min(season), yr2 = max(season), .groups = "drop") %>%
    mutate(label = paste0("Age ", age, " (", yr1, "-", yr2, ")"))

  # Join labels back to long data
  df_plot <- egg %>% left_join(curve_sig, by = c("area", "age", "season")) %>%
    left_join(curve_groups %>% select(area, age, sig, label), by = c("area", "age", "sig")) %>% mutate(value=value/1e+6)

  for(a in unique(df_plot$area)) {
    filename <- filenametopath(rundir,paste("Area",a,"Egg Production.png"))
    plotprep(width=10,height=10,filename=filename,cex=0.9,verbose=FALSE)
    parset(plots=c(1,1))
    p <- df_plot %>%
      filter(area == a) %>%
      ggplot(aes(x = lbin, y = value, colour = label, group = label)) +
      geom_line() +
      labs(x = "Length bin", y = "Egg production (Millions)",
           colour = NULL, title = paste("Area", a)) +
      theme_bw() +
      theme(legend.position = "inside",
            legend.position.inside = c(0.7, 0.3),
            legend.text = element_text(size = 8),
            legend.background = element_rect(fill = alpha("white", 0.7)))
    suppressWarnings(print(p))
    caption <- paste("Egg production curves by area", a,". This is a combination of maturity, multiple spawning and fecundity and used to estimate egg production by area.")
    addplot(filen=filename,rundir=rundir,category="Selectivity_Retention",caption=caption)
  }

  ### Natural Mortality ####
  print("Making Natural Mortality")
  nmort <- findNclean(c('Natural','Mortality'), dat, 2)
  names(nmort) <-c('area','age','year','est')
  nmort %<>% mutate(est=ifelse(est<0,NA,est))
  areas <- unique(nmort$area)
  ages <- unique(nmort$age)
  yLim <- c(0,max(nmort$est,na.rm=T)*2)
  filename <- filenametopath(rundir,paste0("M_Mort.png"))
  plotprep(width=7,height=7,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  nmort %<>% filter(!is.na(est)) %>% mutate(Nme=fleetarea$description[match(area,fleetarea$newarea)])
  print(ggplot(nmort, aes(x=year, y=est, group=age))+
          geom_line()+
          ylim(yLim)+
          labs(y= "Estimate", x='Year')+
          facet_wrap(~Nme)+
          theme(panel.background = element_rect(fill = "white",colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(vjust = 0.5, angle = 45)))
  caption <- "Estimate Natural mortality by area and year (Density dependent mortality)."
  addplot(filen=filename,rundir=rundir,category="Natural_Mortality",caption=caption)

  ## Estimated parameters ####
  print("Making Parameter Diagnostics")
  pars <- findNclean(c('#','Parameter','Par'), dat, 1, char = T)
  nms <- dat[find(c('#','Parameter','Par'), dat, 0),1:15];
  nms <- nms[nms!='#' & nms!='']
  colnames(pars) <- nms[1:ncol(pars)]
  pars %<>% filter(!is.na(Estpar_cnt) | (as.numeric(Link) > 0)) %>%
    mutate(Estimate=round(as.numeric(Estimate),3)) %>%
    select(Parameter,Estimate,SD,Gradient,lwrBound,uprBound,PriorType,PriorMean,PriorSD,Initial,Link)


  ## plot parameters
  ## Parameter distribution plots
  ## Parameter distribution plots
  plot_df <- pars %>%
    mutate(across(c(Estimate, SD, Gradient, lwrBound, uprBound, PriorType, PriorMean, PriorSD, Initial), as.numeric)) %>%
    filter(!is.na(SD), SD > 0, as.numeric(Link) == 0)

  npars_per_page <- 8
  npar_total <- nrow(plot_df)
  if (npar_total > 0) {
    npages <- ceiling(npar_total / npars_per_page)

    for (ipage in 1:npages) {
      idx <- ((ipage - 1) * npars_per_page + 1):min(ipage * npars_per_page, npar_total)
      sub_df <- plot_df[idx, ]

      curve_list <- list()
      for (i in 1:nrow(sub_df)) {
        row <- sub_df[i, ]
        xmin <- row$Estimate - 4 * row$SD
        xmax <- row$Estimate + 4 * row$SD
        if (row$PriorType > 0 && !is.na(row$PriorMean)) {
          xmin <- min(xmin, row$PriorMean - 4 * row$PriorSD)
          xmax <- max(xmax, row$PriorMean + 4 * row$PriorSD)
        }
        xseq <- seq(xmin, xmax, length.out = 200)

        mle_dens <- dnorm(xseq, row$Estimate, row$SD)

        prior_dens <- rep(0, length(xseq))
        if (row$PriorType == 1)
          prior_dens <- dnorm(xseq, row$PriorMean, row$PriorSD)
        if (row$PriorType == 2) {
          shape <- (row$PriorMean / row$PriorSD)^2
          rate <- row$PriorMean / row$PriorSD^2
          prior_dens <- dgamma(pmax(xseq, 0), shape, rate)
        }
        if (row$PriorType == 3) {
          mulog <- log(row$PriorMean) - 0.5 * log(1 + (row$PriorSD / row$PriorMean)^2)
          sdlog <- sqrt(log(1 + (row$PriorSD / row$PriorMean)^2))
          prior_dens <- dlnorm(pmax(xseq, 1e-10), mulog, sdlog)
        }

        curve_list[[i]] <- data.frame(
          Parameter = row$Parameter,
          x = rep(xseq, 2),
          y = c(mle_dens, prior_dens),
          Type = rep(c("max. likelihood", "prior"), each = length(xseq))
        )
      }
      curve_df <- do.call(rbind, curve_list)
      curve_df$Parameter <- factor(curve_df$Parameter, levels = sub_df$Parameter)
      sub_df$Parameter <- factor(sub_df$Parameter, levels = sub_df$Parameter)

      p <- ggplot() +
        geom_line(data = curve_df, aes(x = x, y = y, colour = Type, linewidth = Type)) +
        scale_colour_manual(values = c("max. likelihood" = "blue", "prior" = "black")) +
        scale_linewidth_manual(values = c("max. likelihood" = 0.5, "prior" = 1.2)) +
        geom_vline(data = sub_df, aes(xintercept = Estimate), colour = "blue", linewidth = 0.4) +
        geom_point(data = sub_df, aes(x = Initial, y = 0), colour = "red", shape = 17, size = 3) +
        geom_vline(data = sub_df, aes(xintercept = lwrBound), colour = "orange", linewidth = 0.6) +
        geom_vline(data = sub_df, aes(xintercept = uprBound), colour = "orange", linewidth = 0.6) +
        facet_wrap(~ Parameter, scales = "free", ncol = 2) +
        labs(x = "Parameter value", y = "Density") +
        theme_bw() +
        theme(legend.position = "top", strip.text = element_text(size = 9),
              legend.title = element_blank())

      filename <- filenametopath(rundir, paste0("parameter_distributions_page", ipage, ".png"))
      plotprep(width = 8, height = 10, filename = filename, cex = 0.9, verbose = FALSE)
      parset(plots = c(1, 1))
      suppressWarnings(print(p))
      caption <- paste("Parameter distributions page", ipage, "- prior (black), MLE (blue), initial (red), bounds (orange).")
      addplot(filen = filename, rundir = rundir, category = "Parameters", caption = caption)
    }
  }
  filen <- "Est.Params.csv"
  addtable(intable=pars,filen=filen,rundir=rundir,category="Parameter Table",
           caption="Estimated final parameters and gradients. Link column indicates parameter linking (0 = directly estimated or fixed).")

  ## Correlation matrix heatmap ####
  corfile <- filenametopath(rundir, "CorrelationMatrix.csv")
  if (file.exists(corfile)) {
    cormat <- as.matrix(read.csv(corfile, row.names = 1))
    par_order <- rownames(cormat)
    par_nums <- as.numeric(gsub("\\D+", "", par_order))
    par_order <- par_order[order(par_nums)]
    # Melt for ggplot
    cor_df <- expand.grid(Par1 = rownames(cormat), Par2 = colnames(cormat),
                          stringsAsFactors = FALSE)
    cor_df$r <- as.vector(cormat)

    # Apply numeric ordering
    cor_df$Par1 <- factor(cor_df$Par1, levels = par_order)
    cor_df$Par2 <- factor(cor_df$Par2, levels = par_order)

    # Only plot upper triangle
    cor_df <- cor_df[as.numeric(cor_df$Par1) < as.numeric(cor_df$Par2), ]

    p <- ggplot(cor_df, aes(x = Par1, y = Par2, fill = r)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-1, 1), name = "r") +
      geom_text(data = cor_df[abs(cor_df$r) > 0.85, ],
                aes(label = round(r, 2)), size = 2.5) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6)) +
      labs(x = "", y = "")

    filename <- filenametopath(rundir, "parameter_correlations.png")
    plotprep(width = 10, height = 9, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1, 1))
    suppressWarnings(print(p))
    caption <- "Parameter correlation matrix. Values shown where |r| > 0.85."
    addplot(filen = filename, rundir = rundir, category = "Parameter Table", caption = caption)
  }

  txt5 <- "Built by Simon de Lestang, Andre Punt  and  Klaas Hartmann. Relies on packages developed by Malcolm Haddon."

  runnotes <- matrix(c(txt2,txt2.1,txt3,txt4,txt5), nrow=5)

  endtime <- as.character(Sys.time())

  reportlist <- list(  #these 2 are minimal requirements for the replist
    # though the whole of replist can be NULL if
    starttime=starttime,  # you are feeling lazy.
    endtime=endtime
  )

  ## Copy all files across that made the model
  start <- list.files('..', pattern = 'DAT')
  for(f in 1: length(start))  file.copy(paste0("../",start[f]), paste0(rundir,'/',start[f]))
  ## Copy important output files
  addfiles <- c('Output.RL','SDReport.RL','model final.par')
  for(f in 1: length(addfiles))  file.copy(addfiles[f], paste0(rundir,'/',addfiles[f]))


  make_html(replist=reportlist,rundir=rundir,width=800,openfile=openfile,
            runnotes=runnotes,verbose=FALSE,packagename="makehtml",
            htmlname="IMuLT")

}
