## Update all input files
# Used mainly to change phases and to update parameters with estimates
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#source('maps/maps.R')
source('Modelareas.R')
library(SimeFunc)
library(openxlsx)

## Open up file with all info
wb <- loadWorkbook(file="ModelStructure.xlsx")

#This is the location of the data input files and their associated parameters
dynamics <- readWorkbook(wb,sheet='dynamics')
startseason <- as.numeric(dynamics$value[dynamics$object=='startseason'])
endseason <- as.numeric(dynamics$value[dynamics$object=='endseason'])
projectseason <- as.numeric(dynamics$value[dynamics$object=='projectedseason'])
projectcatch <- as.numeric(dynamics$value[dynamics$object=='projectedcatch'])
burnin <- as.numeric(dynamics$value[dynamics$object=='burnin'])
ages <- as.numeric(dynamics$value[dynamics$object=='ages'])
sexs <- 0:as.numeric(dynamics$value[dynamics$object=='sexs'])
areas <- readWorkbook(wb,sheet='area', startRow = 2)
times <- readWorkbook(wb,sheet='times', startRow = 2)
fleets <- readWorkbook(wb,sheet='fleetcode')
effic <- readWorkbook(wb,sheet='EfficiencyCreep', startRow = 2)
migrate <- readWorkbook(wb,sheet='migrate', startRow = 2)
(zones <- length(unique(areas$newzone)))
area <- areas %>% group_by(newzone) %>% reframe(newarea=unique(newarea)) 
(zoneareas <- split(area$newarea, area$newzone))
lens <- seq(dynamics$value[dynamics$object=='lblwr'],dynamics$value[dynamics$object=='lbupr'],dynamics$value[dynamics$object=='lbgap'])+1
#gauge <- readWorkbook(wb,sheet='legalchanges')
hgrad <- read.csv('high grading for new model.csv') %>% filter(!is.na(year))

## Create a new folder for the model if one does not exist
(files <- list.files(path="../", pattern = 'AgeRun'))
(nfile <- paste(length(unique(area$newarea)),'Area',ages,'AgeRun',substr(startseason,3,4),"_",substr(endseason,3,4),sep=''))
if(!(nfile%in%files)) {
  dir.create(paste('../',nfile,sep='')) ;  dir.create(paste('../',nfile,"/Output",sep=''))
}
(fls <- nfile)
(floc <- paste(dirname(getwd()),fls,sep='/'))  # location of data files

# Make all files from Scratch - only works if you have access to raw data
getwd()
f <- 2
files <- list.files(pattern ='Make')
for(f in 1:length(files)){
  print(paste("Building: ", files[f]))
  source(files[f])  
  }

# Update phases with those defined in "Phases / parameters.xlsx" file - 
#  "Phases / parameters.xlsx" is to be stored in the same folder as this file
source("UpdatePhase.R")          #- Set the phases

#Update parameters with those stored in model1.par which is a file produced when the model is run 
source("UpdateParameters.R")     #- Get previous parameter estimates from model1.par

