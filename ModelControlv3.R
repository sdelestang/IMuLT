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

# Compile and load model
compile("Model1.cpp")#, flags="-Wno-ignored-attributes")
dyn.load(dynlib("Model1"))

# Choose Model run
fls <- csemod(list.files(pattern = 'Run')) 
setwd(paste(getwd(), "\\",fls,sep=""))

## Update starting parameters to last run 
UpdatePars()

# New boolean to write more stuff out
FullOutput <- F

## Load all data
LoadData()

# Set the projection flag (0=n0;1=Yes)
Data$DoProject = 0;

## Load all of the parameters
LoadPars()

## Solve the model - adjust the max iterations for all phases (1st number), last phase (2nd number)
SolveModel(500,3000)

## Solve model again using previous parameters but only short phase to limit parameters and speed up report eval
UpdatePars()
SolveModel(500,3000)

## Run Reports (including s.d. which wont work unless a big run has been done and convergence reached) 
RunReport()
