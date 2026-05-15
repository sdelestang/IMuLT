# IMuLT - Integrated Model using Length Transition

A package to load, control and run the stock assessment model IMuLT for Crustaceans.

## Installation

IMuLT requires several packages, including some only available from GitHub:

Install directly from GitHub using devtools:
```r
# Install pak if you don't have it
install.packages("pak")

# Install a CRAN dependancy
install.packages("TMB")

# Install GitHub dependancies
pak::pak("https://github.com/haddonm/makehtml")
pak::pak("https://github.com/haddonm/codeutils")
pak::pak("https://github.com/haddonm/hplot")

# Install IMuLT (This can take a while (~1 minute) as the model must be compiled)
pak::pak("sdelestang/IMuLT")
```

## Required Dependencies

The following packages will be automatically installed:
- TMB
- readxl
- reshape2
- svDialogs
- stats4
- dplyr
- magrittr
- tidyr
- ggplot2
- openxlsx
- tcltk

## First Time Setup

After installing, copy the required "ModelStructure.xlsx"" template to your working directory:
```r
library(IMuLT)
copy_model_structure()  # Copies to current directory
```

Alternatively, get the file path directly:
```r
get_model_structure_path()
```

## Usage
```r
library(IMuLT)

# Set working directory to your model folder
setwd("C:/Users/YourName/Crust_Model/2026")

# Update data in ModelStructure.xlsx, then build input files
BuildInputFiles()

# Choose your model run
choose_model()

# Run the model
LoadData()
LoadPars()

FitModel(500, 3000, report = TRUE)
MakeDiagReport(is95 = TRUE)

## Re-run model after examining diagnotics ##

# Update length composition weightings 
UpdateLFWeights("Yes")

# If needed update parameters to estimates from previous run
UpdatePars("Yes") # or "No"

# Change phases if needed
AdjustPhase()

# After running multiple models including diagnostics to compare model runs.
compare_legal_biomass()

# To Produce a word document of a certain run  
WriteWord()
```

For quick reference anytime:
```r
imult_usage()  # Displays step-by-step instructions
```

## Authors

Simon de Lestang (DPIRD, Western Australia)
Andre Punt (University of Washington & CSIRO)
Klaas Hartmann (University of Tasmania)
