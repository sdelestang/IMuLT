# IMuLT - Integrated Model using Length Transition

A package to load, control and run the stock assessment model IMuLT for Crustaceans.

## Installation

Install directly from GitHub using devtools:
```r
# Install devtools if you don't have it
install.packages("devtools")

# Install IMuLT
devtools::install_github("sdelestang/IMuLT")
```

## Required Dependencies

The following packages will be automatically installed:
- TMB
- readxl
- reshape2
- svDialogs
- stats4
- makehtml
- hplot
- dplyr
- magrittr
- tidyr
- ggplot2
- openxlsx

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

# Set working directory to your lobster model folder
setwd("C:/Users/YourName/Lobster Model/2026")

# Update data in ModelStructure.xlsx, then build input files
BuildInputFiles()

# Choose your model run
choose_model()

# Update parameters and length composition weightings if needed
UpdatePars("No")
UpdateLFWeights("No")

# Run the model
FullOutput <- FALSE
LoadData()
Data$DoProject <- 0
LoadPars()
SolveModelNew(500, 2000, report = TRUE)
MakeDiagReport(is95 = FALSE)
```

For quick reference anytime:
```r
imult_usage()  # Displays step-by-step instructions
```


## Author

Simon de Lestang (DPIRD, Western Australia)
