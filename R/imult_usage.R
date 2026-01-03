#' Display IMuLT Usage Instructions
#'
#' Prints step-by-step instructions for running the IMuLT model
#' @export
imult_usage <- function() {
  cat("
=================================================================
IMuLT - Integrated Model using Length Transition
Usage Instructions
=================================================================

FIRST TIME SETUP:
-----------------
1. Set working directory to your lobster model folder:
   setwd('C:/Users/YourName/Lobster Model/2026')

2. Copy the model structure template:
   copy_model_structure()

BASIC WORKFLOW:
---------------
1. Update data in 'ModelStructure.xlsx'

2. Build model input files from Excel:
   BuildInputFiles()

3. Choose your model run:
   choose_model()

4. Update parameters and length composition weightings (if needed):
   UpdatePars('No')
   UpdateLFWeights('No')

5. Set output options and run model:
   FullOutput <- FALSE
   LoadData()
   Data$DoProject <- 0
   LoadPars()
   SolveModelNew(500, 2000, report = TRUE)
   MakeDiagReport(is95 = FALSE)

HELPER FUNCTIONS:
-----------------
- get_model_structure_path() : Get path to ModelStructure.xlsx
- copy_model_structure()     : Copy template to current directory
- imult_usage()              : Display these instructions

=================================================================
\n")
}
