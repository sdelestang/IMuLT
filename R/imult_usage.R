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

4. Set output options and run model:
   LoadData()
   LoadPars()
   SolveModelNew(500, 2000, report = TRUE, nRestarts = TRUE, newtonSteps=5)
   MakeDiagReport(is95 = TRUE)

5. Update length composition weightings
   UpdateLFWeights('Yes')

6. If needed update parameters to estimates from previous run
   UpdatePars('Yes') # or 'No'

7. Change phases if needed
   AdjustPhase()

8. After running multiple models including diagnostics compare model runs.
   compare_legal_biomass()

HELPER FUNCTIONS:
-----------------
- get_model_structure_path() : Get path to ModelStructure.xlsx
- copy_model_structure()     : Copy template to current directory
- imult_usage()              : Display these instructions

=================================================================
\n")
}
