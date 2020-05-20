# TimeVaryingData_LTRCforests
Analysis codes for "Ensemble Methods for Right-Censored Survival Data with Time-Varying Covariates"
1.The **pkg** folder contains the R package [LTRCforests].
2.The **analysis** folder provides analysis code in the paper:
  - The subfolder **data** contains functions to create simulated dataset with time-varying covariates, and the function to obtain the real dataset.
  - The subfolder **utils** contains the source functions used to compute the integrated L_2 difference. 
  - The subfolder **analysis** contains the functions to reproduce the analysis:
    - properties_mtry.R: code to evaluate performance of LTRCforests with different _mtry_'s values and _mtry_ tuned by Out-of-Bag procedure.
    - properties_default_vs_proposed.R: code to compare performance of LTRC forests with default parameter settings and proposed parameter settings.
    - comparison.R: code to evaluate performance comparison for the four methods, the Cox model, LTRCCF, LTRCRSF and TSF (all forests trained with proposed parameter settings). 
    - IBSCV.R: code to choose methods by using IBS-based 10-fold CV, and compare the results produced by the selection rule with the best method. 
  - The subfolder **supplemental** contains the source functions used to reproduce the same analysis for right-censored survival data with time-invariant covariates.
    - The subfolder **data** contains functions to create simulated dataset with time-invariant covariates.
    - The subfolder **utils** contains the source functions used to compute the integrated L_2 difference. 
    - The subfolder **analysis** contains the functions to reproduce the analysis.
