# TimeVaryingData_LTRCforests
Analysis codes for "Ensemble Methods for Right-Censored Survival Data with Time-Varying Covariates"

1.The **pkg** folder contains the R package [LTRCforests].

2.The **analysis** folder provides analysis code in the paper:
  - The subfolder **data** contains 
    - functions to create simulated dataset with time-varying covariates
    - functions to create simulated dataset with time-invariant covariates
    - a function to obtain the real dataset.
  - The subfolder **codes** contains the functions to reproduce the analysis:
    - simulations_tvary.R -- codes to reproduce results for simulated datasets with time-varying covariates
    - simulations_tfixed.R -- codes to reproduce results for simulated datasets with time-invariant covariates
    - plot_and_tables_tvary.R -- codes to reproduce plots and tables for simulated datasets with time-varying covariates
    - plot_and_tables_tfixed.R -- codes to reproduce plots and tables for simulated datasets with time-invariant covariates
    - realsetPBC.R -- analysis of real dataset (including functions to reproduce plots)
    
    In particular, simulations_tvary.R and simulations_tvary.R provide results
    - to compare performance of LTRC forests with default parameter settings and proposed parameter settings.
    - to evaluate performance comparison for the four methods, the Cox model, LTRCCF, LTRCRSF and TSF (all forests trained with proposed parameter settings)
    - to choose methods by using IBS-based 10-fold CV, and compare the results produced by the selection rule with the best method. 
  - The subfolder **utils** contains the source functions used to perform the analysis in the folder **codes**, including the functions to compute the integrated L_2 difference. 
  
