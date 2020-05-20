# TimeVaryingData_LTRCforests
Analysis codes for "Ensemble Methods for Right-Censored Survival Data with Time-Varying Covariates"
- The **pkg** folder contains the R package [LTRCforests].
- The **analysis** folder provides analysis code in the paper:
  - The subfolder **data** contains functions to create simulated dataset.
  - The subfolder **utils** contains the source functions used to compute the integrated L_2 difference. 
  - comparison.R: code to evaluate performance comparison for the four methods, IC Cox, IC ctree, IC cforest with default parameter settings and IC cforest with _mtry_ tuned and _minsplit_, _minbucket_, _minprob_ set by "15%-Default-6% Rule".
  - properties_mtry.R: code to evaluate performance of LTRCforests with different _mtry_'s values and _mtry_ tuned by Out-of-Bag procedure.
  - properties_default_vs_proposed.R: code to compare performance of LTRC forests with default parameter settings and proposed parameter settings.
  - The subfolder **supplemental** contains the source functions used to reproduce the same analysis for right-censored survival data with time-invariant covariates.
