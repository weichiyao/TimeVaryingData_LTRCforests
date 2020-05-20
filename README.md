# TimeVaryingData_LTRCforests
Analysis codes for "Ensemble Methods for Right-Censored Survival-Time Data with Time-Varying Covariates"
- The **pkg** folder contains the R package [LTRCforests].
- The **analysis** folder provides analysis code in the paper:
  - The subfolder **data** contains functions to create simulated dataset.
  - The subfolder **util** contains the source functions used to compute the L_2 errors. 
  - comparison.R: code to evaluate performance comparison for the four methods, IC Cox, IC ctree, IC cforest with default parameter settings and IC cforest with _mtry_ tuned and _minsplit_, _minbucket_, _minprob_ set by "15%-Default-6% Rule".
  - properties_of_iccf_mtry.R: code to evaluate performance of IC cforest with different _mtry_'s values and _mtry_ tuned by Out-of-Bag procedure.
  - properties_of_iccf_params.R: code to evaluate performance of IC cforest with different _minsplit_, _minbucket_, _minprob_ settings.
