# LTRCforests for Time-Varying Covariate Data 

In the paper "**Ensemble Methods for Survival Data with Time-Varying Covariates**", we generalize the conditional inference and relative risk forests to allow time-varying covariates and propose two forest algorithms `CIF-TV` and `RSF-TV`. The proposed methods by design can handle survival data with *all* combinations of left-truncation and right-censoring in the survival outcome, and with both time-invariant and time-varying covariates. For this matter, we call the methodology as `LTRC CIF` and `LTRC RRF` for general speaking. 

The [pkg](./pkg/) folder contains the R package [LTRCforests](https://cran.r-project.org/web/packages/LTRCforests/LTRCforests.pdf), available on CRAN. 

We here provide analysis codes for "Ensemble Methods for Survival Data with Time-Varying Covariates", as well as the analysis codes for time-invariant covariate data. 

## Analysis Codes 
We provide analysis codes for "Ensemble Methods for Survival Data with Time-Varying Covariates", as well as the analysis codes for time-invariant covariate data.

The [analysis](./analysis/) folder provides analysis codes in the paper:
  - The subfolder [data](./analysis/data/) contains 
    - functions to create simulated dataset with time-varying covariates
    - functions to create simulated LTRC dataset with time-invariant covariates
    - a function to obtain the real dataset.
  - The subfolder [codes](./analysis/codes/) contains the functions to reproduce the analysis:
    - [simulations_tvary.R](./analysis/codes/simulations_tvary.R) -- codes to reproduce results for simulated datasets with *time-varying* covariates
    - [simulations_tindepLTRC.R](./analysis/codes/simulations_tindepLTRC.R) -- codes to reproduce results for simulated LTRC datasets with *time-invariant* covariates
    - [plot_and_tables_tvary.R](./analysis/codes/plot_and_tables_tvary.R) -- codes to reproduce plots and tables for simulated datasets with time-varying covariates
    - [plot_and_tables_tindepLTRC.R](./analysis/codes/plot_and_tables_tindepLTRC.R) -- codes to reproduce plots and tables for simulated LTRC datasets with time-invariant covariates
    - [realsetPBC.R](./analysis/codes/realsetPBC.R) -- analysis of real dataset (including functions to reproduce plots)
    
    In particular, [simulations_tvary.R](./analysis/codes/simulations_tvary.R) and [simulations_tindepLTRC.R](./analysis/codes/simulations_tindepLTRC.R) provide results
    - to compare performance of LTRC forests with default parameter settings and proposed parameter settings.
    - to evaluate performance comparison for the four methods, the `Cox` model, `CIF`, `RRF` and `TSF` (all forests trained with proposed parameter settings)
    - to choose methods by using IBS-based 10-fold CV, and compare the results produced by the selection rule with the best method. 
  - The subfolder [utils](./analysis/utils/) contains the source functions used to perform the analysis in the folder [codes](./analysis/codes/), including the functions to compute the integrated L2 difference. 


## LTRCforests for Time-invariant covariate data
There are certainly many situations in which only time-invariant (baseline) covariate information is available, and understanding the properties of different methods in that situation is important. In fact, our developed methodology and algorithms allow for estimation using the proposed forests for (left-truncated) right-censored data with time-invariant covariates. While "Ensemble Methods for Survival Data with Time-Varying Covariates" mainly focuses on the analysis of the methodology applied on time-varying covariate data, we here also provide a similar analysis for applying the `LTRC CIF` and `LTRC RRF` on left-truncated right-censored time-invariant covariate data. 
  
The same data-driven guidance for tuning the parameters or selecting a modeling method also applies to the time-invariant covariates case (for both left-truncated right-censored survival data and right-censored survival data), which implies its broad effectiveness regardless of additional left-truncation and regardless of the presence of time-varying effects. 

See below the boxplots of integrated L2 difference for performance comparison. Datasets are generated with time-invariant covariates, left-truncated right-censored survival times following a Weibull-Increasing distribution. The first row shows results for the number of subjects `N=100`, second row for `N=300`, third row for `N=500`, bottom row for `N=1000`; the first column shows results for linear survival relationship, second column for nonlinear, the third column for interaction. In each plot, `LTRC CIF(P)`--LTRC CIF with proposed parameter settings; `LTRC RRF(P)`--LTRC RRF with proposed parameter settings; `LTRC TSF(P)`--LTRC TSF with proposed parameter settings;  Opt--Best method; `IBSCV`--Method chosen by IBS-based 10-fold CV. 

<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/470eddcbf6b3ea3f48132da0bebf73513c6c10bc/analysis/figures/boxplots_RC.pdf" width="350" title="Boxplots of integrated L2 difference for performance comparison on time-invariant covariate data">
