# LTRCforests for Time-Varying Covariate Data 

In the paper "**Ensemble Methods for Survival Function Estimation with Time-Varying Covariates**", we generalize the conditional inference and relative risk forests to allow time-varying covariates and propose two forest algorithms `CIF-TV` and `RSF-TV`. The proposed methods by design can handle survival data with *all* combinations of left-truncation and right-censoring in the survival outcome, and with both time-invariant and time-varying covariates. For this matter, we call the methodology as `LTRC CIF` and `LTRC RRF` for general speaking. 

The [pkg](./pkg/) folder contains the R package [LTRCforests](https://cran.r-project.org/web/packages/LTRCforests/LTRCforests.pdf), available on CRAN. 

We here provide analysis codes for the paper, as well as the analysis codes for time-invariant covariate data in [analysis](./analysis/) folder. Please also find [Supplemental Material](./analysis/doc/supplemental_material.pdf) in the subfolder [doc](./analysis/doc/).

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


## LTRCforests for Time-Varying Covariate Data
Main analysis for applying the methodology on time-varying covariate data have been provided in the paper "Ensemble Methods for Survival Data with Time-Varying Covariates". Here we provide some more detailed information as supplemental material.

### "Out-of-bag" observation-based `mtry` tuning algorithm
The values of `mtry` can be fine-tuned using the "out-of-bag" observation. The simulation results have shown that it can greatly improve the forest performance over the default setting. See the following figure for the performance comparisons using `CIF-TV` for different values of `mtry` vs. the optimal one (`Opt`) vs. the one tuned by the tuning algorithm (`Tuned`).

<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/6a08a7f03d55caca8f9446d55ddee244c1e138f7/analysis/figures/LTRC_time-varying/mtryCIF_WI_c1.png?raw=true" width="900" />  

See [Figure](./analysis/figures/LTRC_time-varying/mtryRRF_WI_c1.png) for the similar results using `RRF-TV` and [Figure](./analysis/figures/LTRC_time-varying/mtryTSF_WI_c1.png) for `TSF-TV`.

### `ntree` in the ensembles
Throughout the experiments, we use `ntree=100L` for all forest ensembles. It has been recommended that a random forest should have a number of trees between 64 and 128 trees (see [Lecture Notes in Computer Science](https://www.researchgate.net/publication/230766603_How_Many_Trees_in_a_Random_Forest)). It is true that generally more trees will result in better accuracy. However, more trees also means higher computational cost, and after a certain number of trees, the improvement is negligible. See the following figure for performance comparisons for different numbers of trees built in the forest methods. 
 
<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/5db207d2d92e9364c78028dd1e847f994d70dd9f/analysis/figures/LTRC_time-varying/ntree_WI_PH_20var_c1.png?raw=true" width="800" />


### Bootstrap pseudo-subject observations vs. subject observations
In forest-like algorithms, bootstrapped samples are typically used to construct each individual tree to increase independence between these base learners. For time-varying covariate data, we have considered two different ways to bootstrap the observations:
 - *Bootstrapping pseudo-subjects*. Namely, it is to bootstrap "independent" observations as the first step of any forest algorithm; this is because all pseudo-subjects are treated as independent observations in the recursive partitioning process; 
 - *Bootstrapping subjects*. It keeps all of the pseudo-subjects for each subject in the bootstrap sample. 

Simulations have shown that the two different bootstrapping mechanisms do not result in fundamentally different levels of performance:
<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/0d1944df06fb41940b984a8d04072a1fd4f69710/analysis/figures/LTRC_time-varying/Bootstrap_pseudosubject_vs_subject.png?raw=true" width="1000" />

## LTRCforests for Time-Invariant Covariate Data
 "Ensemble Methods for Survival Data with Time-Varying Covariates" mainly focuses on the analysis of the methodology applied on time-varying covariate data. There are certainly many situations in which only time-invariant (baseline) covariate information is available, and understanding the properties of different methods in that situation is important. In fact, our developed methodology and algorithms allow for estimation using the proposed forests for (left-truncated) right-censored data with time-invariant covariates. 
  
In fact, the same data-driven guidance for tuning the parameters or selecting a modeling method also applies to the time-invariant covariates case (for both left-truncated right-censored survival data and right-censored survival data). 

### How `mtry` affects the performance and how `mtry`-tuning algorithm performs
The following figures show how `LTRC CIF` performs with different values of `mtry` under the PH setting and non-PH setting, respectively. The datasets are generated with survival times following a Weibull-Increasing distribution, light (right-)censoring rate. This implies its broad effectiveness regardless of additional left-truncation and regardless of the presence of time-varying effects.

<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/86819cc4b1cdc5656ab9ecc2357890728c299366/analysis/figures/LTRC_time-invariant/mtryCIF_PH_20var_WI_c1_LTRC.png?raw=true" width="800" />

**Figure 2.1**. Integrated L2 difference of LTRC CIF with different mtry values distribution under the PH setting.

<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/86819cc4b1cdc5656ab9ecc2357890728c299366/analysis/figures/LTRC_time-invariant/mtryCIF_nonPH_20var_WI_c1_LTRC.png?raw=true" width="800" />

**Figure 2.2**. Integrated L2 difference of LTRC CIF with different mtry values distribution under the non-PH setting.

 
### Using IBS-CV to choose among different methods
See below the boxplots of integrated L2 difference for performance comparison. Datasets are generated with time-invariant covariates, left-truncated right-censored survival times following a Weibull-Increasing distribution. The first row shows results for the number of subjects `N=100`, second row for `N=300`, third row for `N=500`, bottom row for `N=1000`; the first column shows results for linear survival relationship, second column for nonlinear, the third column for interaction. In each plot, `LTRC CIF(P)`--LTRC CIF with proposed parameter settings; `LTRC RRF(P)`--LTRC RRF with proposed parameter settings; `LTRC TSF(P)`--LTRC TSF with proposed parameter settings;  Opt--Best method; `IBSCV`--Method chosen by IBS-based 10-fold CV. 

<img src="https://github.com/ElainaYao/TimeVaryingData_LTRCforests/blob/86819cc4b1cdc5656ab9ecc2357890728c299366/analysis/figures/LTRC_time-invariant/boxplots_LTRC.png?raw=true" width="1200" />

**Figure 3.3**. Boxplots of integrated L2 difference for performance comparison on time-invariant covariate data.
