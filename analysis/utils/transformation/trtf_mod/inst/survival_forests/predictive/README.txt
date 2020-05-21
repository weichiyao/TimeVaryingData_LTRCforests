This is an instruction on running experiments on simulated data with predictive effect.

To generate the data run tr_data_simulation.R. As a result, in the same directory, where
tr_data_simulation.R is located, a new file called simulated_data.rda will be created.
This file will contain the generated data. Approximately 115 Mb of disk space is required.

All predictive methods with the corresponding predict functions and
modified weighting scheme are realized in tr_competitors.R.

To learn all predictive models on the generated data and estimate their log-likelihood
run tr_empeval.R. Runnig this file a folder results_tf will be created, where for each
generated learning-validation pair of datasets a transformation
survival forest model with general splitting score will be saved in the corresponding
.rda file. The log-likelihood estimates will be saved to tr_results_empeval.rda in the
same directory where tr_empeval.R is located. Approximately 54 Kb of disk space is required.
