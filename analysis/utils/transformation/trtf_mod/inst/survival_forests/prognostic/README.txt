This is an instruction on running experiments on simulated data with prognostic effect.

To generate the data run data_simulation.R. As a result, in the same directory, where
data_simulation.R is located, a new file called simulated_data.rda will be created.
This file will contain the generated data. Approximately 122 Mb of disk space is required.

All forest-based prognostic methods with the corresponding predict functions and
modified weighting scheme are realized in competitors.R.

To learn all prognostic models on the generated data and estimate their log-likelihood
run empeval.R. Runnig this file a folder results_tf will be created, where for each
generated learning-validation pair of datasets a transformation
survival forest model with general splitting score will be saved in the corresponding
.rda file. The log-likelihood estimates will be saved to results_empeval.rda in the
same directory where empeval.R is located. Approximately 80 Kb of disk space is required.
