# Slam900_Proteomic_Clock
Functions that together create a cross-sectional or longitudinal aging clock.

**Overview**
To create an aging clock using cross-sectional data, linear regression modeling is used via the glmnet package for feature selection and the caret package (using 'lm' as the method) for training. 
To create an aging clock using longitudinal data, fixed-effect modeling is used via the lme4 package for both feature selection and training. 

## make_clock
This function is designed to easily create a full cloick. It thus makes assumptions as to many of the configuration details. These details are ones which I kept consistent in MY analysis.
| Parameter | Type       | Description|
|-----------|------------|------------|
| train_dt  | data.table | datatable to train the clock on. |
| test_dt   | data.table | dt to test the clock on. If left as NULL, function assumes the train dt IS the test dt. |
| response  | string     | str naming the response variable (the 'y' in our regressions). As with all similar parameters, this should correspond to a col in the train_dt. |
| predictors| list       | list of col names containing predictors |
| covariates| list       | list of col names holding covariable variables (these along with the predictors are our 'x' in our regressions). |
|id | string containing name of column to use for ID. Only use if running mixed effects modeling. |
| prediction_method | string | method used to train prediction model. Choose linear_regression or mixed_effect. |
| title     | string     | name to use for output plots |
| max_selected| numeric  | num of features we want to select. Setting this to Inf will result in a rank deficiency error for most clocks. |
| rep_seeds   | list     | list of seeds to use to create each clock. The number of seeds will be the number of clocks produced. |
| scale       | boolean  | whether or not to scale your continuous covariates (such as age). |

**Return**
List:
| Name | Description | 
|------|-------------|
| "clocks" | list of all clocks. |
| "metrics_plot" | plot of clock metrics. |
| "selected_features" | list of "selected_samples": all predictors selected for the clocks; "selected_info": info for those predictors |
