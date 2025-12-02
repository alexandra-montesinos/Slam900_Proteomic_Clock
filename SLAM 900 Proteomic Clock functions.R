library("glmnet") #v4.1.8
library("dplyr") #v1.1.4
library("lme4") #v1.1.36
library("lmerTest") #v3.1.3
library("Metrics") #v0.1.4
library("ggplot2") #v3.5.1
library("corrplot") #v0.95

source("../../General Use Code/misc_funcs.R")
source("../../General Use Code/clean_data.R")


###The functions in this file together create a proteomic clock. It is recommended you use the make_clock() function to do so, as it
# will run everything else for you--though it has slightly less customization.
# It is also recommended you refer to 'SLAM data prep.R' to see how your dataset should be prepared.

#I should note that these functions are designed to work with things that *aren't* SLAM, but the name of the file just kinda stuck
# and it's too late to go back now. Or at least I don't want to bother.


### FUNCTION TO REMOVE SAMPLES NOT PRESENT FOR ALL VALS OF A CATEGORICAL VARIABLE
# @param datatable   == datatable containing samples and categorical variables.
# @param vars_list   == list of strings naming each column holding categorical variables.
# @param ignore_cols == list of strings naming each column holding non-samples which we want to ignore.
#                       Necessary so we don't remove rows where one of the ignored columns has a value that
#                       isn't present for all values, which isn't important.
#
# @return            == datable passed in, but filtered
remove_incomplete_samples <- function(datatable, vars_list, ignore_cols, threshold = 0) {
  #Save columns we want to ignore to put them back in later
  cols_ignored <- datatable[, ..ignore_cols]
  #Remove columns we want to ignore
  datatable <- datatable[, !colnames(datatable) %in% ignore_cols, with = FALSE]
  
  #Removes all rows which contain NAs for our important categorical variables
  datatable <- datatable[complete.cases(datatable[, vars_list, with = FALSE]), ]
  
  for(var in vars_list) {
    for(val in unique(datatable[[var]])) {
      datatable <- datatable[, .SD, .SDcols = colSums(!is.na(datatable[get(var) == val])) > threshold]
    }
  }
  
  #Append ignored columns back to front of datatable
  datatable <- cbind(cols_ignored, datatable, fill = TRUE)
  #Get rid of that annoying fill column they throw in there
  datatable <- datatable %>% dplyr::select(-fill)
  
  return(datatable)
}

### CHECKS R^2 LEVEL
# @param actual   == list of actual values 
# @param precited == list of predicted values
#
# @return = numeric, resquared value. Closer to 1 is better.
find_rSquared <- function(actual, predicted) {
  return(1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2))
}

### SPLITS GIVEN DATATABLE INTO TRAINING AND TEST SETS
# @param datatable         == datatable
# @param training_percent == percentage of data to use in training set
# @param splitOn_col       == name of the grouping column. Only include if you want to get an even split of GROUPS rather than ALL SAMPLES.
# @param testSet     == the test set to be returned. Only include one if you DON'T want this function to modify the test set.
#
# @return list == list, trainingSet = training set, testSet = test set
train_test_split <- function(dt, training_percent, splitOn_col = NULL, testSet = NULL) {
  #If we have a column to split on it means the rows are grouped somehow, so we don't want to split completely randomly.
  if(!is.null(splitOn_col)) {
    groups <- unique(dt[[splitOn_col]]) #Grab unique groups
    train_groups <- sample(groups, size = training_percent * length(groups)) #Sample portion of GROUPS rather than total samples
    
    trainingSet <- dt[get(splitOn_col) %in% train_groups]
    
    if(is.null(testSet)) {
      testSet  <- dt[!get(splitOn_col) %in% train_groups]
    }
  } else {
    train_indices <- sample(nrow(dt), size = training_percent * nrow(dt))
    trainingSet <- dt[train_indices]
    
    if(is.null(testSet)) {
      testSet <- dt[-train_indices]
    }
  }
  
  return(list("trainingSet" = trainingSet, "testSet" = testSet))
}

### FUNCTION TO SET UP DATATABLE VARIABLES FOR REGRESSION
# Automatically finds and sets the types for common regression variables.
# @param dt
# @param scale == Whether or not to scale continuous variables
set_regression_vars <- function(dt, scale = FALSE) {
  var_table <- data.table(new_col = rep(NA, nrow(dt)))
  
  for(col in names(dt)) {
    col_l <- tolower(col) #Makes checking names more flexible
    #Checking for binary/multinomial variables
    if(col_l == "sex" || col_l == "strain" || col_l ==  "cohort" || col_l == "euthanized" || col_l == "visit" ||
       grepl('(^|[ _])id($|[ _])', col_l, ignore.case = TRUE)) {
      dt[[col]] <- as.factor(dt[[col]]) 
      var_table[[col]] <- dt[[col]]
      #Checking for continuous variables
    } else if(col_l == "age" || col_l == "weight" || 
              grepl('(^|[_])percent($|[_])', col_l, ignore.case = TRUE) || grepl('%', col_l, ignore.case = TRUE)) {
      if(scale) {
        print(dt[[col]])
        print("Scaling")
        dt[[col]] <- as.numeric(scale(as.numeric(dt[[col]])))
        print(dt[[col]])
      } else {
        dt[[col]] <- as.numeric(dt[[col]])
      }
      var_table[[col]] <- dt[[col]]
    }
  }
  
  var_table <- var_table[, -1]
  return(list("dt" = dt, "table" = var_table))
}

### FUNCTION TO RUN CV.GLMNET()
# @param response         == data column containing response variable data.
# @param predictor_matrix == matrix of all covariates & predictors.
#
# @param covariates       == list of strings naming each column holding covariate variables or data columns containing them.
# @param alpha       == alpha to use in glmnet. 0 uses a Ridge model, 1 uses a LASSO model, any number in between uses Elastic Net.
#                       LASSO shrinks some coefs to exactly 0, effectively performing feature selection. When used with data with high colinearity,
#                         it tends to drop features arbitrarily, so may end up with innacurate results.
#                       Ridge shrinks them *towards* zero. Avoids large coeficients, reducing model complexity & avoiding overfitting.
#                       Elastic Net combines the methods. It should perform better than LASSO, but will be more generous in its selection.
# @param nfolds      == number of folds for glmnet to use.
# @param lambda_type == type of lambda to choose for output coefficients datatable.
#                       min is considered more predictive but has potential for overfitting.
#                       1se is the standard choice as it's 'best'.
# @param max_selected == max number of features to select. default: number of predictors passed in + 1, which is the default used for glmnet.
#
# @return            == list of 1: datatable of coefficients; 2: cv.glmnet output object;
#                       3 & 4: predictor matrix & response table, mostly used for further predictions without having to rewrite code.
run_cvGlmnet <- function(response, predictor_matrix, covariates,
                         alpha = 1, nfolds = 10, lambda_type = "lambda.1se", max_selected = ncol(predictor_matrix) + 1) {
  #Ensure our max is not higher than the glmnet default
  if(max_selected > ncol(predictor_matrix) + 1) {
    max_selected = ncol(predictor_matrix) + 1
  }
  print("Running run_cvGlmnet()...")
  ### RUN CV.GLMNET
  print("Running cv.glmnet()...")
  glm <- cv.glmnet(predictor_matrix, response, alpha = alpha, nfolds = nfolds, dfmax = max_selected)
  
  coefs <- coef(glm, s = lambda_type)
  coefs_dt <- data.table(Feature = rownames(coefs), Coefficient = as.numeric(coefs))
  coefs_dt <- subset(coefs_dt, !(Feature %in% covariates)) #Remove covariate predictors as they are not relevant
  coefs_dt <- coefs_dt[coefs_dt$Coefficient != 0, ] #Remove 0-val coefficients
  coefs_dt <- coefs_dt[Feature != "(Intercept)"] #Remove that pesky intercept

  return(list("coefs" = coefs_dt, "glm" = glm))
}

### MAKE_PROTEOMICCLOCK
# Prepares a datatable for use in the proteomic clock functions
# @param dt          == datatable containing response variable, predictors, covariates, and ID columns (latter 2 optional). Other columns will be removed
# @param response          == string naming the response variable
# @param predictors        == list of predictor names
# @param covariates        == list of strings naming each column holding covariate variables
# @param id                == string containing name of column to use for ID OR col itself. Only use if running mixed effects modeling (lme).
# @param scale             == whether or not to scale continuous variables
#
# @return            == prepared dt
prep_clock_data <- function(dt, response, predictors, covariates = NULL, id = NULL, scale = FALSE) {
  #Grab just predictors in dt
  predictors_dt <- dt[, names(dt) %in% predictors, with = FALSE]
  #Scale them!
  predictors_dt_scaled <- as.data.table(scale(predictors_dt))
  #Get rid of predictors with too many NAs
  predictors_dt_imputed <- predictors_dt_scaled[, .SD, .SDcols = which(colMeans(is.na(predictors_dt_scaled)) < 0.1)]
  predictors_dt_imputed <- impute_all_cols(predictors_dt_imputed)
  #Bind together only relevant cols
  dt_complete <- dt[, names(dt) %in% c(response, id, covariates), with = FALSE]
  dt_complete <- cbind(dt_complete, predictors_dt_imputed)
  
  #Remove incomplete samples
  dt_filtered <- remove_incomplete_samples(dt_complete, c(response, covariates), ignore_cols = id, threshold = 0)
  
  #Prep var types
  dt_filtered <- set_regression_vars(dt_filtered, scale = scale)$dt
  
  return(dt_filtered)
}

### SELECT FEATURES FOR CLOCK
# @param dt                == datatable to perform feature selection on. Must contain all response, predictor, and covariate columns.
# @param response          == string naming the response variable (the 'y' in our regressions) 
# @param predictors        == list of relevant predictors
# @param covariates        == list of strings naming each column holding covariable variables (these along with the predictors are our 'x' in our regressions)
# @param prediction_method == method used to train prediction model.
# @param id                == string containing name of column to use for ID. Only use if running mixed effects modeling.
# @param alpha             == alpha to use in glmnet. 0 uses a Ridge model, 1 uses a LASSO model, any number in between uses Elastic Net.
#                               LASSO shrinks some coefs to exactly 0, effectively performing feature selection. When used with data with high colinearity,
#                                it tends to drop features arbitrarily, so may end up with innacurate results.
#                               Ridge shrinks them *towards* zero. Avoids large coeficients, reducing model complexity & avoiding overfitting.
#                               Elastic Net combines the methods. It should perform better than LASSO, but will be more generous in its selection.
# @param nfolds      == number of folds for glmnet to use. (optional; default: 10)
# @param lambda_type == type of lambda to choose for output coefficients datatable.
#                       "lambda.min" is considered more predictive but has potential for overfitting.
#                       "lambda.1se" is the standard choice as it's 'best'. (optional; default: "lambda.1se")
# @param max_selected == Max number of features to select.
# @param seed         == seed to use for randomness. (optional; default: 123)
#
# @return            == list of 1, "selected_dt": dt containing all selected predictors
#                               2, "predictors_info_dt": dt containing all selected predictors and info about them, for reference
feature_select <- function(dt, response, predictors, covariates, prediction_method, id = NULL, alpha = 1, nfolds = 10, 
                           lambda_type = "lambda.1se", max_selected = Inf, seed = 123) {
  set.seed(seed)
  
  if(prediction_method == "glmnet" || prediction_method == "linear_regression") {
    print("Feature selection via glmnet()")
    
    ### PREP PREDICTORS MATRIX
    response_list <- dt[[response]]
    predictor_matrix <- as.matrix(dt[, names(dt) %in% predictors, with = FALSE])
    
    glm <- run_cvGlmnet(response_list, predictor_matrix, covariates, alpha = alpha, nfolds = nfolds, 
                        lambda_type = lambda_type, max_selected = max_selected)

    ### FINDING SELECTED PREDICTORS
    print(paste0("Selected this many: ", length(glm$coefs$Feature)))
    selected_predictors <- glm$coefs$Feature
    #Remove covariates as they are not predictors. Don't train on those!
    selected_predictors <- selected_predictors[!(selected_predictors %in% covariates)]
    #Subsets predictors datatable to only those selected
    selected_dt <- dt[, ..selected_predictors]
    #Append to the response column to create our very good datatable what for splitting
    selected_dt <- bind_cols(dt[[response]], selected_dt)
    colnames(selected_dt)[1] <- response #Retain the name of our response variable column
    
    selected_predictors_dt <- glm
    
  } else if (prediction_method == "mixed_effect") {
    print("Feature selection via lmer()")
    
    #Re-split our training set from its predictors
    training_vars       <- dt %>% dplyr::select(all_of(c(response, id, covariates)))
    training_predictors <- dt[, names(dt) %in% predictors, with = FALSE]
    
    selected_predictors_dt <- data.table(Predictor = character(), Beta = numeric(), FDR = numeric())
    #Feature selection! We run this per protein
    for(i in names(training_predictors)) {
      trainingSet_single <- suppressMessages(bind_cols(training_vars, training_predictors[[i]], .name_repair = "minimal"))
      names(trainingSet_single)[length(trainingSet_single)] <- i
      
      form <- as.formula(paste(response, "~ . -", id, "+ (1 |", id, ")"))
      # trained_model <- suppressMessages(lmer(form, data = trainingSet_single))
      trained_model <- tryCatch(
        {
          suppressMessages(
            lmer(form, data = trainingSet_single)
          )
        },
        error = function(e) {
          message("Skipping ", i, " due to error: ", e$message)
          return(NULL)
        }
      )
      if(is.null(trained_model)) {
        next
      }

      fdr <- p.adjust(summary(trained_model)$coefficients[i, "Pr(>|t|)"], method = "fdr") #Adjust pvals to fdr vals
      #If fdr is insignificant or error, throw out this predictor
      if(fdr >= 0.05) {
        next
      }
      beta <- fixef(trained_model)[[i]]
      # View(fixef(trained_model))
      # print(coef(summary(trained_model)))#[1, "Estimate"]
      # beta <- summary(trained_model)[[10]][5]
      
      #Bind to datatable
      selected_predictors_dt <- rbind(selected_predictors_dt, as.data.frame(list(Predictor = i, Beta = beta, FDR = fdr)))
    }
    
    #Select only top X proteins
    if(max_selected != Inf) {
      selected_predictors_dt <- selected_predictors_dt %>% arrange(Beta) %>% slice_head(n = max_selected)
    }
    
    # top_x <- dt[order(pvalue)][1:x] 
    print(paste0("Selected this many: ", nrow(selected_predictors_dt)))

    #Subset & rebind our selected predictors with our response & ID columns (remember that we don't want to use covariates in our actual model)
    selected_predictors <- training_predictors[, names(training_predictors) %in% selected_predictors_dt$Predictor, with = FALSE]
    selected_dt <- bind_cols(
      setNames(list(dt[[response]]), response),
      setNames(list(dt[[id]]), id),
      selected_predictors
    )
  }
    
  return(list("selected_dt" = selected_dt, "predictors_info_dt" = selected_predictors_dt))
}

### TRAIN A CLOCK
# @param train_dt          == datatable to train clock on.
# @param response          == string naming the response variable (the 'y' in our regressions) 
#                               Elastic Net combines the methods. It should perform better than LASSO, but will be more generous in its selection.
# @param prediction_method == method used to train prediction model.
# @param id                == string containing name of column to use for ID OR col itself. Only use if running mixed effects modeling (lme).
# @param nfolds      == number of folds for glmnet to use. (optional; default: 10)
# @param lambda_type == type of lambda to choose for output coefficients datatable.
#                       "lambda.min" is considered more predictive but has potential for overfitting.
#                       "lambda.1se" is the standard choice as it's 'best'. (optional; default: "lambda.1se")
# @param resampling_method == resampling method used in training clock. (optional; default: cv)
# @param seed       == seed to use for randomness. (optional; default: 123)
#
# @return            == trained clock model, ready for use in predictions
train_clock <- function(trainingSet, response, prediction_method, id = NULL, 
                        nfolds = 10, lambda_type = "lambda.1se", resampling_method = "cv", seed = 123) {
  set.seed(seed)
  
  ### RUN CARET::TRAIN
  if(prediction_method == "glmnet") {
    print("Training clock with caret::train() using glmnet...")
    
    trained_model <- caret::train(
      form = as.formula(paste(response, "~ .")), data = trainingSet,
      method = "glmnet",
      trControl = trainControl(method = resampling_method, number = nfolds),
      tuneGrid = expand.grid(alpha = c(0.5, 1), lambda = glm$glm[[lambda_type]]) #Tune alpha
    )
  } else if(prediction_method == "linear_regression") {
    print("Training clock with caret::train() using lm...")
    
    trained_model <- caret::train(
      form = as.formula(paste(response, "~ .")), data = trainingSet, #Same x & y as feature selection, but without covariates
      method = "lm",
      trControl = trainControl(method = resampling_method, number = nfolds),
      metric = "Rsquared"
    )
  } else if(prediction_method == "mixed_effect") {
    print("Training clock with lme4::lmer()...")
    
    #Does this line actually need to exist? We're pretty sure no -- trainingSet isn't modified in this function
    # testSet <- testSet[, intersect(names(testSet), names(trainingSet)), with = FALSE]
    
    #Run it on our selected predictors datatable
    trained_model <- lmer(as.formula(paste(response, "~ . -", id, "+ (1 |", id, ")")), data = trainingSet)
  }

  return(trained_model)
}

### PREDICT DATASET USING CLOCK
# @param clock == trained model to use for prediction, for example the one generated by train_clock()
# @param testSet == datatable to predict using the clock
# @param response == name of response variable
clock_predict <- function(clock, testSet, response) {
  print("Predicting...")
  #Grab just the name of the response col, if it's not that already
  if(type(response) != "character") {
    response <- response
  }
  
  predictions <- predict(clock, newdata = testSet, allow.new.levels = TRUE)
  
  predictions_table <- data.table(
    Actual    = testSet[[response]],
    Predicted = predictions
  )
  
  metrics <- list("Rsquared" = find_rSquared(predictions_table$Actual, predictions_table$Predicted),
                  "Mean_Sq_Error" = mse(predictions_table$Actual, predictions_table$Predicted),
                  "Mean_Abs_Error" = mae(predictions_table$Actual, predictions_table$Predicted))
  
  return(list("predictions_table" = predictions_table, "metrics" = metrics))
}

### PLOTS CLOCK
# @param predictions == dt containing Actual and Predicted values, such as those generated by clock_pred
plot_clock <- function(predictions, title, scaled = FALSE) {
  print("Plotting predictions...")
  
  plot <- ggplot(predictions$predictions_table, aes(x = Actual, y = Predicted)) +
    geom_point(size = 2.4, alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", color = "red") +
    # coord_fixed(ratio = 1) +
    labs(title = title,
         subtitle = paste0("Adjusted Rsquared: ", round(predictions$metrics$Rsquared, 3), 
                           ", MSE: ", round(predictions$metrics$Mean_Sq_Error, 3), 
                           ", MAE: ", round(predictions$metrics$Mean_Abs_Error, 3))) +
    good_theme_big +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  actual_values <- sort(unique(predictions$predictions_table$Actual))
  predicted_values <- sort(unique(floor(predictions$predictions_table$Predicted)))
  
  if(scaled) { #This section might be broken
    plot <- plot + scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-3, 3)) +
            scale_y_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-3, 3)) +
            xlab("Actual Z-Age") +
            ylab("Predicted Z-Age")
  } else {
    axis_min <- min(c(actual_values, predicted_values)) - 1
    axis_max <- max(c(actual_values, predicted_values)) + 1
    plot <- plot + 
      scale_x_continuous(
        breaks = actual_values,
        limits = c(axis_min, axis_max)
      ) +
      scale_y_continuous(
        breaks = actual_values,
        limits = c(axis_min, axis_max)
      ) +
      coord_equal() +
      xlab("Actual Age") +
      ylab("Predicted Age")
      
  }
  
  return(plot)
}

plot_clock_metrics <- function(metrics, title) {
  print("Plotting metrics...")
  
  metrics_plot <- ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
    geom_point(size = 2, alpha = 0.6, shape = 21, color = "black") + 
    stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.7, width = 0.6) + 
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.3, size = 0.8) + 
    scale_fill_manual(values = c("R-Squared" = "#FFC107", "MSE" = "#648FFF", "MAE" = "#DC267F")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0)) +
    good_theme_big +
    theme(legend.position = "none") +
    labs(title = title, x = NULL, y = "Value")
  
  return(metrics_plot)
}






### CREATE FULL CLOCK & PREDICT ON GIVEN SET
#This function is designed to easily create a full cloick. It thus makes assumptions as to many of the configuration details.
# These details are ones which I kept consistent in MY analysis.
# make_clock 
# @param train_dt == datatable to train the clock on.
# @param test_dt  == datatable to test the clock on. If left as NULL, function assumes the train dt IS the test dt.
# @param response          == string naming the response variable (the 'y' in our regressions). As with all similar parameters, this should correspond to a col in the train_dt
# @param predictors        == list of col names containing predictors
# @param covariates        == list of col names holding covariable variables (these along with the predictors are our 'x' in our regressions)
# @param id                == string containing name of column to use for ID. Only use if running mixed effects modeling.
# @param prediction_method == method used to train prediction model. Choose linear_regression or mixed_effect
# @param title           == name to use for output plots
# @param max_selected    == num of features we want to select. Setting this to Inf will result in a rank deficiency error for most clocks.
# @param rep_seeds             == list of seeds to use to create each clock. The number of seeds will be the number of clocks produced.
# @param scale                 == Whether or not to scale your continuous covariates (such as age)

# @return == list of 1, "clocks": list of all clocks; 2, "metrics_plot": plot of clock metrics; 
#                    3: "selected_features": list of 3.1: "selected_samples": all predictors selected for the clocks; 3.2: "selected_info": info for those predictors
make_clock <- function(train_dt, test_dt = NULL, response, predictors, covariates, id,
                       prediction_method, title = "", max_selected = Inf, rep_seeds = c(), scale = FALSE) {
  set.seed(rep_seeds[1])
  
  # Check passed prediction method
  if(startsWith("linear_regression", prediction_method) || startsWith("lr", prediction_method) || startsWith("lm", prediction_method)) {
    prediction_method <- "linear_regression"
  } else if(startsWith("mixed_effect", prediction_method) || prediction_method == "mixed_fx" || startsWith("mem", prediction_method) || prediction_method == "lme") {
    prediction_method <- "mixed_effect"
  } else {
    stop(paste(prediction_method, "is not a valid prediction type. Please enter 'linear_regression' or 'mixed_effect'."))
  }
  
  #If a covariate only has one value, it doesn't need to be accounted for in our model.
  #This will most often happen because a subset of the entire dataset for one value for a particular covariate is being passed in.
  for(covariate in covariates) {
    if(length(unique(train_dt[[covariate]])) < 2) {
      covariates <- covariates[covariates != covariate]
    }
  }
  
  # Check if we are testing on an independent dataset
  if(!is.null(test_dt) && !identical(train_dt, test_dt)) {
    print("Predicting independent dataset")
    #If so, keep only predictors also in test set
    predictors <- predictors[predictors %in% names(test_dt)]
    
    test_dt <- prep_clock_data(test_dt, response, predictors, id, scale = scale)
    #We may have dropped some predictors. Update the list!
    predictors <- predictors[predictors %in% names(test_dt)]
  }
  #Prep data
  train_dt <- prep_clock_data(train_dt, response, predictors, covariates, id, scale = scale)
  
  #Feature selection
  selected_features <- feature_select(train_dt, response, predictors, covariates, prediction_method, id, max_selected = max_selected, seed = rep_seeds[1])
  selected_dt <- selected_features$selected_dt
  predictors_info <- selected_features$predictors_info_dt
  
  ### CREATE CLOCKS
  clocks <- list()
  for(seed in rep_seeds) {
    set.seed(seed)
    #If 
    if(prediction_method == "linear_regression") {
      train_test_sets <- train_test_split(selected_dt, 0.8, testSet = test_dt)
    } else if (prediction_method == "mixed_effect") {
      train_test_sets <- train_test_split(selected_dt, 0.8, splitOn_col = id, testSet = test_dt)
    } else {
      stop(paste(prediction_method, "this prediction method doesnt exist which shouldve been caught by the code earlier. This really shouldn't be possible! code needs to be tweaked."))
    }
    
    trained_model <- train_clock(train_test_sets$trainingSet, response, id = id, prediction_method = prediction_method, seed = seed)
    train_test_sets$testSet <- train_test_sets$testSet[, intersect(names(train_test_sets$testSet), names(train_test_sets$trainingSet)), with = FALSE]
    
    ### PREDICTING TEST DATA W/ MODEL BUILT ON TRAINING DATA
    predictions <- clock_predict(trained_model, train_test_sets$testSet, response)
    plot <- plot_clock(predictions, title = title)
    
    prediction <- list("predictions" = predictions, "trainingSet" = train_test_sets$trainingSet,
                       "testSet" = train_test_sets$testSet, "plot" = plot)
    clock <- list("model" = trained_model, prediction = prediction)
    print("Clock complete!")
    
    clocks <- c(clocks, list(clock))
  }
  
  print(paste("here is split: ", selected_dt$split))
  
  ### CREATE PLOT FOR AVERAGE METRICS OF CLOCK PREDICTIBILITY
  rsq_list <- c()
  mse_list <- c()
  mae_list <- c()
  
  for(clock in clocks) {
    rsq_list <- c(rsq_list, clock$prediction$predictions$metrics$Rsquared)
    mse_list <- c(mse_list, clock$prediction$predictions$metrics$Mean_Sq_Error)
    mae_list <- c(mae_list, clock$prediction$predictions$metrics$Mean_Abs_Error)
  }
  metrics <- data.frame(
    Metric = rep(c("R-Squared", "MSE", "MAE"), each = length(rsq_list)),
    Value = c(rsq_list, mse_list, mae_list)
  )
  metrics_plot <- plot_clock_metrics(metrics, title)

  #Get rid of non-predictor cols from selected_dt for ease of understanding outside of this function (ncol should equal actual number of predictors)
  selected_dt <- selected_dt[, intersect(predictors, names(selected_dt)), with = FALSE]

  return(list("clocks" = clocks, "metrics_plot" = metrics_plot, 
              "selected_features" = list("selected_samples" = selected_dt, "selected_info" = predictors_info)))
}