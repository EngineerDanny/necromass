library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(fuser)
library(R6)
library(paradox)
library(checkmate)
library(glmnet)
library(ggplot2)

task.dt <- data.table::fread("/projects/genomic-ml/da2343/ml_project_1/data/microbe_ds/necromass_11_15.csv")
taxa_columns <- setdiff(names(task.dt), "Group_ID")
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  # task.list[[task_id]] <- mlr3::TaskRegr$new(
  # task_id, task.dt, target = col_name
  # )$set_col_roles("Group_ID", c("subset", "stratum"))
  
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=col_name
  )
  
  reg.task$col_roles$subset <- "Group_ID"
  # reg.task$col_roles$group <- "Group_ID"
  reg.task$col_roles$stratum <- "Group_ID"
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  #reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id))
  
  #task.list[[task_id]] <- mlr3::TaskRegr$new(
  #  task_id, task.dt, target=task_id
  #)$set_col_roles("Samples",c("subset", "stratum"))
  
  
  task.list[[task_id]] <- reg.task
}

LearnerRegrFuser <- R6Class("LearnerRegrFuser",
                            inherit = LearnerRegr,
                            public = list( 
                              initialize = function() {
                                ps = ps(
                                  lambda = p_dbl(0, default = 0.1, tags = "train"),
                                  gamma = p_dbl(0, default = 0.01, tags = "train"),
                                  tol = p_dbl(lower = 0, upper = 1, default = 1e-3, tags = "train"),
                                  num.it = p_int(default = 2000, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = T, tags = "train")
                                )
                                ps$values = list(lambda = 0.1, gamma = 0.01, 
                                                 tol = 1e-3, num.it = 2000,
                                                 intercept = FALSE, scaling = T)
                                super$initialize(
                                  id = "regr.fuser",
                                  param_set = ps,
                                  feature_types = c("logical", "integer", "numeric"),
                                  label = "Fuser",
                                  packages = c("mlr3learners", "fuser")
                                )
                              },
                              # Helper function for thresholding that can be used in both train and predict
                              threshold_values = function(x, tol, print_info = FALSE) {
                                if(is.null(x)) return(NULL)
                                x <- as.matrix(x)
                                x[is.na(x)] <- 0
                                decimals <- -floor(log10(tol))
                                x <- round(x, decimals)
                                x[abs(x) <= tol] <- 0
                                return(x)
                              }
                            ),
                            private = list(
                              .train = function(task) {
                                subset_ID <- task$col_roles$subset
                                X_train <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                y_train <- as.matrix(task$data(cols = task$target_names))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                pv <- self$param_set$get_values(tags = "train")
                                
                                k <- as.numeric(length(unique(group_ind)))
                                glmnet_model <- mlr3learners::LearnerRegrCVGlmnet$new()$train(task)$model
                                
                                print(pv$gamma)
                                tryCatch({
                                  # Use model parameters from pv instead of hardcoded values
                                  fuser_params <- fuser::fusedLassoProximal(X_train, y_train, group_ind, 
                                                                            lambda = pv$lambda, 
                                                                            G = matrix(1, k, k), 
                                                                            gamma = pv$gamma,
                                                                            num.it = 20000,
                                                                            #intercept = T,
                                                                            #tol = 1e-6,
                                                                            c.flag = TRUE,
                                                                            scaling = pv$scaling)
                                  
                                  # Extract and threshold beta and intercept
                                  beta.estimate <- head(fuser_params, -1)
                                  intercept <- tail(fuser_params, 1)
                                  
                                  # Process beta estimates
                                  if(!is.null(beta.estimate)) {
                                    beta.estimate <- as.matrix(beta.estimate)
                                    colnames(beta.estimate) <- colnames(fuser_params)
                                    
                                    # Apply thresholding during training
                                    beta.estimate <- self$threshold_values(beta.estimate, pv$tol)
                                    intercept <- self$threshold_values(intercept, pv$tol)
                                  }
                                  
                                  self$model <- list(
                                    beta = beta.estimate,
                                    intercept = intercept,
                                    glmnet_model = glmnet_model,
                                    model_type = "fuser",
                                    formula = task$formula(),
                                    data = task$data(),
                                    pv = pv,
                                    groups = group_ind
                                  )
                                  
                                }, error = function(e) {
                                  print("error")
                                  print(e)
                                  self$model <- list(
                                    glmnet_model = glmnet_model,
                                    model_type = "glmnet",
                                    formula = task$formula(),
                                    data = task$data(),
                                    pv = pv,
                                    groups = group_ind
                                  )
                                })
                                
                                self$model
                              },
                              .predict = function(task) {
                                subset_ID <- task$col_roles$subset
                                X_test <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                new_X_test = as.matrix(task$data(cols = task$feature_names))
                                X = apply(X_test, 2, as.numeric)
                                
                                beta = self$model$beta
                                intercept = self$model$intercept
                                
                                # No need to threshold again since values are already thresholded during training
                                model_type <- self$model$model_type
                                train_group_ind <- self$model$groups
                                y.predict <- rep(NA, nrow(X))
                                
                                group_names <- colnames(beta)[unique(group_ind)]
                                
                                for (name in group_names) {
                                  group_rows <- which(group_ind == name)
                                  if (length(group_rows) > 0) {
                                    X.group <- X[group_rows, , drop = FALSE]
                                    y.predict[group_rows] <- as.numeric(X.group %*% beta[, name] + intercept[, name])
                                  } else {
                                    warning(paste("No rows found for group", name))
                                  }
                                }
                                
                                if (any(is.na(y.predict))) {
                                  print("Used glmnet.")
                                  y.predict.new <- as.vector(predict(self$model$glmnet_model, newx = new_X_test))
                                } else {
                                  print("Used normal")
                                  y.predict.new <- y.predict
                                }
                                
                                list(response = y.predict.new)
                              }
                            )
)

set.seed(42)
#mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
#mycv$param_set$values$folds=5
#for(task in task.list){
#  mycv$instantiate(task)
#}
#subtrain.valid.cv <- mlr3::ResamplingCV$new()
#subtrain.valid.cv$param_set$values$folds <- 2

# First create the resampling
mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=5
mycv$instantiate(task.list[["TaxaTomentella"]])
all_iterations <- which(mycv$instance$iteration.dt$train.subsets == "all")
custom_all_cv <- ResamplingCustom$new()
# Get all train and test sets first
all_train_sets <- list()
all_test_sets <- list()
for (i in seq_along(all_iterations)) {
  iter_idx <- all_iterations[i]
  all_train_sets[[i]] <- mycv$train_set(iter_idx)
  all_test_sets[[i]] <- mycv$test_set(iter_idx)
}
# Instantiate with the extracted sets
custom_all_cv$instantiate(
  task = task.list[["TaxaTomentella"]], 
  train_sets = all_train_sets,
  test_sets = all_test_sets
)

fuser.learner =  LearnerRegrFuser$new()
fuser.learner$param_set$values$gamma <- paradox::to_tune(0.01, 1, log=TRUE)
grid.search.5 <- mlr3tuning::TunerGridSearch$new()
grid.search.5$param_set$values$resolution <- 5
fuser.learner.tuned = mlr3tuning::auto_tuner(
  tuner = grid.search.5,
  learner = fuser.learner,
  resampling = custom_all_cv,
  measure = mlr3::msr("regr.mse"))
tuning_result <- fuser.learner.tuned$train(task.list[["TaxaTomentella"]])
gamma <- exp( fuser.learner.tuned$tuning_result$gamma ) 
fuser.learner_fixed <- LearnerRegrFuser$new()
fuser.learner_fixed$param_set$values$gamma <- gamma
tuning_results <- fuser.learner.tuned$tuning_instance$archive$data
tuning_results$gamma_actual <- exp(tuning_results$gamma)
tuning_results <- tuning_results[order(tuning_results$gamma_actual)]
tuning_results

# Now plot the model complexity curve
get_training_error <- function(gamma_value, task) {
  learner <- LearnerRegrFuser$new()
  learner$param_set$values$gamma <- gamma_value
  learner$train(task)
  pred <- learner$predict(task)
  return(pred$score(mlr3::msr("regr.mse")))
}
# Apply this function to each gamma value
tuning_results$train_mse <- sapply(tuning_results$gamma_actual, 
                                   function(g) get_training_error(g, task.list[["TaxaTomentella"]]))


min_val_idx <- which.min(tuning_results$regr.mse)
min_gamma <- tuning_results$gamma_actual[min_val_idx]
min_mse <- tuning_results$regr.mse[min_val_idx]

p <- ggplot(tuning_results) +
  geom_line(aes(x = gamma_actual, y = regr.mse, color = "Validation"), linewidth = 1) +
  geom_point(aes(x = gamma_actual, y = regr.mse, color = "Validation"), size = 3) +
  geom_line(aes(x = gamma_actual, y = train_mse, color = "Subtrain"), linewidth = 1) +
  geom_point(aes(x = gamma_actual, y = train_mse, color = "Subtrain"), size = 3) +
  geom_vline(xintercept = min_gamma, linetype = "dashed", color = "darkgrey") +
  annotate("text", x = min_gamma, y = max(tuning_results$regr.mse, tuning_results$train_mse) * 1, 
           label = paste("Optimal γ =", round(min_gamma, 3)), 
           hjust = -0.1, color = "darkgrey") +
  scale_x_continuous(
    trans = scales::compose_trans(scales::log10_trans(), scales::reverse_trans()),
    breaks = c(0.01, 0.1, 1.0),
    labels = c("0.01", "0.1", "1.00")
  ) +
  labs(title = "TaxaTomentella",
       x = "Gamma (γ)",
       y = "Mean Squared Error",
       color = "Loss") +
  theme_minimal() +
  scale_color_manual(values = c("Validation" = "blue", "Subtrain" = "red"))
ggsave("TaxaTomentella_model_complexity_curve_2.png", p, width = 6, height = 5, dpi = 300)

# Then set up your learners with fallback
glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
glmnet_learner$param_set$values$alpha <- 1
glmnet_learner$fallback <- mlr3::LearnerRegrFeatureless$new()
glmnet_learner$encapsulate <- c(train = "evaluate", predict = "evaluate")

reg.learner.list <- list(
  glmnet_learner,
  mlr3::LearnerRegrFeatureless$new(),
  #fuser.learner.tuned
  fuser.learner_fixed
  #LearnerRegrFuser$new()
)

## For debugging
debug_cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
debug_cv$param_set$values$folds=5
(debug.grid <- mlr3::benchmark_grid(
  task.list["TaxaTomentella"],
  reg.learner.list,
  debug_cv))
debug.result <- mlr3::benchmark(debug.grid)
debug.score.dt <- mlr3resampling::score(debug.result)
filtered_results <- debug.score.dt[debug.score.dt$train.subsets != "other", ]
aggregate_results <- filtered_results[, .(
  mean_mse = mean(regr.mse, na.rm = TRUE),
  sd_mse = sd(regr.mse, na.rm = TRUE),
  n_iterations = .N
), by = .(learner_id, train.subsets)]

print(aggregate_results)
