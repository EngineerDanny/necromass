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

set.seed(42)

dataname <- "HMPv13"
dataname <- "HMPv35"
dataname <- "TwinsUK"
dataname <- "MovingPictures"
dataname <- "qa10394"
dataname <- "necromass"

task.dt <- data.table::fread("/projects/genomic-ml/da2343/ml_project_1/data/microbe_ds/MovingPictures_11_15.csv")
taxa_columns <- setdiff(names(task.dt), "Group_ID")
new_column_names <- paste0("Taxa", taxa_columns)
setnames(task.dt, old = taxa_columns, new = new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=col_name
  )
  reg.task$col_roles$subset <- "Group_ID"
  reg.task$col_roles$stratum <- "Group_ID"
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  task.list[[task_id]] <- reg.task
}

LearnerRegrFuser <- R6Class("LearnerRegrFuser",
                            inherit = LearnerRegr,
                            public = list( 
                              # Add a class variable to store gamma lookup table
                              gamma_lookup = NULL,
                              
                              initialize = function() {
                                ps = ps(
                                  lambda = p_dbl(0, default = 0.1, tags = "train"),
                                  gamma = p_dbl(0, default = 0.01, tags = "train"),
                                  tol = p_dbl(lower = 0, upper = 1, default = 1e-3, tags = "train"),
                                  num.it = p_int(default = 2000, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = T, tags = "train"),
                                  # Add a parameter to control if we should use stored gammas
                                  use_stored_gammas = p_lgl(default = TRUE, tags = "train")
                                )
                                ps$values = list(lambda = 0.1, gamma = 0.01, 
                                                 tol = 1e-3, num.it = 2000,
                                                 intercept = FALSE, scaling = T,
                                                 use_stored_gammas = TRUE)
                                super$initialize(
                                  id = "regr.fuser",
                                  param_set = ps,
                                  feature_types = c("logical", "integer", "numeric"),
                                  label = "Fuser",
                                  packages = c("mlr3learners", "fuser")
                                )
                                
                                # Initialize gamma lookup table if it exists
                                self$load_gamma_lookup()
                              },
                              
                              # Method to load gamma values from file
                              load_gamma_lookup = function() {
                                gamma_file <- "/projects/genomic-ml/da2343/necromass/movingpictures_task_gammas.rds"
                                if (file.exists(gamma_file)) {
                                  self$gamma_lookup <- readRDS(gamma_file)
                                  return(TRUE)
                                }
                                return(FALSE)
                              },
                              
                              # Method to get gamma for a specific task
                              get_gamma_for_task = function(task_id) {
                                if (is.null(self$gamma_lookup)) {
                                  return(NULL)
                                }
                                
                                if (task_id %in% names(self$gamma_lookup)) {
                                  return(self$gamma_lookup[[task_id]])
                                }
                                
                                return(NULL)
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
                                # Check if we should use stored gamma for this task
                                pv <- self$param_set$get_values(tags = "train")
                                
                                # If use_stored_gammas is TRUE, try to get the gamma from our lookup
                                if (pv$use_stored_gammas) {
                                  task_id <- task$id
                                  stored_gamma <- self$get_gamma_for_task(task_id)
                                  
                                  # If we found a stored gamma for this task, use it
                                  if (!is.null(stored_gamma)) {
                                    pv$gamma <- stored_gamma
                                    cat("Using stored gamma value", stored_gamma, "for task", task_id, "\n")
                                  }
                                }
                                
                                subset_ID <- task$col_roles$subset
                                X_train <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                y_train <- as.matrix(task$data(cols = task$target_names))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                
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
                              
                              # No changes needed to .predict method
                              .predict = function(task) {
                                # Existing predict code remains the same
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



# Add these lines at the beginning of your script
library(future)
library(mlr3)

# Set up parallel processing using all available cores
future::plan("multisession", workers = parallel::detectCores() -1)

# Optional: Reduce logging overhead for better performance
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

# Then keep your original code exactly as is
task_gammas <- list()
for (task_id in names(task.list)) {
  task_fuser_learner <- LearnerRegrFuser$new()
  task_fuser_learner$param_set$values$gamma <- paradox::to_tune(0.01, 1, log=TRUE)
  subtrain.valid.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
  subtrain.valid.cv$param_set$values$folds=5
  subtrain.valid.cv$instantiate(task.list[[task_id]])
  all_iterations <- which(subtrain.valid.cv$instance$iteration.dt$train.subsets == "all")
  custom_all_cv <- ResamplingCustom$new()
  all_train_sets <- list()
  all_test_sets <- list()
  for (i in seq_along(all_iterations)) {
    iter_idx <- all_iterations[i]
    all_train_sets[[i]] <- subtrain.valid.cv$train_set(iter_idx)
    all_test_sets[[i]] <- subtrain.valid.cv$test_set(iter_idx)
  }
  # Instantiate with the extracted sets
  custom_all_cv$instantiate(
    task = task.list[[task_id]], 
    train_sets = all_train_sets,
    test_sets = all_test_sets
  )
  grid.search <- mlr3tuning::TunerGridSearch$new()
  grid.search$param_set$values$resolution <- 3
  
  # Create tuning version to find best parameters
  fuser.learner.tuned <- mlr3tuning::auto_tuner(
    tuner = grid.search,
    learner = task_fuser_learner,
    resampling = custom_all_cv,
    measure = mlr3::msr("regr.mse")
  )
  
  tuning_result <- fuser.learner.tuned$train(task.list[[task_id]])
  best_params <- fuser.learner.tuned$tuning_result
  
  # Store only the gamma value, not the entire learner
  task_gammas[[task_id]] <- exp(best_params$gamma)
  cat("Tuned model for task:", task_id, "- Best gamma:", exp(best_params$gamma), "\n")
}
# Save just the gamma values to a file
saveRDS(task_gammas, "qa10394_task_gammas.rds")




# Then set up your learners with fallback
glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
glmnet_learner$param_set$values$alpha <- 1
glmnet_learner$fallback <- mlr3::LearnerRegrFeatureless$new()
glmnet_learner$encapsulate <- c(train = "evaluate", predict = "evaluate")

reg.learner.list <- list(
  glmnet_learner,
  mlr3::LearnerRegrFeatureless$new(),
  LearnerRegrFuser$new()
)

## For debugging
debug_cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
debug_cv$param_set$values$folds=5
(debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa868817"],
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

## For parallel real tests
# First create the resampling
mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=5
future::plan("sequential")
(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- "movingpictures_29_04"
reg <- batchtools::loadRegistry(reg.dir)
unlink(reg.dir, recursive=TRUE)
reg = batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  packages = "mlr3verse"
)
mlr3batchmark::batchmark(
  reg.bench.grid, store_models = TRUE, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 60*480,#seconds
  memory = 1024,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=TRUE), reg=reg)


batchtools::getStatus(reg=reg)
jobs.after <- batchtools::getJobTable(reg=reg)
table(jobs.after$error)
jobs.after[!is.na(error), .(error, task_id=sapply(prob.pars, "[[", "task_id"))][25:26]

ids <- jobs.after[is.na(error), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)

save(bmr, file="/projects/genomic-ml/da2343/necromass/movingpictures_29_04-bmr.RData")
