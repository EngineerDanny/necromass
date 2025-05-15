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
library(Matrix)

set.seed(42)

source("/projects/genomic-ml/da2343/necromass/l2_fusion.R") 

dataname <- "HMPv13"
dataname <- "HMPv35"
dataname <- "TwinsUK"
dataname <- "MovingPictures"
dataname <- "qa10394"
dataname <- "necromass"

task.dt <- data.table::fread( paste0("/projects/genomic-ml/da2343/ml_project_1/data/microbe_ds/", dataname, "_11_15.csv") )
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
                              initialize = function() {
                                ps = ps(
                                  lambda = p_dbl(0, default = 0.01, tags = "train"),
                                  gamma = p_dbl(0, default = 0.1, tags = "train"),
                                  tol = p_dbl(lower = 0, upper = 1, default = 1e-06, tags = "train"),
                                  num.it = p_int(default = 2000, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = TRUE, tags = "train")
                                )
                                ps$values = list(lambda = 0.01, gamma = 0.1, 
                                                 tol = 1e-06, num.it = 2000,
                                                 intercept = FALSE, scaling = TRUE)
                                super$initialize(
                                  id = "regr.fuser",
                                  param_set = ps,
                                  feature_types = c("logical", "integer", "numeric"),
                                  label = "Fuser",
                                  packages = c("mlr3learners", "fuser")
                                )
                              }
                            ),
                            private = list(
                              .train = function(task) {
                                # Check if we should use stored gamma for this task
                                pv <- self$param_set$get_values(tags = "train")
                                subset_ID <- task$col_roles$subset
                                X_train <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                y_train <- as.matrix(task$data(cols = task$target_names))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                
                                if(length(unique(group_ind)) <= 1) {
                                  glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
                                  glmnet_learner$param_set$values$grouped = FALSE
                                  glmnet_model <- glmnet_learner$train(task)$model
                                  
                                  self$model <- list(
                                    glmnet_model = glmnet_model,
                                    model_type = "glmnet",
                                    formula = task$formula(),
                                    data = task$data(),
                                    pv = pv,
                                    groups = group_ind,
                                    single_group = TRUE  # Flag for prediction
                                  )
                                  
                                  return(self$model)
                                }
                                
                                
                                beta.estimate <- fusedL2DescentGLMNet(X_train, y_train, group_ind,
                                                                      lambda = pv$lambda, 
                                                                      gamma = pv$gamma,
                                                                      scaling = pv$scaling)
                                intercept <- matrix(0, nrow=1, ncol=ncol(beta.estimate))
                                colnames(intercept) <- colnames(beta.estimate)
                                self$model <- list(
                                  beta = beta.estimate,
                                  intercept = intercept,
                                  glmnet_model = NULL,
                                  model_type = "fuser",
                                  formula = task$formula(),
                                  data = task$data(),
                                  pv = pv,
                                  groups = group_ind,
                                  single_group = FALSE  # Flag for prediction
                                )
                             
                                
                                self$model
                              },
                              # No changes needed to .predict method
                              .predict = function(task) {
                                subset_ID <- task$col_roles$subset
                                X_test <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                X_test <- apply(X_test, 2, as.numeric)
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                new_X_test = as.matrix(task$data(cols = task$feature_names))
                                
                                # Convert group_ind to vector if it's a matrix
                                group_ind_vec <- as.vector(group_ind)
                                train_groups_vec <- as.vector(self$model$groups)
                                
                                
                                # Check if we used glmnet due to single group
                                if(self$model$single_group) {
                                  y.predict <- as.vector(predict(self$model$glmnet_model, newx = new_X_test))
                                  return(list(response = y.predict))
                                }
                                
                                # Use predictFusedL2
                                y.predict <- predictFusedL2(
                                  beta.mat = self$model$beta,
                                  newX = X_test,
                                  newGroups = group_ind_vec,
                                  groups = train_groups_vec
                                )
                                
                                if (any(is.na(y.predict))) {
                                  print("Used glmnet.")
                                  y.predict <- as.vector(predict(self$model$glmnet_model, newx = new_X_test))
                                } else {
                                  print("Used normal")
                                }
                                
                                list(response = y.predict)
                              }
                              
                            )
)


# Add these lines at the beginning of your script
library(future)
library(mlr3)
#future::plan("multisession", workers = parallel::detectCores() - 1)
withr::with_seed(42, {
  subtrain.valid.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
  subtrain.valid.cv$param_set$values$folds=2
  subtrain.valid.cv$param_set$values$subsets = "A"
  subtrain.valid.cv$instantiate(task.list[[task_id]])
  grid.search <- mlr3tuning::TunerBatchGridSearch$new()
  grid.search$param_set$values$resolution <- 5
})

task_gammas <- list()
for (task_id in names(task.list)) {
  task_fuser_learner <- LearnerRegrFuser$new()
  #task_fuser_learner$param_set$values$gamma <- paradox::to_tune(0.01, 1, log=TRUE)
  task_fuser_learner$param_set$values$gamma <- paradox::to_tune(levels = c(0.001, 0.01, 0.1, 1, 10))
  #task_fuser_learner$param_set$values$lambda <- paradox::to_tune(0.01, 1, log=TRUE)
  # Create tuning version to find best parameters
  fuser.learner.tuned <- mlr3tuning::auto_tuner(
    tuner = grid.search,
    learner = task_fuser_learner,
    resampling = subtrain.valid.cv,
    measure = mlr3::msr("regr.mse")
  )
  print("got here")
  tuning_result <- fuser.learner.tuned$train(task.list[[task_id]])
  best_params <- fuser.learner.tuned$tuning_result
  #task_gammas[[task_id]] <- exp(best_params$gamma)
  task_gammas[[task_id]] <- best_params$gamma
  #print(tuning_result$tuning_result
  cat("Tuned model for task:", task_id, "- Best gamma:", best_params$gamma, "\n")
  #cat("Tuned model for task:", task_id, "- Best lambda:", exp(best_params$lambda), "\n")
  break
}


best_params
# Save just the gamma values to a file
saveRDS(task_gammas, "qa10394_task_gammas_1.rds")


## For debugging
withr::with_seed(42, {
  subtrain.valid.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
  subtrain.valid.cv$param_set$values$folds=2
  subtrain.valid.cv$param_set$values$subsets = "A"
  #subtrain.valid.cv$instantiate(task.list[[task_id]])
  grid.search <- mlr3tuning::TunerBatchGridSearch$new()
  grid.search$param_set$values$resolution <- 5
  
  task_fuser_learner <- LearnerRegrFuser$new()
  task_fuser_learner$param_set$values$gamma <- paradox::to_tune(levels = c(0.001, 0.01, 0.1, 1))
  #task_fuser_learner$param_set$values$lambda <- paradox::to_tune(0.01, 1, log=TRUE)
  fuser.learner.tuned <- mlr3tuning::auto_tuner(
    tuner = grid.search,
    learner = task_fuser_learner,
    resampling = subtrain.valid.cv,
    measure = mlr3::msr("regr.mse")
  )
  glmnet_learner <- mlr3learners::LearnerRegrCVGlmnet$new()
  glmnet_learner$param_set$values$nfolds <- 3
  glmnet_learner$param_set$values$grouped <- T
  glmnet_learner$encapsulate(
    method   = "evaluate",                
    fallback = mlr3::LearnerRegrFeatureless$new()
  )
  fuser_learner <- LearnerRegrFuser$new()
  fuser_learner$encapsulate(
    method   = "evaluate",                
    fallback = mlr3::LearnerRegrFeatureless$new()
  )
  reg.learner.list <- list(
    glmnet_learner,
    #fuser_learner,
    fuser.learner.tuned,
    mlr3::LearnerRegrFeatureless$new()
  )
  debug_cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
  debug_cv$param_set$values$folds=5
  debug_cv$param_set$values$subsets = "A"
  #debug_task <- task.list[["TaxaOvicillium"]]
  #debug_cv$instantiate(debug_task)
  (debug.grid <- mlr3::benchmark_grid(
    task.list[["TaxaRhizobium"]],
    reg.learner.list,
    debug_cv))
  debug.result <- mlr3::benchmark(debug.grid, callbacks = list())
  debug.score.dt <- mlr3resampling::score(debug.result)
  aggregate_results <- debug.score.dt[, .(
    mean_mse = mean(regr.mse, na.rm = TRUE),
    sd_mse = sd(regr.mse, na.rm = TRUE),
    n_iterations = .N
  ), by = .(learner_id, train.subsets)]
  
  print(aggregate_results)
})



## For real
set.seed(42)
mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=5
mycv$param_set$values$subsets = "SA"
future::plan("sequential")
(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- paste0(dataname, "_11_05")
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

debug.score.dt <- mlr3resampling::score(bmr)
aggregate_results <- debug.score.dt[, .(
  mean_mse = mean(regr.mse),
  sd_mse = sd(regr.mse),
  n_iterations = .N
), by = .(learner_id, train.subsets)]
print(aggregate_results)


save(bmr, file="/projects/genomic-ml/da2343/necromass/MovingPictures_11_05-bmr.RData")
