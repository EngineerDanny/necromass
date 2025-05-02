library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(fuser)
library(R6)
library(paradox)
library(checkmate)
library(glmnet)

#task.dt <- data.table::fread("/projects/genomic-ml/da2343/ml_project_1/data/microbe_ds/HMPv35_11_15.csv")
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
  #reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Group_ID"))
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id))
  
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
                                  lambda = p_dbl(0, default = 0.01, tags = "train"),
                                  gamma = p_dbl(0, default = 0.01, tags = "train"),
                                  tol = p_dbl(lower = 0, upper = 1, default = 1e-3, tags = "train"),
                                  num.it = p_int(default = 2000, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = T, tags = "train")
                                )
                                ps$values = list(lambda = 0.01, gamma = 0.01, 
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
                                
                                if(print_info) {
                                  print("Before thresholding:")
                                  print(head(x))
                                  print(paste("Tolerance value:", tol))
                                }
                                
                                # Handle NAs
                                x[is.na(x)] <- 0
                                
                                # Round to fixed decimals based on tolerance
                                decimals <- -floor(log10(tol))
                                x <- round(x, decimals)
                                
                                # Zero out small values strictly
                                x[abs(x) <= tol] <- 0
                                
                                if(print_info) {
                                  print("After thresholding:")
                                  print(head(x))
                                }
                                
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
                                
                                tryCatch({
                                  # Use model parameters from pv instead of hardcoded values
                                  fuser_params <- fuser::fusedLassoProximal(X_train, y_train, group_ind, 
                                                                            lambda = 0.1, 
                                                                            G = matrix(1, k, k), 
                                                                            gamma = 0.01,
                                                                            num.it = 20000,
                                                                            intercept = T,
                                                                            tol = 1e-6,
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


mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=5
#for(task in task.list){
#  mycv$instantiate(task)
#}

#fuser.learner =  LearnerRegrFuser$new()
#fuser.learner$param_set$values$lambda <- paradox::to_tune(0.001, 1, log=TRUE)
#fuser.learner$param_set$values$gamma <- paradox::to_tune(0.001, 1, log=TRUE)
#fuser.learner$param_set$values$tol <- paradox::to_tune(1e-10, 1e-2, log=TRUE)

subtrain.valid.cv <- mlr3::ResamplingCV$new()
subtrain.valid.cv$param_set$values$folds <- 2

#grid.search.5 <- mlr3tuning::TunerGridSearch$new()
#grid.search.5$param_set$values$resolution <- 3
#fuser.learner.tuned = mlr3tuning::auto_tuner(
#  tuner = grid.search.5,
#  learner = fuser.learner,
#  resampling = subtrain.valid.cv,
#  measure = mlr3::msr("regr.mse"))
#reg.learner.list <- list(
#  mlr3::LearnerRegrFeatureless$new(), 
#  fuser.learner.tuned )

# Set up random search tuner
#random.search <- mlr3tuning::TunerRandomSearch$new()
#termination_criterion <- mlr3tuning::trm("evals", n_evals = 20)
#fuser.learner.tuned = mlr3tuning::auto_tuner(
#  tuner = random.search,
#  learner = fuser.learner,
#  resampling = subtrain.valid.cv,
#  measure = mlr3::msr("regr.mse"),
#  terminator = termination_criterion
#)
reg.learner.list <- list(
  mlr3learners::LearnerRegrCVGlmnet$new(),
  #mlr3learners::LearnerRegrGlmnet$new(),
  mlr3::LearnerRegrFeatureless$new(), 
  #fuser.learner.tuned
  #LearnerRegrFuserGlmnet$new(),
  LearnerRegrFuser$new()
)

## For debugging
(debug.grid <- mlr3::benchmark_grid(
  task.list["TaxaMycobacterium"],
  reg.learner.list,
  mycv))
debug.result <- mlr3::benchmark(debug.grid)
debug.score.dt <- mlr3resampling::score(debug.result)

future::plan("sequential")
(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- "necromass_29_04"
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

## 1: Error in approx(lambda, seq(lambda), sfrac) : \n  need at least two non-NA values to interpolate
## 2: Error in elnet(xd, is.sparse, y, weights, offset, type.gaussian, alpha,  : \n  y is constant; gaussian glmnet fails at standardization step
##    task_id
## 1: Kaistia
## 2: Kaistia

ids <- jobs.after[is.na(error), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- mlr3resampling::score(bmr)
#tab = as.data.table(bmr)
save(bmr, file="/projects/genomic-ml/da2343/necromass/necromass_29_04-benchmark.RData")

#ids <- jobs.after[is.na(error), job.id]
#bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
#score.dt <- mlr3resampling::score(bmr)
save(bmr, score.dt, file="/projects/genomic-ml/da2343/necromass/necromass_29_04-benchmark.RData")


