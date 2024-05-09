library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(fuser)
library(R6)
library(paradox)
library(checkmate)
library(glmnet)

#task.dt <- data.table::fread("HMPv35_otu_table_sub_log.csv")
task.dt <- data.table::fread("HMPv35_sub_log.csv")
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
                                  tol = p_dbl(lower = 0, upper = 1, default = 3e-3, tags = "train"),
                                  num.it = p_int(default = 1000, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = T, tags = "train")
                                )
                                ps$values = list(lambda = 0.01, gamma = 0.01, 
                                                 tol = 3e-3, num.it = 1000,
                                                 intercept = FALSE, scaling = T)
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
                                # Assuming 'task' is a task object with data, feature names, and group information
                                subset_ID <- task$col_roles$subset
                                X_train <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                y_train <- as.matrix(task$data(cols = task$target_names))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                #print("Task summary:")
                                #str(reg.task$col_roles)
                                #print("group_ind")
                                #print(group_ind)
                                #print(typeof(group_ind))
                                
                                #group_ind <- task$groups$group
                                #group_ind <- task$groups$subset
                                k <- as.numeric(length(unique(group_ind))) # number of groups
                                n.group <- min(table(group_ind)) # minimum samples per group across all groups
                                
                                
                                # Initialize a new vector to store the balanced group indices
                                balanced_group_ind <- integer(0)
                                
                                # Loop through each unique group and select 'n.group' indices for each
                                for (grp in unique(group_ind)) {
                                  grp_indices <- which(group_ind == grp)
                                  balanced_group_ind <- c(balanced_group_ind, sample(grp_indices, n.group))
                                }
                                
                                # Now 'balanced_group_ind' contains an equal number of indices for each group
                                #group_ind <- group_ind[balanced_group_ind]
                                
                                # Print the final group index
                                #print("updated group_ind")
                                #print(group_ind)
                                
                                # Print the final group index
                                #print("balanced_group_ind")
                                #print(balanced_group_ind)
                                
                                pv <- self$param_set$get_values(tags = "train")
                                
                                # Create the matrices with appropriate dimensions
                                
                                #y <- matrix(y_train[balanced_group_ind], nrow = n.group, ncol = k, byrow = TRUE)
                                #X <- X_train[balanced_group_ind, ]
                                
                                # Initialize glmnet_model before tryCatch
                                glmnet_model <- mlr3learners::LearnerRegrCVGlmnet$new()$train(task)$model
                                
                                tryCatch({
                                  # Use fuser::fusedLassoProximal
                                  # Generate block diagonal matrices for L2 fusion approach
                                  # transformed.data = generateBlockDiagonalMatrices(X, y, group_ind, matrix(1, k, k))
                                  # Use L2 fusion to estimate betas (with near-optimal information sharing among groups)
                                  # beta.estimate = fusedL2DescentGLMNet(transformed.data$X, transformed.data$X.fused, 
                                  #                                     transformed.data$Y, group_ind, lambda=pv$lambda,
                                  #                                     gamma=pv$gamma)
                                  # colnames(beta.estimate) = as.character(sort(unique(group_ind)))
                                  # beta.estimate = beta.estimate[,,1]
                                  # colnames(beta.estimate) = as.character(sort(unique(group_ind)))
                                  # print("beta.estimate")
                                  # print(beta.estimate)
                                  
                                  # Use fuser::fusedLassoProximal
                                  fuser_params <- fuser::fusedLassoProximal(X_train, y_train, group_ind, 
                                                                             lambda = 0.01, 
                                                                             G = matrix(1, k, k), 
                                                                             gamma = 0.1,
                                                                             num.it = 20000,
                                                                             intercept = T,
                                                                             scaling = pv$scaling)
                                  beta.estimate <- head(fuser_params, -1)
                                  intercept <- tail(fuser_params, 1)
                                  
                                  # Update self$model with both beta.estimate and glmnet_model
                                  self$model <- list(beta = beta.estimate, 
                                                     intercept = intercept,
                                                     glmnet_model = glmnet_model,
                                                     model_type = "fuser",
                                                     formula = task$formula(),
                                                     data = task$data(),
                                                     pv = pv,
                                                     groups = group_ind)
                                  
                                }, error = function(e) {
                                  print("error")
                                  print(e)
                                  #print("updated group_ind")
                                 # print(group_ind)
                                  # Update self$model with glmnet_model in case of error
                                  self$model <- list(glmnet_model = glmnet_model,
                                                     model_type = "glmnet",
                                                     formula = task$formula(),
                                                     data = task$data(),
                                                     pv = pv,
                                                     groups = group_ind)
                                })
                                

                                self$model
                              },
                              .predict = function(task) {
                                ordered_features = function(task, learner) {
                                  cols = names(learner$state$data_prototype)
                                  task$data(cols = intersect(cols, task$feature_names))
                                }
                                
                                subset_ID <- task$col_roles$subset
                                X_test <- as.matrix(task$data(cols = setdiff(task$feature_names, subset_ID)))
                                group_ind <- as.matrix(task$data(cols = subset_ID))
                                new_X_test = as.matrix(task$data(cols = task$feature_names))
                                X = apply(X_test, 2, as.numeric)
                                #group_ind <- task$groups$group
                               
                                beta = self$model$beta
                                intercept = self$model$intercept
                                model_type <- self$model$model_type
                                train_group_ind <- self$model$groups
                                y.predict <- rep(NA, nrow(X))
                                
                                group_names <- colnames(beta)[unique(group_ind)]
                                #print("group_names")
                                #print(group_names)
                                
                                for (name in group_names) {
                                  group_rows <- which(group_ind == name)
                                  if (length(group_rows) > 0) {
                                    X.group <- X[group_rows, , drop = FALSE]
                                    y.predict[group_rows] <- as.numeric(X.group %*% beta[, name] + intercept[, name])
                                  } else {
                                    warning(paste("No rows found for group", name))
                                  }
                                }
                                
                                
                                # Check for NAs or errors and use glmnet_model if needed
                                if (any(is.na(y.predict))) {
                                  print("Used glmnet.")
                                  #print(train_group_ind)
                                  y.predict.new <- as.vector(predict(self$model$glmnet_model, newx = new_X_test))
                                  #print("NA values detected in predictions. Falling back to glmnet.")
                                }else{
                                  print("Used normal")
                                  #print(train_group_ind)
                                  y.predict.new <- y.predict
                                }
                                #y.predict <- as.vector(predict(self$model$glmnet_model, newx = X_test))
                                # Return the predictions as a numeric vector
                                list(response = y.predict.new)
                              }
                            )
)


mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=5
#for(task in task.list){
#  mycv$instantiate(task)
#}

fuser.learner =  LearnerRegrFuser$new()
fuser.learner$param_set$values$lambda <- paradox::to_tune(0.001, 1, log=TRUE)
fuser.learner$param_set$values$gamma <- paradox::to_tune(0.001, 1, log=TRUE)
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
random.search <- mlr3tuning::TunerRandomSearch$new()
termination_criterion <- mlr3tuning::trm("evals", n_evals = 20)
fuser.learner.tuned = mlr3tuning::auto_tuner(
  tuner = random.search,
  learner = fuser.learner,
  resampling = subtrain.valid.cv,
  measure = mlr3::msr("regr.mse"),
  terminator = termination_criterion
)
reg.learner.list <- list(
  mlr3learners::LearnerRegrCVGlmnet$new(),
  #mlr3learners::LearnerRegrGlmnet$new(),
  mlr3::LearnerRegrFeatureless$new(), 
  #fuser.learner.tuned
  #LearnerRegrFuserGlmnet$new(),
  LearnerRegrFuser$new()
)


(debug.grid <- mlr3::benchmark_grid(
  task.list["Taxa356733"],
  reg.learner.list,
  mycv))
future::plan("sequential")
debug.result <- mlr3::benchmark(debug.grid)
debug.score.dt <- mlr3resampling::score(debug.result)

(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))
reg.dir <- "HMPv35-2024-05-01-952-benchmark-reg"
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
  memory = 512,#megabytes per cpu
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

save(bmr, file="HMPv35-2024-05-01-952-benchmark.RData")
