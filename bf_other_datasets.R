library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(fuser)
library(R6)
library(paradox)

movingPictures_sub_log <- data.table::fread("MovingPictures_sub_log.csv")
taxa_columns <- setdiff(names(movingPictures_sub_log), "Group_ID")
new_column_names <- paste0("Taxa", taxa_columns)
setnames(movingPictures_sub_log, old = taxa_columns, new = new_column_names)
task.list <- list()
for (col_name in new_column_names) {
  task_id <- col_name
  task.list[[task_id]] <- mlr3::TaskRegr$new(
    task_id, movingPictures_sub_log, target = col_name
  )$set_col_roles("Group_ID", c("group", "stratum"))
}


LearnerRegrFuser <- R6Class("LearnerRegrFuser",
                            inherit = LearnerRegr,
                            public = list( 
                              initialize = function() {
                                ps = ps(
                                  lambda = p_dbl(0, default = 0.001, tags = "train"),
                                  gamma = p_dbl(0, default = 0.001, tags = "train"),
                                  tol = p_dbl(lower = 0, upper = 1, default =  9e-5, tags = "train"),
                                  num.it = p_int(default = 500, tags = "train"),
                                  intercept = p_lgl(default = FALSE, tags = "train"),
                                  scaling = p_lgl(default = FALSE, tags = "train")
                                )
                                ps$values = list(lambda = 0.001, gamma = 0.001, 
                                                 tol = 9e-5, num.it = 5000,
                                                 intercept = FALSE, scaling = FALSE)
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
                                X_train <- as.matrix(task$data(cols = task$feature_names))
                                y_train <- as.matrix(task$data(cols = task$target_names))
                                group_ind <- task$groups$group
                                k <- as.numeric(length(unique(group_ind))) # number of groups
                                n.group <- min(table(group_ind)) # minimum samples per group across all groups
                                max_elements <- n.group * k  # Calculate the maximum number of elements to fit the matrix
                                group_ind <- group_ind[1:max_elements]
                                pv <- self$param_set$get_values(tags = "train")
                                # Create the matrix with appropriate dimensions
                                # Assuming y is your original vector and you have n.group and k defined
                                y <- matrix(y_train[1:max_elements], nrow = n.group, ncol = k, byrow = TRUE)
                                # Check if X has more than one column
                                if (ncol(X_train) > 1) {
                                  X <- apply(X_train, 2, as.numeric)
                                  X <- X[1:max_elements, ]
                                } else {
                                  X <- head(X_train, max_elements)
                                }
                                
                                # Initialize glmnet_model before tryCatch
                                glmnet_model <- mlr3learners::LearnerRegrCVGlmnet$new()$train(task)$model
                                
                                tryCatch({
                                  # Use fuser::fusedLassoProximal
                                  beta.estimate <- fuser::fusedLassoProximal(X, y, group_ind, 
                                                                             lambda = pv$lambda, 
                                                                             G = matrix(1, k, k), 
                                                                             gamma = pv$gamma,
                                                                             tol = pv$tol, 
                                                                             num.it = pv$num.it,
                                                                             intercept = FALSE,
                                                                             scaling = pv$scaling) 
                                  # Update self$model with both beta.estimate and glmnet_model
                                  self$model <- list(beta = beta.estimate, 
                                                     glmnet_model = glmnet_model,
                                                     formula = task$formula(),
                                                     data = task$data(),
                                                     pv = pv,
                                                     groups = group_ind)
                                  
                                }, error = function(e) {
                                  # Update self$model with glmnet_model in case of error
                                  self$model <- list(glmnet_model = glmnet_model,
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
                                X_test = as.matrix(task$data(cols = task$feature_names))
                                X = apply(X_test, 2, as.numeric)
                                group_ind <- task$groups$group
                                k <- as.numeric(length(unique(group_ind)))
                                beta = self$model$beta
                                y.predict <- rep(NA, nrow(X))
                                
                                # Attempt to predict using coefficients
                                tryCatch({
                                  for (k.i in 1:k) {
                                    group_rows <- which(group_ind == k.i)
                                    X.group <- X[group_rows, , drop = FALSE]
                                    y.predict[group_rows] <- as.numeric(X.group %*% beta[, k.i])
                                  }
                                }, error = function(e) {
                                  # If an error occurs, print a warning and leave y.predict as NA
                                  warning("Error occurred with beta coefficient prediction: ", e$message)
                                })
                                
                                # Check for NAs or errors and use glmnet_model if needed
                                if (any(is.na(y.predict))) {
                                  warning("NA values detected in predictions. Falling back to glmnet.")
                                  glmnet_model <- self$model$glmnet_model
                                  y.predict <- predict(glmnet_model, newx = X_test)
                                  y.predict <- as.vector(y.predict)
                                }
                                
                                # Return the predictions as a numeric vector
                                list(response = y.predict)
                              }
                            )
)


MyResampling = R6::R6Class("Resampling",
                           public = list(
                             id = NULL,
                             label = NULL,
                             param_set = NULL,
                             instance = NULL,
                             task_hash = NA_character_,
                             task_nrow = NA_integer_,
                             duplicated_ids = NULL,
                             man = NULL,
                             initialize = function(id, param_set = ps(), duplicated_ids = FALSE, label = NA_character_, man = NA_character_) {
                               self$id = checkmate::assert_string(id, min.chars = 1L)
                               self$label = checkmate::assert_string(label, na.ok = TRUE)
                               self$param_set = paradox::assert_param_set(param_set)
                               self$duplicated_ids = checkmate::assert_flag(duplicated_ids)
                               self$man = checkmate::assert_string(man, na.ok = TRUE)
                             },
                             format = function(...) {
                               sprintf("<%s>", class(self)[1L])
                             },
                             print = function(...) {
                               cat(format(self), if (is.null(self$label) || is.na(self$label)) "" else paste0(": ", self$label))
                               cat("\n* Iterations:", self$iters)
                               cat("\n* Instantiated:", self$is_instantiated)
                               cat("\n* Parameters:\n")
                               str(self$param_set$values)
                             },
                             help = function() {
                               self$man
                             },
                             instantiate = function(task) {
                               task = mlr3::assert_task(mlr3::as_task(task))
                               folds = private$.combine(lapply(task$strata$row_id, private$.sample, task = task))
                               id.fold.groups <- folds[task$groups, on="row_id"]
                               uniq.fold.groups <- setkey(unique(id.fold.groups[, .(
                                 test.fold=fold, test.group=group)]))
                               self$instance <- list(
                                 iteration.dt=data.table(train.groups=c("all","same","other"))[
                                   , data.table(uniq.fold.groups), by=train.groups][, iteration := .I],
                                 id.dt=id.fold.groups)
                               for(iteration.i in 1:nrow(self$instance$iteration.dt)){
                                 split.info <- self$instance$iteration.dt[iteration.i]
                                 is.set.group <- list(
                                   test=id.fold.groups[["group"]] == split.info[["test.group"]])
                                 is.set.group[["train"]] <- switch(
                                   split.info[["train.groups"]],
                                   same=is.set.group[["test"]],
                                   other=!is.set.group[["test"]],
                                   all=rep(TRUE, nrow(id.fold.groups)))
                                 is.set.fold <- list(
                                   test=id.fold.groups[["fold"]] == split.info[["test.fold"]])
                                 is.set.fold[["train"]] <- !is.set.fold[["test"]]
                                 for(set.name in names(is.set.fold)){
                                   is.group <- is.set.group[[set.name]]
                                   is.fold <- is.set.fold[[set.name]]
                                   set(
                                     self$instance$iteration.dt,
                                     i=iteration.i,
                                     j=set.name,
                                     value=list(id.fold.groups[is.group & is.fold, row_id]))
                                 }
                               }
                               self$task_hash = task$hash
                               self$task_nrow = task$nrow
                               invisible(self)
                             },
                             train_set = function(i) {
                               self$instance$iteration.dt$train[[i]]
                             },
                             test_set = function(i) {
                               self$instance$iteration.dt$test[[i]]
                             }
                           ),
                           active = list(
                             is_instantiated = function(rhs) {
                               !is.null(self$instance)
                             },
                             hash = function(rhs) {
                               if (!self$is_instantiated) {
                                 return(NA_character_)
                               }
                               mlr3misc::calculate_hash(list(class(self), self$id, self$param_set$values, self$instance))
                             }
                           )
)
MyResamplingCV = R6::R6Class("MyResamplingCV", inherit = MyResampling,
                             public = list(
                               initialize = function() {
                                 ps = paradox::ps(
                                   folds = paradox::p_int(2L, tags = "required")
                                 )
                                 ps$values = list(folds = 10L)
                                 super$initialize(
                                   id = "mycv",
                                   param_set = ps,
                                   label = "Cross-Validation",
                                   man = "TODO")
                               }
                             ),
                             active = list(
                               iters = function(rhs) {
                                 nrow(mycv$instance$iteration.dt)
                               }
                             ),
                             private = list(
                               .sample = function(ids, ...) {
                                 data.table(
                                   row_id = ids,
                                   fold = sample(seq(0, length(ids)-1) %% as.integer(self$param_set$values$folds) + 1L),
                                   key = "fold"
                                 )
                               },
                               .combine = function(instances) {
                                 rbindlist(instances, use.names = TRUE)
                               },
                               deep_clone = function(name, value) {
                                 switch(name,
                                        "instance" = copy(value),
                                        "param_set" = value$clone(deep = TRUE),
                                        value
                                 )
                               }
                             )
)

mycv <- MyResamplingCV$new()
mycv$param_set$values$folds=2
for(task in task.list){
  mycv$instantiate(task)
}

(reg.learner.list <- list(
  mlr3learners::LearnerRegrCVGlmnet$new(),
  mlr3::LearnerRegrFeatureless$new(),
  LearnerRegrFuser$new()
))

#(reg.learner.list <- list(
#  LearnerRegrFuser$new()
#))

(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))
future::plan("sequential")
reg.dir <- "movingpictures-2024-04-18-benchmark-reg"
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
  walltime = 60*5,#seconds
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
save(bmr, file="movingpictures-2024-04-18-benchmark.RData")
