library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(fuser)
library(R6)
library(paradox)

data.list <- list()
cons.csv.vec <- Sys.glob("/projects/genomic-ml/necromass/data-2023-12-22/*.csv")
for(cons.csv in cons.csv.vec){
  king.dt <- nc::capture_first_vec(cons.csv, "/Dec22_", kingdom="[^_]+")
  data.list[[king.dt$kingdom]] <- fread(cons.csv)[order(Necrobag_ID)]
}
str(data.list)
names.list <- lapply(data.list, names)
sapply(data.list, dim)
names.meta <- do.call(intersect, unname(names.list))
cat(names.meta, sep=", ")
do.call(identical, unname(lapply(data.list, "[[", "Necrobag_ID")))
meta.long.list <- list()
for(kingdom in names(data.list)){
  king.full.dt <- data.list[[kingdom]]
  king.meta.dt <- king.full.dt[,names.meta,with=FALSE]
  meta.long.list[[kingdom]] <- data.table(kingdom, king.meta.dt)
}
meta.long <- rbindlist(meta.long.list)
meta.longer <- melt(meta.long, id=c("kingdom", "Necrobag_ID"))
meta.wide <- dcast(meta.longer, Necrobag_ID + variable ~ kingdom)
(question <- meta.wide[bacteria != fungi])

meta.names <- c(
  "Necrobag_ID", "sample_ID", "Domain", "Habitat", "Melanization", 
  "Incubation_time", "Plot", "Comp", "OTU_ID")
necro.dt.list <- list()
for(data.name in names(data.list)){
  meta.and.data <- data.list[[data.name]]
  is.meta <- names(meta.and.data) %in% meta.names
  meta.dt <- meta.and.data[, is.meta, with=FALSE]
  necro.dt.list[[data.name]] <- meta.and.data[, !is.meta, with=FALSE]
}
(necro.dt <- do.call(data.table, necro.dt.list))

necro.tall <- melt(necro.dt, measure.vars=names(necro.dt))

if (F){
library(ggplot2)
ggplot()+
  geom_histogram(aes(
    value),
    data=necro.tall)+
  facet_wrap(~variable)

ggplot()+
  geom_histogram(aes(
    log10(value+1)),
    data=necro.tall)+
  facet_wrap(~variable)
}

log.necro.dt <- log10(necro.dt+1)

meta.dt[, Samples := ifelse(Habitat=="soil", "Habitat=soil", paste0("Melanization=",Melanization))][, table(Samples)]
task.list <- list()
for(out.i in 1:ncol(log.necro.dt)){
  task.dt <- data.table(log.necro.dt, Samples=meta.dt$Samples)
  task_id <- names(log.necro.dt)[[out.i]]
  reg.task <- mlr3::TaskRegr$new(
    task_id, task.dt, target=task_id
  )
  
  reg.task$col_roles$subset <- "Samples"
  #reg.task$col_roles$group <- "Samples"
  reg.task$col_roles$stratum <- "Samples"
  reg.task$col_roles$feature <- setdiff(names(task.dt), c(task_id, "Samples"))
  
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
                                print("got here")
                                print(group_ind)
                                print("got here 2")
                                print(task$strata)
                                
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

mycv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
mycv$param_set$values$folds=2
for(task in task.list){
  mycv$instantiate(task)
}
mycv$instance$iteration.dt

(reg.learner.list <- list(
  mlr3learners::LearnerRegrCVGlmnet$new(),
  mlr3::LearnerRegrFeatureless$new(),
 LearnerRegrFuser$new()
))

(debug.grid <- mlr3::benchmark_grid(
  task.list["bacteria.Bradyrhizobium"],
  reg.learner.list,
  mycv))
future::plan("sequential")
debug.result <- mlr3::benchmark(debug.grid)
debug.result

(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

reg.dir <- "data-2024-04-23-benchmark-reg-2"
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
save(bmr, file="data-2024-04-23-benchmark-2.RData")

# search for predictions with NAs
tab = as.data.table(bmr)
tab[map_lgl(tab$prediction, function(pred) any(is.na(pred$response)))]
tab$prediction
