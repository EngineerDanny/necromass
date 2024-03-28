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

log.necro.dt <- log10(necro.dt+1)

meta.dt[, Samples := ifelse(Habitat=="soil", "Habitat=soil", paste0("Melanization=",Melanization))][, table(Samples)]
task.list <- list()
for(out.i in 1:ncol(log.necro.dt)){
  task.dt <- data.table(log.necro.dt, Samples=meta.dt$Samples)
  task_id <- names(log.necro.dt)[[out.i]]
  task.list[[task_id]] <- mlr3::TaskRegr$new(
    task_id, task.dt, target=task_id
  )$set_col_roles("Samples",c("group","stratum"))
}


LearnerRegrFuser <- R6Class("LearnerRegrFuser",
  inherit = LearnerRegr,
  public = list( 
    initialize = function() {
       ps = ps(
          lambda = p_dbl(0, default = 1, tags = "train"),
          gamma = p_dbl(0, default = 0, tags = "train"),
          tol = p_dbl(lower = 0, upper = 1, default = 0.01, tags = "train"),
          num.it = p_int(default = 500, tags = "train"),
          intercept = p_lgl(default = TRUE, tags = "train"),
          scaling = p_lgl(default = FALSE, tags = "train")
       )
       ps$values = list(lambda = 0.01, gamma = 0.01, 
                        tol = 0.01, num.it = 5000,
                        intercept = TRUE, scaling = FALSE)
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
     # Check if k = 1, then default to glmnet
      if (k == 1) {
        glmnet_model <- glmnet::glmnet(X_train, y_train)
        self$model = list(glmnet_model = glmnet_model,
                          formula = task$formula(),
                          data = task$data(),
                          pv = pv,
                          groups = group_ind)
      } else {
        # Use fuser::fusedLassoProximal
        beta.estimate <- fuser::fusedLassoProximal(X, y, group_ind, 
                                           lambda = pv$lambda, 
                                           G = matrix(1, k, k), 
                                           gamma = pv$gamma,
                                           tol = pv$tol, 
                                           num.it = pv$num.it,
                                           intercept = FALSE,
                                           scaling = pv$scaling) 
        self$model = list(beta = beta.estimate, 
                          formula = task$formula(),
                          data = task$data(),
                          pv = pv,
                          groups = group_ind)
      }

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
            # Extract the coefficients from the learner
            beta <- self$model$beta
            y.predict <- rep(NA, nrow(X))
            for (k.i in 1:k) {
                group_rows <- which(group_ind == k.i)
                X.group <- X[group_rows, , drop = FALSE]
                y.predict[group_rows] <- as.numeric(X.group %*% beta[, k.i])
            }
        }, error = function(e) {
            #warning("Error occurred with beta coefficient prediction. Falling back to glmnet.")
            glmnet_model <- self$model$glmnet_model
            y.predict <- predict(glmnet_model, newx = X_test)
        })
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

(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  mycv))

(debug.grid <- mlr3::benchmark_grid(
  task.list["bacteria.Kaistia"],
  reg.learner.list,
  mycv))
future::plan("sequential")
debug.result <- mlr3::benchmark(debug.grid)

reg.dir <- "data-2024-03-28-benchmark-reg"
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
  walltime = 60*60,#seconds
  memory = 4000,#megabytes per cpu
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
save(bmr, file="data-2024-03-28-benchmark.RData")
