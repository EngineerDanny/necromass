library(data.table)

bfcpt <- fread("bacteria_fungi_conservative_power_transformed.csv")
meta.dt <- fread("bacteria_conservative_r_same.csv")[
  order(Necrobag_ID)
, c("Necrobag_ID", "Sample_ID", "Domain", "Habitat", "Melanization", 
"Incubation_time", "Plot", "Comp")]
same_other_cv <- mlr3resampling::ResamplingSameOtherCV$new()
meta.dt[, Samples := ifelse(Habitat=="soil", "Habitat=soil", paste0("Melanization=",Melanization))][, table(Samples)]
task.list <- list()
for(out.i in 1:ncol(bfcpt)){
  task.dt <- data.table(bfcpt, Samples=meta.dt$Samples)
  task_id <- names(bfcpt)[[out.i]]
  task.list[[task_id]] <- mlr3::TaskRegr$new(
    task_id, task.dt, target=task_id
  )$set_col_roles("Samples",c("group","stratum"))
}

(reg.learner.list <- list(
  mlr3learners::LearnerRegrCVGlmnet$new(),
  mlr3::LearnerRegrFeatureless$new()))

(reg.bench.grid <- mlr3::benchmark_grid(
  task.list,
  reg.learner.list,
  same_other_cv))

(debug.grid <- mlr3::benchmark_grid(
  task.list["Kaistia"],
  reg.learner.list,
  same_other_cv))
future::plan("sequential")
debug.result <- mlr3::benchmark(debug.grid)

if(F){
  future::plan(
    future.batchtools::batchtools_slurm,
    template="~/slurm-future.tmpl",
    resources=list(
      walltime=60*60,#seconds as defined https://mllg.github.io/batchtools/reference/submitJobs
      memory=1000,
      ncpus=1,
      ntasks=1,
      chunks.as.arrayjobs=TRUE))
  ## mlr3::benchmark uses future_map defined in https://github.com/mlr-org/mlr3/blob/545873fbdd7a55be9ca74ef5b264ddaca0de2f8e/R/helper_exec.R#L25 which uses future.apply::future_mapply
  future.apply::future_mapply(function(i)Sys.sleep(60), 1:100)
  (reg.bench.result <- mlr3::benchmark(
    reg.bench.grid, store_models = TRUE))
}

unlink("danny-data-reg", recursive=TRUE)
reg = batchtools::makeExperimentRegistry(
  file.dir = "danny-data-reg",
  seed = 1,
  packages = "mlr3verse"
)
mlr3batchmark::batchmark(
  reg.bench.grid, store_models = TRUE, reg=reg)
job.table <- batchtools::getJobTable(reg=reg)
chunks <- data.frame(job.table, chunk=1)
batchtools::submitJobs(chunks, resources=list(
  walltime = 60*60,#seconds
  memory = 2000,#megabytes per cpu
  ncpus=1,  #>1 for multicore/parallel jobs.
  ntasks=1, #>1 for MPI jobs.
  chunks.as.arrayjobs=TRUE), reg=reg)
batchtools::getStatus(reg=reg)
## Status for 1944 jobs at 2023-12-22 07:10:54:
##   Submitted    : 1944 (100.0%)
##   -- Queued    :    0 (  0.0%)
##   -- Started   : 1944 (100.0%)
##   ---- Running :    0 (  0.0%)
##   ---- Done    : 1873 ( 96.3%)
##   ---- Error   :   71 (  3.7%)
##   ---- Expired :    0 (  0.0%)

table(jobs.after$error)
jobs.after <- batchtools::getJobTable(reg=reg)
jobs.after[!is.na(error), .(error, task_id=sapply(prob.pars, "[[", "task_id"))][25:26]
                                                                                 
## 1:                                            Error in approx(lambda, seq(lambda), sfrac) : \n  need at least two non-NA values to interpolate
## 2: Error in elnet(xd, is.sparse, y, weights, offset, type.gaussian, alpha,  : \n  y is constant; gaussian glmnet fails at standardization step
##    task_id
## 1: Kaistia
## 2: Kaistia


ids <- jobs.after[is.na(error), job.id]
bmr = mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
score.dt <- mlr3resampling::score(bmr)
save(bmr, score.dt, file="data-danny-same-other-cv.RData")
