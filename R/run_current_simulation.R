options(stringsAsFactors=FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R/current_study.R"))
s <- load_current_study(root)
args <- commandArgs(TRUE)
arg <- function(key,default) {
  a <- grep(paste0("^--",key,"="),args,value=TRUE)
  if(length(a)) sub(paste0("^--",key,"="),"",a[1]) else default
}
stage <- match.arg(arg("stage","simulate"),c("simulate","bootstrap","population"))
experiment <- match.arg(arg("experiment","primary"),"primary")
mode <- match.arg(arg("mode","smoke"),c("smoke","full"))
reps <- as.integer(arg("reps",if(mode=="full") "10000" else "5"))
samples <- as.integer(strsplit(arg("samples",if(stage=="bootstrap") "200" else
  if(mode=="full") "200,600,800,1200" else "200"),",",fixed=TRUE)[[1]])
effects <- strsplit(arg("effects",if(stage=="bootstrap") "null" else "null,alternative"),",",fixed=TRUE)[[1]]
cores <- as.integer(arg("cores","1"))
stopifnot(reps>=1L,reps<=10000L,!anyNA(samples),!anyDuplicated(samples),
  all(samples%in%c(200L,600L,800L,1200L)),all(effects%in%c("null","alternative")),
  !anyDuplicated(effects),cores>=1L)
out <- arg("output",file.path(root,"simulation_runs",paste(experiment,stage,mode,sep="_")))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out)
published <- normalizePath(file.path(root,"simulation_results"))
if(out==published || startsWith(out,paste0(published,"/")))
  stop("Published simulation_results are read-only inputs; choose a different output directory.")
p <- s$frozen_parameters[[experiment]]
truth <- s$read_study_truth(experiment)
oracle <- data.table::fread(file.path(root,"simulation_results",experiment,"oracle_if_pvr_by_W.csv"))
if(stage=="population") {
  if(!requireNamespace("survival",quietly=TRUE)) stop("Package survival is required for population diagnostics.")
  settings <- s$frozen_parameters[[paste0(experiment,"_population_settings")]]
  if(mode=="smoke") {
    settings$train_n <- settings$validation_n <- 2000L
    settings$truth_pairs <- 10000L
    settings$truth_chunk <- 2000L
    settings$reference_n <- settings$oracle_train_n <- settings$oracle_validation_n <- 2000L
  }
  s$evaluate_component_prediction(p,settings,out)
  tr <- s$calculate_population_win_truth(p,settings,out)
  s$calculate_oracle_if_pvr(p,tr,settings,out)
} else {
  if(stage=="bootstrap" && (experiment!="primary"||!identical(samples,200L)))
    stop("The retained bootstrap comparison is primary DGP at N=200 only.")
  plan <- s$reproduction_tasks(samples,effects,reps,bootstrap=stage=="bootstrap")
  checkpoint_dir <- file.path(out,"checkpoints"); dir.create(checkpoint_dir,showWarnings=FALSE)
  K <- as.integer(arg("K",if(mode=="full") "2000" else "20"))
  settings <- list(stage=stage,experiment=experiment,reps=reps,samples=samples,effects=effects,K=K,
    params=p,truth=truth,source_streams=plan$streams,bootstrap_streams=plan$bootstrap_streams,
    implementation_hashes=setNames(unname(tools::md5sum(s$source_files)),basename(s$source_files)))
  config_file <- file.path(out,"run_configuration.rds")
  if(file.exists(config_file)) stopifnot(identical(settings,readRDS(config_file))) else saveRDS(settings,config_file)
  ids <- split(seq_len(nrow(plan$tasks)),(seq_len(nrow(plan$tasks))-1L)%/%100L)
  chunk_files <- character(length(ids))
  for(j in seq_along(ids)) {
    idx <- ids[[j]]; f <- file.path(checkpoint_dir,sprintf("chunk_%04d.csv.gz",j))
    chunk_files[j] <- f
    if(file.exists(f)&&file.exists(paste0(f,".md5"))) {
      stopifnot(identical(unname(tools::md5sum(f)),readLines(paste0(f,".md5"))))
      next
    }
    if(stage=="simulate") {
      raw <- s$run_trial_tasks(plan$tasks[idx],plan$streams[idx],p,truth,cores=cores)
    } else {
      stopifnot(K>=2L)
      results <- lapply(idx,function(i) s$bootstrap_trial(plan$tasks[i],plan$streams[[i]],
        plan$bootstrap_streams[[i]],p,truth,K))
      raw <- data.table::rbindlist(lapply(results,"[[","results"))
      data.table::fwrite(data.table::rbindlist(lapply(results,"[[","failures")),
        file.path(checkpoint_dir,sprintf("failures_%04d.csv",j)))
    }
    data.table::fwrite(raw,paste0(f,".tmp.gz"),compress="gzip")
    stopifnot(file.rename(paste0(f,".tmp.gz"),f))
    writeLines(unname(tools::md5sum(f)),paste0(f,".md5"))
    message("Completed chunk ",j," / ",length(ids))
  }
  raw <- data.table::rbindlist(lapply(chunk_files,data.table::fread))
  stopifnot(nrow(raw)==nrow(plan$tasks)*27L,!anyDuplicated(raw[,.(task_id,W,estimand)]))
  if(stage=="simulate") {
    summary <- s$summarize_operating_characteristics(raw,oracle,out)
  } else {
    summary <- raw[,{
      paired <- bootstrap_complete & is.finite(se_boot)
      list(n_attempted=.N,n_complete=sum(paired),ESD=sd(estimate),ASE_IF=mean(se),
        ASE_boot=mean(se_boot[paired]),CP_IF=mean(covered),CP_boot=mean(covered_boot[paired]),
        reject_IF=mean(reject),reject_boot=mean(reject_boot[paired]),
        failed_draws=sum(bootstrap_failed),truth=truth[1])
    },by=.(N,effect,W,estimand)]
    summary[,':='(ASE_ESD_IF=ASE_IF/ESD,ASE_ESD_boot=ASE_boot/ESD)]
    ci <- s$wilson_ci(summary$reject_boot,summary$n_complete)
    summary[,':='(reject_boot_mc_low=ci[,1],reject_boot_mc_high=ci[,2])]
    data.table::fwrite(summary,file.path(out,"bootstrap_summary.csv"))
  }
  print(summary[,.(N,effect,W,estimand,n_attempted)])
}
writeLines(c(paste("Experiment:",experiment),paste("Stage:",stage),paste("Mode:",mode)),
  file.path(out,"completion.txt"))
