# Outcome-generation parameters are frozen, not recalibrated by these entry points.
working_sets <- list(W0=integer(),W1=c(3L,4L),W2=c(1L,2L),W3=1:4,W4=2:4,
  W5=3:5,W6=1:5,W7=1:10,W8=1:15)
working_labels <- c(W0="Unadjusted",W1="Shared",W2="D1-pair (exploratory)",
  W3="D1-prog",W4="D2-prog",W5="D3-prog",W6="Union",W7="Union+5",W8="Union+10")
population_means <- setNames(c(0,.5,0,.35,rep(0,2),.5,0,.3,rep(0,6)),paste0("X",1:15))
binary_indices <- c(2L,4L,7L,9L)
legacy_margins <- list(death_epsilon=0,alsfrs_epsilon=2,nfl_epsilon=log(1.25))
default_dgm_params <- function() frozen_parameters$primary

generate_rct_replication <- function(N,effect,params) {
  generate_hce_outcomes(make_population_draws(N,A=NULL,e=params$e),params,effect)
}

read_study_truth <- function(experiment) {
  data.table::fread(file.path(project_dir,"simulation_results",experiment,"simulation_truth.csv"))
}

reproduction_tasks <- function(sample_sizes,effects,reps,bootstrap=FALSE) {
  file <- file.path(project_dir,"simulation_results",if(bootstrap) "bootstrap" else "primary",
    "random_stream_schedule.csv")
  blocks <- data.table::fread(file)
  blocks <- blocks[N%in%sample_sizes & effect%in%effects]
  stopifnot(nrow(blocks)==length(sample_sizes)*length(effects),reps>=1L,
    all(reps<=blocks$n_reps))
  decode <- function(s) as.integer(strsplit(s,",",fixed=TRUE)[[1]])
  tasks <- list(); streams <- list(); bootstrap_streams <- list()
  for(j in seq_len(nrow(blocks))) {
    r <- blocks[j]
    s <- decode(r$trial_rng_start)
    b <- if(bootstrap) decode(r$bootstrap_rng_start) else NULL
    tasks[[j]] <- data.table::data.table(N=r$N,effect=r$effect,rep=seq_len(reps),
      task_id=r$first_task_id+seq_len(reps)-1L)
    for(i in seq_len(reps)) {
      streams[[length(streams)+1L]] <- s
      s <- parallel::nextRNGStream(s)
      if(bootstrap) {
        bootstrap_streams[[length(bootstrap_streams)+1L]] <- b
        b <- parallel::nextRNGStream(b)
      }
    }
  }
  list(tasks=data.table::rbindlist(tasks),streams=streams,bootstrap_streams=bootstrap_streams)
}

bootstrap_trial <- function(task,trial_stream,bootstrap_stream,params,truth,K) {
  use_stream(trial_stream)
  dat <- generate_rct_replication(task$N,task$effect,params)
  tr <- truth[effect==task$effect]
  original <- run_win_estimators_for_all_W(dat,params,c(NB=tr$NB,LWR=tr$LWR,LWO=tr$LWO))
  ctx <- bootstrap_context(dat,params)
  use_stream(bootstrap_stream)
  counts <- rmultinom(K,ctx$N,rep(1/ctx$N,ctx$N))
  placements <- bootstrap_placements(ctx,counts)
  draws <- array(NA_real_,c(K,length(working_sets),3L),
    dimnames=list(NULL,names(working_sets),c("NB","LWR","LWO")))
  failures <- list()
  for(b in seq_len(K)) {
    fit <- bootstrap_refit(ctx,counts[,b],placements,b)
    draws[b,,] <- fit$estimates
    bad <- which(nzchar(fit$reason))
    if(length(bad)) failures[[length(failures)+1L]] <- data.table::data.table(
      task_id=task$task_id,N=task$N,effect=task$effect,rep=task$rep,
      bootstrap_id=b,W=names(working_sets)[bad],reason=fit$reason[bad],detail=fit$detail[bad])
  }
  out <- data.table::rbindlist(lapply(names(working_sets),function(w) {
    values <- matrix(draws[,w,],K,3L,dimnames=list(NULL,c("NB","LWR","LWO")))
    good <- apply(values,1,function(x) all(is.finite(x)))
    bootstrap_se <- if(all(good)) apply(values,2,sd) else setNames(rep(NA_real_,3),colnames(values))
    r <- data.table::copy(original[W==w])
    r[,':='(N=task$N,effect=task$effect,rep=task$rep,task_id=task$task_id,
      bootstrap_K=K,bootstrap_success=sum(good),bootstrap_failed=K-sum(good),
      bootstrap_complete=all(good),se_boot=unname(bootstrap_se[estimand]))]
    r[,':='(ci_low_boot=estimate-qnorm(.975)*se_boot,ci_high_boot=estimate+qnorm(.975)*se_boot,
      p_value_boot=2*pnorm(-abs(estimate/se_boot)),reject_boot=abs(estimate/se_boot)>qnorm(.975))]
    r[,covered_boot:=ci_low_boot<=truth & truth<=ci_high_boot]
    r
  }))
  list(results=out,failures=if(length(failures)) data.table::rbindlist(failures) else
    data.table::data.table(task_id=integer(),N=integer(),effect=character(),rep=integer(),
      bootstrap_id=integer(),W=character(),reason=character(),detail=character()))
}
