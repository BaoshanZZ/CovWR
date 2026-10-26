options(stringsAsFactors=FALSE)
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R/current_study.R"))
s <- load_current_study(root)
checks <- list()
for(kind in "primary") {
  p <- s$frozen_parameters[[kind]]
  tr <- s$read_study_truth(kind)
  plan <- s$reproduction_tasks(c(200L,600L,800L,1200L),c("null","alternative"),1L)
  for(i in seq_len(nrow(plan$tasks))) {
    task <- plan$tasks[i]
    s$use_stream(plan$streams[[i]])
    dat <- s$generate_rct_replication(task$N,task$effect,p)
    pair <- s$compute_pairwise_quantities(dat,p$e,p)
    reverse <- s$revised_pair_matrices(dat[A==0],dat[A==1],p)
    truth <- tr[effect==task$effect]
    rows <- s$run_win_estimators_for_all_W(dat,p,c(NB=truth$NB,LWR=truth$LWR,LWO=truth$LWO))
    stopifnot(nrow(rows)==27L,all(rows$success),s$check_kernel_partition(pair$matrices),
      all(pair$matrices$win==t(reverse$loss)),pair$allocation_prob==.5,
      all(dat$icu_days>=0 & dat$icu_days<=floor(dat$death_time)),
      all(is.na(dat$biomarker)==(dat$death_event==1)))
    X <- as.matrix(dat[,paste0("X",1:5),with=FALSE])
    fit <- s$fit_projection_qr(X,dat$A,pair$psi,p$e)
    D <- sweep(X,2,colMeans(X),"-")*(dat$A-p$e)
    gram_beta <- solve(crossprod(D),crossprod(D,pair$psi))
    stopifnot(max(abs(fit$beta-gram_beta))<1e-10)
    tau <- c(pair$tau_w,pair$tau_l)-colMeans(fit$adjustment)
    measure <- s$calc_measures(tau[[1]],tau[[2]],p)
    observed <- rows[W=="W6"]
    stopifnot(max(abs(observed$estimate-unname(measure[observed$estimand])))<1e-12)
    for(m in c("NB","LWR","LWO")) {
      g <- s$measure_gradient(m,tau[[1]],tau[[2]],p)
      psi <- as.vector(fit$residual%*%g)
      expected <- sqrt(mean((psi-mean(psi))^2)/nrow(dat))
      stopifnot(abs(observed[estimand==m,se]-expected)<1e-12)
      g0 <- s$measure_gradient(m,pair$tau_w,pair$tau_l,p)
      stopifnot(var(as.vector(fit$residual%*%g0))<=var(as.vector(pair$psi%*%g0))+1e-10)
    }
    checks[[length(checks)+1L]] <- data.table::data.table(experiment=kind,N=task$N,effect=task$effect,
      kernel_partition=TRUE,QR_Gram_equivalent=TRUE,exact_transformations=TRUE,
      covariance_variance_equivalent=TRUE,common_gradient_projection_no_loss=TRUE)
  }
}
plan <- s$reproduction_tasks(200L,"null",1L,bootstrap=TRUE)
p <- s$frozen_parameters$primary
s$use_stream(plan$streams[[1]])
draws <- s$make_population_draws(200L,A=NULL,e=.5)
dat <- s$generate_hce_outcomes(draws,p)
ctx <- s$bootstrap_context(dat,p)
s$use_stream(plan$bootstrap_streams[[1]])
counts <- rmultinom(3L,200L,rep(1/200,200))
placements <- s$bootstrap_placements(ctx,counts)
for(b in 1:3) {
  fit <- s$bootstrap_refit(ctx,counts[,b],placements,b)
  expanded <- data.table::copy(dat[rep(seq_len(nrow(dat)),counts[,b])])
  expanded[,id:=seq_len(.N)]
  rows <- s$run_win_estimators_for_all_W(expanded,p,c(NB=0,LWR=0,LWO=0))
  for(w in names(s$working_sets)) {
    r <- rows[W==w]
    stopifnot(max(abs(r$estimate-fit$estimates[w,r$estimand]))<1e-10)
  }
}
boot <- s$bootstrap_trial(plan$tasks[1],plan$streams[[1]],plan$bootstrap_streams[[1]],p,
  s$read_study_truth("primary"),20L)
stopifnot(nrow(boot$results)==27L,all(boot$results$bootstrap_complete),
  all(is.finite(boot$results$se_boot)),all(boot$results$se_boot>0),nrow(boot$failures)==0L)
print(data.table::rbindlist(checks))
cat("Validated all eight primary N/effect cells, all 27 method-measure rows, and pooled bootstrap multiplicities.\n")
