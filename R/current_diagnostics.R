evaluate_component_prediction <- function(p,settings,out) {
  if(!requireNamespace("survival",quietly=TRUE)) stop("Package survival is required.")
  set_stream_seed(settings$seed_train)
  train <- make_population_draws(settings$train_n)
  set_stream_seed(settings$seed_validation)
  validation <- make_population_draws(settings$validation_n)
  tr <- generate_hce_outcomes(train,p)
  va <- generate_hce_outcomes(validation,p)
  at <- tr$death_event==0; av <- va$death_event==0
  zeta <- staggered_predictors(validation,p)
  full <- data.table::data.table(metric=c("mortality_S","Harrell_C","AUC_S","D2_observed_R2",
    "D3_log_R2","D2_mean_survivors","D3_median_survivors"),value=c(
      mean(va$death_event),survival_metrics(va,zeta[,1]),
      prediction_r2(va$icu_days[av],expected_icu(zeta[av,2],p$sigma2,p$S)),
      prediction_r2(va$log_biomarker[av],zeta[av,3]),
      mean(va$icu_days[av]),median(va$biomarker[av])))
  by_w <- data.table::rbindlist(lapply(names(working_sets),function(w) {
    cols <- names(population_means)[working_sets[[w]]]
    if(!length(cols)) return(data.table::data.table(W=w,p=0L,Harrell_C=.5,AUC_S=.5,D2_observed_R2=0,D3_log_R2=0))
    fit <- survival::coxph(reformulate(cols,response="survival::Surv(death_time,death_event)"),
      data=as.data.frame(tr),ties="breslow",singular.ok=FALSE)
    sm <- survival_metrics(va,as.numeric(predict(fit,newdata=as.data.frame(va),type="lp")))
    Xt <- cbind(1,as.matrix(tr[at,..cols])); Xv <- cbind(1,as.matrix(va[av,..cols]))
    b2 <- lm.fit(Xt,tr$icu_days[at])$coefficients
    b3 <- lm.fit(Xt,tr$log_biomarker[at])$coefficients
    stopifnot(!anyNA(c(b2,b3)))
    data.table::data.table(W=w,p=length(cols),Harrell_C=sm[1],AUC_S=sm[2],
      D2_observed_R2=prediction_r2(va$icu_days[av],as.vector(Xv%*%b2)),
      D3_log_R2=prediction_r2(va$log_biomarker[av],as.vector(Xv%*%b3)))
  }))
  save_csv(full,"component_calibration_full",out)
  save_csv(by_w,"component_calibration_by_W",out)
  invisible(list(full=full,by_W=by_w))
}
