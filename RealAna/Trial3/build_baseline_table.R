options(stringsAsFactors=FALSE)
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root <- dirname(normalizePath(entry))
dat <- as.data.frame(readRDS(file.path(root,"trial3.rds")))
audit <- read.csv(file.path(root,"output/trial3_tableB_covariates.csv"))
balance <- read.csv(file.path(root,"output/trial3_baseline_balance.csv"))
dictionary <- read.csv(file.path(root,"trial3_dictionary.csv"))
selected <- audit$Variable[audit$Selected]
expected <- c("X_hospital_0m","X_age_over_90_0m","X_sex_0m","X_icu_0m",
              "X_elx_score_0m","X_diabetes_0m","X_chf_0m")
stopifnot(setequal(selected,expected),nrow(dat)==6030L,
  !anyNA(dat$Treatment),setequal(as.character(dat$Treatment),c("Alert","Usual Care")),
  all(vapply(dat[selected],function(x) !anyNA(x),logical(1))))
outcomes <- c("YS_death_event","YS_dialysis_event","YS_aki_progression_event")
stopifnot(all(vapply(dat[outcomes],function(x) !anyNA(x),logical(1))),
  sum(is.na(dat$X_age_0m))==507L,sum(is.na(dat$X_creat_0m))==2L,
  all(dictionary$variable_role[match(selected,dictionary$variable_name)]=="Baseline covariate"))
A <- as.character(dat$Treatment)=="Alert"
stopifnot(sum(A)==3059L,sum(!A)==2971L)
rows <- list()
add <- function(variable,label,level=NULL) {
  binary <- !is.null(level)
  x <- if(binary) as.numeric(as.character(dat[[variable]])==level) else dat[[variable]]
  t <- x[A]; c <- x[!A]
  denom <- sqrt(((length(t)-1)*var(t)+(length(c)-1)*var(c))/(length(t)+length(c)-2))
  smd <- (mean(t)-mean(c))/denom
  feature <- if(binary) paste0(variable,level) else variable
  old <- balance[balance$Feature==feature,]
  if(nrow(old)) stopifnot(nrow(old)==1L,
    abs(old$Intervention_mean-mean(t))<1e-12,
    abs(old$Control_mean-mean(c))<1e-12,
    abs(old$Standardized_mean_difference-smd)<1e-12)
  rows[[length(rows)+1L]] <<- data.frame(
    Variable=variable,Characteristic=label,Level=if(binary) level else "",
    Statistic=if(binary) "n (%)" else "mean (SD)",
    Alert_N=length(t),Usual_care_N=length(c),
    Alert_count=if(binary) sum(t) else NA_real_,Usual_care_count=if(binary) sum(c) else NA_real_,
    Alert_mean=mean(t),Usual_care_mean=mean(c),Alert_SD=sd(t),Usual_care_SD=sd(c),
    Alert_display=if(binary) sprintf("%d (%.2f)",sum(t),100*mean(t)) else sprintf("%.2f (%.2f)",mean(t),sd(t)),
    Usual_care_display=if(binary) sprintf("%d (%.2f)",sum(c),100*mean(c)) else sprintf("%.2f (%.2f)",mean(c),sd(c)),
    Signed_SMD=smd,Absolute_SMD=abs(smd),Missing=0L)
}
for(site in 1:6) add("X_hospital_0m",paste("Hospital",site),as.character(site))
add("X_age_over_90_0m","Age over 90 years","1")
add("X_sex_0m","Sex (code 1)","1")
add("X_icu_0m","ICU status at randomization","1")
add("X_elx_score_0m","Elixhauser comorbidity score")
add("X_diabetes_0m","Diabetes","1")
add("X_chf_0m","Congestive heart failure","1")
result <- do.call(rbind,rows)
stopifnot(nrow(result)==12L,
  sum(result$Alert_count[result$Variable=="X_hospital_0m"])==sum(A),
  sum(result$Usual_care_count[result$Variable=="X_hospital_0m"])==sum(!A))
write.csv(result,file.path(root,"output/trial3_baseline_characteristics.csv"),row.names=FALSE,na="")
text <- c("% Requires booktabs. Sex code 1 is deliberately not assigned a male/female label.",
  r"(\begin{table*}[t])",r"(\centering\small)",
  r"(\caption{Baseline characteristics of participants in the AKI alert trial analysis.})",
  r"(\label{tab:trial3_baseline})",r"(\setlength{\tabcolsep}{6pt})",
  r"(\renewcommand{\arraystretch}{1.08})",r"(\begin{tabular}{@{}lrrr@{}})",r"(\toprule)",
  r"(Characteristic & Alert ($n=3059$) & Usual care ($n=2971$) & Absolute SMD \\)",
  r"(\midrule)",r"(Hospital, $n$ (\%) & & & \\)")
for(i in seq_len(nrow(result))) {
  r <- result[i,]
  label <- if(r$Variable=="X_hospital_0m") paste0("\\quad ",r$Characteristic) else r$Characteristic
  label <- paste0(label,if(r$Variable=="X_hospital_0m") "" else
    if(r$Statistic=="n (%)") r"(, $n$ (\%))" else ", mean (SD)")
  if(i==7L) text <- c(text,r"(\addlinespace[3pt])")
  text <- c(text,sprintf(r"(%s & %s & %s & %.3f \\)",label,r$Alert_display,r$Usual_care_display,r$Absolute_SMD))
}
text <- c(text,r"(\bottomrule\end{tabular})",r"(\par\smallskip)",
  r"(\begin{minipage}{0.92\textwidth}\footnotesize)",
  r"(Values are $n$ (\%) or mean (SD). SMD is the absolute difference in arm-specific means divided by the pooled within-arm SD, $\{[(n_1-1)s_1^2+(n_0-1)s_0^2]/(n_1+n_0-2)\}^{1/2}$.)",
  r"(Hospital SMDs refer to separate category indicators, not a single overall hospital SMD. Hospital numbers are the identifiers in the de-identified dataset.)",
  r"(The available dictionary labels sex only as 0/1; code 1 is reported without assuming a male/female mapping.)",
  r"(All variables displayed and all three 14-day component outcomes were complete in the analysis dataset. Continuous age and baseline serum creatinine had 507 and two missing values, respectively, and were not included in the primary adjustment set; the fully observed age-over-90 indicator was used instead of continuous age.)",
  r"(No imputation or participant exclusion was performed. ICU, intensive care unit.)",
  r"(\end{minipage})",r"(\end{table*})")
writeLines(text,file.path(root,"trial3_baseline_characteristics.tex"))
print(result[,c("Characteristic","Alert_display","Usual_care_display","Absolute_SMD")],row.names=FALSE)
cat("Baseline characteristics table generated; existing treatment-effect results unchanged.\n")
