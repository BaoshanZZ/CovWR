options(stringsAsFactors = FALSE, scipen = 999)

local_project_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    candidate <- sub("^--file=", "", file_arg)
    if (candidate != "-" && file.exists(candidate)) {
      return(dirname(dirname(normalizePath(candidate))))
    }
  }
  if (dir.exists("WRCov/R")) return(normalizePath("WRCov"))
  if (dir.exists("R")) return(normalizePath("."))
  stop("Run from the WRCov project directory or its parent.", call. = FALSE)
}

parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 1L) sub(paste0("^--", name, "="), "", hit) else default
}

wilson_interval <- function(p, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  data.frame(lower = pmax(0, center - half_width), upper = pmin(1, center + half_width))
}

escape_latex <- function(x) {
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}

format_fixed <- function(x, digits) formatC(x, format = "f", digits = digits)
format_percent_ci <- function(p, low, high) {
  sprintf(
    "\\shortstack[c]{%.2f\\\\(%.2f, %.2f)}",
    100 * p, 100 * low, 100 * high
  )
}
format_method <- function(method) {
  ifelse(method == "UNADJ", "$\\bar\\eta_m$", "$\\hat\\eta_m^{\\mathrm{lin}}$")
}
format_pvr <- function(method, pvr_percent) {
  ifelse(method == "UNADJ", "--", format_fixed(pvr_percent, 1))
}
format_bias <- function(bias) {
  scaled <- 1000 * bias
  scaled[abs(scaled) < 0.005] <- 0
  format_fixed(scaled, 2)
}
format_configuration <- function(label) {
  parts <- strsplit(label, ": ", fixed = TRUE)[[1L]]
  if (length(parts) == 2L) {
    return(paste0(
      "\\multirow[c]{2}{*}{\\shortstack[c]{",
      escape_latex(parts[[1L]]), "\\\\", escape_latex(parts[[2L]]), "}}"
    ))
  }
  paste0("\\multirow[c]{2}{*}{", escape_latex(label), "}")
}

project_dir <- local_project_dir()
input_dir <- parse_arg(
  "input-dir",
  file.path(project_dir, "simulation_results")
)
output_dir <- parse_arg("output-dir", file.path(input_dir, "main_report"))
requested_estimands <- trimws(strsplit(parse_arg("estimands", "LWR"), ",", fixed = TRUE)[[1L]])
valid_estimands <- c("NB", "LWR", "LWO")
unknown_estimands <- setdiff(requested_estimands, valid_estimands)
if (length(unknown_estimands) > 0L) {
  stop("Unknown estimand(s): ", paste(unknown_estimands, collapse = ", "), call. = FALSE)
}
requested_estimands <- valid_estimands[valid_estimands %in% requested_estimands]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
table_suffix <- ""
label_suffix <- ""
status_comment <- "% FINAL: generated from the fixed-e, centered no-intercept LADJ simulation output."

summary_path <- file.path(input_dir, "simulation_summary.csv")
grid_path <- file.path(input_dir, "scenario_grid.csv")
if (!file.exists(summary_path) || !file.exists(grid_path)) {
  stop("Missing simulation_summary.csv or scenario_grid.csv in input directory.", call. = FALSE)
}

summary <- utils::read.csv(summary_path, check.names = FALSE)
grid <- unique(utils::read.csv(grid_path, check.names = FALSE)[c("scenario_id", "scenario_change")])
grid_full <- unique(utils::read.csv(grid_path, check.names = FALSE))
reference_n <- unique(grid_full$N[grid_full$scenario_change == "baseline"])
comparator_n <- unique(grid_full$N[grepl("^N_", grid_full$scenario_change)])
status_note <- paste0(
  "The reference configuration used $N=", reference_n[[1L]],
  "$ and the smaller-sample configuration used $N=", comparator_n[[1L]],
  "$. Alternative truths were evaluated using $10^7$ independently generated treatment--control pairs."
)
dat <- merge(summary, grid, by = "scenario_id", all.x = TRUE, sort = FALSE)
dat <- dat[
  dat$estimand %in% requested_estimands & dat$method %in% c("UNADJ", "LADJ"),
  ,
  drop = FALSE
]
expected_rows <- 36L * length(requested_estimands)
if (nrow(dat) != expected_rows) {
  stop(
    "Expected ", expected_rows, " rows (9 configurations x 2 effects x 2 methods x ",
    length(requested_estimands), " estimands); found ", nrow(dat), call. = FALSE
  )
}

scenario_map <- data.frame(
  scenario_change = c(
    "baseline", "N_200", "N_400", "N_500", "allocation_1_to_2", "weak_prognosis",
    "strong_prognosis", "heterogeneous_effect", "endpoint_specific_prognosis",
    "two_noise_covariates", "omit_X3"
  ),
  scenario_order = c(1L, 2L, 2L, 2L, 3:9),
  scenario_label = c(
    "S1: Reference", "S2: $N=200$", "S2: $N=400$", "S2: $N=500$",
    "S3: Allocation 1:2", "S4: Weak prognosis", "S5: Strong prognosis",
    "S6: Heterogeneity", "S7: Endpoint-specific", "S8: Add noise",
    "S9: Omit key"
  ),
  stringsAsFactors = FALSE
)
dat <- merge(dat, scenario_map, by = "scenario_change", all.x = TRUE, sort = FALSE)
if (anyNA(dat$scenario_order)) stop("Unmapped scenario_change value.", call. = FALSE)
dat$method_order <- match(dat$method, c("UNADJ", "LADJ"))
dat$estimand_order <- match(dat$estimand, requested_estimands)
dat <- dat[
  order(dat$estimand_order, dat$effect, dat$scenario_order, dat$method_order),
  ,
  drop = FALSE
]

dat$ASE_ESD <- dat$ASE / dat$ESD
dat$PVR_percent <- 100 * dat$PVR
cp_ci <- wilson_interval(dat$CP, dat$n_rep)
reject_ci <- wilson_interval(dat$reject_rate, dat$n_rep)
dat$CP_low <- cp_ci$lower
dat$CP_high <- cp_ci$upper
dat$reject_low <- reject_ci$lower
dat$reject_high <- reject_ci$upper

make_null_table <- function(x, estimand_key, estimand_title) {
  rows <- character(nrow(x))
  for (i in seq_len(nrow(x))) {
    configuration <- if (i %% 2L == 1L) {
      format_configuration(x$scenario_label[[i]])
    } else {
      ""
    }
    rows[[i]] <- paste(
      configuration,
      format_method(x$method[[i]]),
      format_bias(x$bias[[i]]),
      format_fixed(x$ESD[[i]], 4),
      format_fixed(x$ASE[[i]], 4),
      format_fixed(x$ASE_ESD[[i]], 3),
      format_pvr(x$method[[i]], x$PVR_percent[[i]]),
      format_percent_ci(x$CP[[i]], x$CP_low[[i]], x$CP_high[[i]]),
      format_percent_ci(x$reject_rate[[i]], x$reject_low[[i]], x$reject_high[[i]]),
      sep = " & "
    )
    rows[[i]] <- paste0(rows[[i]], " \\\\")
    if (i %% 2L == 0L && i < nrow(x)) {
      rows[[i]] <- paste0(
        rows[[i]],
        "\n\\addlinespace[1pt]",
        "\n\\cmidrule(lr){1-9}",
        "\n\\addlinespace[1pt]"
      )
    }
  }
  c(
    status_comment,
    "% Requires \\usepackage{multirow} in the manuscript preamble.",
    "\\begin{table*}[t]",
    "\\centering",
    "\\scriptsize",
    paste0("\\caption{Null operating characteristics for the ", estimand_title, ".}"),
    paste0("\\label{tab:sim_", tolower(estimand_key), "_null", label_suffix, "}"),
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{@{}ccccccccc@{}}",
    "\\toprule",
    "\\multirow[c]{2}{*}{Configuration} & \\multirow[c]{2}{*}{Method} & \\multicolumn{1}{c}{Point Est.} & \\multicolumn{4}{c}{Precision and efficiency} & \\multicolumn{2}{c}{Inference, \\%} \\\\",
    "\\cmidrule(lr){3-3} \\cmidrule(lr){4-7} \\cmidrule(lr){8-9}",
    " & & Bias ($\\times10^{-3}$) & ESD & ASE & ASE/ESD & PVR, \\% & \\shortstack[c]{CProb\\\\(95\\% MC CI)} & \\shortstack[c]{Type I error\\\\(95\\% MC CI)} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{minipage}{0.98\\textwidth}",
    "\\footnotesize $\\bar\\eta_m$, unadjusted estimator; $\\hat\\eta_m^{\\mathrm{lin}}$, linearly adjusted estimator; Bias is displayed in units of $10^{-3}$; ESD, empirical standard deviation; ASE, average analytical standard error; PVR, percent variance reduction relative to the unadjusted estimator computed from ESD$^2$; CProb, empirical coverage probability of the nominal 95\\% Wald interval; MC CI, Wilson Monte Carlo confidence interval, shown in parentheses below each percentage. A dash denotes the unadjusted variance reference. Results use 20,000 Monte Carlo replications.",
    "Scenarios S1--S9 are defined in Table~\\ref{tab:sim_scenarios}.",
    status_note,
    "\\end{minipage}",
    "\\end{table*}"
  )
}

make_alternative_table <- function(x, estimand_key, estimand_title) {
  rows <- character(nrow(x))
  for (i in seq_len(nrow(x))) {
    configuration <- if (i %% 2L == 1L) {
      format_configuration(x$scenario_label[[i]])
    } else {
      ""
    }
    rows[[i]] <- paste(
      configuration,
      format_method(x$method[[i]]),
      format_bias(x$bias[[i]]),
      format_fixed(x$ESD[[i]], 4),
      format_fixed(x$ASE[[i]], 4),
      format_fixed(x$ASE_ESD[[i]], 3),
      format_pvr(x$method[[i]], x$PVR_percent[[i]]),
      format_percent_ci(x$CP[[i]], x$CP_low[[i]], x$CP_high[[i]]),
      format_percent_ci(x$reject_rate[[i]], x$reject_low[[i]], x$reject_high[[i]]),
      sep = " & "
    )
    rows[[i]] <- paste0(rows[[i]], " \\\\")
    if (i %% 2L == 0L && i < nrow(x)) {
      rows[[i]] <- paste0(
        rows[[i]],
        "\n\\addlinespace[1pt]",
        "\n\\cmidrule(lr){1-9}",
        "\n\\addlinespace[1pt]"
      )
    }
  }
  c(
    status_comment,
    "% Requires \\usepackage{multirow} in the manuscript preamble.",
    "\\begin{table*}[t]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{@{}ccccccccc@{}}",
    "\\toprule",
    "\\multirow[c]{2}{*}{Configuration} & \\multirow[c]{2}{*}{Method} & \\multicolumn{1}{c}{Point Est.} & \\multicolumn{4}{c}{Precision and efficiency} & \\multicolumn{2}{c}{Inference, \\%} \\\\",
    "\\cmidrule(lr){3-3} \\cmidrule(lr){4-7} \\cmidrule(lr){8-9}",
    " & & Bias ($\\times10^{-3}$) & ESD & ASE & ASE/ESD & PVR, \\% & \\shortstack[c]{CProb\\\\(95\\% MC CI)} & \\shortstack[c]{Power\\\\(95\\% MC CI)} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    paste0("\\caption{Alternative operating characteristics for the ", estimand_title, ".}"),
    paste0("\\label{tab:sim_", tolower(estimand_key), "_alternative", label_suffix, "}"),
    "\\begin{minipage}{0.98\\textwidth}",
    "\\footnotesize $\\bar\\eta_m$, unadjusted estimator; $\\hat\\eta_m^{\\mathrm{lin}}$, linearly adjusted estimator; Bias is displayed in units of $10^{-3}$; ESD, empirical standard deviation; ASE, average analytical standard error; PVR, percent variance reduction relative to the unadjusted estimator computed from ESD$^2$; CProb, empirical coverage probability; MC CI, Wilson Monte Carlo confidence interval, shown in parentheses below each percentage. A dash denotes the unadjusted variance reference. Results use 20,000 Monte Carlo replications.",
    "Scenarios S1--S9 are defined in Table~\\ref{tab:sim_scenarios}.",
    status_note,
    "\\end{minipage}",
    "\\end{table*}"
  )
}

estimand_titles <- c(
  NB = "net benefit (NB)",
  LWR = "log win ratio (LWR)",
  LWO = "log win odds (LWO)"
)
written_paths <- character()
combined_lines <- character()

for (estimand_key in requested_estimands) {
  estimand_dat <- dat[dat$estimand == estimand_key, , drop = FALSE]
  estimand_title <- unname(estimand_titles[[estimand_key]])
  null_lines <- make_null_table(
    estimand_dat[estimand_dat$effect == "null", , drop = FALSE],
    estimand_key,
    estimand_title
  )
  alternative_lines <- make_alternative_table(
    estimand_dat[estimand_dat$effect == "alternative", , drop = FALSE],
    estimand_key,
    estimand_title
  )

  file_stem <- tolower(estimand_key)
  null_path <- file.path(
    output_dir,
    paste0(file_stem, "_unadj_ladj_null_results", table_suffix, ".tex")
  )
  alternative_path <- file.path(
    output_dir,
    paste0(file_stem, "_unadj_ladj_alternative_results", table_suffix, ".tex")
  )
  writeLines(null_lines, null_path)
  writeLines(alternative_lines, alternative_path)
  written_paths <- c(written_paths, null_path, alternative_path)

  estimand_block <- c(null_lines, "", "\\clearpage", "", alternative_lines)
  if (length(combined_lines) > 0L) {
    combined_lines <- c(combined_lines, "", "\\clearpage", "")
  }
  combined_lines <- c(combined_lines, estimand_block)
}

combined_stem <- paste(tolower(requested_estimands), collapse = "_")
combined_path <- file.path(
  output_dir,
  paste0(combined_stem, "_unadj_ladj_results_tables", table_suffix, ".tex")
)
writeLines(combined_lines, combined_path)
written_paths <- c(written_paths, combined_path)

cat("Wrote:\n", paste(normalizePath(written_paths), collapse = "\n"), "\n", sep = "")
