required_plot_packages <- function() {
  c("data.table", "ggplot2", "scales")
}

check_plot_packages <- function() {
  missing <- required_plot_packages()[!vapply(required_plot_packages(), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required plotting package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

scenario_label <- function(N,
                           e,
                           prognosis,
                           outcome_shape,
                           effect_heterogeneity,
                           endpoint_prognosis_pattern = rep("shared", length(N))) {
  out <- rep("Baseline", length(N))
  out[N == 500L] <- "N = 500"
  out[abs(e - 2 / 3) < 1e-8] <- "e = 2/3"
  out[outcome_shape == "nonlinear"] <- "Nonlinear shape"
  out[prognosis == "strong"] <- "Strong prognosis"
  out[prognosis == "weak"] <- "Weak prognosis"
  out[effect_heterogeneity == "heterogeneous"] <- "Heterogeneous effect"
  out[endpoint_prognosis_pattern == "endpoint_specific"] <- "Endpoint-specific prognosis"
  out
}

prepare_plot_data <- function(summary_file) {
  dat <- data.table::fread(summary_file)
  dat[, scenario_label := scenario_label(
    N = N,
    e = e,
    prognosis = prognosis,
    outcome_shape = outcome_shape,
    effect_heterogeneity = effect_heterogeneity,
    endpoint_prognosis_pattern = if ("endpoint_prognosis_pattern" %in% names(dat)) {
      endpoint_prognosis_pattern
    } else {
      rep("shared", .N)
    }
  )]

  scenario_levels <- c(
    "Baseline",
    "N = 500",
    "e = 2/3",
    "Nonlinear shape",
    "Strong prognosis",
    "Weak prognosis",
    "Heterogeneous effect",
    "Endpoint-specific prognosis"
  )
  dat[, scenario_label := factor(scenario_label, levels = scenario_levels)]
  dat[, method := factor(method, levels = c("UNADJ", "LADJ", "FADJ"))]
  dat[, estimand := factor(estimand, levels = c("NB", "LWR", "LWO"))]
  dat[, effect := factor(effect, levels = c("null", "alternative"), labels = c("Null", "Alternative"))]
  dat[]
}

report_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      plot.title.position = "plot"
    )
}

method_colors <- c(
  "UNADJ" = "#5B6770",
  "LADJ" = "#1F77B4",
  "FADJ" = "#D95F02"
)

save_report_plot <- function(plot, output_dir, stem, width = 12, height = 7) {
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
}

plot_reject_rate <- function(dat) {
  alternative_ref <- data.frame(
    effect = factor("Alternative", levels = levels(dat$effect)),
    yintercept = 0.90
  )

  ggplot2::ggplot(dat, ggplot2::aes(scenario_label, reject_rate, color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    ggplot2::geom_hline(
      data = alternative_ref,
      ggplot2::aes(yintercept = yintercept),
      linetype = "dashed",
      color = "#B2182B",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(effect ~ estimand, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(
      title = "Rejection Rate by Scenario",
      x = NULL,
      y = "Reject rate",
      color = NULL
    ) +
    report_theme()
}

plot_pvr <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(scenario_label, PVR, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
    ggplot2::facet_grid(effect ~ estimand) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "Variance Reduction Relative to UNADJ",
      x = NULL,
      y = "PVR",
      fill = NULL
    ) +
    report_theme()
}

plot_esd <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(scenario_label, ESD, color = method, group = method)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(effect ~ estimand, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(
      title = "Empirical Standard Deviation by Scenario",
      x = NULL,
      y = "ESD",
      color = NULL
    ) +
    report_theme()
}

plot_coverage <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(scenario_label, CP, color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(effect ~ estimand) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0.90, 1.00)) +
    ggplot2::labs(
      title = "95% Wald Interval Coverage by Scenario",
      x = NULL,
      y = "Coverage",
      color = NULL
    ) +
    report_theme()
}

plot_bias <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(scenario_label, bias, color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(effect ~ estimand, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(
      title = "Empirical Bias by Scenario",
      x = NULL,
      y = "Bias",
      color = NULL
    ) +
    report_theme()
}

plot_alternative_power <- function(dat) {
  alt <- dat[effect == "Alternative"]
  ggplot2::ggplot(alt, ggplot2::aes(scenario_label, reject_rate, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
    ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
    ggplot2::facet_wrap(~ estimand) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(
      title = "Power Under the Alternative",
      x = NULL,
      y = "Power",
      fill = NULL
    ) +
    report_theme()
}

make_simulation_report_figures <- function(input_dir = "WRCov/simulation_output_one_factor_zeta",
                                           output_dir = NULL,
                                           estimand_filter = NULL) {
  check_plot_packages()
  if (!is.null(estimand_filter)) {
    estimand_filter <- match.arg(estimand_filter, c("NB", "LWR", "LWO"))
  }
  if (is.null(output_dir)) {
    output_dir <- if (is.null(estimand_filter)) {
      file.path(input_dir, "figures")
    } else {
      file.path(input_dir, paste0("figures_", tolower(estimand_filter)))
    }
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_file <- file.path(input_dir, "simulation_summary.csv")
  dat <- prepare_plot_data(summary_file)
  if (!is.null(estimand_filter)) {
    dat <- dat[estimand == estimand_filter]
    dat[, estimand := droplevels(estimand)]
  }

  plots <- list(
    reject_rate_by_scenario = plot_reject_rate(dat),
    pvr_by_scenario = plot_pvr(dat),
    esd_by_scenario = plot_esd(dat),
    coverage_by_scenario = plot_coverage(dat),
    bias_by_scenario = plot_bias(dat),
    power_alternative = plot_alternative_power(dat)
  )

  for (nm in names(plots)) {
    save_report_plot(plots[[nm]], output_dir, nm)
  }

  data.table::fwrite(dat, file.path(output_dir, "plotting_data_with_labels.csv"))
  invisible(list(data = dat, plots = plots, output_dir = output_dir))
}

running_as_direct_plot_script <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    return(FALSE)
  }
  script_file <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  this_file <- normalizePath("WRCov/R/plot_simulation_results.R", mustWork = FALSE)
  identical(script_file, this_file)
}

if (running_as_direct_plot_script()) {
  args <- commandArgs(trailingOnly = TRUE)
  input_dir_arg <- grep("^--input-dir=", args, value = TRUE)
  output_dir_arg <- grep("^--output-dir=", args, value = TRUE)
  estimand_arg <- grep("^--estimand=", args, value = TRUE)
  input_dir <- if (length(input_dir_arg) == 1L) {
    sub("^--input-dir=", "", input_dir_arg)
  } else {
    "WRCov/simulation_output_one_factor_zeta"
  }
  output_dir <- if (length(output_dir_arg) == 1L) {
    sub("^--output-dir=", "", output_dir_arg)
  } else {
    NULL
  }
  estimand_filter <- if (length(estimand_arg) == 1L) {
    sub("^--estimand=", "", estimand_arg)
  } else {
    NULL
  }
  res <- make_simulation_report_figures(
    input_dir = input_dir,
    output_dir = output_dir,
    estimand_filter = estimand_filter
  )
  cat("Figures written to:", res$output_dir, "\n")
}
