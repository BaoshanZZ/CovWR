required <- c("data.table", "ggplot2", "scales", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)

parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 1L) sub(paste0("^--", name, "="), "", hit) else default
}
wilson <- function(p, n) {
  z <- stats::qnorm(0.975)
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  list(low = pmax(0, center - half), high = pmin(1, center + half))
}

input_dir <- parse_arg("input-dir", "simulation_results")
output_dir <- parse_arg("output-dir", file.path(input_dir, "main_report"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary <- data.table::fread(file.path(input_dir, "simulation_summary.csv"))
grid <- unique(data.table::fread(file.path(input_dir, "scenario_grid.csv"))[, .(scenario_id, scenario_change)])
dat <- merge(summary[estimand == "LWR" & method %in% c("UNADJ", "LADJ")], grid, by = "scenario_id")
scenario_map <- data.table::data.table(
  scenario_change = c(
    "baseline", "N_400", "N_500", "allocation_1_to_2", "weak_prognosis",
    "strong_prognosis", "heterogeneous_effect", "endpoint_specific_prognosis",
    "two_noise_covariates", "omit_X3"
  ),
  order = c(1L, 2L, 2L, 3:9),
  code = c("S1", "S2", "S2", paste0("S", 3:9)),
  label = c(
    "Ref.", "N=400", "N=500", "Alloc.\n1:2", "Weak\nprog.",
    "Strong\nprog.", "Effect\nhet.", "Endpoint",
    "Add\nnoise", "Omit\nkey"
  )
)
dat <- merge(dat, scenario_map, by = "scenario_change", all.x = TRUE)
if (anyNA(dat$order)) stop("Unmapped scenario.", call. = FALSE)
dat[, label := factor(label, levels = scenario_map$label)]
dat[, code := factor(code, levels = paste0("S", 1:9))]
dat[, method := factor(method, levels = c("UNADJ", "LADJ"))]
ci <- wilson(dat$reject_rate, dat$n_rep)
dat[, `:=`(mc_low = ci$low, mc_high = ci$high)]
coverage_ci <- wilson(dat$CP, dat$n_rep)
dat[, `:=`(
  coverage_mc_low = coverage_ci$low,
  coverage_mc_high = coverage_ci$high,
  bias_mc_low = bias - stats::qnorm(0.975) * MC_SE_mean,
  bias_mc_high = bias + stats::qnorm(0.975) * MC_SE_mean
)]

colors <- c(UNADJ = "#5B6770", LADJ = "#1F77B4")
dodge <- ggplot2::position_dodge(width = 0.42)
theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 0, hjust = 0.5, vjust = 1, size = 9, face = "bold"
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    plot.title.position = "plot"
  )

type1 <- ggplot2::ggplot(
  dat[effect == "null"],
  ggplot2::aes(label, reject_rate, color = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mc_low, ymax = mc_high),
    width = 0.14, linewidth = 0.65, position = dodge
  ) +
  ggplot2::geom_point(size = 2.2, position = dodge) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    limits = c(0.0415, 0.057)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Type I error (95% Monte Carlo CI)",
    title = "C. Type I error"
  ) + theme

alt <- dat[effect == "alternative"]
power_wide <- data.table::dcast(alt, label ~ method, value.var = "reject_rate")
power_limits <- range(c(alt$mc_low, alt$mc_high), finite = TRUE)
power_limits <- c(
  max(0, floor((power_limits[[1L]] - 0.03) * 20) / 20),
  min(1, ceiling((max(power_limits[[2L]], 0.90) + 0.02) * 20) / 20)
)
power <- ggplot2::ggplot(
  alt,
  ggplot2::aes(label, reject_rate, color = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed", color = "#B2182B") +
  ggplot2::geom_segment(
    data = power_wide,
    ggplot2::aes(x = label, xend = label, y = UNADJ, yend = LADJ),
    inherit.aes = FALSE, color = "grey65", linewidth = 0.7
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mc_low, ymax = mc_high),
    width = 0.14, linewidth = 0.6, position = dodge
  ) +
  ggplot2::geom_point(size = 2.1, position = dodge) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), limits = power_limits
  ) +
  ggplot2::labs(title = "A. Power", x = NULL, y = "Power (95% Monte Carlo CI)") + theme

pvr <- ggplot2::ggplot(
  alt[method == "LADJ"],
  ggplot2::aes(label, PVR, group = method)
) +
  ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  ggplot2::geom_segment(
    ggplot2::aes(x = label, xend = label, y = 0, yend = PVR),
    linewidth = 0.7, color = "grey70"
  ) +
  ggplot2::geom_point(size = 2.2, color = colors[["LADJ"]]) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), limits = c(0, 0.41)
  ) +
  ggplot2::labs(title = "B. Variance reduction", x = NULL, y = "PVR relative to UNADJ") +
  theme + ggplot2::theme(legend.position = "none")

bias_plot <- ggplot2::ggplot(
  alt,
  ggplot2::aes(label, bias, color = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = bias_mc_low, ymax = bias_mc_high),
    width = 0.14, linewidth = 0.6, position = dodge
  ) +
  ggplot2::geom_point(size = 2.1, position = dodge) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::labs(
    title = "A. Bias under the alternative", x = NULL,
    y = "Bias (95% Monte Carlo interval)"
  ) + theme

coverage_plot <- ggplot2::ggplot(
  alt,
  ggplot2::aes(label, CP, color = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = coverage_mc_low, ymax = coverage_mc_high),
    width = 0.14, linewidth = 0.6, position = dodge
  ) +
  ggplot2::geom_point(size = 2.1, position = dodge) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1), limits = c(0.94, 0.96)
  ) +
  ggplot2::labs(
    title = "B. Coverage under the alternative", x = NULL,
    y = "Coverage (95% Monte Carlo CI)"
  ) + theme

compact_x_theme <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(
    angle = 32, hjust = 1, vjust = 1, size = 8, face = "bold"
  )
)
bias_plot <- bias_plot + compact_x_theme
coverage_plot <- coverage_plot + compact_x_theme
type1_validity <- type1 + compact_x_theme

ggplot2::ggsave(file.path(output_dir, "lwr_type1_error.pdf"), type1, width = 10.5, height = 4.7)
ggplot2::ggsave(file.path(output_dir, "lwr_type1_error.png"), type1, width = 10.5, height = 4.7, dpi = 300)

draw_two_panel <- function(filename, png = FALSE) {
  if (png) {
    ragg::agg_png(filename, width = 13, height = 5.2, units = "in", res = 300)
  } else {
    grDevices::pdf(filename, width = 13, height = 5.2)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  layout <- grid::grid.layout(1, 2)
  grid::pushViewport(grid::viewport(layout = layout))
  print(power, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pvr, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}
draw_two_panel(file.path(output_dir, "lwr_alternative_power_pvr.pdf"))
draw_two_panel(file.path(output_dir, "lwr_alternative_power_pvr.png"), png = TRUE)

draw_validity_panels <- function(filename, png = FALSE) {
  if (png) {
    ragg::agg_png(filename, width = 16.5, height = 5.2, units = "in", res = 300)
  } else {
    grDevices::pdf(filename, width = 16.5, height = 5.2)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  layout <- grid::grid.layout(1, 3)
  grid::pushViewport(grid::viewport(layout = layout))
  print(bias_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(coverage_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(type1_validity, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
}
draw_validity_panels(file.path(output_dir, "lwr_validity_bias_coverage_type1.pdf"))
draw_validity_panels(file.path(output_dir, "lwr_validity_bias_coverage_type1.png"), png = TRUE)
data.table::fwrite(dat, file.path(output_dir, "lwr_plotting_data.csv"))
cat("Final LWR figures written to:", normalizePath(output_dir), "\n")
