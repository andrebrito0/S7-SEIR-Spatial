#================================================================#
#  model_results.R  —  Visualise Simulation Results
#
#  Loads a saved workspace from run_scenarios.R and produces all
#  result plots used in the presentation.
#
#  Usage
#  ─────────────────────────────────────────────────────────────
#  Set sim_session to the .RData filename produced by run_scenarios.R,
#  then source this file from the repository root:
#    source("model_results.R")
#
#  The workspace image contains every object defined in config.R
#  and run_scenarios.R (parameters, ensemble data frames, summaries).
#  No additional sourcing is required.
#
#  Plots produced
#  ─────────────────────────────────────────────────────────────
#  Plot 1  — Global epidemic ribbon per scenario
#  Plot 2  — Per-patch mean curves, faceted by scenario (2 × 2)
#  Plot 3  — Final attack rate by scenario × patch (boxplot)
#  Plot 4  — Peak timing by scenario × patch (boxplot)
#  Plot 5  — Two-panel: epidemic curves + hospitalisation data
#  Plot 6  — Combined epidemic curves / hospitalisation /
#             temperature, with calendar-dated x-axis
#            (uncomment pdf() / dev.off() to export as PDF)
#================================================================#

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#----------------------------------------------------------------#
# § 0  LOAD WORKSPACE
#----------------------------------------------------------------#

rm(list = ls())

sim_session <- "sim_session_20260602.RData"   # <- update date as needed

source("S7_SEIR_spatial_core.R")              # S7 class definitions must be
# present before load() restores
# disease_spatial objects
load(file.path("results", sim_session))

# Colour palettes — defined here so they are available for all plots
# without being stored in the workspace image
scenario_cols <- c(
  Baseline    = "steelblue",
  Mobility    = "darkorange",
  Temperature = "darkgreen",
  Combined    = "firebrick"
)

patch_cols <- setNames(
  c("#378ADD", "#1D9E75", "#D85A30", "#7F77DD",
    "#BA7517", "#D4537E", "#888780"),
  PATCH_NAMES
)

dir.create(DIR_FIGURES, showWarnings = FALSE)


#================================================================#
# § 1  PLOT 1 — GLOBAL EPIDEMIC RIBBON
#
#  Mean active cases across all patches with 10–90% stochastic
#  band, one line per scenario.
#================================================================#

# png(file.path(DIR_FIGURES, "plot1_global_ribbon.png"), width = 900, height = 550)
par(mar = c(4, 4, 3, 2))

plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, max(global_rep$I_global) * 1.05),
     xlab = "Day", ylab = "Active cases (all patches)",
     main = "Global epidemic \u2014 mean \u00b1 10\u201390% by scenario")

for (sc in levels(ensemble_all$scenario)) {
  col <- scenario_cols[sc]
  d   <- global_ribbon[global_ribbon$scenario == sc, ]
  polygon(c(d$day, rev(d$day)), c(d$lo, rev(d$hi)),
          col = adjustcolor(col, 0.15), border = NA)
  lines(d$day, d$mean, col = col, lwd = 2.5)
}

legend("topright",
       legend = names(scenario_cols),
       col    = scenario_cols,
       lty = 1, lwd = 2.5, bty = "n")
# dev.off()


#================================================================#
# § 2  PLOT 2 — PER-PATCH CURVES FACETED BY SCENARIO
#
#  2 × 2 grid; one panel per scenario showing mean curves for
#  all patches.
#================================================================#

# png(file.path(DIR_FIGURES, "plot2_patch_curves_by_scenario.png"), width = 1100, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

for (sc in levels(ensemble_all$scenario)) {
  plot(NULL,
       xlim = c(1, N_DAYS),
       ylim = c(0, max(patch_ribbon$hi) * 1.05),
       xlab = "Day", ylab = "Active cases",
       main = paste("Patch curves \u2014", sc))
  
  for (p in PATCH_NAMES) {
    d <- patch_ribbon[patch_ribbon$scenario == sc & patch_ribbon$patch == p, ]
    polygon(c(d$day, rev(d$day)), c(d$lo, rev(d$hi)),
            col = adjustcolor(patch_cols[p], 0.15), border = NA)
    lines(d$day, d$mean, col = patch_cols[p], lwd = 2)
  }
  
  legend("topright",
         legend = PATCH_NAMES,
         col    = patch_cols,
         lty = 1, lwd = 2, bty = "n", cex = 0.75)
}

par(mfrow = c(1, 1))
# dev.off()


#================================================================#
# § 3  PLOT 3 — FINAL ATTACK RATE BY SCENARIO × PATCH
#================================================================#

# png(file.path(DIR_FIGURES, "plot3_attack_rate.png"), width = 1000, height = 550)
par(mar = c(5, 4, 3, 2))

boxplot(attack_rate * 100 ~ scenario + patch, data = epi_chars,
        col    = rep(scenario_cols, length(PATCH_NAMES)),
        border = "gray35",
        las = 2, xlab = "", ylab = "Final attack rate (%)",
        main = "Attack rate by scenario and patch")

stripchart(attack_rate * 100 ~ scenario + patch, data = epi_chars,
           method = "jitter", jitter = 0.15,
           pch = 16, cex = 0.7,
           col = adjustcolor("black", 0.4),
           vertical = TRUE, add = TRUE)
# dev.off()


#================================================================#
# § 4  PLOT 4 — PEAK TIMING BY SCENARIO × PATCH
#================================================================#

# png(file.path(DIR_FIGURES, "plot4_peak_timing.png"), width = 1000, height = 550)
par(mar = c(5, 4, 3, 2))

boxplot(peak_day ~ scenario + patch, data = epi_chars,
        col    = rep(scenario_cols, length(PATCH_NAMES)),
        border = "gray35",
        las = 2, xlab = "", ylab = "Day of peak incidence",
        main = "Peak timing by scenario and patch")

stripchart(peak_day ~ scenario + patch, data = epi_chars,
           method = "jitter", jitter = 0.15,
           pch = 16, cex = 0.7,
           col = adjustcolor("black", 0.4),
           vertical = TRUE, add = TRUE)
# dev.off()


#================================================================#
# § 5  HOSPITALISATION DATA PREPARATION
#
#  Shared across § 6 and § 7.  Builds hosp_plot (7-day moving
#  average of respiratory emergency visits in ARS LVT) and
#  temp_mean_df (cross-patch mean temperature with 7-day MA).
#================================================================#

hosp_plot <- hosp_gripe %>%
  filter(
    type == "n_episdios_urgencia_infecao",
    ars  == "ARS Lisboa e Vale do Tejo"
  ) %>%
  left_join(df_date, by = "date") %>%
  filter(!is.na(day)) %>%
  arrange(day) %>%
  mutate(ma_value = zoo::rollapply(value, 7, mean, align = "right", fill = NA))

# temp_matrix rows = patches, columns = days; rowMeans gives daily cross-patch mean
temp_mean_df <- data.frame(
  day    = seq_len(N_DAYS),
  date   = df_date$date,
  t_mean = rowMeans(t(temp_matrix))
) %>%
  mutate(t_ma = zoo::rollapply(t_mean, 7, mean, align = "right", fill = NA))

# Calendar axis tick positions and labels (used in §§ 6–7)
date_at     <- seq(SIM_START, SIM_END, by = "month")
date_labels <- format(date_at, "%b\n%Y")
day_at      <- as.numeric(date_at - SIM_START) + 1


#================================================================#
# § 6  PLOT 5 — TWO-PANEL: EPIDEMIC CURVES + HOSPITALISATIONS
#================================================================#

par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# Panel 1 — Simulated epidemic curves
plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, 600),
     xlab = "Day", ylab = "Active cases (all patches)",
     main = "Global epidemic \u2014 mean \u00b1 10\u201390% by scenario")

for (sc in levels(ensemble_all$scenario)) {
  col <- scenario_cols[sc]
  d   <- global_ribbon[global_ribbon$scenario == sc, ]
  polygon(c(d$day, rev(d$day)), c(d$lo, rev(d$hi)),
          col = adjustcolor(col, 0.15), border = NA)
  lines(d$day, d$mean, col = col, lwd = 2.5)
}

legend("topright",
       legend = names(scenario_cols),
       col    = scenario_cols,
       lty = 1, lwd = 2.5, bty = "n")

# Panel 2 — Observed respiratory emergency visits
plot(hosp_plot$day, hosp_plot$ma_value,
     type = "l", lwd = 2, col = "black",
     xlim = c(1, N_DAYS),
     ylim = c(0, max(hosp_plot$ma_value, na.rm = TRUE) * 1.05),
     xlab = "Day", ylab = "Emergency visits (7-day MA)",
     main = "Respiratory emergency visits \u2014 ARS Lisboa e Vale do Tejo")

par(mfrow = c(1, 1))


#================================================================#
# § 7  PLOT 6 — COMBINED EPIDEMIC CURVES PRESENTATION FIGURE
#
#  Layout: rows with heights 2 (epidemic curves) : 1 (hosp + temp)
#  Calendar-dated x-axis on both panels.
#  Temperature is scaled to the hospitalisation axis with a right-
#  hand secondary axis; T_mid is marked with a dashed reference line.
#
#  Uncomment pdf() / dev.off() to export for the Beamer slides.
#================================================================#

# pdf(file.path(DIR_FIGURES, "combined_curve_plot.pdf"),
#     width = 13, height = 8, bg = "transparent")

layout(matrix(1:2, nrow = 2), heights = c(2, 1))
par(mar = c(3, 4, 3, 2), oma = c(2, 0, 0, 0))

# ── Panel 1: Epidemic curves ──────────────────────────────────────
plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, 1500),
     xaxt = "n",
     xlab = "", ylab = "Active cases (all patches)",
     main = "Global epidemic \u2014 mean \u00b1 10\u201390% by scenario")
axis(1, at = day_at, labels = date_labels, cex.axis = 0.75)
grid()

for (sc in levels(ensemble_all$scenario)) {
  col <- scenario_cols[sc]
  d   <- global_ribbon[global_ribbon$scenario == sc, ]
  polygon(c(d$day, rev(d$day)), c(d$lo, rev(d$hi)),
          col = adjustcolor(col, 0.15), border = NA)
  lines(d$day, d$mean, col = col, lwd = 2.5)
}

legend("topright",
       legend = names(scenario_cols),
       col    = scenario_cols,
       lty = 1, lwd = 2.5, bty = "n", cex = 0.8)

# ── Panel 2: Hospitalisations (left axis) + Temperature (right axis)
# Map temperature onto the hospitalisation axis using a linear
# scaling so both series are visible in the same panel.
hosp_max   <- max(hosp_plot$ma_value, na.rm = TRUE) * 1.05
temp_min   <- min(temp_mean_df$t_ma,  na.rm = TRUE) - 1
temp_max   <- max(temp_mean_df$t_ma,  na.rm = TRUE) + 1
temp_range <- temp_max - temp_min

temp_to_hosp <- function(t) (t - temp_min) / temp_range * hosp_max
hosp_to_temp <- function(h) h / hosp_max * temp_range + temp_min

plot(hosp_plot$day, hosp_plot$ma_value,
     type = "l", lwd = 2, col = "black",
     xlim = c(1, N_DAYS),
     ylim = c(0, hosp_max),
     xaxt = "n",
     xlab = "", ylab = "Emergency visits (7-day MA)",
     main = "Respiratory emergency visits and mean temperature")
axis(1, at = day_at, labels = date_labels, cex.axis = 0.75)

# Overlay temperature on scaled hospitalisation axis
lines(temp_mean_df$day, temp_to_hosp(temp_mean_df$t_ma),
      lwd = 2, col = "#C96A6A")

# Reference line and label at T_mid
abline(h = temp_to_hosp(T_MID), lty = 2, col = "grey50")
text(x = N_DAYS, y = temp_to_hosp(T_MID) + hosp_max * 0.02,
     labels = sprintf("T_mid = %g\u00b0C", T_MID),
     col = "grey40", adj = 1, cex = 0.75)

# Secondary (right) axis showing actual temperature values
axis(4,
     at       = temp_to_hosp(seq(ceiling(temp_min), floor(temp_max), by = 2)),
     labels   = seq(ceiling(temp_min), floor(temp_max), by = 2),
     col.axis = "#C96A6A", col = "#C96A6A", cex.axis = 0.8)
mtext("Temperature (\u00b0C)", side = 4, line = 3, col = "#C96A6A", cex = 0.85)

legend("topright",
       legend = c("Emergency visits", "Mean temperature"),
       col    = c("black", "#C96A6A"),
       lty = 1, lwd = 2, bty = "n", cex = 0.8)
grid()

# Reset layout
layout(1)
par(oma = c(0, 0, 0, 0))

# dev.off()
