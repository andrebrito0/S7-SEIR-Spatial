#----------------------------------------------------------------#
#  S7-SEIR Spatial: 4-Scenario Ensemble Run
#
#  Scenarios:
#    1. Baseline     — no mobility coupling, no temperature scaling
#    2. Mobility     — commuting matrix couples patches, no temperature
#    3. Temperature  — local temperature scales beta, no mobility
#    4. Combined     — both mobility and temperature active
#
#  Usage:
#    Rscript run_scenarios.R
#    — or —
#    source("run_scenarios.R")   # from an interactive R session at repo root
#
#  Outputs written to results/ :
#    ensemble_baseline.rds, ensemble_mobility.rds,
#    ensemble_temperature.rds, ensemble_combined.rds,
#    ensemble_4scenarios.rds, epi_summary_4scenarios.csv,
#    processed_inputs.rds,
#    plot1_global_ribbon.png … plot4_peak_timing.png
#----------------------------------------------------------------#

# ── 0. Config & core model ────────────────────────────────────────
# config.R sets all shared parameters, installs missing packages,
# and provides check_data_file() and here::here() path helpers.
source("config.R")
source("S7_SEIR_spatial_core.R")

library(dplyr)
library(tidyr)
library(zoo)
library(reshape2)

#================================================================#
# 1. DATA LOADING & VALIDATION
#================================================================#

# Fail with a clear message if any required data file is absent,
# rather than a cryptic "cannot open connection" error.
check_data_file(PATH_TEMP_ARS)
check_data_file(PATH_OD_COMPLETE)
check_data_file(PATH_POP_INFO)

temp_ars    <- readRDS(PATH_TEMP_ARS)
od_complete <- readRDS(PATH_OD_COMPLETE)
pop_info    <- readRDS(PATH_POP_INFO)

#================================================================#
# 2. DATA PROCESSING
#================================================================#

# ── 2a. Temperature matrix ────────────────────────────────────────
sim_time <- as.numeric(SIM_END - SIM_START) + 1

temp_wide <- temp_ars %>%
  filter(municipio %in% PATCH_MUNICIPIO_MAP,
         date >= SIM_START, date <= SIM_END) %>%
  select(date, municipio, t_mean) %>%
  pivot_wider(names_from = municipio, values_from = t_mean) %>%
  arrange(date) %>%
  select(date, all_of(unname(PATCH_MUNICIPIO_MAP)))

temp_matrix <- t(as.matrix(na.approx(select(temp_wide, -date), rule = 2)))
rownames(temp_matrix) <- PATCH_NAMES

# Hard assertions: catch data problems before simulation
stopifnot(
  "temp_matrix row count must equal number of patches" =
    nrow(temp_matrix) == length(PATCH_NAMES),
  "temp_matrix column count must equal simulation length" =
    ncol(temp_matrix) == sim_time,
  "temp_matrix must not contain NAs after interpolation" =
    !anyNA(temp_matrix)
)

# ── 2b. Commuting matrix (Φ) ──────────────────────────────────────
# Rows = Origin, Cols = Destination
# Row-normalised: Φ[i,j] = probability of resident i visiting patch j
od_filtered <- od_complete %>%
  filter(Origin      %in% PATCH_MUNICIPIO_MAP,
         Destination %in% PATCH_MUNICIPIO_MAP)

node_order <- sort(unique(od_filtered$Origin))

mob_matrix <- od_filtered %>%
  dcast(Origin ~ Destination, value.var = "count", fill = 0) %>%
  arrange(factor(Origin, levels = node_order)) %>%
  select(all_of(node_order)) %>%
  as.matrix()

commuting_matrix <- mob_matrix / rowSums(mob_matrix)
commuting_matrix[is.nan(commuting_matrix)] <- 0   # safety for empty rows

# ── 2c. Save processed inputs ─────────────────────────────────────
# Saved so plots.R can load them directly — no duplicated processing.
dir.create(DIR_RESULTS, showWarnings = FALSE)
saveRDS(
  list(temp_matrix = temp_matrix, commuting_matrix = commuting_matrix),
  file.path(DIR_RESULTS, "processed_inputs.rds")
)

#================================================================#
# 3. DISEASE & SCENARIO OBJECTS
#================================================================#

influenza <- disease_spatial(
  name      = "Influenza (H3N2, 2017-18)",
  beta_core = BETA_CORE,   # ~0.16 contribution from 2 household contacts
  beta_ext  = BETA_EXT,    # ~0.19 contribution from 8 community contacts
  sigma     = SIGMA,       # 1/1.9 days latent period
  gamma     = GAMMA,       # 1/4.1 days infectious period
  use_temp  = FALSE
)
summary(influenza)
plot_beta_R0_profile(influenza)

influenza_temp <- disease_spatial(
  name      = "Influenza (H3N2, 2017-18)",
  beta_core = 0.16,   # ~0.16 contribution from 2 household contacts
  beta_ext  = 0.025,    # ~0.19 contribution from 8 community contacts
  sigma     = SIGMA,        # 1/1.9 days latent period
  gamma     = GAMMA,        # 1/4.1 days infectious period
  y_min     = Y_MIN,
  delta_y   = DELTA_Y,
  k_temp    = K_TEMP,
  T_mid     = T_MID,
  use_temp  = TRUE
)
summary(influenza_temp)
plot_beta_R0_profile(influenza_temp)

# Phi_mob  : census-derived commuting matrix — enables spatial coupling
# Phi_id   : identity matrix  — no cross-patch movement (baseline/temp scenarios)
Phi_mob <- commuting_matrix
Phi_id  <- diag(length(PATCH_NAMES))

N_DAYS  <- sim_time
N_total <- sum(PATCH_SIZES)

#================================================================#
# 4. HELPER FUNCTIONS
#================================================================#

#' extract_history
#'
#' Pull per-patch daily SEIR history from a finished metapopulation run
#' and annotate with patch metadata, rep ID, scenario label, and seed.
#'
#' @param mpop           A metapopulation S7 object (post-simulation)
#' @param rep_id         Integer replicate index
#' @param scenario_label Character label for this scenario
#' @param patch_names    Character vector of patch names (same order as patches)
#' @param patch_sizes    Integer vector of patch sizes (same order as patches)
#' @param seed_used      The RNG seed used for this replicate
#' @return A data frame with one row per patch × day
extract_history <- function(mpop, rep_id, scenario_label,
                            patch_names, patch_sizes, seed_used) {
  lapply(seq_along(mpop@patches), function(pi) {
    h            <- mpop@patches[[pi]]@history
    h$day        <- seq_len(nrow(h))
    h$patch      <- patch_names[pi]
    h$N          <- patch_sizes[pi]
    h$rep        <- rep_id
    h$scenario   <- scenario_label
    h$seed       <- seed_used      # trace which seed produced this replicate
    h
  }) |> bind_rows()
}

#' run_scenario_ensemble
#'
#' Run N_REPS stochastic replicates for one scenario configuration.
#'
#' @param label    Character label ("Baseline", "Mobility", etc.)
#' @param use_mob  Logical: use census commuting matrix (TRUE) or identity (FALSE)
#' @param use_temp Logical: use sigmoid temperature scaling (TRUE) or κ≡1 (FALSE)
#' @param n_reps   Number of stochastic replicates
#' @param n_days   Simulation length in days
#' @return A long data frame of all replicates, with seed column
run_scenario_ensemble <- function(label, use_mob, use_temp, n_reps, n_days, disease_obj) {
  
  Phi_arg  <- if (use_mob)  Phi_mob else Phi_id
  temp_arg <- if (use_temp) temp_matrix else NULL
  
  cat(sprintf("\n--- Scenario: %s ---\n", label))
  
  results <- vector("list", n_reps)
  
  for (rep in seq_len(n_reps)) {
    
    # Set seed BEFORE any random call, including create_metapopulation().
    # Record the seed in the output so every row of results is traceable.
    seed_used <- rep * 100L
    set.seed(seed_used)
    
    mpop <- create_metapopulation(
      patch_sizes = PATCH_SIZES,
      patch_names = PATCH_NAMES,
      n_infected  = N_INFECTED_SEED,
      k_core_pop  = K_CORE_POP,
      k_ext_pop   = K_EXT_POP,
      disease_obj = disease_obj,
      Phi         = Phi_arg,
      temp_series = temp_arg
    )
    
    for (day in seq_len(n_days)) mpop <- step_metapop(mpop)
    
    # Pass patch_names and patch_sizes explicitly — no implicit globals
    results[[rep]] <- extract_history(
      mpop,
      rep_id         = rep,
      scenario_label = label,
      patch_names    = PATCH_NAMES,
      patch_sizes    = PATCH_SIZES,
      seed_used      = seed_used
    )
    
    total_R <- sum(vapply(mpop@patches, \(p) tail(p@history$R, 1), integer(1)))
    cat(sprintf("  Rep %2d (seed %d) | Final recovered: %d (%.1f%%)\n",
                rep, seed_used, total_R, 100 * total_R / N_total))
  }
  
  bind_rows(results)
}

#================================================================#
# 5. RUN ALL 4 SCENARIOS
#================================================================#
if(TRUE){
  ensemble_baseline    <- run_scenario_ensemble("Baseline",
                                                use_mob  = FALSE, use_temp = FALSE,
                                                n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza)
  
  ensemble_mobility    <- run_scenario_ensemble("Mobility",
                                                use_mob  = TRUE,  use_temp = FALSE,
                                                n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza)
  
  ensemble_temperature <- run_scenario_ensemble("Temperature",
                                                use_mob  = FALSE, use_temp = TRUE,
                                                n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza_temp)
  
  ensemble_combined    <- run_scenario_ensemble("Combined",
                                                use_mob  = TRUE,  use_temp = TRUE,
                                                n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza_temp)
  
  # Combine all scenarios into one long data frame with ordered factor
  ensemble_all <- bind_rows(
    ensemble_baseline,
    ensemble_mobility,
    ensemble_temperature,
    ensemble_combined
  ) %>%
    mutate(scenario = factor(scenario,
                             levels = c("Baseline", "Mobility",
                                        "Temperature", "Combined")))
}

#================================================================#
# 6A. UPLOAD
#================================================================#
if(FALSE){
  ensemble_all         <- readRDS(file.path(DIR_RESULTS, "ensemble_all.rds"))
  ensemble_baseline    <- readRDS(file.path(DIR_RESULTS, "ensemble_baseline.rds"))
  ensemble_mobility    <- readRDS(file.path(DIR_RESULTS, "ensemble_mobility.rds"))
  ensemble_temperature <- readRDS(file.path(DIR_RESULTS, "ensemble_temperature.rds"))
  ensemble_combined    <- readRDS(file.path(DIR_RESULTS, "ensemble_combined.rds"))
}


#================================================================#
# 6. DERIVED SUMMARIES
#================================================================#

# Global (all patches) daily I per scenario × rep
global_rep <- ensemble_all %>%
  group_by(scenario, rep, day) %>%
  summarise(I_global = sum(I), .groups = "drop")

# Ribbon: mean + 10–90% quantile band across reps
global_ribbon <- global_rep %>%
  group_by(scenario, day) %>%
  summarise(
    mean   = mean(I_global),
    median = median(I_global),
    lo     = quantile(I_global, 0.10),
    hi     = quantile(I_global, 0.90),
    .groups = "drop"
  )

# Per-patch ribbon per scenario
patch_ribbon <- ensemble_all %>%
  group_by(scenario, patch, day) %>%
  summarise(
    mean = mean(I),
    lo   = quantile(I, 0.10),
    hi   = quantile(I, 0.90),
    .groups = "drop"
  )

# Epidemic characteristics per scenario × rep × patch
epi_chars <- ensemble_all %>%
  group_by(scenario, rep, patch, N) %>%
  summarise(
    attack_rate = last(R) / first(N),
    peak_I      = max(I),
    peak_day    = which.max(I),
    .groups     = "drop"
  )

# Patch-level summary table per scenario
epi_summary <- epi_chars %>%
  group_by(scenario, patch) %>%
  summarise(
    AR_mean_pct   = round(mean(attack_rate) * 100, 1),
    AR_sd_pct     = round(sd(attack_rate)   * 100, 1),
    peak_I_mean   = round(mean(peak_I), 0),
    peak_day_mean = round(mean(peak_day), 1),
    .groups = "drop"
  )

cat("\n=== Ensemble summary across scenarios ===\n")
print(epi_summary)

#================================================================#
# 7. SAVE RESULTS
#================================================================#

# One RDS per scenario (convenient for loading individually)
saveRDS(left_join(ensemble_all, df_date)       , file.path(DIR_RESULTS, "ensemble_all.rds"))
saveRDS(left_join(ensemble_baseline, df_date)  , file.path(DIR_RESULTS, "ensemble_baseline.rds"))
saveRDS(left_join(ensemble_mobility, df_date)  , file.path(DIR_RESULTS, "ensemble_mobility.rds"))
saveRDS(left_join(ensemble_temperature,df_date), file.path(DIR_RESULTS, "ensemble_temperature.rds"))
saveRDS(left_join(ensemble_combined, df_date)  , file.path(DIR_RESULTS, "ensemble_combined.rds"))

# Combined long-format ensemble + all summaries in one file
saveRDS(list(
  ensemble_all  = ensemble_all,
  global_rep    = global_rep,
  global_ribbon = global_ribbon,
  patch_ribbon  = patch_ribbon,
  epi_chars     = epi_chars,
  epi_summary   = epi_summary
), file.path(DIR_RESULTS, "ensemble_4scenarios.rds"))

write.csv(epi_summary,
          file.path(DIR_RESULTS, "epi_summary_4scenarios.csv"),
          row.names = FALSE)

cat("\nResults saved to", DIR_RESULTS, "\n")

#================================================================#
# 8. PLOTS
#================================================================#

dir.create(DIR_RESULTS, showWarnings = FALSE)

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

# ── Plot 1: Global epidemic ribbon per scenario ───────────────────
#png(file.path(DIR_RESULTS, "plot1_global_ribbon.png"), width = 900, height = 550)
par(mar = c(4, 4, 3, 2))

plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, max(global_rep$I_global) * 1.05),
     xlab = "Day", ylab = "Active cases (all patches)",
     main = "Global epidemic — mean \u00b1 10\u201390% by scenario")

for (sc in levels(ensemble_all$scenario)) {
  col <- scenario_cols[sc]
  d   <- global_ribbon[global_ribbon$scenario == sc, ]
  polygon(c(d$day, rev(d$day)), c(d$lo, rev(d$hi)),
          col = adjustcolor(col, 0.15), border = NA)
  lines(d$day, d$mean, col = col, lwd = 2.5)
}

legend("topright", legend = names(scenario_cols),
       col = scenario_cols, lty = 1, lwd = 2.5, bty = "n")
#dev.off()

# ── Plot 2: Per-patch mean curves, faceted by scenario ────────────
# png(file.path(DIR_RESULTS, "plot2_patch_curves_by_scenario.png"), width = 1100, height = 900)
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
  
  legend("topright", legend = PATCH_NAMES,
         col = patch_cols, lty = 1, lwd = 2, bty = "n", cex = 0.75)
}

par(mfrow = c(1, 1))
# dev.off()

# ── Plot 3: Attack rate by scenario × patch ───────────────────────
# png(file.path(DIR_RESULTS, "plot3_attack_rate.png"), width = 1000, height = 550)
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

# ── Plot 4: Peak timing by scenario × patch ───────────────────────
# png(file.path(DIR_RESULTS, "plot4_peak_timing.png"), width = 1000, height = 550)
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


# --- Prep hosp reference curve ---
df_date <- data.frame(
  date = seq(as.Date("2017-08-01"), by = "day", length.out = N_DAYS),
  day  = seq_len(N_DAYS)
)

hosp_plot <- hosp_gripe %>%
  filter(
    type == "n_episdios_urgencia_infecao",
    ars  == "ARS Lisboa e Vale do Tejo"
  ) %>%
  left_join(df_date, by = "date") %>%
  filter(!is.na(day)) %>%
  arrange(day) %>%
  mutate(ma_value = zoo::rollapply(value, 7, mean, align = "right", fill = NA))

# --- Two-panel layout ---
par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# --- Panel 1: Epidemic curves ---
plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, 600),
     xlab = "Day", ylab = "Active cases (all patches)",
     main = "Global epidemic — mean ± 10–90% by scenario")

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

# --- Panel 2: Hospitalizations ---
plot(hosp_plot$day, hosp_plot$ma_value,
     type = "l", lwd = 2, col = "black",
     xlim = c(1, N_DAYS),
     ylim = c(0, max(hosp_plot$ma_value, na.rm = TRUE) * 1.05),
     xlab = "Day", ylab = "Emergency visits (7-day MA)",
     main = "Respiratory emergency visits — ARS Lisboa e Vale do Tejo")

par(mfrow = c(1, 1))


# --- Prep dates ---
df_date <- data.frame(
  date = seq(as.Date("2017-08-01"), by = "day", length.out = N_DAYS),
  day  = seq_len(N_DAYS)
)

# --- Prep hosp ---
hosp_plot <- hosp_gripe %>%
  filter(
    type == "n_episdios_urgencia_infecao",
    ars  == "ARS Lisboa e Vale do Tejo"
  ) %>%
  left_join(df_date, by = "date") %>%
  filter(!is.na(day)) %>%
  arrange(day) %>%
  mutate(ma_value = zoo::rollapply(value, 7, mean, align = "right", fill = NA))

# --- Prep temperature (mean across patches) ---
temp_mean_df <- data.frame(
  day      = seq_len(N_DAYS),
  date     = df_date$date,
  t_mean   = rowMeans(t(temp_matrix))  # temp_matrix is [n_patches x N_DAYS]
) %>%
  mutate(t_ma = zoo::rollapply(t_mean, 7, mean, align = "right", fill = NA))

# --- Date axis helper ---
date_at     <- seq(as.Date("2017-08-01"), as.Date("2018-06-01"), by = "month")
date_labels <- format(date_at, "%b\n%Y")
day_at      <- as.numeric(date_at - as.Date("2017-08-01")) + 1

# --- Three-panel layout ---
par(mfrow = c(3, 1), mar = c(3, 4, 3, 2), oma = c(2, 0, 0, 0))

# --- Panel 1: Epidemic curves ---
plot(NULL,
     xlim = c(1, N_DAYS),
     ylim = c(0, 600),
     xaxt = "n",
     xlab = "", ylab = "Active cases (all patches)",
     main = "Global epidemic — mean ± 10–90% by scenario")

axis(1, at = day_at, labels = date_labels, cex.axis = 0.75)

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

# --- Panel 2: Hospitalizations ---
plot(hosp_plot$day, hosp_plot$ma_value,
     type = "l", lwd = 2, col = "black",
     xlim = c(1, N_DAYS),
     ylim = c(0, max(hosp_plot$ma_value, na.rm = TRUE) * 1.05),
     xaxt = "n",
     xlab = "", ylab = "Emergency visits (7-day MA)",
     main = "Respiratory emergency visits — ARS Lisboa e Vale do Tejo")

axis(1, at = day_at, labels = date_labels, cex.axis = 0.75)

# --- Panel 3: Temperature ---
plot(temp_mean_df$day, temp_mean_df$t_ma,
     type = "l", lwd = 2, col = "#C96A6A",
     xlim = c(1, N_DAYS),
     ylim = range(temp_mean_df$t_ma, na.rm = TRUE) + c(-1, 1),
     xaxt = "n",
     xlab = "", ylab = "Temperature (°C)",
     main = "Mean daily temperature across patches (7-day MA)")

axis(1, at = day_at, labels = date_labels, cex.axis = 0.75)

abline(h = 10, lty = 2, col = "grey50")  # T_mid reference line
text(x = N_DAYS, y = 10.4, labels = "T_mid = 10°C",
     col = "grey40", adj = 1, cex = 0.75)

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
