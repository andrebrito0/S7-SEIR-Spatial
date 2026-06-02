#================================================================#
#  run_scenarios.R  —  Four-Scenario Ensemble Run
#
#  Runs N_REPS stochastic replicates for each of four scenarios
#  that independently toggle mobility and temperature coupling,
#  then derives epidemic summaries and saves all results.
#
#  Scenarios
#  ─────────────────────────────────────────────────────────────
#  Baseline    — identity Φ,  κ ≡ 1   (homogeneous reference)
#  Mobility    — census Φ,    κ ≡ 1   (spatial spread only)
#  Temperature — identity Φ,  κ(T)    (seasonal forcing only)
#  Combined    — census Φ,    κ(T)    (full coupled model)
#
#  Usage
#  ─────────────────────────────────────────────────────────────
#  From an interactive session at the repository root:
#    source("run_scenarios.R")
#  From the command line:
#    Rscript run_scenarios.R
#
#  Prerequisites
#  ─────────────────────────────────────────────────────────────
#  All data files listed in config.R § 5 must be present in data/.
#  See data/README_data.md for download instructions.
#
#  Outputs  →  results/
#  ─────────────────────────────────────────────────────────────
#  processed_inputs.rds        — temp_matrix + commuting_matrix
#  ensemble_baseline.rds       — long data frame, Baseline reps
#  ensemble_mobility.rds       — long data frame, Mobility reps
#  ensemble_temperature.rds    — long data frame, Temperature reps
#  ensemble_combined.rds       — long data frame, Combined reps
#  ensemble_4scenarios.rds     — all four scenarios in one list
#  epi_summary_4scenarios.csv  — attack rate & peak timing table
#  sim_session_YYYYMMDD.RData  — full workspace image
#================================================================#


#----------------------------------------------------------------#
# § 0  DEPENDENCIES
#----------------------------------------------------------------#

source("config.R")
source("S7_SEIR_spatial_core.R")

library(dplyr)
library(tidyr)
library(zoo)
library(reshape2)


#================================================================#
# § 1  DATA LOADING & VALIDATION
#================================================================#

# Fail with a clear message if any required file is absent
check_data_file(PATH_TEMP_ARS)
check_data_file(PATH_OD_COMPLETE)
check_data_file(PATH_POP_INFO)
check_data_file(PATH_HEALTHCARE)

temp_ars    <- readRDS(PATH_TEMP_ARS)
od_complete <- readRDS(PATH_OD_COMPLETE)
pop_info    <- readRDS(PATH_POP_INFO)
hosp_gripe  <- readRDS(PATH_HEALTHCARE)


#================================================================#
# § 2  DATA PROCESSING
#================================================================#

# ── 2a. Temperature matrix  [n_patches × sim_time] ───────────────
# Filter to simulation window and selected municipalities, pivot to
# wide format (one column per municipality), then transpose so that
# rows = patches and columns = days, matching step_metapop()'s
# expected layout.  na.approx fills any isolated missing days.
temp_wide <- temp_ars %>%
  filter(municipio %in% PATCH_MUNICIPIO_MAP,
         date >= SIM_START, date <= SIM_END) %>%
  select(date, municipio, t_mean) %>%
  pivot_wider(names_from = municipio, values_from = t_mean) %>%
  arrange(date) %>%
  select(date, all_of(unname(PATCH_MUNICIPIO_MAP)))

temp_matrix <- t(as.matrix(na.approx(select(temp_wide, -date), rule = 2)))
rownames(temp_matrix) <- PATCH_NAMES

# Hard assertions — catch data problems before the simulation starts
stopifnot(
  "temp_matrix row count must equal number of patches" =
    nrow(temp_matrix) == length(PATCH_NAMES),
  "temp_matrix column count must equal simulation length" =
    ncol(temp_matrix) == sim_time,
  "temp_matrix must not contain NAs after interpolation" =
    !anyNA(temp_matrix)
)

# ── 2b. Commuting matrix Φ  [n_patches × n_patches] ─────────────
# Rows = origin patch, columns = destination patch.
# Row-normalised: Φ[i,j] = fraction of time residents of patch i
# spend in patch j on a typical day.  Rows sum to 1 by construction.
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
commuting_matrix[is.nan(commuting_matrix)] <- 0   # guard: empty-row patches

# ── 2c. Save processed inputs ─────────────────────────────────────
# Persisted so plots.R can load them without re-processing raw data.
dir.create(DIR_RESULTS, showWarnings = FALSE)
saveRDS(
  list(temp_matrix = temp_matrix, commuting_matrix = commuting_matrix),
  file.path(DIR_RESULTS, "processed_inputs.rds")
)


#================================================================#
# § 3  DISEASE & SCENARIO OBJECTS
#================================================================#

# Non-temperature scenarios: κ ≡ 1 (use_temp = FALSE)
influenza <- disease_spatial(
  name      = "Influenza (H3N2, 2017-18)",
  beta_core = BETA_CORE,
  beta_ext  = BETA_EXT,
  sigma     = SIGMA,
  gamma     = GAMMA,
  use_temp  = FALSE
)
summary(influenza)
plot_beta_R0_profile(influenza)

# Temperature-coupled scenarios: κ(T) active (use_temp = TRUE).
# beta_core / beta_ext are recalibrated upward — see config.R § 2.
influenza_temp <- disease_spatial(
  name      = "Influenza (H3N2, 2017-18)",
  beta_core = BETA_CORE_TEMP,
  beta_ext  = BETA_EXT_TEMP,
  sigma     = SIGMA,
  gamma     = GAMMA,
  y_min     = Y_MIN,
  delta_y   = DELTA_Y,
  k_temp    = K_TEMP,
  T_mid     = T_MID,
  use_temp  = TRUE
)
summary(influenza_temp)
plot_beta_R0_profile(influenza_temp)

# Mobility matrix variants
Phi_mob <- commuting_matrix          # Census OD — enables spatial coupling
Phi_id  <- diag(length(PATCH_NAMES)) # Identity  — no cross-patch movement

N_DAYS  <- sim_time
N_total <- sum(PATCH_SIZES)


#================================================================#
# § 4  HELPER FUNCTIONS
#================================================================#

#' extract_history
#'
#' Pull per-patch daily SEIR history from a completed metapopulation
#' run and annotate each row with patch metadata, replicate ID,
#' scenario label, and the RNG seed used.
#'
#' @param mpop           A metapopulation S7 object (post-simulation).
#' @param rep_id         Integer replicate index.
#' @param scenario_label Character label for this scenario.
#' @param patch_names    Character vector of patch names (same order as patches).
#' @param patch_sizes    Integer vector of patch sizes (same order as patches).
#' @param seed_used      The RNG seed set before this replicate.
#' @return A data frame with one row per patch × day.
extract_history <- function(mpop, rep_id, scenario_label,
                            patch_names, patch_sizes, seed_used) {
  lapply(seq_along(mpop@patches), function(pi) {
    h          <- mpop@patches[[pi]]@history
    h$day      <- seq_len(nrow(h))
    h$patch    <- patch_names[pi]
    h$N        <- patch_sizes[pi]
    h$rep      <- rep_id
    h$scenario <- scenario_label
    h$seed     <- seed_used
    h
  }) |> bind_rows()
}


#' run_scenario_ensemble
#'
#' Run N_REPS stochastic replicates for one scenario configuration.
#' Seeds are set to rep × 100 before every random call, including
#' population construction, so results are exactly reproducible.
#'
#' @param label       Character label ("Baseline", "Mobility", etc.).
#' @param use_mob     Logical: use census Φ (TRUE) or identity Φ (FALSE).
#' @param use_temp    Logical: use sigmoid κ(T) (TRUE) or κ ≡ 1 (FALSE).
#' @param n_reps      Number of stochastic replicates.
#' @param n_days      Simulation length in days.
#' @param disease_obj A disease_spatial object.
#' @return A long data frame of all replicates, including a seed column.
run_scenario_ensemble <- function(label, use_mob, use_temp,
                                  n_reps, n_days, disease_obj) {
  
  Phi_arg  <- if (use_mob)  Phi_mob else Phi_id
  temp_arg <- if (use_temp) temp_matrix else NULL
  
  cat(sprintf("\n--- Scenario: %s ---\n", label))
  results <- vector("list", n_reps)
  
  for (rep in seq_len(n_reps)) {
    
    # Seed is set BEFORE any random call, including create_metapopulation(),
    # so population structure and infection draws are fully reproducible.
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
    
    results[[rep]] <- extract_history(
      mpop,
      rep_id         = rep,
      scenario_label = label,
      patch_names    = PATCH_NAMES,
      patch_sizes    = PATCH_SIZES,
      seed_used      = seed_used
    )
    
    total_R <- sum(vapply(mpop@patches, \(p) tail(p@history$R, 1), integer(1)))
    cat(sprintf("  Rep %2d (seed %4d) | Final recovered: %d (%.1f%%)\n",
                rep, seed_used, total_R, 100 * total_R / N_total))
  }
  
  bind_rows(results)
}


#================================================================#
# § 5  RUN ALL FOUR SCENARIOS
#================================================================#

ensemble_baseline <- run_scenario_ensemble(
  "Baseline",
  use_mob = FALSE, use_temp = FALSE,
  n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza
)

ensemble_mobility <- run_scenario_ensemble(
  "Mobility",
  use_mob = TRUE, use_temp = FALSE,
  n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza
)

ensemble_temperature <- run_scenario_ensemble(
  "Temperature",
  use_mob = FALSE, use_temp = TRUE,
  n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza_temp
)

ensemble_combined <- run_scenario_ensemble(
  "Combined",
  use_mob = TRUE, use_temp = TRUE,
  n_reps = N_REPS, n_days = N_DAYS, disease_obj = influenza_temp
)

# Combine into one long data frame with an ordered scenario factor
# and a calendar date column derived from config.R's df_date
ensemble_all <- bind_rows(
  ensemble_baseline,
  ensemble_mobility,
  ensemble_temperature,
  ensemble_combined
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Baseline", "Mobility",
                                      "Temperature", "Combined"))) %>%
  mutate(date = SIM_START + day - 1L)


#================================================================#
# § 6  DERIVED SUMMARIES
#================================================================#

# Global (all patches) daily infectious count per scenario × rep
global_rep <- ensemble_all %>%
  group_by(scenario, rep, day) %>%
  summarise(I_global = sum(I), .groups = "drop")

# Ribbon summary: mean + 10–90% interquantile range across reps
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

# Patch-level summary table per scenario (mean ± SD across reps)
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
# § 7  SAVE RESULTS
#================================================================#

# Individual scenario data frames
saveRDS(ensemble_baseline,    file.path(DIR_RESULTS, "ensemble_baseline.rds"))
saveRDS(ensemble_mobility,    file.path(DIR_RESULTS, "ensemble_mobility.rds"))
saveRDS(ensemble_temperature, file.path(DIR_RESULTS, "ensemble_temperature.rds"))
saveRDS(ensemble_combined,    file.path(DIR_RESULTS, "ensemble_combined.rds"))

# Combined list with derived summaries
saveRDS(
  list(
    ensemble_all   = ensemble_all,
    global_rep     = global_rep,
    global_ribbon  = global_ribbon,
    patch_ribbon   = patch_ribbon,
    epi_chars      = epi_chars,
    epi_summary    = epi_summary
  ),
  file.path(DIR_RESULTS, "ensemble_4scenarios.rds")
)

# Summary CSV for quick inspection outside R
write.csv(epi_summary,
          file.path(DIR_RESULTS, "epi_summary_4scenarios.csv"),
          row.names = FALSE)

# Full workspace image — loaded by model_results.R
save_date <- format(Sys.Date(), format = "%Y%m%d")
save.image(file.path(DIR_RESULTS, paste0("sim_session_", save_date, ".RData")))

cat("\nAll results saved to", DIR_RESULTS, "\n")
