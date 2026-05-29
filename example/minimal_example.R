#================================================================#
#  Minimal Worked Example — S7-SEIR Spatial
#
#  Self-contained: uses synthetic data only.
#  No external data files required.
#  Expected runtime: < 1 minute
#
#  What this script demonstrates:
#    1. Constructing a disease_spatial object
#    2. Building a 3-patch metapopulation with a mobility matrix
#    3. Running a 120-day simulation
#    4. Visualising SEIR curves and a prevalence heatmap
#    5. Comparing baseline vs. temperature-coupled scenarios
#================================================================#

source("S7_SEIR_spatial_core.R")

set.seed(42)

# ── 1. Define the pathogen ────────────────────────────────────────
# Parameters calibrated loosely to influenza (H3N2)
influenza_example <- disease_spatial(
  name      = "Example Influenza",
  beta_core = 0.08,    # Per-contact household transmission risk
  beta_ext  = 0.024,   # Per-contact community transmission risk
  sigma     = 0.526,   # 1 / 1.9 days incubation
  gamma     = 0.244,   # 1 / 4.1 days infectious
  y_min     = 0.20,    # Minimum transmission scaling (warm weather)
  delta_y   = 0.80,    # Range of scaling
  k_temp    = 0.30,    # Sigmoid steepness
  T_mid     = 10.0,     # Inflection temperature (°C)
  use_temp  = TRUE
)

summary(influenza_example)

# ── 2. Visualise the temperature coupling curve ───────────────────
par(mfrow = c(1, 1))
plot_temp_scaling_curve(influenza_example)

# ── 3. Synthetic temperature series (120 days, autumn → spring) ──
# Simulates a seasonal temperature cycle dropping into winter
n_days   <- 120
n_patches <- 3

# Three patches with slightly different temperatures (e.g. urban vs rural)
day_seq   <- seq_len(n_days)
temp_base <- 15 - 10 * sin(pi * day_seq / n_days)   # 15°C → 5°C → 15°C arc

temp_synthetic <- rbind(
  temp_base,           # Patch 1: baseline
  temp_base - 2,       # Patch 2: 2°C cooler (e.g. higher elevation)
  temp_base + 1        # Patch 3: 1°C warmer (e.g. coastal)
)

# ── 4. Mobility matrix ─────────────────────────────────────────────
# Residents mostly stay local; some commuting between patches 1 & 2
Phi_example <- matrix(c(
  0.80, 0.15, 0.05,   # Patch 1: 80% local, 15% to patch 2, 5% to patch 3
  0.10, 0.85, 0.05,   # Patch 2: 10% to patch 1, 85% local
  0.05, 0.05, 0.90    # Patch 3: mostly local
), nrow = 3, byrow = TRUE)

# Verify row-stochastic
stopifnot(all(abs(rowSums(Phi_example) - 1) < 1e-10))

# ── 5. Build and run: COMBINED scenario ──────────────────────────
cat("\n--- Running: Combined (mobility + temperature) ---\n")
mpop_combined <- create_metapopulation(
  patch_sizes = c(500L, 400L, 300L),
  patch_names = c("City", "Suburb", "Rural"),
  n_infected  = c(5L, 1L, 0L),     # Seed only in City
  k_core_pop  = 3L,
  k_ext_pop   = 8L,
  disease_obj = influenza_example,
  Phi         = Phi_example,
  temp_series = temp_synthetic
)

for (day in seq_len(n_days)) mpop_combined <- step_metapop(mpop_combined)

summary(mpop_combined)

# ── 6. Build and run: BASELINE scenario (no mobility, no temp) ───
cat("\n--- Running: Baseline ---\n")
set.seed(42)
mpop_baseline <- create_metapopulation(
  patch_sizes = c(500L, 400L, 300L),
  patch_names = c("City", "Suburb", "Rural"),
  n_infected  = c(5L, 1L, 0L),
  k_core_pop  = 3L,
  k_ext_pop   = 8L,
  disease_obj = influenza_example
  # Phi defaults to identity; temp_series defaults to T_mid → κ = constant
)

for (day in seq_len(n_days)) mpop_baseline <- step_metapop(mpop_baseline)

# ── 7. Plot SEIR curves for each scenario ────────────────────────
cat("\nPlotting SEIR curves...\n")

par(mfrow = c(1, 2))

# Combined
h_c <- do.call(rbind, lapply(mpop_combined@patches, \(p) {
  h <- p@history; h$patch <- p@name; h
}))
plot(NULL,
     xlim = c(1, n_days), ylim = c(0, 500),
     xlab = "Day", ylab = "Agents",
     main = "Combined (mobility + temperature)")
patch_colors <- c(City = "steelblue", Suburb = "darkorange", Rural = "darkgreen")
for (p in c("City", "Suburb", "Rural")) {
  d <- h_c[h_c$patch == p, ]
  lines(d$day, d$I, col = patch_colors[p], lwd = 2)
}
legend("topright", legend = names(patch_colors),
       col = patch_colors, lty = 1, lwd = 2, bty = "n")

# Baseline
h_b <- do.call(rbind, lapply(mpop_baseline@patches, \(p) {
  h <- p@history; h$patch <- p@name; h
}))
plot(NULL,
     xlim = c(1, n_days), ylim = c(0, 500),
     xlab = "Day", ylab = "Agents",
     main = "Baseline (no mobility, κ ≡ 1)")
for (p in c("City", "Suburb", "Rural")) {
  d <- h_b[h_b$patch == p, ]
  lines(d$day, d$I, col = patch_colors[p], lwd = 2)
}
legend("topright", legend = names(patch_colors),
       col = patch_colors, lty = 1, lwd = 2, bty = "n")

par(mfrow = c(1, 1))

# ── 8. Prevalence heatmap ─────────────────────────────────────────
cat("\nPlotting prevalence heatmap (Combined scenario)...\n")
plot_prevalence_heatmap(mpop_combined)

# ── 9. Attack rate comparison ─────────────────────────────────────
cat("\n=== Final Attack Rates ===\n")
for (scenario_name in c("Combined", "Baseline")) {
  mpop <- if (scenario_name == "Combined") mpop_combined else mpop_baseline
  cat(sprintf("\n%s:\n", scenario_name))
  for (p in mpop@patches) {
    ar <- tail(p@history$R, 1) / p@pop_size * 100
    cat(sprintf("  %-8s: %.1f%%\n", p@name, ar))
  }
}

