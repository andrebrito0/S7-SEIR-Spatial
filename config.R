#================================================================#
#  config.R — Single source of truth for all shared parameters
#
#  SOURCE THIS FILE at the top of run_scenarios.R and plots.R.
#  Never define these values in more than one place.
#
#  To change any parameter for a sensitivity analysis, edit here
#  and re-run both scripts — everything stays in sync.
#================================================================#

# ── R & package version guards ────────────────────────────────────
stopifnot(
  "R >= 4.3.0 is required" = getRversion() >= "4.3.0",
  "S7 >= 0.2.0 is required" = packageVersion("S7") >= "0.2.0"
)

# ── Required packages (install if missing) ────────────────────────
required_pkgs <- c(
  "S7", "dplyr", "tidyr", "ggplot2", "zoo", "lubridate",
  "reshape2", "readxl", "sf", "gganimate", "transformr",
  "gifski", "magick", "here"
)
missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

library(here)

# ── Simulation time window ────────────────────────────────────────
SIM_START  <- as.Date("2017-08-01")   # Start of simulation period
SIM_END    <- as.Date("2018-06-01")   # End of simulation period

df_date <- tibble(date = seq.Date(SIM_START, SIM_END, by = '1 day'))
df_date <- df_date %>% mutate(day = row_number())

# ── Plotting / data time window ───────────────────────────────────
# Used in plots.R to filter healthcare data before COVID disruption
PLOT_START  <- as.Date("2017-09-01")
PLOT_END    <- as.Date("2018-09-01")
DATE_BREAK  <- as.Date("2020-03-01")  # Cut-off: exclude COVID period

# ── Disease parameters (Influenza H3N2, 2017-18) ──────────────────
# Sources: Cowling et al. 2010 (serial interval);
#          Carrat et al. 2008 (infectious period)
BETA_CORE = 0.15   # ~0.16 contribution from 2 household contacts
BETA_EXT  = 0.02    # ~0.19 contribution from 8 community contacts
SIGMA     <- 0.5   # Rate E→I  (1 / 2 days incubation)
GAMMA     <- 0.2   # Rate I→R  (1 / 5 days infectious)

# ── Temperature sigmoid parameters ───────────────────────────────
Y_MIN   <- 0.40   # Floor: minimum κ at high temperature
DELTA_Y <- 1 - Y_MIN   # Range: κ ceiling = Y_MIN + DELTA_Y = 1.0
K_TEMP  <- 0.90   # Steepness (positive → colder = higher κ)
T_MID   <- 12.5   # Inflection temperature (°C)

# ── Patch definitions ─────────────────────────────────────────────
PATCH_NAMES <- c("City_A", "City_B", "City_C", "City_D",
                 "City_E", "City_F", "City_G")

PATCH_SIZES <- c(3000L, 2000, 1500L, 1200L, 800L, 2000L, 700L)

# Maps patch labels to municipality names in the temperature/OD data
PATCH_MUNICIPIO_MAP <- c(
  City_A = "Lisboa",
  City_B = "Santarém",
  City_C = "Setúbal",
  City_D = "Sintra",
  City_E = "Torres Vedras",
  City_F = "Cascais",
  City_G = "Mafra"
)

# ── Agent-level contact structure ────────────────────────────────
K_CORE_POP <- 3L   # Mean household (core group) size
K_EXT_POP  <- 8L   # Mean daily external contacts

# ── Initial infection seeds (one integer per patch) ──────────────
N_INFECTED_SEED <- c(5L, 5L, 1L, 1L, 1L, 1L, 1L)

# ── Ensemble settings ─────────────────────────────────────────────
N_REPS <- 10   # Number of stochastic replicates per scenario

# ── Data file paths (all relative to repo root via here::here()) ──
PATH_TEMP_ARS    <- here("data", "temp_ars.rds")
PATH_OD_COMPLETE <- here("data", "od_complete.rds")
PATH_POP_INFO    <- here("data", "pop_info.rds")
PATH_HEALTHCARE  <- here("data", "healthcare_behaviour.rds")
PATH_TEMP_PER_ARS <- here("data", "temp_per_ARS.rds")
PATH_CAOP        <- here("data", "Continente_CAOP2024_1.gpkg")
PATH_ARS_EXCEL   <- here("data", "ARS_Portugal_Continental_Estrutura.xlsx")

# ── Output directories ────────────────────────────────────────────
DIR_RESULTS <- here("results")
DIR_FIGURES <- here("figures")
DIR_FRAMES  <- here("frames")

# ── Helper: check a required data file exists, stop informatively ─
check_data_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required data file not found:\n  ", path, "\n",
      "See data/README_data.md for download and preparation instructions.\n",
      call. = FALSE
    )
  }
}

