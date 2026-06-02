#================================================================#
#  config.R  —  Parameters
#
#  SOURCE this file at the top of every analysis script:
#    source("config.R")
#
#  Never redefine these values elsewhere. For a sensitivity
#  analysis, edit here and re-run all downstream scripts.
#
#  Contents
#  ────────────────────────────────────────────────────────────
#  § 0   R & package version guards
#  § 1   Simulation time window
#  § 2   Disease parameters
#  § 3   Temperature sigmoid parameters
#  § 4   Spatial & demographic parameters
#  § 5   File paths & output directories
#  § 6   Utility helpers
#================================================================#


#----------------------------------------------------------------#
# § 0  R & PACKAGE VERSION GUARDS
#----------------------------------------------------------------#

stopifnot(
  "R >= 4.3.0 is required" = getRversion() >= "4.3.0",
  "S7 >= 0.2.0 is required" = packageVersion("S7") >= "0.2.0"
)

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

invisible(lapply(required_pkgs, library, character.only = TRUE))
message("All packages loaded successfully.")

library(here)


#----------------------------------------------------------------#
# § 1  SIMULATION TIME WINDOW
#----------------------------------------------------------------#

SIM_START <- as.Date("2017-08-01")   # First simulated day (day 1)
SIM_END   <- as.Date("2018-06-01")   # Last simulated day  (day sim_time)

# Integer number of simulated days — used for array sizing & assertions
sim_time <- as.integer(SIM_END - SIM_START) + 1L   # 305 days

# Day-to-date lookup table: joined onto ensemble output in run_scenarios.R
df_date <- tibble(date = seq.Date(SIM_START, SIM_END, by = "1 day")) %>%
  mutate(day = row_number())

# Plotting window (healthcare data, excluding COVID disruption)
PLOT_START <- as.Date("2017-09-01")
PLOT_END   <- as.Date("2018-09-01")
DATE_BREAK <- as.Date("2020-03-01")   # Exclude all data on or after this date


#----------------------------------------------------------------#
# § 2  DISEASE PARAMETERS  (Influenza H3N2, 2017–18 season)
#
#  References:
#    Incubation / infectious periods — Carrat et al. (2008)
#    Serial interval               — Cowling et al. (2010)
#
#  Two β calibrations are required because the sigmoid κ(T)
#  suppresses mean transmissibility below 1 over the season.
#  Temperature-coupled scenarios use raised base rates so that
#  the seasonal-mean attack rate matches the non-coupled baseline.
#----------------------------------------------------------------#

# Non-temperature scenarios (Baseline, Mobility) — κ ≡ 1
BETA_CORE      <- 0.10    # Transmission risk per infectious household contact
BETA_EXT       <- 0.021   # Transmission risk per infectious community contact

# Temperature-coupled scenarios (Temperature, Combined) — κ(T) active
BETA_CORE_TEMP <- 0.13
BETA_EXT_TEMP  <- 0.033

SIGMA <- 0.526   # Progression rate E→I  (mean incubation = 1.9 days)
GAMMA <- 0.244   # Recovery rate    I→R  (mean infectious period = 4.1 days)


#----------------------------------------------------------------#
# § 3  TEMPERATURE SIGMOID PARAMETERS
#
#  Thermal response curve:
#    κ(T) = y_min + Δy / (1 + exp(k_temp · (T − T_mid)))
#
#  κ approaches 1.0 as T falls below T_mid (cold → more transmission)
#  κ approaches y_min as T rises above T_mid (warm → less transmission)
#----------------------------------------------------------------#

Y_MIN   <- 0.40          # Floor: minimum κ at high temperature
DELTA_Y <- 1 - Y_MIN    # Range: κ ceiling = Y_MIN + DELTA_Y = 1.0
K_TEMP  <- 0.50          # Steepness (positive → colder = higher κ)
T_MID   <- 13            # Inflection temperature (°C)


#----------------------------------------------------------------#
# § 4  SPATIAL & DEMOGRAPHIC PARAMETERS
#----------------------------------------------------------------#

# Patch labels used throughout — order must match PATCH_SIZES and
# PATCH_MUNICIPIO_MAP
PATCH_NAMES <- c("City_A", "City_B", "City_C", "City_D",
                 "City_E", "City_F", "City_G")

# Census population size for each patch (must sum to N_total = 11 200)
PATCH_SIZES <- c(3000L, 2000L, 1500L, 1200L, 800L, 2000L, 700L)

# Maps patch labels → municipality names used in the temperature and OD data
PATCH_MUNICIPIO_MAP <- c(
  City_A = "Lisboa",
  City_B = "Santarém",
  City_C = "Setúbal",
  City_D = "Sintra",
  City_E = "Torres Vedras",
  City_F = "Cascais",
  City_G = "Mafra"
)

# Contact structure
K_CORE_POP <- 3L   # Target mean household (core group) size
K_EXT_POP  <- 8L   # Target mean daily external (community) contacts

# Initial infectious seeds per patch at t = 0 (length must equal PATCH_NAMES)
N_INFECTED_SEED <- c(5L, 5L, 1L, 1L, 1L, 1L, 1L)

# Ensemble size — set to 15 for the final BioInference run;
# use 1 during development for fast iteration
N_REPS <- 15L


#----------------------------------------------------------------#
# § 5  FILE PATHS & OUTPUT DIRECTORIES
#----------------------------------------------------------------#

# Input data (not included in repo — see data/README_data.md)
PATH_TEMP_ARS     <- here("data", "temp_ars.rds")
PATH_OD_COMPLETE  <- here("data", "od_complete.rds")
PATH_POP_INFO     <- here("data", "pop_info.rds")
PATH_HEALTHCARE   <- here("data", "healthcare_behaviour.rds")
PATH_TEMP_PER_ARS <- here("data", "temp_per_ARS.rds")
PATH_CAOP         <- here("data", "Continente_CAOP2024_1.gpkg")
PATH_ARS_EXCEL    <- here("data", "ARS_Portugal_Continental_Estrutura.xlsx")

# Output directories (created on first run by each script)
DIR_RESULTS <- here("results")
DIR_FIGURES <- here("figures")
DIR_FRAMES  <- here("frames")


#----------------------------------------------------------------#
# § 6  UTILITY HELPERS
#----------------------------------------------------------------#

#' check_data_file
#'
#' Stop with an informative message if a required data file is absent.
#' Called at the top of run_scenarios.R and plots.R so missing files
#' are reported clearly rather than producing a cryptic read error.
#'
#' @param path Character path to the file to check.
check_data_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required data file not found:\n  ", path, "\n",
      "See data/README_data.md for download and preparation instructions.\n",
      call. = FALSE
    )
  }
}
