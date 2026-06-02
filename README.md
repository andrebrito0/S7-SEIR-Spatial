 # S7-SEIR Spatial: Temperature-Dependent Disease Dynamics in a Metapopulation

> **Research question:** How does coupling temperature-dependent transmission with human mobility alter influenza dynamics in a SEIR framework?

This repository contains the full software implementation for the presentation at **BioInference 2026**. The model built with R's S7 object-oriented system, and it's a stochastic agent-based SEIR metapopulation model that simulates seasonal influenza transmission across municipalities in the Lisbon region, coupling temperature-driven transmissibility with inter-municipal commuting flows to isolate each mechanism's effect on epidemic dynamics.

---

## Overview

Respiratory infections such as influenza exhibit strong seasonal patterns driven by a combination of host resistance, virus survival, and behavioural changes. Rather than acting as direct causes, environmental factors — particularly temperature — act as **modulators of transmission intensity**.

This project implements a spatial SEIR agent-based model (ABM) that couples two mechanisms simultaneously:

- **Temperature-dependent transmission** — a logistic sigmoid maps daily mean temperature to a transmission scaling factor $\kappa(T)$, capturing the increased transmissibility of influenza at lower temperatures
- **Human mobility** — an origin–destination matrix $\theta$ derived from the 2021 Portuguese Census routes external contacts across municipalities, allowing infection to spread spatially

Four scenarios are compared in a stochastic ensemble to isolate each mechanism's contribution to outbreak dynamics.

---

## Repository Structure

```
.
├── config.R                  # Single source of truth: all parameters, paths, version guards
├── S7_SEIR_spatial_core.R    # Core model: S7 class definitions, simulation engine, plotting
├── run_scenarios.R           # Four-scenario ensemble run — main analysis script
├── plots.R                   # Data visualisation: seasonal dynamics, mobility & temp animations
├── example/
│   └── minimal_example.R    # Self-contained tutorial — no external data required
├── data/
│   └── README_data.md       # Data sources and download instructions
├── results/                  # Created on first run — simulation outputs
├── figures/                  # Created on first run — figures and animations
├── LICENSE
└── README.md
```

---

## Quickstart

### Requirements

- **R ≥ 4.3.0** — download from [r-project.org](https://www.r-project.org/)
- All packages are installed automatically when you source `config.R`

### Run the self-contained example (no data required, < 1 minute)

Clone the repo and run the minimal example directly — no data downloads needed:

```r
source("example/minimal_example.R")
```

This runs a complete 3-patch, 120-day simulation on synthetic data and produces SEIR curves, a prevalence heatmap, and an attack-rate comparison between the Baseline and Combined scenarios.

### Run the full analysis

The full analysis requires external data files (see [Data](#data) below). Once the `data/` folder is populated:

```r
# Step 1: run scenarios and save results
source("run_scenarios.R")

# Step 2: produce all figures (set RENDER_ANIMATIONS <- FALSE to skip GIFs)
source("plots.R")
```

All outputs are written to `results/` and `figures/`. The working directory must be the repository root; `here::here()` is used throughout so paths resolve correctly regardless of how R is launched.

---

## Model Description

### Compartments

Each agent occupies one of four states: **S**usceptible → **E**xposed → **I**nfectious → **R**ecovered. The metapopulation consists of **P = 7 patches** (municipalities), each with its own population of agents, temperature time series, and history of daily SEIR counts.

### Force of Infection

The daily infection risk for a susceptible agent in home patch *i* combines two pathways:

**Household (core) risk** — contacts with infected household members; not temperature-modulated (indoor, climate-controlled environments are assumed):

$$P(\text{core infection}) = 1 - (1 - \beta_{\text{core}})^{n_{\text{infected household members}}}$$

**Community (external) risk** — contacts distributed across destination patches *j* according to the mobility matrix $\theta$; temperature-modulated in each destination:

$$P(\text{ext. context}) = \prod_j \left(1 - \kappa(T_j) \cdot \beta_{\text{ext}} \cdot \pi_j \right)^{k_{ij}}$$

where $k_{ij} \sim \text{Poisson}(k_{\text{ext}} \cdot \theta_{ij})$ and $\pi_j$ is the effective infectious prevalence in patch *j* weighted by the inflow from all origins.

### Temperature Coupling

A logistic sigmoid maps daily mean temperature to a scaling factor $\kappa \in [y_{\min}, y_{\min} + \Delta y]$:

$$\kappa(T) = y_{\min} + \frac{\Delta y}{1 + e^{\,k_T (T - T_{\text{mid}})}}$$

Low temperature conditions $(T < T_{\text{mid}})$ push ${\kappa \rightarrow 1}$ (maximum transmission); Higher temperature conditions suppress it toward $y_{\min} = 0.20$. Parameters are set in `config.R` and documented there.

### Mobility Matrix

Φ is derived from pendular (commuting) movement flows in the **2021 Portuguese Census**. It is row-normalised so that $\theta_{ij}$ gives the proportion of time residents of patch *i* spend in patch *j*, with $\sum_j \theta_{ij} = 1$ for all *i*.

### Disease Parameters (Influenza H3N2, 2017–18)

| Parameter | Symbol | Value | Interpretation |
|---|---|---|---|
| Household transmission risk | β_core | 0.15 | Per infectious household contact |
| Community transmission risk | β_ext | 0.015 | Per infectious external contact |
| Incubation rate | σ | 0.526 day⁻¹ | Mean incubation 1.9 days |
| Recovery rate | γ | 0.244 day⁻¹ | Mean infectious period 4.1 days |
| κ floor | y_min | 0.20 | Minimum scaling at high temperature |
| κ range | Δy | 0.80 | Scaling range (max κ = 1.0) |
| Sigmoid steepness | k_T | 0.30 | Positive → colder = higher κ |
| Sigmoid inflection | T_mid | 10 °C | Temperature of half-maximal effect |

### Simulation Design

Fifteen stochastic replicates are run per scenario. Each replicate uses a fixed seed (`rep × 100`) set before any random call, including population construction, so results are exactly reproducible. The seed is stored as a column in every output data frame.

| Scenario | Mobility Φ | Temperature κ(T) |
|---|---|---|
| **Baseline** | Identity (no cross-patch movement) | κ ≡ 1 (no temperature effect) |
| **Mobility** | Census OD matrix | κ ≡ 1 |
| **Temperature** | Identity | Sigmoid κ(T) |
| **Combined** | Census OD matrix | Sigmoid κ(T) |

---

## S7 Class Architecture

The model is implemented using R's **S7** OOP system. The class hierarchy is:

```
disease
└── disease_spatial          adds sigmoid κ(T) parameters

agent                        individual with SEIR status + contact lists

population                   single-patch agent pool (used by step_sim)

patch                        spatial unit storing aggregate SEIR counts + history

metapopulation               multi-patch system: patches + flat agent pool + Φ + temp_series
```

Key functions:

| Function | Description |
|---|---|
| `create_metapopulation()` | Build a P-patch system from scratch |
| `step_metapop()` | Advance the metapopulation by one day (spatial FOI) |
| `temp_scaling()` | Evaluate κ(T) for a given temperature and disease object |
| `summary(metapopulation)` | Print current SEIR counts per patch |
| `plot_metapop_curves()` | SEIR time series, one panel per patch |
| `plot_prevalence_heatmap()` | Patch × day prevalence intensity map |
| `plot_temp_scaling_curve()` | Visualise the sigmoid κ(T) function |

---

## Data

The full analysis uses three external datasets. They are not included in the repository due to size and licensing. Place the processed files in `data/` before running `run_scenarios.R` or `plots.R`. Full instructions are in [`data/README_data.md`](data/README_data.md).

| File | Source | Description |
|---|---|---|
| `data/temp_ars.rds` | [Copernicus E-OBS](https://cds.climate.copernicus.eu/) | Daily mean temperature per municipality, 2017–2018 |
| `data/od_complete.rds` | [INE Census 2021](https://www.ine.pt/) | Origin–destination pendular mobility flows |
| `data/pop_info.rds` | [INE Census 2021](https://www.ine.pt/) | Municipality population sizes |
| `data/healthcare_behaviour.rds` | [DGS / SNS](https://transparencia.sns.gov.pt/) | Respiratory emergency visits per health region |
| `data/temp_per_ARS.rds` | [Copernicus E-OBS](https://cds.climate.copernicus.eu/) | Temperature aggregated to ARS level (for plots) |
| `data/ARS_Portugal_Continental_Estrutura.xlsx` | [SNS](https://www.sns.gov.pt/) | Municipality–ARS mapping |
| `data/lvt_map_proj` | [CAOP](https://www.dgterritorio.gov.pt/atividades/cartografia/cartografia-tematica/caop) | Polygon Projection |

> **No data needed for the minimal example.** `example/minimal_example.R` is fully self-contained and uses synthetic temperature and mobility data.

---

## Outputs

### From `run_scenarios.R`

| File | Description |
|---|---|
| `results/ensemble_baseline.rds` | Long-format data frame: all reps, Baseline scenario |
| `results/ensemble_mobility.rds` | Long-format data frame: all reps, Mobility scenario |
| `results/ensemble_temperature.rds` | Long-format data frame: all reps, Temperature scenario |
| `results/ensemble_combined.rds` | Long-format data frame: all reps, Combined scenario |
| `results/ensemble_4scenarios.rds` | All scenarios + summaries in one list |
| `results/epi_summary_4scenarios.csv` | Attack rate and peak timing summary table |
| `results/processed_inputs.rds` | Processed temperature matrix and commuting matrix (reused by `plots.R`) |
| `results/plot1_global_ribbon.png` | Global epidemic curves with 10–90% stochastic band |
| `results/plot2_patch_curves_by_scenario.png` | Per-patch curves faceted by scenario |
| `results/plot3_attack_rate.png` | Final attack rate distributions by scenario and patch |
| `results/plot4_peak_timing.png` | Day of peak incidence by scenario and patch |

### From `plots.R`

| File | Description |
|---|---|
| `figures/seasonal_plot.pdf` | 7-day moving average of respiratory emergency visits with temperature overlay |
| `figures/od_flows_temporal.gif` | Animated daily commuting flows across LVT municipalities |
| `figures/temp_coupling.gif` | Animated draw of the sigmoid κ(T) curve |
| `frames/od/frame-N.png` | Individual frames extracted for LaTeX `\animategraphics` |
| `frames/temp/frame-N.png` | Individual frames extracted for LaTeX `\animategraphics` |

---

## Reproducing the Analysis

The exact results presented at BioInference 2026 can be reproduced by:

1. Cloning this repository
2. Installing R ≥ 4.3.0
3. Placing the data files in `data/` (see [`data/README_data.md`](data/README_data.md))
4. Running, **in order**, from the repository root:

```r
source("run_scenarios.R")   # ~30–60 min depending on hardware
source("plots.R")           # ~10–15 min including GIF rendering
```

Stochastic replicates use seeds `rep × 100` (100, 200, …, 1500), set before any random call including population initialisation. The `seed` column in every output data frame records which seed produced each replicate.

To reproduce a single replicate exactly:

```r
source("config.R")
source("S7_SEIR_spatial_core.R")

set.seed(100)   # replicate 1
mpop <- create_metapopulation(
  patch_sizes = PATCH_SIZES,
  patch_names = PATCH_NAMES,
  n_infected  = N_INFECTED_SEED,
  k_core_pop  = K_CORE_POP,
  k_ext_pop   = K_EXT_POP,
  disease_obj = disease_spatial(
    name = "Influenza (H3N2, 2017-18)",
    beta_core = BETA_CORE, beta_ext = BETA_EXT,
    sigma = SIGMA, gamma = GAMMA,
    y_min = Y_MIN, delta_y = DELTA_Y,
    k_temp = K_TEMP, T_mid = T_MID
  )
)
```

---

## Package Dependencies

All packages are installed automatically by `config.R` if not already present. For an exact snapshot of the package versions used to produce the presented results, see `renv.lock`.

| Package | Role |
|---|---|
| `S7` ≥ 0.2.0 | S7 OOP system — core model class definitions |
| `dplyr`, `tidyr` | Data wrangling |
| `ggplot2` | Static plots |
| `zoo` | Rolling averages and NA interpolation |
| `lubridate` | Date handling |
| `reshape2` | OD matrix casting (`dcast`) |
| `readxl` | ARS–municipality mapping spreadsheet |
| `sf` | Spatial data (municipality boundaries) |
| `gganimate`, `transformr`, `gifski` | GIF animation rendering |
| `magick` | GIF frame extraction for LaTeX |
| `here` | Reproducible relative file paths |

---

## Citation

If you use this software, please cite:

```
Brito, A. (2026). S7-SEIR Spatial: Temperature-dependent influenza dynamics
in a metapopulation [Software]. Presented at BioInference 2026.
GitHub: https://github.com/andrebrito0/s7-seir-spatial
```

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for full terms. You are free to use, modify, and distribute this software with attribution.

---

## Author

**André Brito**
Centre for Mathematics and Applications (NOVA FCT)
📧 anm.brito@campus.fct.unl.pt
🐙 [github.com/andrebrito0](https://github.com/andrebrito0)

*This work was supported by national funds through FCT – Fundação para a Ciência e a Tecnologia under projects UIDB/00297/2020, UIDP/00297/2020 (Center for Mathematics and Applications), and 2024.00664.BDANA.*
