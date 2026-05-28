#================================================================#
#  S7-SEIR  —  Unified Agent-Based Epidemic Model
#
#  STRUCTURE
#  ─────────────────────────────────────────────────────────────
#  SECTION 0  Version & dependency checks
#
#  SECTION 1  Class definitions
#               1a. disease          – base pathogen parameters
#               1b. disease_spatial  – extends disease with sigmoid
#                                      temperature response
#               1c. agent            – individual with SEIR status
#               1d. population       – single-patch agent pool
#               1e. patch            – spatial unit (SEIR counts)
#               1f. metapopulation   – multi-patch system
#
#  SECTION 2  Core methods  (single-patch model)
#               update_agent        – E→I, I→R transitions
#               assign_core_group   – household structure
#               assign_infected     – seed infections
#               step_sim            – advance population one day
#
#  SECTION 3  Constructors
#               create_population   – build a single-patch population
#               create_metapopulation – build a multi-patch system
#
#  SECTION 4  Spatial step
#               step_metapop        – advance metapopulation one day
#                                     (mobility + temperature FOI)
#
#  SECTION 5  Summary methods
#               summary(disease)
#               summary(disease_spatial)
#               summary(population)
#               summary(metapopulation)
#
#  SECTION 6  Plotting helpers
#               plot_population_structure
#               plot_disease_profile
#               plot_temp_scaling_curve
#               plot_metapop_curves
#               plot_prevalence_heatmap
#
#  SECTION 7  Example usage  (wrapped in if(FALSE) — safe to source)
#
#  Temperature FOI
#  ─────────────────────────────────────────────────────────────
#  Transmission is modulated by a sigmoid function of daily
#  mean temperature, capturing the monotone increase seen in
#  cold-adapted respiratory pathogens:
#
#    κ(T) = y_min + Δy / (1 + exp(k_temp · (T − T_mid)))
#
#    y_min  : floor scaling at high temperature  (>0, e.g. 0.2)
#    Δy     : range  (ceiling = y_min + Δy,      e.g. 0.8 → max 1.0)
#    k_temp : steepness; positive → colder = higher κ
#    T_mid  : inflection temperature (°C)
#
#  Core (household) contacts are not temperature-modulated —
#  indoor household exposure is assumed to be climate-controlled.
#================================================================#


#----------------------------------------------------------------#
# SECTION 0  VERSION & DEPENDENCY CHECKS
#----------------------------------------------------------------#

# Fail immediately with a clear message if the R version or key
# packages are too old, rather than producing cryptic errors later.
if (getRversion() < "4.3.0") {
  stop(
    "R >= 4.3.0 is required (current: ", getRversion(), ").\n",
    "Download the latest R from https://www.r-project.org/\n",
    call. = FALSE
  )
}

required_pkgs <- c("S7")
missing_pkgs  <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

if (packageVersion("S7") < "0.2.0") {
  stop(
    "S7 >= 0.2.0 is required (current: ", packageVersion("S7"), ").\n",
    "Run: install.packages('S7')\n",
    call. = FALSE
  )
}

library(S7)


#----------------------------------------------------------------#
# SECTION 1  CLASS DEFINITIONS
#----------------------------------------------------------------#

# ── 1a. Base disease ─────────────────────────────────────────────
disease <- new_class(
  "disease",
  properties = list(
    name      = class_character,
    beta_core = class_numeric,   # Transmission risk per infectious core contact
    beta_ext  = class_numeric,   # Transmission risk per infectious external contact
    sigma     = class_numeric,   # Rate E→I  (1 / mean incubation days)
    gamma     = class_numeric    # Rate I→R  (1 / mean infectious days)
  )
)

# ── 1b. Spatially-aware disease — sigmoid temperature response ───
#
#  Inherits all base disease properties and adds four parameters
#  that define the sigmoid thermal response curve κ(T).
#  All downstream spatial code expects disease_spatial; the base
#  disease class is left untouched so single-patch code is unaffected.

disease_spatial <- new_class(
  "disease_spatial",
  parent     = disease,
  properties = list(
    y_min   = class_numeric,     # Minimum κ at high temperature
    delta_y = class_numeric,     # Range of κ  (max κ = y_min + delta_y)
    k_temp  = class_numeric,     # Steepness; positive → colder = more transmission
    T_mid   = class_numeric,      # Inflection temperature (°C)
    use_temp = class_logical
  )
)

#' temp_scaling  —  evaluate the sigmoid thermal kernel κ(T)
#'
#' @param T_today  Scalar daily mean temperature (°C)
#' @param dis      A disease_spatial object
#' @return         Scalar in [y_min, y_min + delta_y]
temp_scaling <- function(T_today, dis) {
  if (!dis@use_temp) return(1.0)
  dis@y_min + dis@delta_y / (1 + exp(dis@k_temp * (T_today - dis@T_mid)))
}

# ── 1c. Agent ────────────────────────────────────────────────────
agent <- new_class(
  "agent",
  properties = list(
    id            = class_integer,
    age           = class_integer,
    status        = class_character,   # "S", "E", "I", or "R"
    timer         = class_numeric,     # Days remaining in current state
    core_contacts = class_integer,     # Global agent IDs of household members
    ext_contacts  = class_integer      # Mean daily external contacts (λ)
  )
)

# ── 1d. Population (single patch) ────────────────────────────────
population <- new_class(
  "population",
  properties = list(
    agents     = class_list,
    disease    = disease,
    k_core_pop = class_integer,
    k_ext_pop  = class_integer,
    history    = class_data.frame
  )
)

# ── 1e. Patch ────────────────────────────────────────────────────
#  Stores aggregate SEIR counts for one spatial unit.
#  Agents are not duplicated per patch; they live in a flat list
#  on the metapopulation and are mapped via patch_ids.

patch <- new_class(
  "patch",
  properties = list(
    id       = class_integer,
    name     = class_character,
    pop_size = class_integer,
    S        = class_integer,
    E        = class_integer,
    I        = class_integer,
    R        = class_integer,
    history  = class_data.frame   # One row per simulated day
  )
)

#' patch_prevalence  —  infectious fraction of patch residents
patch_prevalence <- function(p) {
  total <- p@S + p@E + p@I + p@R
  if (total == 0L) return(0)
  p@I / total
}

# ── 1f. Metapopulation ───────────────────────────────────────────
metapopulation <- new_class(
  "metapopulation",
  properties = list(
    patches     = class_list,          # list of patch objects
    agents      = class_list,          # flat list of ALL agents across patches
    patch_ids   = class_integer,       # length = n_agents; home patch per agent
    disease     = disease_spatial,
    # Mobility matrix: Phi[i,j] = fraction of patch-i residents whose
    # external contacts occur in patch j on a typical day.
    # Row-stochastic: rowSums(Phi) == 1.
    Phi         = class_any,           # matrix [n_patches × n_patches]
    # Daily mean temperatures: matrix [n_patches × n_days].
    # Column index is clamped at ncol if simulation runs longer than the series.
    temp_series = class_any,           # matrix [n_patches × n_days]
    day         = class_integer,       # Current day counter (starts at 1)
    k_core_pop  = class_integer,
    k_ext_pop   = class_integer,
    history     = class_data.frame     # Aggregated SEIR across all patches
  )
)


#----------------------------------------------------------------#
# SECTION 2  CORE METHODS  (single-patch model)
#----------------------------------------------------------------#

# ── Biological state transitions ─────────────────────────────────
update_agent <- new_generic("update_agent", "agent")

method(update_agent, agent) <- function(agent, pop, dt = 1) {
  dis <- pop@disease
  
  if (agent@status %in% c("E", "I")) {
    agent@timer <- agent@timer - dt
  }
  
  if (agent@status == "E" && agent@timer <= 0) {
    agent@status <- "I"
    agent@timer  <- rexp(1, rate = dis@gamma)
  } else if (agent@status == "I" && agent@timer <= 0) {
    agent@status <- "R"
    agent@timer  <- Inf
  }
  agent
}

# ── Household structure assignment ───────────────────────────────
assign_core_group <- new_generic("assign_core_group", "pop")

method(assign_core_group, population) <- function(pop) {
  indices   <- sample(seq_along(pop@agents))
  mean_size <- pop@k_core_pop
  
  while (length(indices) > 0) {
    size        <- rpois(1, lambda = mean_size - 1) + 1
    actual_size <- min(size, length(indices))
    group_idx   <- indices[1:actual_size]
    
    if (length(group_idx) > 1) {
      for (idx in group_idx) {
        pop@agents[[idx]]@core_contacts <- as.integer(setdiff(group_idx, idx))
      }
    }
    indices <- indices[-(1:actual_size)]
  }
  pop
}

# ── Infection seeding ─────────────────────────────────────────────
assign_infected <- new_generic("assign_infected", "pop")

method(assign_infected, population) <- function(pop, n_infected) {
  idx <- sample(seq_along(pop@agents), n_infected)
  for (i in idx) {
    pop@agents[[i]]@status <- "I"
    pop@agents[[i]]@timer  <- rexp(1, rate = pop@disease@gamma)
  }
  pop
}

# ── Single-patch simulation step ──────────────────────────────────
step_sim <- new_generic("step_sim", "pop")

method(step_sim, population) <- function(pop) {
  dis          <- pop@disease
  agents       <- pop@agents
  all_statuses <- vapply(agents, \(a) a@status, character(1))
  prevalence   <- sum(all_statuses == "I") / length(all_statuses)
  
  # Record current state before updating
  counts      <- as.data.frame(
    t(as.matrix(table(factor(all_statuses, levels = c("S", "E", "I", "R")))))
  )
  pop@history <- rbind(pop@history, counts)
  
  pop@agents <- lapply(seq_along(agents), function(i) {
    a <- agents[[i]]
    
    if (a@status == "S") {
      # Core (household) risk
      n_core_I  <- sum(all_statuses[a@core_contacts] == "I")
      risk_core <- 1 - (1 - dis@beta_core)^n_core_I
      
      # External stochastic risk
      daily_k  <- rpois(1, lambda = a@ext_contacts)
      risk_ext <- 1 - (1 - (prevalence * dis@beta_ext))^daily_k
      
      if (runif(1) < 1 - (1 - risk_core) * (1 - risk_ext)) {
        a@status <- "E"
        a@timer  <- rexp(1, rate = dis@sigma)
      }
    }
    update_agent(a, pop)
  })
  pop
}


#----------------------------------------------------------------#
# SECTION 3  CONSTRUCTORS
#----------------------------------------------------------------#

# ── Single-patch population ───────────────────────────────────────
#' create_population
#'
#' @param n           Total number of agents
#' @param n_infected  Number of agents seeded as infectious (I)
#' @param k_core_pop  Target mean household size
#' @param k_ext_pop   Target mean daily external contacts
#' @param disease_obj A disease or disease_spatial object
#' @return A population S7 object ready for step_sim()

create_population <- function(n, n_infected, k_core_pop, k_ext_pop, disease_obj) {
  agents <- lapply(seq_len(n), \(i) agent(
    id           = as.integer(i),
    age          = as.integer(runif(1, 0, 80)),
    status       = "S",
    timer        = 0,
    ext_contacts = as.integer(rgamma(1, shape = 2, rate = 2 / k_ext_pop))
  ))
  
  pop <- population(
    agents     = agents,
    disease    = disease_obj,
    k_core_pop = as.integer(k_core_pop),
    k_ext_pop  = as.integer(k_ext_pop),
    history    = data.frame()
  )
  
  pop |> assign_core_group() |> assign_infected(n_infected)
}


# ── Metapopulation ────────────────────────────────────────────────
#' create_metapopulation
#'
#' Builds a multi-patch system by constructing one sub-population per
#' patch via create_population(), re-indexing agent IDs globally, and
#' assembling the metapopulation object.
#'
#' @param patch_sizes   Integer vector [n_patches] — residents per patch
#' @param patch_names   Character vector of patch labels (auto-generated if NULL)
#' @param n_infected    Integer vector [n_patches] of seeds, OR a single integer
#'                      applied to patch 1 only (all others start clean)
#' @param k_core_pop    Mean household size (applied uniformly)
#' @param k_ext_pop     Mean daily external contacts (applied uniformly)
#' @param disease_obj   A disease_spatial object
#' @param Phi           Row-stochastic mobility matrix [n_patches × n_patches];
#'                      defaults to identity matrix (no cross-patch movement)
#' @param temp_series   Matrix [n_patches × n_days] of daily mean temperatures;
#'                      a plain vector is replicated across all patches;
#'                      defaults to T_mid everywhere (κ = midpoint of sigmoid)
#' @return A metapopulation S7 object ready for step_metapop()

create_metapopulation <- function(patch_sizes,
                                  patch_names  = NULL,
                                  n_infected   = 5L,
                                  k_core_pop   = 3L,
                                  k_ext_pop    = 8L,
                                  disease_obj,
                                  Phi          = NULL,
                                  temp_series  = NULL) {
  
  n_patches <- length(patch_sizes)
  
  # Default patch labels
  if (is.null(patch_names)) patch_names <- paste0("Patch_", seq_len(n_patches))
  
  # Default mobility: identity (everyone stays in their home patch)
  if (is.null(Phi)) Phi <- diag(n_patches)
  stopifnot(nrow(Phi) == n_patches, ncol(Phi) == n_patches)
  row_sums              <- rowSums(Phi)
  row_sums[row_sums == 0] <- 1
  Phi                   <- Phi / row_sums    # enforce row-stochastic
  
  # Default temperature: T_mid everywhere → κ = y_min + Δy/2 (midpoint)
  if (is.null(temp_series)) {
    temp_series <- matrix(disease_obj@T_mid, nrow = n_patches, ncol = 1)
  }
  if (is.vector(temp_series)) {
    temp_series <- matrix(rep(temp_series, each = n_patches), nrow = n_patches)
  }
  
  # Seed vector
  if (length(n_infected) == 1L) {
    seeds <- c(as.integer(n_infected), rep(0L, n_patches - 1L))
  } else {
    seeds <- as.integer(n_infected)
  }
  stopifnot(length(seeds) == n_patches)
  
  # Build agents patch by patch, concatenating into a global flat list
  all_agents       <- list()
  all_patch_ids    <- integer(0)
  global_id_offset <- 0L
  patch_list       <- vector("list", n_patches)
  
  for (p in seq_len(n_patches)) {
    n_p <- as.integer(patch_sizes[p])
    
    sub_pop <- create_population(
      n           = n_p,
      n_infected  = seeds[p],
      k_core_pop  = k_core_pop,
      k_ext_pop   = k_ext_pop,
      disease_obj = disease_obj
    )
    
    # Re-index IDs to be globally unique
    for (i in seq_along(sub_pop@agents)) {
      sub_pop@agents[[i]]@id <- as.integer(global_id_offset + i)
    }
    
    all_agents    <- c(all_agents, sub_pop@agents)
    all_patch_ids <- c(all_patch_ids, rep(as.integer(p), n_p))
    global_id_offset <- global_id_offset + n_p
    
    statuses <- vapply(sub_pop@agents, \(a) a@status, character(1))
    patch_list[[p]] <- patch(
      id       = as.integer(p),
      name     = patch_names[p],
      pop_size = n_p,
      S        = as.integer(sum(statuses == "S")),
      E        = as.integer(sum(statuses == "E")),
      I        = as.integer(sum(statuses == "I")),
      R        = as.integer(sum(statuses == "R")),
      history  = data.frame()
    )
  }
  
  metapopulation(
    patches     = patch_list,
    agents      = all_agents,
    patch_ids   = all_patch_ids,
    disease     = disease_obj,
    Phi         = Phi,
    temp_series = temp_series,
    day         = 1L,
    k_core_pop  = as.integer(k_core_pop),
    k_ext_pop   = as.integer(k_ext_pop),
    history     = data.frame()
  )
}


#----------------------------------------------------------------#
# SECTION 4  SPATIAL STEP
#----------------------------------------------------------------#

#' step_metapop  —  advance the metapopulation by one day
#'
#' Each day proceeds in five sub-steps:
#'
#'  (a) Temperature scaling
#'      κ_j = temp_scaling(T_j(d), disease)  for each patch j.
#'      Uses the sigmoid function; κ rises as temperature falls.
#'
#'  (b) Effective prevalence in each destination patch
#'      People physically present in patch j on day d come from all
#'      origin patches i in proportion to Phi[i,j].  The infectious
#'      fraction among all people in j is:
#'
#'        π_j = Σᵢ(Φ[i,j] · Iᵢ) / Σᵢ(Φ[i,j] · Nᵢ)
#'
#'  (c) Agent-level infection (susceptibles only)
#'      Each susceptible in home patch r distributes their external
#'      contacts across destination patches:
#'        k_ij ~ Poisson(ext_contacts · Φ[r,j])
#'      Escape probability from patch j:
#'        P_escape_j = (1 − κ_j · β_ext · π_j)^k_ij
#'      Combined survival across all destinations (independence):
#'        P_escape = Π_j P_escape_j
#'      Core (household) risk is local and not temperature-modulated.
#'
#'  (d) State transitions  E→I, I→R via update_agent()
#'
#'  (e) Refresh patch SEIR counts and append to patch histories.

step_metapop <- new_generic("step_metapop", "mpop")

method(step_metapop, metapopulation) <- function(mpop) {
  
  dis       <- mpop@disease
  agents    <- mpop@agents
  patch_ids <- mpop@patch_ids
  Phi       <- mpop@Phi
  n_patches <- length(mpop@patches)
  day       <- mpop@day
  
  # ── (a) Temperature scaling per destination patch ─────────────
  temp_col <- min(day, ncol(mpop@temp_series))   # clamp to series length
  kappa    <- vapply(seq_len(n_patches), function(j) {
    temp_scaling(mpop@temp_series[j, temp_col], dis)
  }, numeric(1))
  
  # ── (b) Effective infectious prevalence in each patch ─────────
  all_statuses <- vapply(agents, \(a) a@status, character(1))
  
  I_vec <- vapply(seq_len(n_patches), \(p)
                  as.integer(sum(all_statuses[patch_ids == p] == "I")), integer(1))
  N_vec <- vapply(seq_len(n_patches), \(p)
                  sum(patch_ids == p), integer(1))
  
  pi_dest <- numeric(n_patches)
  for (j in seq_len(n_patches)) {
    denom       <- sum(Phi[, j] * N_vec)
    pi_dest[j]  <- if (denom > 0) sum(Phi[, j] * I_vec) / denom else 0
  }
  
  # ── (c & d) Agent-level infection + state transitions ─────────
  # A minimal population shell is constructed once to satisfy
  # update_agent()'s signature (it only reads @disease, @k_core_pop,
  # @k_ext_pop — never the agent list — so an empty agents slot is safe).
  pop_shell <- population(
    agents     = list(),
    disease    = dis,
    k_core_pop = mpop@k_core_pop,
    k_ext_pop  = mpop@k_ext_pop,
    history    = data.frame()
  )
  
  mpop@agents <- lapply(seq_along(agents), function(idx) {
    a  <- agents[[idx]]
    ri <- patch_ids[[idx]]
    
    if (a@status == "S") {
      
      # Core risk — local, no temperature effect
      n_core_I  <- sum(all_statuses[a@core_contacts] == "I")
      risk_core <- 1 - (1 - dis@beta_core)^n_core_I
      
      # External risk — distributed across destination patches
      survival_ext <- 1.0
      for (j in seq_len(n_patches)) {
        phi_ij <- Phi[ri, j]
        if (phi_ij < 1e-9) next
        k_ij  <- rpois(1, lambda = a@ext_contacts * phi_ij)
        if (k_ij == 0L) next
        p_esc <- (1 - kappa[j] * dis@beta_ext * pi_dest[j])^k_ij
        survival_ext <- survival_ext * max(p_esc, 0)
      }
      risk_ext <- 1 - survival_ext
      
      # Exposure event
      if (runif(1) < 1 - (1 - risk_core) * (1 - risk_ext)) {
        a@status <- "E"
        a@timer  <- rexp(1, rate = dis@sigma)
      }
    }
    
    update_agent(a, pop_shell)
  })
  
  # ── (e) Refresh patch counts and histories ─────────────────────
  new_statuses <- vapply(mpop@agents, \(a) a@status, character(1))
  
  for (p in seq_len(n_patches)) {
    mask <- patch_ids == p
    s_p  <- new_statuses[mask]
    S_p  <- as.integer(sum(s_p == "S"))
    E_p  <- as.integer(sum(s_p == "E"))
    I_p  <- as.integer(sum(s_p == "I"))
    R_p  <- as.integer(sum(s_p == "R"))
    
    mpop@patches[[p]]@history <- rbind(
      mpop@patches[[p]]@history,
      data.frame(day = day, patch = p,
                 patch_name = mpop@patches[[p]]@name,
                 S = S_p, E = E_p, I = I_p, R = R_p)
    )
    mpop@patches[[p]]@S <- S_p
    mpop@patches[[p]]@E <- E_p
    mpop@patches[[p]]@I <- I_p
    mpop@patches[[p]]@R <- R_p
  }
  
  mpop@history <- rbind(mpop@history, data.frame(
    day = day,
    S   = as.integer(sum(new_statuses == "S")),
    E   = as.integer(sum(new_statuses == "E")),
    I   = as.integer(sum(new_statuses == "I")),
    R   = as.integer(sum(new_statuses == "R"))
  ))
  mpop@day <- day + 1L
  
  mpop
}


#----------------------------------------------------------------#
# SECTION 5  SUMMARY METHODS
#----------------------------------------------------------------#

method(summary, disease) <- function(object) {
  cat(sprintf("--- Disease Profile: %s ---\n", object@name))
  cat("Biological timings:\n")
  cat(sprintf("  Avg incubation (E→I) : %.1f days\n", 1 / object@sigma))
  cat(sprintf("  Avg infectious period: %.1f days\n", 1 / object@gamma))
  cat("Transmission risk (β):\n")
  cat(sprintf("  Core group (close)   : %.1f%%\n", object@beta_core * 100))
  cat(sprintf("  External (casual)    : %.1f%%\n", object@beta_ext  * 100))
  cat("-----------------------------------\n")
  invisible(list(incubation = 1 / object@sigma, infectious = 1 / object@gamma))
}

method(summary, disease_spatial) <- function(object) {
  cat(sprintf("--- Disease Profile: %s ---\n", object@name))
  cat("Biological timings:\n")
  cat(sprintf("  Avg incubation (E→I) : %.1f days\n", 1 / object@sigma))
  cat(sprintf("  Avg infectious period: %.1f days\n", 1 / object@gamma))
  cat("Transmission risk (β):\n")
  cat(sprintf("  Core group (close)   : %.1f%%\n", object@beta_core * 100))
  cat(sprintf("  External (casual)    : %.1f%%\n", object@beta_ext  * 100))
  if (object@use_temp) {
    cat("Sigmoid temperature response  κ(T) = y_min + Δy / (1 + exp(k·(T−T_mid))):\n")
    cat(sprintf("  T_mid   : %.1f °C  (inflection point)\n", object@T_mid))
    cat(sprintf("  k_temp  : %.2f    (steepness)\n",          object@k_temp))
    cat(sprintf("  y_min   : %.2f    (floor κ)\n",            object@y_min))
    cat(sprintf("  Δy      : %.2f    (range; ceiling = %.2f)\n",
                object@delta_y, object@y_min + object@delta_y))
  } else {
    cat("Temperature coupling         : DISABLED (κ ≡ 1.0)\n")
  }
  cat("-----------------------------------\n")
  invisible(list(incubation = 1 / object@sigma, infectious = 1 / object@gamma))
}

method(summary, population) <- function(object) {
  n_total    <- length(object@agents)
  ages       <- vapply(object@agents, \(a) a@age,                        integer(1))
  hh_sizes   <- vapply(object@agents, \(a) length(a@core_contacts) + 1L, numeric(1))
  k_ext_vals <- vapply(object@agents, \(a) a@ext_contacts,               integer(1))
  
  ss_threshold <- 2 * object@k_ext_pop
  ss_indices   <- which(k_ext_vals >= ss_threshold)
  n_ss         <- length(ss_indices)
  
  cat("--- S7-SEIR Population Structural Summary ---\n")
  cat(sprintf("Total agents             : %d\n", n_total))
  cat(sprintf("Potential super-spreaders: %d (%.1f%%)\n",
              n_ss, (n_ss / n_total) * 100))
  cat(sprintf("  Threshold              : > %d ext. contacts/day\n", ss_threshold))
  cat(sprintf("  Max individual k_ext   : %d\n", max(k_ext_vals)))
  cat("Social structure:\n")
  cat(sprintf("  Avg household size     : %.2f\n", mean(hh_sizes)))
  cat(sprintf("  Avg community contacts : %.2f\n", mean(k_ext_vals)))
  cat("---------------------------------------------\n")
  invisible(list(n_ss = n_ss, ss_indices = ss_indices))
}

method(summary, metapopulation) <- function(object) {
  n_patches <- length(object@patches)
  n_agents  <- length(object@agents)
  cat("=== S7-SEIR Spatial Metapopulation ===\n")
  cat(sprintf("  Patches  : %d\n", n_patches))
  cat(sprintf("  Agents   : %d\n", n_agents))
  cat(sprintf("  Disease  : %s\n", object@disease@name))
  cat(sprintf("  Days run : %d\n", object@day - 1L))
  cat("\n--- Current Patch Status ---\n")
  cat(sprintf("  %-20s  %6s  %6s  %6s  %6s  %8s\n",
              "Patch", "S", "E", "I", "R", "Prev(%)"))
  for (p in object@patches) {
    prev <- if (p@pop_size > 0) round(p@I / p@pop_size * 100, 2) else 0
    cat(sprintf("  %-20s  %6d  %6d  %6d  %6d  %7.2f%%\n",
                p@name, p@S, p@E, p@I, p@R, prev))
  }
  cat("=======================================\n")
  invisible(object)
}


#----------------------------------------------------------------#
# SECTION 6  PLOTTING HELPERS
#----------------------------------------------------------------#

#' plot_population_structure
#'
#' Four-panel diagnostic: age distribution, household sizes,
#' external contact heterogeneity, and sociality by age.
#'
#' @param pop A population S7 object (after assign_core_group())
plot_population_structure <- function(pop) {
  ages       <- vapply(pop@agents, \(a) a@age,                        integer(1))
  hh_sizes   <- vapply(pop@agents, \(a) length(a@core_contacts) + 1L, numeric(1))
  k_ext_vals <- vapply(pop@agents, \(a) a@ext_contacts,               integer(1))
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  hist(ages, breaks = 15, col = "darkseagreen3", border = "white",
       main = "Demographic profile", xlab = "Age", ylab = "Agents")
  
  barplot(table(hh_sizes), col = "steelblue", border = "white",
          main = "Core group (household) sizes",
          xlab = "Members", ylab = "Frequency")
  
  hist(k_ext_vals, breaks = 20, col = "plum3", border = "white",
       main = "External sociality (k_ext)",
       xlab = "Mean daily contacts", ylab = "Frequency")
  abline(v = pop@k_ext_pop, col = "red", lwd = 2, lty = 2)
  
  plot(ages, k_ext_vals, pch = 16, col = rgb(0.1, 0.1, 0.1, 0.2),
       main = "Sociality by age", xlab = "Age", ylab = "Ext. contact rate")
  
  par(mfrow = c(1, 1))
}


#' plot_disease_profile
#'
#' Two panels: transition time densities and theoretical R0 vs. sociality.
#'
#' @param dis A disease or disease_spatial object
plot_disease_profile <- function(dis) {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  
  days_seq <- seq(0, 30, length.out = 200)
  matplot(days_seq,
          cbind(dexp(days_seq, rate = dis@sigma),
                dexp(days_seq, rate = dis@gamma)),
          type = "l", lty = 1, lwd = 2, col = c("orange", "red"),
          main = "Transition time densities",
          xlab = "Days", ylab = "Density")
  legend("topright", legend = c("E→I", "I→R"),
         col = c("orange", "red"), lty = 1, bty = "n")
  grid()
  
  k_range  <- 1:20
  r0_vals  <- ((2 * dis@beta_core) + (k_range * dis@beta_ext)) * (1 / dis@gamma)
  matplot(c(1, k_range), c(0, r0_vals), type = "b", pch = 16,
          col = "darkorchid",
          main = expression(bold("Theoretical ") ~ R[0] ~ bold(" vs. sociality")),
          xlab = "Avg. daily external contacts (k)", ylab = "R0",
          xaxt = "n")
  axis(1, at = k_range)
  abline(v = k_range, col = "lightgray", lty = "dotted")
  abline(h = seq(0, max(r0_vals) + 1, by = 0.5), col = "lightgray", lty = "dotted")
  abline(h = 1, col = "red", lty = 2, lwd = 2)
  text(min(k_range), 1.15,
       labels = bquote(bold("Epidemic threshold ") ~ R[0] == 1),
       col = "red", pos = 4, font = 2, cex = 0.7)
  
  par(mfrow = c(1, 1))
}


#' plot_beta_R0_profile
#'
#' Two-panel diagnostic for a disease_spatial object showing how
#' temperature drives the effective transmission rate β(T) and the
#' implied basic reproduction number R₀(T) across a range of
#' daily mean temperatures.
#'
#' R₀ is approximated as the expected secondary cases generated by
#' one infectious agent over the infectious period, combining core
#' (household) and external contacts:
#'
#'   R₀(T) = [β_core · k_core  +  κ(T) · β_ext · k_ext] / γ
#'
#' where k_core and k_ext are the mean contact counts supplied
#' (defaulting to the values used in the simulation).
#'
#' @param dis       A disease_spatial object
#' @param k_core    Mean household contacts per agent (default 2,
#'                  i.e. household size 3 minus self)
#' @param k_ext     Mean daily external contacts (default 8)
#' @param T_range   Numeric vector c(T_min, T_max) in °C (default -5 to 30)
#' @param n_points  Resolution of the temperature grid (default 300)
#' @param show_Tmid Logical — draw a vertical reference at T_mid (default TRUE)

plot_beta_R0_profile <- function(dis,
                                 k_core    = 2L,
                                 k_ext     = 8L,
                                 T_range   = c(5, 30),
                                 n_points  = 300,
                                 show_Tmid = TRUE) {
  
  use_temp <- inherits(dis, "disease_spatial") && dis@use_temp
  
  T_seq <- seq(T_range[1], T_range[2], length.out = n_points)
  kappa <- if (use_temp) {
    vapply(T_seq, \(T) temp_scaling(T, dis), numeric(1))
  } else {
    rep(1, n_points)
  }
  
  beta_core_contrib <- dis@beta_core * k_core
  beta_ext_contrib  <- kappa * dis@beta_ext * k_ext
  beta_eff          <- (beta_core_contrib + beta_ext_contrib) / (k_core + k_ext)
  R0                <- (beta_core_contrib + beta_ext_contrib) / dis@gamma
  
  # ── Floor / ceiling references ───────────────────────────────────
  if (use_temp) {
    kappa_floor   <- dis@y_min
    kappa_ceiling <- dis@y_min + dis@delta_y
  } else {
    kappa_floor   <- 1
    kappa_ceiling <- 1
  }
  beta_floor   <- (beta_core_contrib + kappa_floor   * dis@beta_ext * k_ext) / (k_core + k_ext)
  beta_ceiling <- (beta_core_contrib + kappa_ceiling * dis@beta_ext * k_ext) / (k_core + k_ext)
  R0_floor     <- beta_floor   * (k_core + k_ext) / dis@gamma
  R0_ceiling   <- beta_ceiling * (k_core + k_ext) / dis@gamma
  
  # ── Layout ──────────────────────────────────────────────────────
  par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 2), oma = c(0, 0, 3, 0))
  
  # ── Panel 1: β_eff(T) ───────────────────────────────────────────
  plot(T_seq, beta_eff,
       type = "l", lwd = 2.5, col = "tomato",
       ylim = c(min(beta_eff)*0.9, max(beta_eff) * 1.1),
       xlab = "Daily mean temperature (°C)",
       ylab = expression(beta[eff] ~ "(effective transmission rate)"),
       main = expression(bold("Effective") ~ beta ~ bold("vs. Temperature")),
       las  = 1)
  
  if (use_temp) {
    abline(h = beta_floor,   lty = 3, col = "steelblue", lwd = 1.5)
    abline(h = beta_ceiling, lty = 3, col = "firebrick",  lwd = 1.5)
    text(T_range[2], beta_floor   - diff(range(beta_eff)) * 0.03,
         sprintf("floor = %.4f", beta_floor),     pos = 2, col = "steelblue", cex = 0.8)
    text(T_range[2], beta_ceiling + diff(range(beta_eff)) * 0.03,
         sprintf("ceiling = %.4f", beta_ceiling), pos = 2, col = "firebrick",  cex = 0.8)
    
    if (show_Tmid) {
      abline(v = dis@T_mid, lty = 2, col = "grey40", lwd = 1.2)
      text(dis@T_mid, max(beta_eff) * 0.05,
           sprintf("T_mid = %.1f°C", dis@T_mid),
           pos = 4, col = "grey30", cex = 0.8)
    }
  } else {
    abline(h = beta_floor, lty = 3, col = "steelblue", lwd = 1.5)
    text(T_range[2], beta_floor * 1.05,
         sprintf("constant = %.4f", beta_floor),
         pos = 2, col = "steelblue", cex = 0.8)
  }
  
  grid(col = "grey90", lty = 1)
  
  # ── Panel 2: R₀(T) ──────────────────────────────────────────────
  R0_threshold <- 1
  
  plot(T_seq, R0,
       type = "n",
       ylim = c(min(R0)*0.9, max(R0) * 1.1),
       xlab = "Daily mean temperature (°C)",
       ylab = expression(R[0] ~ "(basic reproduction number)"),
       main = expression(bold(R[0]) ~ bold("vs. Temperature")),
       las  = 1)
  
  polygon(c(T_seq, rev(T_seq)),
          c(pmax(R0, R0_threshold), rev(rep(R0_threshold, n_points))),
          col = adjustcolor("firebrick", 0.15), border = NA)
  polygon(c(T_seq, rev(T_seq)),
          c(pmin(R0, R0_threshold), rev(rep(R0_threshold, n_points))),
          col = adjustcolor("steelblue", 0.15), border = NA)
  
  lines(T_seq, R0, lwd = 2.5, col = "darkorchid")
  
  abline(h = R0_threshold, lty = 2, col = "black", lwd = 1.8)
  text(T_range[1], 1.06,
       expression(bold("Epidemic threshold   ") ~ R[0] == 1),
       pos = 4, col = "black", cex = 0.8)
  
  if(use_temp){
    abline(h = R0_floor,   lty = 3, col = "steelblue", lwd = 1.5)
    abline(h = R0_ceiling, lty = 3, col = "firebrick",  lwd = 1.5)
    text(T_range[2], R0_floor - diff(range(R0) + 1e-9) * 0.03,
         sprintf("R0 floor = %.2f",   R0_floor),   pos = 2, col = "steelblue", cex = 0.8)
    text(T_range[2], R0_ceiling + diff(range(R0) + 1e-9) * 0.03,
         sprintf("R0 ceiling = %.2f", R0_ceiling), pos = 2, col = "firebrick",  cex = 0.8)
  } else{
    abline(h = R0_floor,   lty = 3, col = "steelblue", lwd = 1.5)
    text(T_range[2], (R0_floor   + diff(range(R0) + 1e-9)) * 1.05,
         sprintf("R0 = %.2f",   R0_floor),   pos = 2, col = "steelblue", cex = 0.8)
  }
  
  if (use_temp && show_Tmid) {
    abline(v = dis@T_mid, lty = 2, col = "grey40", lwd = 1.2)
    text(dis@T_mid, max(R0) * 0.05,
         sprintf("T_mid = %.1f°C", dis@T_mid),
         pos = 4, col = "grey30", cex = 0.8)
  }
  
  grid(col = "grey90", lty = 1)
  
  # ── Shared title ─────────────────────────────────────────────────
  mtext(sprintf("Transmission profile — %s  (k_core = %d, k_ext = %d)",
                dis@name, k_core, k_ext),
        outer = TRUE, cex = 1.1, font = 2, line = 1)
  
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  invisible(list(T = T_seq, kappa = kappa, beta_eff = beta_eff, R0 = R0))
}

#' plot_temp_scaling_curve
#'
#' Visualise the sigmoid thermal response κ(T) for a disease_spatial object.
#' Marks the floor (y_min), ceiling (y_min + delta_y), and inflection (T_mid).
#'
#' @param dis     A disease_spatial object
#' @param T_range Numeric vector of length 2: temperature axis limits (°C)
plot_temp_scaling_curve <- function(dis, T_range = c(5, 30)) {
  T_seq <- seq(T_range[1], T_range[2], length.out = 300)
  kappa <- vapply(T_seq, \(T) temp_scaling(T, dis), numeric(1))
  
  plot(T_seq, kappa, type = "l", lwd = 2, col = "tomato",
       ylim = c(0, 1),
       main = sprintf("Sigmoid thermal response: %s", dis@name),
       xlab = "Temperature (°C)",
       ylab = expression(kappa ~ "(transmission scaling)"))
  
  abline(v = dis@T_mid, lty = 2, col = "grey40")
  text(dis@T_mid, 0.04, labels = sprintf("T_mid = %.1f°C", dis@T_mid),
       pos = 4, col = "grey30", cex = 0.85)
  
  abline(h = dis@y_min,               lty = 3, col = "steelblue", lwd = 1.5)
  abline(h = dis@y_min + dis@delta_y, lty = 3, col = "tomato",    lwd = 1.5)
  
  text(T_range[2], dis@y_min + 0.02,
       labels = sprintf("floor = %.2f", dis@y_min),
       pos = 2, col = "steelblue", cex = 0.8)
  text(T_range[2], dis@y_min + dis@delta_y - 0.03,
       labels = sprintf("ceiling = %.2f", dis@y_min + dis@delta_y),
       pos = 2, col = "tomato", cex = 0.8)
  
  grid()
}

#' plot_metapop_curves
#'
#' One SEIR panel per patch plus an optional aggregate panel.
#'
#' @param mpop      A metapopulation S7 object (after simulation)
#' @param aggregate Logical: add an extra panel summing across all patches
plot_metapop_curves <- function(mpop, aggregate = TRUE) {
  n_patches <- length(mpop@patches)
  n_panels  <- n_patches + as.integer(aggregate)
  n_cols    <- min(3L, n_panels)
  n_rows    <- ceiling(n_panels / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(3, 4, 3, 1))
  cols <- c(S = "steelblue", E = "orange", I = "red", R = "darkgreen")
  
  for (p in mpop@patches) {
    h <- p@history
    if (nrow(h) == 0) next
    matplot(h$day, h[, c("S", "E", "I", "R")],
            type = "l", lty = 1, lwd = 2, col = cols,
            main = p@name, xlab = "Day", ylab = "Agents",
            ylim = c(0, p@pop_size))
    legend("topright", legend = names(cols), col = cols,
           lty = 1, lwd = 2, bty = "n", cex = 0.7)
  }
  
  if (aggregate && nrow(mpop@history) > 0) {
    h       <- mpop@history
    N_total <- sum(vapply(mpop@patches, \(p) p@pop_size, integer(1)))
    matplot(h$day, h[, c("S", "E", "I", "R")],
            type = "l", lty = 1, lwd = 2, col = cols,
            main = "Aggregate (all patches)", xlab = "Day", ylab = "Agents",
            ylim = c(0, N_total))
    legend("topright", legend = names(cols), col = cols,
           lty = 1, lwd = 2, bty = "n", cex = 0.7)
  }
  par(mfrow = c(1, 1))
}


#' plot_prevalence_heatmap
#'
#' Patches × days intensity map of daily prevalence (%).
#'
#' @param mpop A metapopulation S7 object (after simulation)
plot_prevalence_heatmap <- function(mpop) {
  n_patches <- length(mpop@patches)
  n_days    <- nrow(mpop@patches[[1]]@history)
  if (n_days == 0) { message("No history to plot."); return(invisible(NULL)) }
  
  prev_mat    <- matrix(NA_real_, nrow = n_patches, ncol = n_days)
  patch_names <- character(n_patches)
  
  for (p in seq_len(n_patches)) {
    h              <- mpop@patches[[p]]@history
    patch_names[p] <- mpop@patches[[p]]@name
    prev_mat[p, ]  <- (h$I / mpop@patches[[p]]@pop_size) * 100
  }
  
  image(t(prev_mat),
        x    = seq_len(n_days),
        y    = seq_len(n_patches),
        col  = hcl.colors(64, "YlOrRd", rev = FALSE),
        xlab = "Day", ylab = "",
        main = "Prevalence (%) by patch and day",
        yaxt = "n")
  axis(2, at = seq_len(n_patches), labels = patch_names, las = 2, cex.axis = 0.8)
}