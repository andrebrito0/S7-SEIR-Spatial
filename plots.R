#----------------------------------------------------------------#
#  plots.R — Data Visualisation
#
#  Produces all figures for the presentation:
#    § 1  Seasonal dynamics plot      → figures/seasonal_plot.pdf
#    § 2  Spatial data preparation    (shared objects for §§ 3–5)
#    § 3  Temperature data            (shared objects for § 5)
#    § 4  Map preparation             (shared objects for § 5)
#    § 5  OD flow animation           → figures/od_flows_temporal.gif
#    § 6  Temperature coupling anim.  → figures/temp_coupling.gif
#    § 7  GIF frame extraction        → frames/od/, frames/temp/
#
#  Usage:
#    Run run_scenarios.R first to generate results/processed_inputs.rds.
#    Then: source("plots.R")   # from repo root
#    — or —  Rscript plots.R
#
#  Note: Sections 5–7 (animations) are slow (~5–10 min).
#        Set RENDER_ANIMATIONS <- FALSE to skip them.
#----------------------------------------------------------------#

RENDER_ANIMATIONS <- TRUE   # Set FALSE to skip GIF rendering (§§ 5–7)

# ── 0. Config & dependencies ──────────────────────────────────────
source("config.R")

library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(tidyr)
library(reshape2)
library(readxl)
library(sf)

if (RENDER_ANIMATIONS) {
  library(gganimate)
  library(transformr)
  library(gifski)
  library(magick)
}

# Create output directories if they do not exist
dir.create(DIR_FIGURES, showWarnings = FALSE)
dir.create(DIR_RESULTS, showWarnings = FALSE)

#================================================================#
# 1. SEASONAL DYNAMICS PLOT
#================================================================#

# --- 1.1 Load & validate data ---
check_data_file(PATH_HEALTHCARE)
check_data_file(PATH_TEMP_PER_ARS)

hosp_gripe   <- readRDS(PATH_HEALTHCARE)
temp_per_ARS <- readRDS(PATH_TEMP_PER_ARS) %>% rename(ars = ARS)

# --- 1.2 Build analysis data frame ---
df_hosp <- left_join(hosp_gripe, temp_per_ARS, by = c("date", "ars")) %>%
  filter(
    date < DATE_BREAK,
    type %in% c("n_consultas_gripe", "n_episdios_urgencia_infecao"),
    ars  == "ARS Lisboa e Vale do Tejo"
  ) %>%
  arrange(date) %>%
  group_by(type) %>%
  mutate(
    ma_value = rollapply(value,  7, mean, align = "right", fill = NA),
    temp_ma  = rollapply(t_mean, 7, mean, align = "right", fill = NA),
    season   = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                          ~ "Autumn"
    )
  ) %>%
  ungroup()

# --- 1.3 Season background rectangles ---
season_rects <- hosp_gripe %>%
  filter(date < DATE_BREAK) %>%
  mutate(
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                          ~ "Autumn"
    ),
    season_year = if_else(month(date) == 12, year(date) + 1, year(date))
  ) %>%
  distinct(date, season, season_year) %>%
  group_by(season, season_year) %>%
  summarise(xmin = min(date), xmax = max(date), .groups = "drop")

# --- 1.4 Dual-axis scaling parameters ---
y_mean    <- mean(df_hosp$ma_value, na.rm = TRUE)
y_sd      <- sd(df_hosp$ma_value,   na.rm = TRUE)
t_ma_mean <- mean(df_hosp$temp_ma,  na.rm = TRUE)
t_sd      <- sd(df_hosp$temp_ma,    na.rm = TRUE)

# --- 1.5 Plot ---
p_seasonal <- df_hosp %>%
  filter(type == "n_episdios_urgencia_infecao") %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data = season_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = season),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  geom_line(aes(y = ma_value), color = "black", linewidth = 0.8) +
  geom_line(
    aes(y = (temp_ma - t_ma_mean) / t_sd * y_sd + y_mean),
    color = "#C96A6A", linewidth = 0.4
  ) +
  facet_grid(
    type ~ .,
    scales   = "free_y",
    labeller = as_labeller(c(
      n_episdios_urgencia_infecao = "Respiratory Emergency Visits"
    ))
  ) +
  scale_y_continuous(
    name     = "Daily count",
    sec.axis = sec_axis(
      trans = ~ (. - y_mean) / y_sd * t_sd + t_ma_mean,
      name  = "Mean temperature"
    )
  ) +
  labs(fill = "Season", x = "Date") +
  theme_bw(base_size = 16) +
  theme(
    legend.position    = "bottom",
    axis.title.y.right = element_text(color = "#C96A6A"),
    axis.text.y.right  = element_text(color = "#C96A6A"),
    panel.background   = element_rect(fill = "#FDFEFE", color = NA),
    plot.background    = element_rect(fill = "#FDFEFE", color = NA),
    panel.grid.major   = element_line(color = "grey80"),
    panel.grid.minor   = element_blank(),
    legend.background  = element_rect(fill = "#FDFEFE", color = NA),
    legend.key         = element_rect(fill = "#FDFEFE", color = NA)
  )

# --- 1.6 Export ---
ggsave(
  filename = file.path(DIR_FIGURES, "seasonal_plot.pdf"),
  plot     = p_seasonal,
  width    = 16, height = 4, units = "in",
  device   = cairo_pdf
)
cat("Saved: seasonal_plot.pdf\n")

#================================================================#
# 2. MOBILITY & SPATIAL DATA
#    Load processed_inputs.rds saved by run_scenarios.R.
#    If not found, recompute from raw data to keep plots.R
#    self-contained when run_scenarios.R has not been run yet.
#================================================================#

processed_inputs_path <- file.path(DIR_RESULTS, "processed_inputs.rds")

if (file.exists(processed_inputs_path)) {
  # Fast path: reuse commuting matrix already built by run_scenarios.R
  cat("Loading pre-processed inputs from run_scenarios.R...\n")
  inputs           <- readRDS(processed_inputs_path)
  commuting_matrix <- inputs$commuting_matrix
  
} else {
  # Fallback: recompute from raw OD data
  cat("processed_inputs.rds not found — computing commuting matrix from raw data.\n")
  cat("(Run run_scenarios.R first to avoid this step.)\n")
  
  check_data_file(PATH_OD_COMPLETE)
  check_data_file(PATH_POP_INFO)
  
  od_complete <- readRDS(PATH_OD_COMPLETE)
  node_order  <- sort(unique(od_complete$Origin))
  
  mob_matrix <- od_complete %>%
    dcast(Origin ~ Destination, value.var = "count", fill = 0) %>%
    arrange(factor(Origin, levels = node_order)) %>%
    select(all_of(node_order)) %>%
    as.matrix()
  
  commuting_matrix <- mob_matrix / rowSums(mob_matrix)
  commuting_matrix[is.nan(commuting_matrix)] <- 0
}

# Population data frame (needed for centroid sizing in §5)
check_data_file(PATH_POP_INFO)
pop_info   <- readRDS(PATH_POP_INFO)
node_order <- sort(unique(colnames(commuting_matrix)))

pop_vector <- pop_info %>%
  arrange(factor(Origin, levels = node_order)) %>%
  pull(total_pop)
names(pop_vector) <- node_order

pop_df <- data.frame(
  municipio  = names(pop_vector),
  population = as.numeric(pop_vector)
)

#================================================================#
# 3. TEMPERATURE DATA (for animation in §6)
#================================================================#

check_data_file(PATH_TEMP_ARS)
temp_ars <- readRDS(PATH_TEMP_ARS)

temp_wide <- temp_ars %>%
  select(date, municipio, t_mean) %>%
  filter(date >= PLOT_START & date <= PLOT_END) %>%
  pivot_wider(names_from = municipio, values_from = t_mean) %>%
  arrange(date) %>%
  select(date, all_of(node_order))

temp_matrix_plots <- as.matrix(na.approx(select(temp_wide, -date), rule = 2))

#================================================================#
# 4. SPATIAL DATA — MAP PREPARATION
#================================================================#

check_data_file(PATH_CAOP)
check_data_file(PATH_ARS_EXCEL)

tabela_ars            <- read_excel(PATH_ARS_EXCEL) %>% rename(municipio = Município)
limites_freguesias_sf <- st_read(PATH_CAOP, quiet = TRUE)

limites_municipios_sf <- limites_freguesias_sf %>%
  mutate(caop_mun = substr(dtmnfr, 1, 4)) %>%
  group_by(caop_mun) %>%
  summarise(
    municipio = first(municipio),
    distrito  = first(distrito_ilha),
    .groups   = "drop"
  )

# Name mismatch check — print any ARS entries not matched in CAOP
unmatched <- tabela_ars$municipio[
  !(tabela_ars$municipio %in% limites_municipios_sf$municipio)
]
if (length(unmatched) > 0) {
  warning("ARS entries not found in CAOP limits (check name spelling):\n  ",
          paste(unmatched, collapse = "\n  "))
}

ars_sf_final <- left_join(limites_municipios_sf, tabela_ars, by = "municipio") %>%
  filter(!is.na(ARS))

lvt_map      <- ars_sf_final %>% filter(ARS == "ARS Lisboa e Vale do Tejo")
lvt_map_proj <- st_transform(lvt_map, 3763)   # ETRS89 / Portugal TM06

#================================================================#
# 5. OD FLOW ANIMATION
#================================================================#

if (RENDER_ANIMATIONS) {
  
  # --- 5.1 Centroids with population ---
  centroids_pop <- lvt_map_proj %>%
    st_centroid() %>%
    left_join(pop_df, by = "municipio") %>%
    mutate(pop_scaled = sqrt(population / max(population, na.rm = TRUE)))
  
  centroids_coords <- centroids_pop %>%
    mutate(
      lon = st_coordinates(geom)[, 1],
      lat = st_coordinates(geom)[, 2]
    ) %>%
    st_drop_geometry() %>%
    select(municipio, lon, lat)
  
  # --- 5.2 Build flow edge table ---
  # Busiest routes depart earliest; lightest routes trickle in later
  check_data_file(PATH_OD_COMPLETE)
  od_complete <- readRDS(PATH_OD_COMPLETE)
  
  flow_df <- od_complete %>%
    filter(Origin != Destination) %>%
    left_join(centroids_coords, by = c("Origin"      = "municipio")) %>%
    rename(lon_o = lon, lat_o = lat) %>%
    left_join(centroids_coords, by = c("Destination" = "municipio")) %>%
    rename(lon_d = lon, lat_d = lat) %>%
    filter(!is.na(lon_o), !is.na(lon_d)) %>%
    mutate(
      flow_weight = count / max(count),
      flow_rank   = rank(-count, ties.method = "random") / n(),
      depart_time = flow_rank
    )
  
  # --- 5.3 Build pulsing frames ---
  n_frames_od <- 120        # total animation frames
  n_pulse     <- 12         # frames per pulse cycle
  fade_start  <- 0.75       # fraction of day at which flows begin fading
  
  flow_frames <- lapply(seq_len(n_frames_od), function(f) {
    t         <- f / n_frames_od
    pulse_val <- 0.5 + 0.5 * sin(2 * pi * f / n_pulse)
    decay     <- if (t > fade_start) {
      1 - ((t - fade_start) / (1 - fade_start))^2
    } else 1
    
    flow_df %>%
      filter(depart_time <= t) %>%
      mutate(
        anim_frame  = f,
        time_active = t - depart_time,
        maturity    = pmin(time_active / 0.1, 1),
        envelope    = maturity * decay,
        arrow_width = flow_weight * (0.6 + 0.4 * pulse_val) * envelope,
        arrow_alpha = (0.2 + 0.6 * flow_weight) * (0.7 + 0.3 * pulse_val) * envelope
      )
  }) %>% bind_rows()
  
  # --- 5.4 Clock labels (06:00 → 22:00) ---
  time_labels <- data.frame(
    anim_frame = seq_len(n_frames_od),
    clock      = format(
      as.POSIXct("1970-01-01 06:00:00") +
        (seq_len(n_frames_od) / n_frames_od) * 16 * 3600,
      "%H:%M"
    )
  )
  flow_frames <- flow_frames %>% left_join(time_labels, by = "anim_frame")
  
  # --- 5.5 Base map ---
  base_map <- ggplot() +
    geom_sf(
      data      = lvt_map_proj,
      aes(fill  = as.numeric(factor(municipio))),
      colour    = "white",
      linewidth = 0.4
    ) +
    scale_fill_gradientn(
      colours = c("#b8d4e8", "#90c4d4", "#96dbc0", "#a8b8e0", "#c8b4d8"),
      guide   = "none"
    ) +
    theme_void(base_size = 13) +
    theme(
      legend.position   = "none",
      plot.title        = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle     = element_text(size = 11, hjust = 0.5, colour = "grey40"),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA)
    )
  
  # --- 5.6 Build animation ---
  flow_anim <- base_map +
    geom_curve(
      data      = flow_frames,
      aes(x = lon_o, y = lat_o, xend = lon_d, yend = lat_d,
          linewidth = flow_weight * 0.3,
          alpha     = envelope * 0.12,
          group     = interaction(Origin, Destination, anim_frame)),
      colour    = "#1a6faf",
      curvature = 0.25,
      arrow     = arrow(length = unit(0.010, "npc"), type = "closed")
    ) +
    geom_curve(
      data      = flow_frames,
      aes(x = lon_o, y = lat_o, xend = lon_d, yend = lat_d,
          linewidth = arrow_width,
          alpha     = arrow_alpha,
          group     = interaction(Origin, Destination, anim_frame)),
      colour    = "#1a6faf",
      curvature = 0.25,
      arrow     = arrow(length = unit(0.010, "npc"), type = "closed")
    ) +
    scale_linewidth_continuous(range = c(0.1, 3.0), guide = "none") +
    scale_alpha_continuous(range  = c(0.0, 0.95),  guide = "none") +
    transition_manual(anim_frame)
  
  # --- 5.7 Render ---
  animate(
    flow_anim,
    nframes  = n_frames_od,
    fps      = 12,
    width    = 900,
    height   = 900,
    bg       = "transparent",
    renderer = gifski_renderer(file.path(DIR_FIGURES, "od_flows_temporal.gif"))
  )
  cat("Saved: od_flows_temporal.gif\n")
  
} # end if (RENDER_ANIMATIONS)

#================================================================#
# 6. TEMPERATURE COUPLING ANIMATION
#================================================================#

if (RENDER_ANIMATIONS) {
  
  # Parameters come from config.R — no local re-definition
  # n_frames_kappa: frames built; rendered frames must match
  n_frames_kappa <- 300           # both built AND rendered at this value
  T_full         <- seq(-5, 30, length.out = 300)
  T_reveals      <- seq(-5, 30, length.out = n_frames_kappa)
  
  df_kappa <- do.call(rbind, lapply(seq_along(T_reveals), function(f) {
    T_sub <- T_full[T_full <= T_reveals[f]]
    data.frame(
      T     = T_sub,
      kappa = Y_MIN + DELTA_Y / (1 + exp(K_TEMP * (T_sub - T_MID))),
      frame = f
    )
  }))
  
  df_dot <- do.call(rbind, lapply(seq_along(T_reveals), function(f) {
    T_tip <- T_reveals[f]
    data.frame(
      T     = T_tip,
      kappa = Y_MIN + DELTA_Y / (1 + exp(K_TEMP * (T_tip - T_MID))),
      frame = f
    )
  }))
  
  p_kappa <- ggplot(df_kappa, aes(x = T, y = kappa)) +
    geom_area(fill = "steelblue", alpha = 0.15) +
    geom_line(colour = "steelblue", linewidth = 1.2) +
    geom_point(data = df_dot, colour = "steelblue", size = 3) +
    geom_hline(yintercept = Y_MIN,           linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = Y_MIN + DELTA_Y, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = T_MID, linetype = "dashed", colour = "red", alpha = 0.7) +
    annotate("text", x = T_MID + 0.5, y = 0.05,
             label = sprintf("T[mid]==%g*'\u00b0C'", T_MID), parse = TRUE,
             colour = "red", hjust = 0, size = 3.5) +
    annotate("text", x = -4.5, y = Y_MIN + 0.03,
             label = sprintf("y[min]==%.2f", Y_MIN), parse = TRUE,
             colour = "grey40", hjust = 0, size = 3) +
    annotate("text", x = -4.5, y = Y_MIN + DELTA_Y + 0.03,
             label = sprintf("y[min]+Delta*y==%.2f", Y_MIN + DELTA_Y), parse = TRUE,
             colour = "grey40", hjust = 0, size = 3) +
    scale_x_continuous(limits = c(-5, 30), expand = 0) +
    scale_y_continuous(limits = c(0, 1.1),  expand = 0) +
    labs(
      x     = "Daily Temperature T (\u00b0C)",
      y     = "\u03ba(T)",
      title = "Temperature Coupling Function \u03ba(T)"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA)
    ) +
    transition_manual(frame)
  
  # nframes matches n_frames_kappa — no wasted computation
  animate(
    p_kappa,
    nframes  = n_frames_kappa,
    fps      = 30,
    width    = 800,
    height   = 400,
    bg       = "transparent",
    renderer = gifski_renderer(file.path(DIR_FIGURES, "temp_coupling.gif"))
  )
  cat("Saved: temp_coupling.gif\n")
  
} # end if (RENDER_ANIMATIONS)

#================================================================#
# 7. GIF FRAME EXTRACTION (for LaTeX animategraphics)
#================================================================#

if (RENDER_ANIMATIONS) {
  
  od_gif_path   <- file.path(DIR_FIGURES, "od_flows_temporal.gif")
  temp_gif_path <- file.path(DIR_FIGURES, "temp_coupling.gif")
  
  # OD flows frames
  od_frames_dir <- file.path(DIR_FRAMES, "od")
  dir.create(od_frames_dir, recursive = TRUE, showWarnings = FALSE)
  od_gif <- image_read(od_gif_path)
  for (i in seq_along(od_gif)) {
    image_write(od_gif[i],
                path   = file.path(od_frames_dir, sprintf("frame-%d.png", i - 1)),
                format = "png")
  }
  cat(sprintf("Extracted %d OD flow frames to %s\n", length(od_gif), od_frames_dir))
  
  # Temperature coupling frames
  temp_frames_dir <- file.path(DIR_FRAMES, "temp")
  dir.create(temp_frames_dir, recursive = TRUE, showWarnings = FALSE)
  temp_gif <- image_read(temp_gif_path)
  for (i in seq_along(temp_gif)) {
    image_write(temp_gif[i],
                path   = file.path(temp_frames_dir, sprintf("frame-%d.png", i - 1)),
                format = "png")
  }
  cat(sprintf("Extracted %d temp coupling frames to %s\n",
              length(temp_gif), temp_frames_dir))
  
} # end if (RENDER_ANIMATIONS)

cat("\nplots.R complete.\n")

