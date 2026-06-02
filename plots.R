#================================================================#
#  plots.R  —  Presentation Figures & Animations
#
#  Produces all figures used in the Beamer slides.
#  Run run_scenarios.R first so that results/processed_inputs.rds
#  exists; otherwise the commuting matrix is recomputed from raw data.
#
#  Usage
#  ─────────────────────────────────────────────────────────────
#  From an interactive session at the repository root:
#    source("plots.R")
#  From the command line:
#    Rscript plots.R
#
#  Toggle animations
#  ─────────────────────────────────────────────────────────────
#  Set RENDER_ANIMATIONS <- TRUE to produce GIF outputs (§§ 5–7).
#  This adds ~5–10 minutes of rendering time.
#
#  Outputs  →  figures/
#  ─────────────────────────────────────────────────────────────
#  seasonal_plot.pdf         — §1  Seasonal respiratory visits
#  temp_coupling.pdf         — §6  Sigmoid κ(T) curve
#  od_flows_temporal.gif     — §5  OD flow animation (optional)
#
#  Outputs  →  frames/
#  ─────────────────────────────────────────────────────────────
#  od/frame-N.png            — §7  Individual GIF frames for LaTeX
#                                   \animategraphics
#================================================================#

RENDER_ANIMATIONS <- FALSE   # Set TRUE to render GIF animations (§§ 5–7)

#----------------------------------------------------------------#
# § 0  DEPENDENCIES
#----------------------------------------------------------------#

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

dir.create(DIR_FIGURES, showWarnings = FALSE)
dir.create(DIR_RESULTS, showWarnings = FALSE)


#================================================================#
# § 1  SEASONAL DYNAMICS PLOT
#
#  7-day moving average of respiratory emergency visits in ARS LVT
#  with mean temperature overlaid on a secondary axis.  Season
#  backgrounds derived from month of year.
#  Output: figures/seasonal_plot.pdf
#================================================================#

# ── 1.1  Load data ───────────────────────────────────────────────
check_data_file(PATH_HEALTHCARE)
check_data_file(PATH_TEMP_PER_ARS)

hosp_gripe   <- readRDS(PATH_HEALTHCARE)
temp_per_ARS <- readRDS(PATH_TEMP_PER_ARS) %>% rename(ars = ARS)

# ── 1.2  Build analysis data frame ───────────────────────────────
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

# ── 1.3  Season background rectangles ────────────────────────────
season_rects <- hosp_gripe %>%
  filter(date < DATE_BREAK) %>%
  mutate(
    season      = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                          ~ "Autumn"
    ),
    season_year = if_else(month(date) == 12, year(date) + 1L, year(date))
  ) %>%
  distinct(date, season, season_year) %>%
  group_by(season, season_year) %>%
  summarise(xmin = min(date), xmax = max(date), .groups = "drop")

# ── 1.4  Dual-axis scaling parameters ────────────────────────────
# Standardise temperature to the same scale as the visit count so
# both series can share the same panel with a labelled secondary axis.
y_mean    <- mean(df_hosp$ma_value, na.rm = TRUE)
y_sd      <- sd(df_hosp$ma_value,   na.rm = TRUE)
t_ma_mean <- mean(df_hosp$temp_ma,  na.rm = TRUE)
t_sd      <- sd(df_hosp$temp_ma,    na.rm = TRUE)

# ── 1.5  Build & export plot ─────────────────────────────────────
p_seasonal <- df_hosp %>%
  filter(type == "n_episdios_urgencia_infecao") %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data        = season_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = season),
    inherit.aes = FALSE,
    alpha       = 0.08
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

ggsave(
  filename = file.path(DIR_FIGURES, "seasonal_plot.pdf"),
  plot     = p_seasonal,
  width    = 16, height = 4, units = "in",
  device   = cairo_pdf
)
cat("Saved: seasonal_plot.pdf\n")


#================================================================#
# § 2  MOBILITY & SPATIAL DATA
#
#  Load the commuting matrix from processed_inputs.rds if available
#  (written by run_scenarios.R); otherwise recompute from raw data.
#  Also loads population sizes used for centroid sizing in § 5.
#================================================================#

processed_inputs_path <- file.path(DIR_RESULTS, "processed_inputs.rds")

if (file.exists(processed_inputs_path)) {
  cat("Loading pre-processed inputs from run_scenarios.R...\n")
  inputs           <- readRDS(processed_inputs_path)
  commuting_matrix <- inputs$commuting_matrix
  
} else {
  # Fallback: recompute directly from raw OD data
  cat("processed_inputs.rds not found — recomputing from raw data.\n")
  cat("(Run run_scenarios.R first to avoid this step.)\n")
  
  check_data_file(PATH_OD_COMPLETE)
  
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

# node_order is the canonical column ordering of the commuting matrix
node_order <- colnames(commuting_matrix)

# Population data frame (used for centroid sizing in § 5)
check_data_file(PATH_POP_INFO)
pop_info <- readRDS(PATH_POP_INFO)

pop_vector <- pop_info %>%
  arrange(factor(Origin, levels = node_order)) %>%
  pull(total_pop)
names(pop_vector) <- node_order

pop_df <- data.frame(
  municipio  = names(pop_vector),
  population = as.numeric(pop_vector)
)


#================================================================#
# § 3  TEMPERATURE DATA
#
#  Wide temperature data frame for the plotting window; used in
#  § 5 (animation colouring) if extended to overlay temperature.
#================================================================#

check_data_file(PATH_TEMP_ARS)
temp_ars <- readRDS(PATH_TEMP_ARS)

temp_wide <- temp_ars %>%
  select(date, municipio, t_mean) %>%
  filter(date >= PLOT_START, date <= PLOT_END) %>%
  pivot_wider(names_from = municipio, values_from = t_mean) %>%
  arrange(date) %>%
  select(date, all_of(node_order))


#================================================================#
# § 4  SPATIAL DATA — MAP PREPARATION
#
#  Loads the pre-projected LVT municipality polygons.
#  The original CAOP → ARS join is commented out because the raw
#  CAOP file is too large to include in the repository.
#  Re-run the commented block once to regenerate lvt_map_proj.rds.
#================================================================#

# Re-generation block (run once if data/lvt_map_proj.rds is missing):
# ─────────────────────────────────────────────────────────────────
# check_data_file(PATH_CAOP)
# check_data_file(PATH_ARS_EXCEL)
#
# tabela_ars            <- read_excel(PATH_ARS_EXCEL) %>% rename(municipio = Município)
# limites_freguesias_sf <- st_read(PATH_CAOP, quiet = TRUE)
#
# limites_municipios_sf <- limites_freguesias_sf %>%
#   mutate(caop_mun = substr(dtmnfr, 1, 4)) %>%
#   group_by(caop_mun) %>%
#   summarise(
#     municipio = first(municipio),
#     distrito  = first(distrito_ilha),
#     .groups   = "drop"
#   )
#
# unmatched <- tabela_ars$municipio[
#   !(tabela_ars$municipio %in% limites_municipios_sf$municipio)
# ]
# if (length(unmatched) > 0)
#   warning("ARS entries not matched in CAOP: ", paste(unmatched, collapse = ", "))
#
# ars_sf_final <- left_join(limites_municipios_sf, tabela_ars, by = "municipio") %>%
#   filter(!is.na(ARS))
#
# lvt_map_proj <- ars_sf_final %>%
#   filter(ARS == "ARS Lisboa e Vale do Tejo") %>%
#   st_transform(3763)   # ETRS89 / Portugal TM06
#
# saveRDS(lvt_map_proj, "data/lvt_map_proj.rds")
# ─────────────────────────────────────────────────────────────────

lvt_map_proj <- readRDS("data/lvt_map_proj.rds")


#================================================================#
# § 5  OD FLOW ANIMATION
#
#  Animated daily commuting flows across LVT municipalities.
#  Flows are revealed progressively (busiest routes first) and
#  pulse with a sinusoidal envelope to suggest movement.
#  Output: figures/od_flows_temporal.gif
#================================================================#

if (RENDER_ANIMATIONS) {
  
  # ── 5.1  Centroids with population ───────────────────────────
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
  
  # ── 5.2  Flow edge table ──────────────────────────────────────
  # Rank flows by volume so that the busiest routes appear first
  # in the animation; flow_rank drives the depart_time offset.
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
  
  # ── 5.3  Build per-frame animation data ───────────────────────
  n_frames_od <- 120    # total frames
  n_pulse     <- 12     # frames per sinusoidal pulse cycle
  fade_start  <- 0.75   # fraction of the animation at which flows begin fading
  
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
  
  # ── 5.4  Clock labels (06:00 → 22:00) ────────────────────────
  time_labels <- data.frame(
    anim_frame = seq_len(n_frames_od),
    clock      = format(
      as.POSIXct("1970-01-01 06:00:00") +
        (seq_len(n_frames_od) / n_frames_od) * 16 * 3600,
      "%H:%M"
    )
  )
  flow_frames <- flow_frames %>% left_join(time_labels, by = "anim_frame")
  
  # ── 5.5  Base map ─────────────────────────────────────────────
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
      legend.position  = "none",
      plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
  
  # ── 5.6  Construct animation ──────────────────────────────────
  flow_anim <- base_map +
    # Background ghost layer — very faint permanent arc
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
    # Foreground layer — pulsing, weight-scaled arrows
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
  
  # ── 5.7  Render & save ────────────────────────────────────────
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
# § 6  TEMPERATURE COUPLING CURVE
#
#  Static plot of the sigmoid κ(T) function with annotated floor,
#  ceiling, and T_mid reference line.  Parameters read from config.R.
#  Output: figures/temp_coupling.pdf
#================================================================#

T_full <- seq(5, 30, length.out = 300)

df_kappa <- data.frame(
  t     = T_full,
  kappa = Y_MIN + DELTA_Y / (1 + exp(K_TEMP * (T_full - T_MID)))
)

p_kappa <- ggplot(df_kappa, aes(x = t, y = kappa)) +
  geom_area(fill = "steelblue", alpha = 0.15) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  # Reference lines for floor, ceiling, and T_mid
  geom_hline(yintercept = Y_MIN,           linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = Y_MIN + DELTA_Y, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = T_MID,           linetype = "dashed", colour = "red", alpha = 0.7) +
  # Annotations — use T_MID from config so label always matches the line
  annotate("text", x = T_MID + 0.5, y = 0.05,
           label  = sprintf("T[mid]==%g*'\u00b0C'", T_MID), parse = TRUE,
           colour = "red", hjust = 0, size = 6) +
  annotate("text", x = 6.5, y = Y_MIN + 0.03,
           label  = sprintf("y[min]==%.2f", Y_MIN), parse = TRUE,
           colour = "grey40", hjust = 0, size = 6) +
  annotate("text", x = 6.5, y = Y_MIN + DELTA_Y + 0.03,
           label  = sprintf("y[min]+Delta*y==%.2f", Y_MIN + DELTA_Y), parse = TRUE,
           colour = "grey40", hjust = 0, size = 6) +
  scale_x_continuous(limits = c(5, 30), expand = 0) +
  scale_y_continuous(limits = c(0, 1.1), expand = 0) +
  labs(
    x = "Daily Temperature T (\u00b0C)",
    y = "\u03ba(T)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave(
  filename = file.path(DIR_FIGURES, "temp_coupling.pdf"),
  plot     = p_kappa,
  width    = 6, height = 7, units = "in",
  device   = cairo_pdf
)
cat("Saved: temp_coupling.pdf\n")


#================================================================#
# § 7  GIF FRAME EXTRACTION FOR LATEX
#
#  Extracts individual PNG frames from the OD GIF so that the
#  LaTeX \animategraphics command can embed the animation in the
#  Beamer slides without requiring an external GIF viewer.
#  Output: frames/od/frame-N.png  (0-indexed, matching \animategraphics)
#================================================================#

if (RENDER_ANIMATIONS) {
  
  od_gif_path   <- file.path(DIR_FIGURES, "od_flows_temporal.gif")
  od_frames_dir <- file.path(DIR_FRAMES, "od")
  dir.create(od_frames_dir, recursive = TRUE, showWarnings = FALSE)
  
  od_gif <- image_read(od_gif_path)
  for (i in seq_along(od_gif)) {
    image_write(
      od_gif[i],
      path   = file.path(od_frames_dir, sprintf("frame-%d.png", i - 1)),
      format = "png"
    )
  }
  cat(sprintf("Extracted %d OD flow frames to %s\n",
              length(od_gif), od_frames_dir))
  
} # end if (RENDER_ANIMATIONS)

cat("\nplots.R complete.\n")
