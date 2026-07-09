# ============================================================
# SEIRD model for the spread of HPAI virus in NFLD seabirds in 2022
# ============================================================

# ============================================================
# Libraries & global settings
# ============================================================

library(macpan2)
library(tidyverse)
library(stringr)
library(splines)
library(tidyverse)
library(ggtext)

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_figures")
options(scipen = 999)

# ============================================================
# Load data
# ============================================================

dat3 <- readRDS("dat.rds") %>% ungroup()

# ============================================================
# Parameters
# ============================================================

alpha <- 1/7
gamma <- 1/11

patches <- paste0("Patch", 1:7)

muI_species <- c(
  "Northern Gannet" = 0.10,
  "Common Murre"    = 0.08
)

# ============================================================
# Carrying capacities 
# ============================================================

K_list <- list(
  "Northern Gannet" = c(
    Patch1 = 50000, Patch2 = 5000, Patch3 = 3488,
    Patch4 = 21928, Patch5 = 1000, Patch6 = 5000, Patch7 = 1000
  ),
  "Common Murre" = c(
    Patch1 = 20000, Patch2 = 750000, Patch3 = 5000,
    Patch4 = 944518, Patch5 = 5000, Patch6 = 10000, Patch7 = 15000
  )
)

flows <- list(
  mp_per_capita_flow("S", "E", "beta * I / N", "exposure"),
  mp_per_capita_flow("E", "I", "alpha", "infection"),
  mp_per_capita_flow("I", "R", "gamma", "recovery"),
  mp_per_capita_flow("I", "D", "muI", "death")
)

# ============================================================
# Model function
# ============================================================

run_model <- function(patch, species) {
  
  obs_all_fit <- dat3 %>%
    transmute(
      patch_name,
      CommonName = str_replace(bird, "_", " "),
      Date = as.Date(date),
      CumObserved = cumtotal_fill
    ) %>%
    filter(patch_name == patch, CommonName == species) %>%
    arrange(Date) %>%
    mutate(time = row_number())
  
  if (nrow(obs_all_fit) == 0) return(NULL)
  
  T <- max(obs_all_fit$time)
  D0 <- first(obs_all_fit$CumObserved)
  
  obsdat <- obs_all_fit %>%
    transmute(
      time,
      matrix = "D",
      value = CumObserved
    )
  
  muI <- muI_species[[species]]
  K <- K_list[[species]][[patch]]
  
  # ============================================================
  # Model spec
  # ============================================================
  
  spec <- mp_tmb_model_spec(
    
    before = list(
      D ~ D0,
      S ~ N - E - I - R - D
    ),
    
    during = flows,
    
    default = list(
      beta = 0.2,
      alpha = alpha,
      gamma = gamma,
      muI = muI,
      N = K,
      E = 20,
      I = 5,
      R = 0,
      D0 = D0
    )
  )
  
  # ============================================================
  # Time-varying beta
  # ============================================================
  
  beta_basis <- ns(1:T, df = 6)
  
  spec <- mp_tmb_insert_glm_timevar(
    spec,
    "beta",
    beta_basis,
    rep(0, ncol(beta_basis)),
    link_function = mp_log
  )
  
  # ============================================================
  # Calibration
  # ============================================================
  
  cal <- mp_tmb_calibrator(
    
    spec = spec |> mp_rk4(),
    
    data = obsdat,
    
    time = mp_sim_bounds(1, T),
    
    traj = list(
      D = mp_pois()
    ),
    
    default = list(
      muI = muI,
      N = K
    ),
    
    par = list(
      time_var_beta = mp_norm(0, 0.5),
      I = mp_norm(5, 3),
      E = mp_norm(20, 10)
    ),
    
    outputs = c("beta", "D")
  )
  
  mp_optimize(cal)
  
  fit <- mp_trajectory_sd(cal, conf.int = TRUE)
  
  fit %>%
    mutate(
      patch = patch,
      species = species
    )
}

# ============================================================
# Run all patches
# ============================================================

all_fits <- list()

for (patch in patches) {
  for (species in names(muI_species)) {
    
    message("Fitting: ", patch, " - ", species)
    
    res <- try(run_model(patch, species), silent = TRUE)
    
    if (!inherits(res, "try-error") && !is.null(res)) {
      all_fits[[paste(patch, species, sep = "_")]] <- res
    }
  }
}

# ============================================================
# Combine output
# ============================================================

fitted_data <- bind_rows(all_fits)

saveRDS(fitted_data, "fitted_data.rds")

# ============================================================
# Plotting
# ============================================================

# ============================================================
# Load data
# ============================================================

dat3 <- readRDS("dat.rds")
fitted_data <- readRDS("fitted_data.rds")

# ============================================================
# Constants
# ============================================================

species_colors <- c(
  "Northern Gannet" = "red",
  "Common Murre"    = "blue"
)

custom_patch_names <- c(
  "Patch1" = "Cape St. Mary's",
  "Patch2" = "Witless Bay",
  "Patch3" = "Baccalieu Island",
  "Patch4" = "Funk Island",
  "Patch5" = "Hare Bay",
  "Patch6" = "West Coast",
  "Patch7" = "Lawn Islands"
)

month_breaks <- function(x) {
  seq(
    from = as.Date(format(min(x), "%Y-%m-01")),
    to   = as.Date(format(max(x), "%Y-%m-01")),
    by   = "1 month"
  )
}

beta_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 20),
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_blank()
  )

# ============================================================
# Observed data
# ============================================================

obs_all_fit <- dat3 %>%
  transmute(
    patch = patch_name,
    CommonName = str_replace(bird, "_", " "),
    Date = as.Date(date),
    CumObserved = cumtotal_fill
  ) %>%
  filter(CommonName %in% names(species_colors)) %>%
  arrange(patch, CommonName, Date) %>%
  group_by(patch, CommonName) %>%
  mutate(time = row_number()) %>%
  ungroup()

obs_plot <- obs_all_fit %>%
  rename(species = CommonName) %>%
  mutate(
    patch = factor(
      patch,
      levels = names(custom_patch_names),
      labels = custom_patch_names
    )
  )

# ============================================================
# Time lookup
# ============================================================

time_lookup <- obs_all_fit %>%
  transmute(
    patch,
    species = CommonName,
    time,
    Date
  )

time_lookup0 <- obs_all_fit %>%
  group_by(patch, CommonName) %>%
  slice(1) %>%
  transmute(
    patch,
    species = CommonName,
    time = 0,
    Date
  )

time_lookup <- bind_rows(time_lookup0, time_lookup)

# ============================================================
# Death model output
# ============================================================

fit_D <- fitted_data %>%
  filter(matrix == "D") %>%
  mutate(
    patch = as.character(patch),
    species = as.character(species),
    time = as.integer(time)
  )

D0_table <- obs_all_fit %>%
  group_by(patch, CommonName) %>%
  summarise(D0 = first(CumObserved), .groups = "drop") %>%
  transmute(
    patch,
    species = CommonName,
    D0
  )

death_plot_data <- fit_D %>%
  left_join(time_lookup, by = c("patch", "species", "time")) %>%
  left_join(D0_table, by = c("patch", "species")) %>%
  arrange(patch, species, time) %>%
  group_by(patch, species) %>%
  mutate(
    value = replace(value, 1, D0),
    conf.low = replace(conf.low, 1, D0),
    conf.high = replace(conf.high, 1, D0),
    value_pos = pmax(value, 0),
    ci_low = pmax(conf.low, 0),
    ci_high = pmax(conf.high, 0),
    patch = factor(
      patch,
      levels = names(custom_patch_names),
      labels = custom_patch_names
    )
  ) %>%
  ungroup()

# ============================================================
# Death plot: Fig 3.4
# ============================================================

p_deaths <- ggplot(death_plot_data) +
  
  geom_ribbon(
    aes(Date, ymin = ci_low + 1, ymax = ci_high + 1, fill = species),
    alpha = 0.25
  ) +
  
  geom_line(
    aes(Date, value_pos + 1, colour = species, group = interaction(species, patch)),
    linewidth = 1.3
  ) +
  
  geom_point(
    data = obs_plot,
    aes(Date, CumObserved + 1, colour = species),
    size = 3
  ) +
  
  facet_wrap(~patch, scales = "free", labeller = label_wrap_gen()) +
  
  scale_colour_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  
  scale_y_log10() +
  scale_x_date(breaks = month_breaks, date_labels = "%b") +
  
  labs(
    x = NULL,
    y = "Cumulative mortality",
    colour = NULL,
    fill = NULL
  ) +
  
  theme_bw(base_size = 20) +
  
  theme(
    legend.position = c(0.50, 0.25),
    legend.justification = c(0.5, 1),
    legend.text = element_text(size = 18),      
    legend.key.size = unit(1.2, "cm"),          
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold")
  )

ggsave(
  "Fig_3.4.png",
  p_deaths,
  width = 16,
  height = 10,
  dpi = 300
)

# ============================================================
# Beta model output
# ============================================================

beta_data <- fitted_data %>%
  filter(matrix == "beta") %>%
  mutate(
    patch = as.character(patch),
    species = as.character(species),
    time = as.integer(time)
  )

beta_plot_data <- beta_data %>%
  left_join(time_lookup, by = c("patch", "species", "time")) %>%
  mutate(
    value_pos = pmax(value, 0),
    ci_low    = pmax(conf.low, 0),
    ci_high   = pmax(conf.high, 0),
    patch = factor(patch, levels = names(custom_patch_names),
                   labels = custom_patch_names)
  )

# ============================================================
# Beta timeseries: Fig 3.5
# ============================================================

july_start <- min(
  beta_plot_data$Date[format(beta_plot_data$Date, "%m-%d") >= "07-01"],
  na.rm = TRUE
)

p_beta <-
  ggplot(beta_plot_data) +
  
  geom_ribbon(
    aes(
      x = Date,
      ymin = ci_low,
      ymax = ci_high,
      fill = patch
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_line(
    aes(
      x = Date,
      y = value_pos,
      colour = patch,
      group = interaction(species, patch)
    ),
    linewidth = 1
  ) +
  
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  
  coord_cartesian(
    xlim = c(july_start, max(beta_plot_data$Date, na.rm = TRUE))
  ) +
  
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  
  labs(y = "Transmission rate <i>β(t)</i>") +
  
  beta_theme

ggsave(
  "Fig_3.5.png",
  p_beta,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# Beta distribution by species: Fig 3.6A
# ============================================================

beta_dist <- beta_plot_data %>%
  select(species, patch, beta = value_pos)

p_beta_species <- ggplot(beta_dist,
                         aes(species, beta, fill = species)) +
  
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
  
  scale_fill_manual(values = species_colors) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    y = "Transmission rate <i>β(t)</i>",
    tag = "A"
  ) +
  
  beta_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.tag = element_text(face = "bold")
  )

ggsave("Fig_3.6A.png",
       p_beta_species, width = 7, height = 5, dpi = 300)

# ============================================================
# Beta distribution by patch: Fig 3.6B
# ============================================================

p_beta_patch <- ggplot(beta_dist,
                       aes(species, beta, fill = species)) +
  
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
  
  facet_wrap(~patch, ncol = 3) +
  
  scale_fill_manual(values = species_colors) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(y = "Transmission rate <i>β(t)</i>", tag = "B") +
  
  beta_theme +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = c(0.5, 0.15),
    plot.tag = element_text(face = "bold")
  )

ggsave("Fig_3.6B.png",
       p_beta_patch, width = 12, height = 8, dpi = 300)

# ============================================================
# Extract beta
# ============================================================

beta_data <- fitted_data %>%
  filter(matrix == "beta") %>%
  mutate(
    value = pmax(value, 0)
  )

# ============================================================
# Summary by patch x species
# ============================================================

beta_summary_patch <- beta_data %>%
  group_by(patch, species) %>%
  summarise(
    Mean_Beta   = mean(value, na.rm = TRUE),
    Median_Beta = median(value, na.rm = TRUE),
    SD_Beta     = sd(value, na.rm = TRUE),
    
    Q1_Beta     = quantile(value, 0.25, na.rm = TRUE),
    Q3_Beta     = quantile(value, 0.75, na.rm = TRUE),
    IQR_Beta    = IQR(value, na.rm = TRUE),
    
    Min_Beta    = min(value, na.rm = TRUE),
    Max_Beta    = max(value, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(patch, species)

print(beta_summary_patch)

# ============================================================
# Summary by species x patch
# ============================================================

beta_summary_species_patch <- beta_data %>%
  group_by(species, patch) %>%
  summarise(
    Mean_Beta   = mean(value, na.rm = TRUE),
    Median_Beta = median(value, na.rm = TRUE),
    SD_Beta     = sd(value, na.rm = TRUE),
    
    Q1_Beta     = quantile(value, 0.25, na.rm = TRUE),
    Q3_Beta     = quantile(value, 0.75, na.rm = TRUE),
    IQR_Beta    = IQR(value, na.rm = TRUE),
    
    Min_Beta    = min(value, na.rm = TRUE),
    Max_Beta    = max(value, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(species, patch)

print(beta_summary_species_patch)

# ============================================================
# Overall summary by species 
# ============================================================

beta_summary_overall <- beta_data %>%
  group_by(species) %>%
  summarise(
    Mean_Beta   = mean(value, na.rm = TRUE),
    Median_Beta = median(value, na.rm = TRUE),
    SD_Beta     = sd(value, na.rm = TRUE),
    
    Q1_Beta     = quantile(value, 0.25, na.rm = TRUE),
    Q3_Beta     = quantile(value, 0.75, na.rm = TRUE),
    IQR_Beta    = IQR(value, na.rm = TRUE),
    
    Min_Beta    = min(value, na.rm = TRUE),
    Max_Beta    = max(value, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(species)

print(beta_summary_overall)
