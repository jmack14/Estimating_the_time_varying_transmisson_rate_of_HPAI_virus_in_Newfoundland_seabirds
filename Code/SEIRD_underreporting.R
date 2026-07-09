# ============================================================
# SEIRD model for the spread of HPAI virus in NFLD seabirds in 2022 with underreporting scenarios
# ============================================================

# ============================================================
# Libraries & global settings
# ============================================================

library(macpan2)
library(tidyverse)
library(stringr)
library(splines)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggtext)

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
# Uunderreporting scenarios
# ============================================================

reporting_scenarios <- c(
  low_reporting    = 0.2,
  medium_reporting = 0.5,
  high_reporting   = 0.8
)

# ============================================================
# Fitting loop
# ============================================================

all_fits <- list()

for (scenario in names(reporting_scenarios)) {
  
  rho <- reporting_scenarios[[scenario]]
  
  for (patch in patches) {
    for (species in names(muI_species)) {
      
      message("Fitting: ", scenario, " | ", patch, " - ", species)
      
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
      
      if (nrow(obs_all_fit) == 0) next
      
      T  <- max(obs_all_fit$time)
      D0 <- first(obs_all_fit$CumObserved)
      
      muI <- muI_species[[species]]
      K   <- K_list[[species]][[patch]]
      
      # ========================================================
      # Observation model
      # ========================================================
      
      obsdat <- obs_all_fit %>%
        transmute(
          time,
          matrix = "D_obs",
          value = CumObserved
        )
      
      reporting_map <- list(
        D_obs ~ rho * D
      )
      
      # ========================================================
      # Model spec
      # ========================================================
      
      spec <- mp_tmb_model_spec(
        
        before = list(
          D ~ D0,
          S ~ N - E - I - R - D
        ),
        
        during = c(
          flows,
          reporting_map
        ),
        
        default = list(
          beta  = 0.2,
          alpha = alpha,
          gamma = gamma,
          muI   = muI,
          rho   = rho,
          N     = K,
          E     = 20,
          I     = 5,
          R     = 0,
          D0    = D0
        )
      )
      
      # ========================================================
      # Time-varying beta
      # ========================================================
      
      beta_basis <- ns(1:T, df = 6)
      
      spec <- mp_tmb_insert_glm_timevar(
        spec,
        "beta",
        beta_basis,
        rep(0, ncol(beta_basis)),
        link_function = mp_log
      )
      
      # ========================================================
      # Calibration
      # ========================================================
      
      cal <- mp_tmb_calibrator(
        
        spec = spec |> mp_rk4(),
        
        data = obsdat,
        
        time = mp_sim_bounds(1, T),
        
        traj = list(
          D_obs = mp_pois()
        ),
        
        default = list(
          muI = muI,
          N   = K,
          rho = rho
        ),
        
        par = list(
          time_var_beta = mp_norm(0, 0.5),
          I = mp_norm(5, 3),
          E = mp_norm(20, 10)
        ),
        
        outputs = c("beta", "D", "D_obs")
      )
      
      mp_optimize(cal)
      
      fit <- try(
        mp_trajectory_sd(cal, conf.int = TRUE),
        silent = TRUE
      )
      
      if (inherits(fit, "try-error")) next
      
      all_fits[[paste(scenario, patch, species, sep = "_")]] <-
        fit %>%
        mutate(
          scenario = scenario,
          rho = rho,
          patch = patch,
          species = species
        )
    }
  }
}

# ============================================================
# Combine output
# ============================================================

fitted_data <- bind_rows(all_fits)

saveRDS(fitted_data, "fitted_data_underreporting.rds")

# ============================================================
# Load data
# ============================================================

dat3 <- readRDS("dat.rds")
fitted_data <- readRDS("fitted_data_underreporting.rds")

# ============================================================
# Plotting
# ============================================================

reporting_scenarios <- c(
  low_reporting = 0.2,
  medium_reporting = 0.5,
  high_reporting = 0.8
)

custom_patch_names <- c(
  Patch1 = "Cape St. Mary's",
  Patch2 = "Witless Bay",
  Patch3 = "Baccalieu Island",
  Patch4 = "Funk Island",
  Patch5 = "Hare Bay",
  Patch6 = "West Coast",
  Patch7 = "Lawn Islands"
)

# ============================================================
# Rebuild time lookup
# ============================================================

time_lookup <- dat3 %>%
  transmute(
    patch = patch_name,
    species = str_replace(bird, "_", " "),
    Date = as.Date(date)
  ) %>%
  group_by(patch, species) %>%
  mutate(time = row_number()) %>%
  ungroup() %>%
  distinct(patch, species, time, Date)

# ============================================================
# Extract beta 
# ============================================================

beta_data <- fitted_data %>%
  filter(matrix == "beta") %>%
  mutate(
    patch = factor(
      patch,
      levels = names(custom_patch_names),
      labels = custom_patch_names
    ),
    species = as.character(species),
    time = as.integer(time),
    beta = pmax(value, 0)
  )

beta_data <- beta_data %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("low_reporting", "medium_reporting", "high_reporting"),
      labels = c("Low", "Medium", "High")
    )
  )


# ============================================================
# Global scaling 
# ============================================================

beta_ymax <- quantile(beta_data$beta, 0.95, na.rm = TRUE)

# ============================================================
# Theme
# ============================================================

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
# Beta distribution by species: Fig 3.7A
# ============================================================

p_beta_species <- ggplot(beta_data, aes(scenario, beta, fill = scenario)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~species) +
  scale_fill_manual(
    values = c(
      "Low" = "#1b9e77",
      "Medium" = "#7570b3",
      "High" = "#d95f02"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(
    ylim = c(0, beta_ymax)
  ) +
  labs(
    tag = "A",
    x = NULL,
    y = "Transmission rate <i>β(t)</i>"
  ) +
  beta_theme +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "Fig_3.7A.png",
  p_beta_species,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# Beta distribution by patch: Fig 3.7B
# ============================================================

p_beta_patch <- ggplot(beta_data, aes(scenario, beta, fill = scenario)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~patch, ncol = 3, labeller = label_wrap_gen()) +
  scale_fill_manual(
    values = c(
      "Low" = "#1b9e77",
      "Medium" = "#7570b3",
      "High" = "#d95f02"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(
    ylim = c(0, beta_ymax)
  ) +
  labs(
    tag = "B",
    x = NULL,
    y = "Transmission rate <i>β(t)</i>"
  ) +
  beta_theme +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.5, 0.15),
    legend.justification = c(0.5, 0),
    legend.background = element_blank()
  )

ggsave(
  "Fig_3.7B.png",
  p_beta_patch,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# Beta trajectory data
# ============================================================

beta_plot_data <- beta_data %>%
  left_join(time_lookup, by = c("patch", "species", "time")) %>%
  mutate(
    value_pos = beta,
    ci_low = pmax(conf.low, 0),
    ci_high = pmax(conf.high, 0),
    patch = factor(patch, levels = names(custom_patch_names))
  )

# ============================================================
# Step structure
# ============================================================

beta_step <- beta_plot_data %>%
  arrange(scenario, species, patch, Date) %>%
  group_by(scenario, species, patch) %>%
  mutate(Date_end = lead(Date)) %>%
  mutate(Date_end = if_else(is.na(Date_end), max(Date) + 1, Date_end)) %>%
  ungroup()

# ============================================================
# Beta summaries 
# ============================================================

beta_summary_patch <- beta_data %>%
  group_by(scenario, patch) %>%
  summarise(
    Mean_Beta   = mean(beta, na.rm = TRUE),
    Median_Beta = median(beta, na.rm = TRUE),
    SD_Beta     = sd(beta, na.rm = TRUE),
    Q1_Beta     = quantile(beta, 0.25, na.rm = TRUE),
    Q3_Beta     = quantile(beta, 0.75, na.rm = TRUE),
    IQR_Beta    = IQR(beta, na.rm = TRUE),
    Min_Beta    = min(beta, na.rm = TRUE),
    Max_Beta    = max(beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, patch)

print(beta_summary_patch)

beta_summary_species <- beta_data %>%
  group_by(scenario, species) %>%
  summarise(
    Mean_Beta   = mean(beta, na.rm = TRUE),
    Median_Beta = median(beta, na.rm = TRUE),
    SD_Beta     = sd(beta, na.rm = TRUE),
    Q1_Beta     = quantile(beta, 0.25, na.rm = TRUE),
    Q3_Beta     = quantile(beta, 0.75, na.rm = TRUE),
    IQR_Beta    = IQR(beta, na.rm = TRUE),
    Min_Beta    = min(beta, na.rm = TRUE),
    Max_Beta    = max(beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, species)

print(beta_summary_species)
