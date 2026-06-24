# ============================================================
# SEIRD model for HPAI virus in NFLD seabirds with underreporting scenarios 
# ============================================================

# ============================================================
# Libraries and setup
# ============================================================

library(macpan2)
library(tidyverse)
library(patchwork)
library(zoo)
library(ggtext)
library(stringr)
library(splines)
library(tidyr)
library(ggplot2)
options(scipen = 999)

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_figures")
theme_set(theme_bw())

# ============================================================
# Parameters
# ============================================================

alpha <- 1/7
gamma <- 1/11

patches <- paste0("Patch", 1:7)
i_vec <- 1:7

species_list <- c("Northern Gannet", "Common Murre")

muI_species <- c(
  "Northern Gannet" = 0.10,
  "Common Murre"    = 0.08
)

species_colors <- c(
  "Northern Gannet" = "red",
  "Common Murre"    = "blue"
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
# Carrying capacities
# ============================================================

K_list <- list(
  "Northern Gannet" = c(Patch1=50000, Patch2=5000, Patch3=3488,
                        Patch4=21928, Patch5=1000, Patch6=5000, Patch7=1000),
  "Common Murre"    = c(Patch1=20000, Patch2=750000, Patch3=5000,
                        Patch4=944518, Patch5=5000, Patch6=10000, Patch7=15000)
)

# ============================================================
# Data prep
# ============================================================

obs_all <- read.csv(
  "NFLD_mortalities_7patch_with_backdistributed_gannets_NB.csv",
  stringsAsFactors = FALSE
) %>%
  
  filter(
    patch_name %in% patches,
    CommonName %in% species_list
  ) %>%
  
  mutate(
    TotalObserved = as.numeric(TotalObserved),
    Date = as.Date(DateObserved,
                   tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))
  ) %>%
  
  group_by(patch_name, CommonName) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    time = as.integer(Date - min(Date)) + 1
  ) %>%
  ungroup() %>%
  
  # ------------------------------------------------------------
# Remove early observation outliers
# ------------------------------------------------------------
filter(
  !(patch_name == "Patch5" & CommonName == "Common Murre" & time == 1),
  !(patch_name == "Patch6" & CommonName == "Common Murre" & time == 1),
  !(patch_name == "Patch6" & CommonName == "Common Murre" & time == 18),
  !(patch_name == "Patch7" & CommonName == "Northern Gannet" & time == 1),
  !(patch_name == "Patch1" & CommonName == "Northern Gannet" & time == 1),
  !(patch_name == "Patch2" & CommonName == "Northern Gannet" & time %in% c(1, 5, 13)),
  !(patch_name == "Patch3" & CommonName == "Northern Gannet" & time %in% c(1, 8))
) %>%
  
  # ------------------------------------------------------------
# Aggregate multiple observations per day
# ------------------------------------------------------------
group_by(patch_name, CommonName, Date) %>%
  summarise(
    TotalObserved = sum(TotalObserved, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # ------------------------------------------------------------
# Recompute time + cumulative counts 
# ------------------------------------------------------------
group_by(patch_name, CommonName) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    time = as.integer(Date - min(Date)) + 1,
    CumObserved = cumsum(TotalObserved)
  ) %>%
  ungroup() %>%
  arrange(patch_name, CommonName, Date)
time_lookup <- obs_all %>%
  select(patch_name, CommonName, time, Date) %>%
  distinct()

obs_all_fit <- obs_all
time_steps_fit <- max(obs_all_fit$time)

# ============================================================
# Underreporting scenarios
# ============================================================

reporting_scenarios <- c(
  low_reporting    = 0.2,
  medium_reporting = 0.5,
  high_reporting   = 0.8
)

# ============================================================
# Fitting loop with reporting scenarios
# ============================================================

all_fits <- list()

for (scenario in names(reporting_scenarios)) {
  
  rho <- reporting_scenarios[[scenario]]
  
  for (p in patches) {
    for (sp in species_list) {
      
      dat <- obs_all_fit %>%
        filter(patch_name == p, CommonName == sp)
      
      if (nrow(dat) == 0) next
      
      obsdat <- dat %>%
        mutate(
          patch_id = str_remove(patch_name, "Patch"),
          matrix = paste0("D_obs", patch_id)
        ) %>%
        select(time, matrix, value = CumObserved) %>%
        distinct()
      
      i <- as.integer(str_remove(p, "Patch"))
      
      K_vals <- K_list[[sp]][p]
      
      D0 <- dat %>%
        arrange(time) %>%
        slice(1) %>%
        pull(CumObserved)
      
      flows <- list(
        mp_per_capita_flow(
          paste0("S", i),
          paste0("E", i),
          paste0("beta", i, " * I", i, " / N", i),
          paste0("exposure", i)
        ),
        
        mp_per_capita_flow(
          paste0("E", i),
          paste0("I", i),
          "alpha",
          paste0("infection", i)
        ),
        
        mp_per_capita_flow(
          paste0("I", i),
          paste0("R", i),
          "gamma",
          paste0("recovery", i)
        ),
        
        mp_per_capita_flow(
          paste0("I", i),
          paste0("D", i),
          "muI",
          paste0("death", i)
        )
      )
      
      flows_beta_map <- list(
        as.formula(
          paste0("beta", i, "_thing ~ beta", i)
        )
      )
      
      reporting_map <- list(
        as.formula(
          paste0("D_obs", i, " ~ rho * D", i)
        )
      )
      
      initialize_state <- list(
        as.formula(
          paste0("D", i, " ~ D0")
        ),
        
        as.formula(
          paste0(
            "S", i, " ~ N", i,
            " - E", i,
            " - I", i,
            " - R", i,
            " - D", i
          )
        )
      )
      
      default <- c(
        setNames(0.2, paste0("beta", i)),
        alpha = as.numeric(alpha),
        gamma = as.numeric(gamma),
        muI   = as.numeric(muI_species[[sp]]),
        rho   = as.numeric(rho),
        
        setNames(as.numeric(K_vals), paste0("N", i)),
        
        D0 = as.numeric(D0),
        
        setNames(5, paste0("E", i)),
        setNames(2, paste0("I", i)),
        setNames(0, paste0("R", i)),
        setNames(0, paste0("D", i))
      )
      
      spec <- mp_tmb_model_spec(
        before = initialize_state,
        during = c(
          flows,
          flows_beta_map,
          reporting_map
        ),
        default = default
      )
      
      X <- diag(time_steps_fit)
      
      timevar_spec <- mp_tmb_insert_glm_timevar(
        spec,
        paste0("beta", i),
        X,
        rep(0, time_steps_fit),
        link_function = mp_log
      )
      
      death_traj <- setNames(
        list(mp_pois()),
        paste0("D_obs", i)
      )
      
      cal <- mp_tmb_calibrator(
        spec = timevar_spec |> mp_rk4(),
        data = obsdat,
        time = mp_sim_bounds(1, time_steps_fit),
        traj = death_traj,
        default = default,
        par = setNames(
          list(mp_norm(0, 0.3)),
          paste0("time_var_beta", i)
        ),
        outputs = c(
          paste0("beta", i, "_thing"),
          paste0("D_obs", i)
        )
      )
      
      mp_optimize(cal)
      
      res <- mp_trajectory_sd(
        cal,
        conf.int = TRUE
      ) %>%
        mutate(
          patch = p,
          species = sp,
          scenario = scenario,
          rho = rho
        )
      
      all_fits[[paste(
        scenario,
        p,
        sp,
        sep = "_"
      )]] <- res
      
    }
  }
}

fitted_data <- bind_rows(all_fits)

# ============================================================
# Plotting prep 
# ============================================================

time_lookup_fit <- time_lookup %>%
  rename(
    patch = patch_name,
    species = CommonName
  )

obs_plot <- obs_all_fit %>%
  rename(
    patch = patch_name,
    species = CommonName
  )

plot_data <- fitted_data %>%
  mutate(
    patch = paste0("Patch", str_extract(matrix, "\\d+")),
    variable = if_else(
      str_detect(matrix, "^D_obs\\d+$"),
      "Deaths",
      "Beta"
    ),
    value_pos = pmax(value, 0),
    ci_low = pmax(conf.low, 0),
    ci_high = pmax(conf.high, 0)
  ) %>%
  left_join(
    time_lookup_fit,
    by = c("patch", "species", "time")
  )

# ============================================================
# β-only dataset for distributions
# ============================================================

time_grid_beta <- expand.grid(
  Date = sort(unique(plot_data$Date)),
  species = unique(plot_data$species),
  patch = unique(plot_data$patch),
  scenario = unique(plot_data$scenario)
)

plot_data_beta <- time_grid_beta %>%
  left_join(
    plot_data %>%
      filter(variable == "Beta"),
    by = c("Date", "species", "patch", "scenario")
  ) %>%
  mutate(
    patch_label = factor(
      patch,
      levels = names(custom_patch_names),
      labels = custom_patch_names
    )
  ) %>%
  filter(
    !is.na(Date),
    !is.na(value_pos)
  )

# ============================================================
# Global scaling 
# ============================================================

beta_ymax <- quantile(
  plot_data_beta$ci_high,
  0.95,
  na.rm = TRUE
)

plot_data_beta <- plot_data_beta %>%
  mutate(
    beta_raw      = value_pos,
    beta_analysis = pmax(beta_raw, 0)
  )

# ============================================================
# Distribution dataset
# ============================================================

beta_dist <- plot_data_beta %>%
  transmute(
    species,
    patch_label,
    scenario,
    beta = beta_analysis
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("low_reporting", "medium_reporting", "high_reporting"),
      labels = c("low", "medium", "high")
    )
  )

y_fmt <- scale_y_continuous(
  limits = c(0, beta_ymax),
  expand = expansion(mult = c(0, 0.05))
)

no_x_axis <- scale_x_discrete(
  breaks = NULL,
  labels = NULL
)

# ============================================================
# Fig 3.5A: Species distribution
# ============================================================

p_beta_species <- ggplot(
  beta_dist,
  aes(scenario, beta, fill = scenario)
) +
  
  geom_violin(alpha = 0.4) +
  
  geom_boxplot(
    width = 0.15,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  
  facet_wrap(~species) +
  
  y_fmt +
  no_x_axis +
  
  labs(
    y = "Transmission rate <i>β(t)</i>",
    x = NULL
  ) +
  
  beta_theme

ggsave(
  "Fig_3.5A_underreporting.png",
  p_beta_species,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# Fig 3.5B: Patch distribution
# ============================================================

p_beta_patch <- ggplot(
  beta_dist,
  aes(scenario, beta, fill = scenario)
) +
  
  geom_violin(alpha = 0.4) +
  
  geom_boxplot(
    width = 0.15,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  
  facet_wrap(~patch_label, ncol = 3) +
  
  y_fmt +
  no_x_axis +
  
  labs(
    y = "Transmission rate <i>β(t)</i>",
    x = NULL
  ) +
  
  beta_theme +
  
  theme(
    legend.position = c(0.5, 0.25),   # moved up inside plot
    legend.justification = c(0.5, 1),
    legend.background = element_blank()
  )

ggsave(
  "Fig_3.5B_underreporting.png",
  p_beta_patch,
  width = 12,
  height = 8,
  dpi = 300
)
