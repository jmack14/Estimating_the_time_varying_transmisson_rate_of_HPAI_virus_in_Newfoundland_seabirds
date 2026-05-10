# ============================================================
# SEIRD model for the spread of HPAI virus in NFLD seabirds in 2022
# ============================================================

# ============================================================
# Libraries & global settings
# ============================================================
library(macpan2)
library(tidyverse)
library(patchwork)
library(zoo)
library(ggtext)
library(stringr)
library(splines)

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_figures")
theme_set(theme_bw())

# ============================================================
# Parameters
# ============================================================
alpha     <- 1/7
gamma     <- 1/11
DEATH_CAP <- 100
BURN_IN   <- 7

patches      <- paste0("Patch", 1:7)
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
obs_all <- read.csv("NFLD_mortalities_7patch.csv", stringsAsFactors = FALSE) %>%
  filter(patch_name %in% patches,
         CommonName %in% species_list) %>%
  mutate(
    TotalObserved = pmin(TotalObserved, DEATH_CAP),
    Date = as.Date(DateObserved)
  ) %>%
  group_by(patch_name, CommonName) %>%
  arrange(Date) %>%
  ungroup() %>%
  group_by(patch_name) %>%
  mutate(time = as.integer(Date - min(Date)) + 1) %>%
  ungroup()

time_lookup <- obs_all %>%
  select(patch_name, CommonName, time, Date) %>%
  distinct()

# ---- burn-in removal ----
obs_all_fit   <- obs_all %>% filter(time > BURN_IN)
time_steps_fit <- max(obs_all_fit$time)

# ============================================================
# Observation format
# ============================================================
make_obsdat <- function(dat){
  dat %>%
    mutate(
      patch_id = str_remove(patch_name, "Patch"),
      matrix   = paste0("death", patch_id)
    ) %>%
    select(time, matrix, value = TotalObserved) %>%
    distinct()
}

obsdat_list <- species_list %>%
  set_names() %>%
  map(~ make_obsdat(filter(obs_all_fit, CommonName == .x)))

# ============================================================
# Model specification
# ============================================================
flows_main <- map(1:7, function(i){
  list(
    mp_per_capita_flow(paste0("S", i), paste0("E", i),
                       paste0("beta", i, " * I", i, " / N", i),
                       paste0("exposure", i)),
    mp_per_capita_flow(paste0("E", i), paste0("I", i), "alpha",
                       paste0("infection", i)),
    mp_per_capita_flow(paste0("I", i), paste0("R", i), "gamma",
                       paste0("recovery", i)),
    mp_per_capita_flow(paste0("I", i), paste0("D", i), "muI",
                       paste0("death", i))
  )
}) %>% flatten()

flows_beta_map <- map(1:7, ~ as.formula(paste0("beta", .x, "_thing ~ beta", .x)))
flows <- c(flows_main, flows_beta_map)

default <- c(
  setNames(rep(0.2, 7), paste0("beta", 1:7)),
  alpha = alpha,
  gamma = gamma,
  muI   = 0.1,
  setNames(as.numeric(K_list[["Northern Gannet"]]), paste0("N", 1:7)),
  setNames(rep(5, 7), paste0("E", 1:7)),
  setNames(rep(2, 7), paste0("I", 1:7)),
  setNames(rep(0, 7), paste0("R", 1:7)),
  setNames(rep(0, 7), paste0("D", 1:7))
)

initialize_state <- map(1:7, ~ as.formula(
  paste0("S", .x, " ~ N", .x, " - E", .x, " - I", .x, " - R", .x, " - D", .x)
))

spec <- mp_tmb_model_spec(
  before  = initialize_state,
  during  = flows,
  default = default
)

# ============================================================
# Time-varying beta
# ============================================================
basis_cols <- 15
t_scaled   <- seq_len(time_steps_fit) / time_steps_fit
X <- ns(t_scaled, df = basis_cols, intercept = TRUE)

timevar_spec <- reduce(1:7, function(sp, i){
  mp_tmb_insert_glm_timevar(
    sp,
    paste0("beta", i),
    X,
    rep(0, basis_cols),
    link_function = mp_log
  )
}, .init = spec)

# ============================================================
# Calibration
# ============================================================
fit_species <- function(obsdat, muI_value, species_name){
  
  N_defaults <- setNames(
    as.numeric(K_list[[species_name]]),
    paste0("N", 1:7)
  )
  
  death_traj <- setNames(
    map(unique(obsdat$matrix), ~ mp_pois()),
    unique(obsdat$matrix)
  )
  
  cal <- mp_tmb_calibrator(
    spec = timevar_spec |> mp_rk4(),
    data = obsdat,
    time = mp_sim_bounds(1, time_steps_fit),
    traj = death_traj,
    
    default = c(list(muI = muI_value), as.list(N_defaults)),
    
    par = setNames(
      rep(list(mp_norm(0, 0.3)), 7),
      paste0("time_var_beta", 1:7)
    ),
    
    outputs = c(
      paste0("beta", 1:7, "_thing"),
      paste0("death", 1:7)
    )
  )
  
  mp_optimize(cal)
  mp_trajectory_sd(cal, conf.int = TRUE)
}

# ============================================================
# Fit model 
# ============================================================
fitted_data <- map_dfr(species_list, function(s){
  fit_species(obsdat_list[[s]], muI_species[[s]], s) %>%
    mutate(species = s)
})

# ============================================================
# Plot prep
# ============================================================
time_lookup_fit <- time_lookup %>% filter(time > BURN_IN)

obs_plot <- obs_all_fit %>%
  rename(patch = patch_name,
         species = CommonName)

plot_data <- fitted_data %>%
  mutate(
    patch    = paste0("Patch", str_extract(matrix, "\\d")),
    variable = if_else(str_detect(matrix, "death"), "Deaths", "Beta"),
    value_pos = pmax(value, 0),
    ci_low    = pmax(conf.low, 0),
    ci_high   = pmax(conf.high, 0)
  ) %>%
  left_join(time_lookup_fit,
            by = c("patch" = "patch_name",
                   "species" = "CommonName",
                   "time" = "time"))

# ============================================================
# Themes
# ============================================================
deaths_theme <- theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    strip.text   = element_text(size = 20, face = "bold"),
    legend.title = element_blank()
  )

beta_theme <- theme_bw(base_size = 18) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 20),
    strip.text   = element_text(size = 18, face = "bold"),
    legend.title = element_blank()
  )

# ============================================================
# Deaths plot
# ============================================================
p_deaths <- ggplot() +
  geom_ribbon(
    data = filter(plot_data, variable == "Deaths"),
    aes(Date,
        ymin = ci_low + 1,
        ymax = pmin(ci_high, quantile(ci_high, 0.95, na.rm = TRUE)) + 1,
        fill = species),
    alpha = 0.25
  ) +
  geom_line(
    data = filter(plot_data, variable == "Deaths"),
    aes(Date, value_pos + 1, color = species),
    linewidth = 1.3
  ) +
  geom_point(
    data = obs_plot,
    aes(Date, TotalObserved + 1, color = species),
    size = 3
  ) +
  facet_wrap(~patch, scales = "free",
             labeller = labeller(patch = custom_patch_names)) +
  scale_y_log10() +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(y = "Daily deaths") +
  deaths_theme

# ============================================================
# Beta plot 
# ============================================================
plot_data_beta <- plot_data %>%
  filter(variable == "Beta") %>%
  mutate(
    patch_label = factor(patch,
                         levels = names(custom_patch_names),
                         labels = custom_patch_names)
  )

ylim_df <- summarize(plot_data_beta,
                     ymax = quantile(ci_high, 0.95, na.rm = TRUE))

july_start <- as.Date(
  paste0(format(min(plot_data_beta$Date, na.rm = TRUE), "%Y"), "-07-01")
)

p_beta <- ggplot(plot_data_beta, aes(Date, value_pos)) +
  geom_ribbon(aes(ymin = ci_low,
                  ymax = pmin(ci_high, ylim_df$ymax),
                  fill = patch_label),
              alpha = 0.2) +
  geom_line(aes(color = patch_label)) +
  facet_wrap(~species, scales = "free_y", ncol = 1) +
  coord_cartesian(
    xlim = c(july_start, max(plot_data_beta$Date, na.rm = TRUE)),
    ylim = c(0, ylim_df$ymax)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "Transmission rate <i>β(t)</i>") +
  beta_theme

# ============================================================
# Beta distributions
# ============================================================
beta_dist <- plot_data_beta %>%
  filter(Date >= july_start) %>%
  mutate(
    beta = pmax(value_pos, 0),
    beta = pmin(beta, ylim_df$ymax)
  ) %>%
  select(species, patch_label, beta)

y_fmt <- scale_y_continuous(
  limits = c(-0.05 * ylim_df$ymax, ylim_df$ymax),
  expand = expansion(mult = c(0, 0.05))
)

tag_theme_A <- theme(
  plot.tag = element_text(size = 20, face = "bold"),
  plot.tag.position = c(0.01, 1.01),
  plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
)

tag_theme_B <- theme(
  plot.tag = element_text(size = 20, face = "bold"),
  plot.tag.position = c(0.01, 0.98),
  plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
)

# ---- A: species ----
p_beta_species <- ggplot(beta_dist,
                         aes(species, beta, fill = species)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = species_colors) +
  y_fmt +
  labs(y = "Transmission rate <i>β(t)</i>", tag = "A") +
  beta_theme +
  theme(legend.position = "none") +
  tag_theme_A +
  coord_cartesian(clip = "off")

# ---- B: patch ----
p_beta_patch <- ggplot(beta_dist,
                       aes(species, beta, fill = species)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~patch_label, nrow = 2) +
  scale_fill_manual(values = species_colors) +
  y_fmt +
  labs(y = "Transmission rate <i>β(t)</i>", tag = "B") +
  beta_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.76, 0.2),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.75, "cm"),
    legend.spacing.y = unit(0.25, "cm"),
    
    panel.spacing = unit(1.2, "lines"),
    aspect.ratio = 1
  ) +
  tag_theme_B +
  coord_cartesian(clip = "off")

# ============================================================
# Beta summary for each species 
# ============================================================

beta_summary_species <- beta_dist %>%
  group_by(species) %>%
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    q1 = quantile(beta, 0.25, na.rm = TRUE),
    median = median(beta, na.rm = TRUE),
    q3 = quantile(beta, 0.75, na.rm = TRUE),
    IQR = IQR(beta, na.rm = TRUE)
  )

beta_summary_species

# ============================================================
# Save outputs
# ============================================================

ggsave(
  filename = "SEIRD_7patch_deaths.png",
  plot     = p_deaths,
  width    = 16,
  height   = 10,
  dpi      = 300
)

ggsave(
  filename = "SEIRD_7patch_beta.png",
  plot     = p_beta,
  width    = 10,
  height   = 6,
  dpi      = 300
)

ggsave(
  filename = "beta_species_plot.png",
  plot     = p_beta_species,
  width    = 8,
  height   = 5,
  dpi      = 300
)

ggsave(
  filename = "beta_patch_plot.png",
  plot     = p_beta_patch,
  width    = 12,
  height   = 8,
  dpi      = 300
)
