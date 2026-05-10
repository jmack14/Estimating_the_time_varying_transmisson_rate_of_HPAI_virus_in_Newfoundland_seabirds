# ============================================================
# SEIRD model for the spread of HPAI virus in NFLD seabirds with inter-patch connectivity
# ============================================================

# ============================================================
# Libraries
# ============================================================
library(macpan2)
library(tidyverse)
library(patchwork)
library(zoo)
library(ggtext)
library(stringr)
library(splines)
library(sf)
library(geosphere)
library(reshape2)
library(viridis)

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_figures")
theme_set(theme_bw())

# ============================================================
# Parameters
# ============================================================
alpha     <- 1/7
gamma     <- 1/11
DEATH_CAP <- 100
BURN_IN   <- 7
epsilon   <- 0.01

patches <- paste0("Patch", 1:7)
n_patches <- length(patches)

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
  "Northern Gannet" = c(30000, 5000, 2000, 12000, 1000, 5000, 1000),
  "Common Murre"    = c(20000, 155000, 5000, 792000, 5000, 10000, 15000)
)
names(K_list[[1]]) <- names(K_list[[2]]) <- patches

# ============================================================
# Data prep
# ============================================================
obs_all <- read.csv("NFLD_mortalities_7patch.csv") %>%
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

obs_all_fit <- obs_all %>% filter(time > BURN_IN)
time_steps_fit <- max(obs_all_fit$time)

# ============================================================
# Observation formatting
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

obsdat_list <- set_names(
  map(species_list, ~ make_obsdat(filter(obs_all_fit, CommonName == .x))),
  species_list
)

# ============================================================
# Connectivity matrix
# ============================================================
patch_coords <- my_kml_nl %>%
  st_transform(4326) %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  ) %>%
  arrange(patch_name)

coords <- cbind(patch_coords$lon, patch_coords$lat)

dist_mat <- geosphere::distm(coords) / 1000

M <- 1 / (dist_mat^2 + 1e-6)
diag(M) <- 0
M <- M / apply(M, 1, max)

# ============================================================
# Model specification
# ============================================================
flows_main <- map(1:n_patches, function(i){
  
  infect_term <- paste0(
    "I", i,
    " + epsilon * (",
    paste0("(", M[i, -i], " * I", setdiff(1:n_patches, i), ")", collapse = " + "),
    ")"
  )
  
  list(
    mp_per_capita_flow(
      paste0("S", i), paste0("E", i),
      paste0("beta", i, " * (", infect_term, ") / N", i),
      paste0("exposure", i)
    ),
    mp_per_capita_flow(paste0("E", i), paste0("I", i), "alpha",
                       paste0("infection", i)),
    mp_per_capita_flow(paste0("I", i), paste0("R", i), "gamma",
                       paste0("recovery", i)),
    mp_per_capita_flow(paste0("I", i), paste0("D", i), "muI",
                       paste0("death", i))
  )
}) %>% flatten()

flows_beta_map <- map(1:n_patches,
                      ~ as.formula(paste0("beta", .x, "_thing ~ beta", .x)))

flows <- c(flows_main, flows_beta_map)

# ============================================================
# Defaults
# ============================================================
default <- c(
  setNames(rep(0.2, n_patches), paste0("beta", 1:n_patches)),
  alpha = alpha,
  gamma = gamma,
  muI   = 0.1,
  epsilon = epsilon,
  setNames(K_list[[1]], paste0("N", 1:n_patches)),
  setNames(rep(5, n_patches), paste0("E", 1:n_patches)),
  setNames(rep(2, n_patches), paste0("I", 1:n_patches)),
  setNames(rep(0, n_patches), paste0("R", 1:n_patches)),
  setNames(rep(0, n_patches), paste0("D", 1:n_patches))
)

initialize_state <- map(1:n_patches, ~ as.formula(
  paste0("S", .x, " ~ N", .x, " - E", .x, " - I", .x, " - R", .x, " - D", .x)
))

spec <- mp_tmb_model_spec(
  before = initialize_state,
  during = flows,
  default = default
)

# ============================================================
# Time-varying beta
# ============================================================
basis_cols <- 15
t_scaled <- seq_len(time_steps_fit) / time_steps_fit
X <- ns(t_scaled, df = basis_cols, intercept = TRUE)

timevar_spec <- reduce(1:n_patches, function(sp, i){
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
  
  N_defaults <- setNames(K_list[[species_name]], paste0("N", 1:n_patches))
  
  death_traj <- setNames(map(unique(obsdat$matrix), ~ mp_pois()),
                         unique(obsdat$matrix))
  
  cal <- mp_tmb_calibrator(
    spec = timevar_spec |> mp_rk4(),
    data = obsdat,
    time = mp_sim_bounds(1, time_steps_fit),
    traj = death_traj,
    
    default = c(
      list(muI = muI_value, epsilon = epsilon),
      as.list(N_defaults)
    ),
    
    par = setNames(rep(list(mp_norm(0, 0.3)), n_patches),
                   paste0("time_var_beta", 1:n_patches)),
    
    outputs = c(
      paste0("beta", 1:n_patches, "_thing"),
      paste0("death", 1:n_patches)
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
# Connectivity heatmap
# ============================================================

# 1. Long format
M_df <- melt(M)
colnames(M_df) <- c("ReceiverPatch", "SourcePatch", "Weight")

# ensure numeric indexing safety
M_df$ReceiverPatch <- as.integer(M_df$ReceiverPatch)
M_df$SourcePatch   <- as.integer(M_df$SourcePatch)

# 2. Geographic ordering
patch_names_geo <- c(
  "West Coast",
  "Hare Bay",
  "Funk Island",
  "Baccalieu Island",
  "Witless Bay",
  "Cape St. Mary's",
  "Lawn Islands"
)

original_names <- unname(custom_patch_names)

M_df$ReceiverPatch <- original_names[M_df$ReceiverPatch]
M_df$SourcePatch   <- original_names[M_df$SourcePatch]

M_df$ReceiverPatch <- factor(M_df$ReceiverPatch, levels = patch_names_geo)
M_df$SourcePatch   <- factor(M_df$SourcePatch,   levels = patch_names_geo)

# 3. Plot
conn_plot <- ggplot(M_df, aes(x = SourcePatch, y = ReceiverPatch, fill = Weight)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Connectivity") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )

ggsave("Intercolony_Connectivity_Heatmap.png",
       conn_plot, width = 8, height = 6, dpi = 300)

conn_plot