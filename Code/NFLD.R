# ============================================================
# Spatial distribution of HPAI-associated seabird mortality in Newfoundland (2022)
# ============================================================

# ============================================================
# Libraries
# ============================================================
library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(gganimate)
library(gifski)
library(transformr)
library(cowplot)
library(cli)
library(purrr)
library(ggrepel)
library(grid)

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_figures")
sf::sf_use_s2(FALSE)

# ============================================================
# Base map (Newfoundland & Labrador)
# ============================================================
map_nl <- ne_states(country = "canada", returnclass = "sf") %>%
  filter(name == "Newfoundland and Labrador")

bbox_nl_sfc <- st_as_sfc(st_bbox(map_nl))

# ============================================================
# Mortality data
# ============================================================
mortality_data_path <- "C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_code/Data S1. Reported mortalities and morbidities in eastern Canada 2022.csv"

filtered_data <- read.csv(mortality_data_path, stringsAsFactors = FALSE) %>%
  filter(
    Prov == "NL",
    Double_Count_Analysis_Dataset == "retain",
    CommonName %in% c("Northern Gannet","Common Murre")
  ) %>%
  mutate(
    DateObserved = as.Date(DateObserved, "%m/%d/%Y"),
    Species = factor(CommonName,
                     levels = c("Northern Gannet","Common Murre"))
  )

# ============================================================
# IBA polygons
# ============================================================
my_kml_nl <- st_read(
  "C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_code/Important Bird and Biodiversity Areas of Canada  Zones importantes pour la conservation des oiseaux et de la biodiversité du Can.kml",
  quiet = TRUE
) %>%
  st_make_valid() %>%
  filter(grepl("^(NF|LB)", Name)) %>%
  st_intersection(bbox_nl_sfc)

# ---- patch definitions ----
patch_groups <- list(
  Patch1 = c("NF001","NF028"),
  Patch2 = c("NF015","NF024","NF002"),
  Patch3 = c("NF022","NF021","NF003","NF019"),
  Patch4 = c("NF025","NF013","NF004"),
  Patch5 = c("NF010","NF009","NF008"),
  Patch6 = c("NF045"),
  Patch7 = c("NF036","NF032","NF031","NF030")
)

# assign patches
my_kml_nl$patch_name <- map_chr(my_kml_nl$Name, function(nm){
  p <- names(patch_groups)[map_lgl(patch_groups, ~ nm %in% .x)]
  if(length(p) == 0) NA else p
})

# union polygons
my_kml_nl <- my_kml_nl %>%
  filter(!is.na(patch_name)) %>%
  group_by(patch_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(
    patch_id = factor(row_number()),
    color    = hue_pal()(n())
  )

patch_colors <- setNames(my_kml_nl$color, my_kml_nl$patch_name)

# ============================================================
# Seabird ecological reserves
# ============================================================
seabird_sf <- data.frame(
  Name = c("Cape St. Mary's","Witless Bay","Baccalieu Island",
           "Funk Island","Hare Bay","Lawn Islands"),
  lon  = c(-54.20,-52.83,-52.80,-53.18,-55.92,-55.62),
  lat  = c(46.82,47.28,48.14,49.76,51.25,46.87)
) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# west coast label 
west_coast_sf <- my_kml_nl %>%
  filter(patch_name == "Patch6") %>%
  st_point_on_surface() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  transmute(Name = "West Coast",
            lon = X - 0.15,
            lat = Y) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

label_points <- bind_rows(seabird_sf, west_coast_sf)

# ============================================================
# IBA map
# ============================================================
p_iba <- ggplot() +
  geom_sf(data = map_nl, fill = "grey92", color = "black", linewidth = 0.5) +
  geom_sf(data = my_kml_nl, aes(fill = patch_name),
          color = "black", alpha = 0.8, linewidth = 0.3) +
  geom_sf(data = label_points, color = "black", size = 2.5) +
  geom_label_repel(
    data = label_points,
    aes(label = Name, geometry = geometry),
    stat = "sf_coordinates",
    size = 5, fontface = "bold",
    fill = "white", label.size = 0.2,
    box.padding = 0.6, point.padding = 0.5,
    segment.color = NA
  ) +
  coord_sf(xlim = c(-60,-52), ylim = c(46,52), expand = FALSE) +
  scale_fill_manual(values = patch_colors, guide = "none") +
  scale_x_continuous(breaks = seq(-60,-52,2)) +
  scale_y_continuous(breaks = seq(46,52,1)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

ggsave("IBAs.png", plot = p_iba, width = 8, height = 6, dpi = 600, bg = "white")

# ============================================================
# Points → sf + patch assignment
# ============================================================
points_sf <- st_as_sf(filtered_data,
                      coords = c("Long","Lat"),
                      crs = st_crs(my_kml_nl))

idx <- st_nearest_feature(points_sf, my_kml_nl)

points_sf <- points_sf %>%
  mutate(
    patch_id   = my_kml_nl$patch_id[idx],
    patch_name = my_kml_nl$patch_name[idx],
    color      = my_kml_nl$color[idx]
  )

# ---- mortality bins ----
size_values <- c("1–10"=3,"10–50"=5,"50–100"=7,"100–500"=9,"500–1000"=11,"1000+"=13)

points_sf <- points_sf %>%
  mutate(
    MortBin = cut(TotalObserved,
                  breaks = c(1,10,50,100,500,1000,Inf),
                  labels = names(size_values),
                  include.lowest = TRUE)
  ) %>%
  filter(!is.na(MortBin))

# ---- export for model ----
points_model_df <- points_sf %>%
  st_drop_geometry() %>%
  mutate(
    Long = st_coordinates(points_sf)[,1],
    Lat  = st_coordinates(points_sf)[,2]
  ) %>%
  select(DateObserved, CommonName, TotalObserved,
         MortBin, patch_id, patch_name, Long, Lat)

write.csv(points_model_df, "NFLD_mortalities_7patch.csv", row.names = FALSE)

# ============================================================
# Animation (cumulative)
# ============================================================
points_sf <- points_sf %>%
  arrange(DateObserved) %>%
  mutate(anim_id = row_number())

p_anim <- ggplot() +
  geom_sf(data = map_nl, fill = "grey95") +
  geom_sf(data = my_kml_nl, aes(fill = patch_name),
          alpha = 0.35, color = "black") +
  geom_sf(data = points_sf,
          aes(color = patch_name, shape = Species, size = MortBin),
          alpha = 0.85) +
  scale_fill_manual(values = patch_colors, guide = "none") +
  scale_color_manual(values = patch_colors, guide = "none") +
  scale_shape_manual(values = c("Northern Gannet"=16,"Common Murre"=17)) +
  scale_size_manual(values = size_values, name = "Mortalities") +
  scale_x_continuous(breaks = seq(-60,-52,2)) +
  scale_y_continuous(breaks = seq(46,52,1)) +
  coord_sf(xlim = c(-60,-52), ylim = c(46,52), expand = FALSE) +
  theme_bw() +
  transition_manual(frames = DateObserved, cumulative = TRUE) +
  labs(title = "{current_frame}")

animate(p_anim, fps = 5, duration = 15, width = 800, height = 700)
anim_save("NFLD_Mortalities_7patch_Cumulative.gif")