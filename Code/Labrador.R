# ============================================================
# Proposed Labrador patches
# ============================================================

# ============================================================
# Libraries and setup
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(rnaturalearth)

sf::sf_use_s2(FALSE)

# ============================================================
# Base map
# ============================================================
canada_map <- ne_states(country = "canada", returnclass = "sf")
map_nl <- canada_map %>% filter(name == "Newfoundland and Labrador")

# ============================================================
# Read IBA polygons
# ============================================================
my_kml <- st_read(
  "C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Chp_2/Chp2_code/Important Bird and Biodiversity Areas of Canada  Zones importantes pour la conservation des oiseaux et de la biodiversité du Can.kml",
  quiet = TRUE
) %>%
  st_make_valid() %>%
  filter(grepl("^(NF|LB)", Name))

# ============================================================
# Patch definitions
# ============================================================
patch_groups <- list(
  "Patch1" = c("NF001","NF028"),
  "Patch2" = c("NF015","NF024","NF002"),
  "Patch3" = c("NF022","NF021","NF003","NF019"),
  "Patch4" = c("NF025","NF013","NF004"),
  "Patch5" = c("NF010","NF009","NF008"),
  "Patch6" = "NF045",
  "Patch7" = c("NF036","NF032","NF031","NF030"),
  
  "Patch8_SevenIslandsBay"  = c("LB003","LB024"),
  "Patch9_NainCoastline"    = c("LB006","LB005"),
  "Patch10_QuakerHatIsland" = "LB009",
  "Patch11_NE_GroswaterBay" = c("LB012","LB011","LB013","LB025","LB020"),
  "Patch12_GannetIslands"   = c("LB001","LB019","LB027"),
  "Patch13_StPeterBay"      = "LB023",
  "Patch14_PointAmour"      = "LB022"
)

# Assign patch names (UNCHANGED LOGIC)
my_kml$patch_name <- sapply(my_kml$Name, function(nm) {
  patch <- names(patch_groups)[sapply(patch_groups, function(x) nm %in% x)]
  if (length(patch) == 0) NA else patch
})

# ============================================================
# Union patches
# ============================================================
my_kml_union <- my_kml %>%
  filter(!is.na(patch_name)) %>%
  group_by(patch_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

nf_union <- my_kml_union %>% filter(grepl("^Patch[1-7]$", patch_name))
lb_union <- my_kml_union %>% filter(grepl("^Patch(8|9|10|11|12|13|14)", patch_name))

# ============================================================
# Points
# ============================================================

# Labrador (centroids — unchanged method)
lb_points <- st_point_on_surface(lb_union)

# Newfoundland fixed points (unchanged)
nf_points <- data.frame(
  Name = c("Cape St. Mary's", "Witless Bay", "Baccalieu Island",
           "Funk Island", "Hare Bay", "Lawn Islands"),
  lon = c(-54.20, -52.83, -52.80, -53.18, -55.92, -55.62),
  lat = c(46.82, 47.28, 48.14, 49.76, 51.25, 46.87)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# West Coast label (unchanged logic)
west_coast_point <- nf_union %>%
  filter(patch_name == "Patch6") %>%
  st_point_on_surface()

wc_coords <- st_coordinates(west_coast_point)

west_coast_sf <- data.frame(
  Name = "West Coast",
  lon = wc_coords[,1] - 0.15,
  lat = wc_coords[,2]
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# ============================================================
# Labrador patch names
# ============================================================
labrador_name_map <- c(
  "Patch8_SevenIslandsBay"  = "Seven Islands Bay",
  "Patch9_NainCoastline"    = "Nain Coastline",
  "Patch10_QuakerHatIsland" = "Quaker Hat Island",
  "Patch11_NE_GroswaterBay" = "NE Groswater Bay",
  "Patch12_GannetIslands"   = "Gannet Islands",
  "Patch13_StPeterBay"      = "St. Peter Bay",
  "Patch14_PointAmour"      = "Point Amour"
)

lb_labels <- lb_points %>%
  mutate(Name = labrador_name_map[patch_name])

label_points <- bind_rows(nf_points, lb_labels, west_coast_sf)

# ============================================================
# Colors
# ============================================================
patch_colors <- setNames(
  scales::hue_pal()(nrow(my_kml_union)),
  my_kml_union$patch_name
)

# ============================================================
# Plot proposed Labrador patches: Fig 3.4
# ============================================================
p_full <- ggplot() +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  
  geom_sf(data = map_nl,
          fill = "grey92",
          color = "black",
          linewidth = 0.5) +
  
  geom_sf(data = my_kml_union,
          aes(fill = patch_name),
          color = "black",
          alpha = 0.6,
          linewidth = 0.2) +
  
  geom_sf(data = bind_rows(lb_points, nf_points, west_coast_sf),
          color = "black",
          size = 2.8) +
  
  geom_label_repel(
    data = label_points,
    aes(label = Name, geometry = geometry),
    stat = "sf_coordinates",
    size = 4.5,
    fontface = "bold",
    fill = "white",
    color = "black",
    label.size = 0.2,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 3,
    max.iter = 5000,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  coord_sf(
    xlim = c(-68.5, -52),
    ylim = c(46, 60.5),
    expand = FALSE
  ) +
  
  scale_fill_manual(values = patch_colors, guide = "none")

# ============================================================
# Save
# ============================================================
ggsave(
  "Fig_3.4.png",
  p_full,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)