# ============================================================
# Connectivity matrix + heatmap 
# ============================================================

# ============================================================
# Libraries
# ============================================================
library(sf)
library(geosphere)
library(reshape2)
library(viridis)

# ------------------------------------------------------------
# Patch coordinates 
# ------------------------------------------------------------
patch_coords <- my_kml_nl %>%
  st_transform(4326) %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  ) %>%
  arrange(patch_name)

coords <- cbind(patch_coords$lon, patch_coords$lat)

# ------------------------------------------------------------
# Distance-based connectivity
# ------------------------------------------------------------
dist_mat <- geosphere::distm(coords) / 1000  # km

M <- 1 / (dist_mat^2 + 1e-6)
diag(M) <- 0
M <- M / apply(M, 1, max)   # row-normalized connectivity

# ------------------------------------------------------------
# Convert to long format
# ------------------------------------------------------------
M_df <- melt(M)
colnames(M_df) <- c("ReceiverPatch", "SourcePatch", "Weight")

M_df$ReceiverPatch <- as.integer(M_df$ReceiverPatch)
M_df$SourcePatch   <- as.integer(M_df$SourcePatch)

# ------------------------------------------------------------
# Patch name mapping 
# ------------------------------------------------------------
patch_levels <- c(
  "Cape St. Mary's",
  "Witless Bay",
  "Baccalieu Island",
  "Funk Island",
  "Hare Bay",
  "West Coast",
  "Lawn Islands"
)

M_df$ReceiverPatch <- custom_patch_names[M_df$ReceiverPatch]
M_df$SourcePatch   <- custom_patch_names[M_df$SourcePatch]

M_df$ReceiverPatch <- factor(M_df$ReceiverPatch, levels = patch_levels)
M_df$SourcePatch   <- factor(M_df$SourcePatch, levels = patch_levels)

# ------------------------------------------------------------
# Heatmap plot: Fig_3.6
# ------------------------------------------------------------
conn_plot <- ggplot(M_df, aes(x = SourcePatch, y = ReceiverPatch, fill = Weight)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Connectivity") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave("Fig_3.6.png",
       conn_plot, width = 8, height = 6, dpi = 300)