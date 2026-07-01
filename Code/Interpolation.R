# ============================================================
# Data interpolation for SEIRD model 
# ============================================================

# ============================================================
# Libraries 
# ============================================================

library(tidyverse); theme_set(theme_bw())
library(lme4)
library(shellpipes)

# ============================================================
# Load data
# ============================================================

rpcall("dat.Rout dat.R NFLD_mortalities_7patch_with_backdistributed_gannets_NB.csv")

rawdat <- csvRead()

dat <- rawdat %>%
  transmute(
    date = as.Date(DateObserved),
    bird = gsub(" ", "_", CommonName),
    count = TotalObserved,
    patch_name,
    Source
  )

# ============================================================
# Remove outliers
# ============================================================

dat_clean <- dat %>%
  filter(
    !(patch_name == "Patch1" & bird == "Northern_Gannet" & date == as.Date("2022-06-05")),
    !(patch_name == "Patch5" & bird == "Common_Murre" & date == as.Date("2022-05-11")),
    !(patch_name == "Patch6" & bird == "Common_Murre" & date %in% as.Date(c("2022-04-23", "2022-05-10"))),
    !(patch_name == "Patch7" & bird == "Northern_Gannet" & date == as.Date("2022-05-22")),
    
    # NEW: Patch 2 Gannets
    !(patch_name == "Patch2" & bird == "Northern_Gannet" & date %in% as.Date(c("2022-06-02", "2022-06-06", "2022-06-14"))),
    
    # NEW: Patch 3 Gannets
    !(patch_name == "Patch3" & bird == "Northern_Gannet" & date %in% as.Date(c("2022-05-22", "2022-05-29")))
  )

# ============================================================
# Plot deaths
# ============================================================

gg <- ggplot(dat_clean, aes(date, count)) +
  geom_point(aes(color = Source)) +
  facet_wrap(~patch_name, scales = "free")

print(gg)

# ===========================================================
# Aggregation
# ============================================================

dat2 <- dat_clean %>%
  group_by(date, patch_name, bird) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(patch_name, bird, date) %>%
  group_by(patch_name, bird) %>%
  mutate(
    cumtotal = cumsum(total),
    mindate = min(date),
    maxdate = max(date)
  ) %>%
  ungroup()

# ============================================================
# Interpolation
# ============================================================

dat3 <- dat2 %>%
  group_by(patch_name, bird) %>%
  complete(date = seq.Date(min(date), max(date), by = "day")) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    cumtotal = zoo::na.approx(cumtotal, na.rm = FALSE),
    cumtotal_fill = cumtotal
  ) %>%
  ungroup()

# ============================================================
# Plot interpolated deaths
# ============================================================

gg2 <- ggplot(dat3, aes(date, cumtotal, color = bird)) +
  geom_point() +
  geom_line(aes(y = cumtotal_fill)) +
  facet_wrap(~patch_name, ncol = 1)

print(gg2)

# ============================================================
# Save output
# ============================================================

rdsSave(dat3)