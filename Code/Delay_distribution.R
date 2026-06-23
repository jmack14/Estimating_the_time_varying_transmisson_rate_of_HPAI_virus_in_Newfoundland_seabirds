# ============================================================
# Negative binomial delay distribution for Northern Gannet aerial survey counts
# ============================================================

# ============================================================
# Libraries
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# Remove original aerial survey observations
# ============================================================

points_model_df <- points_model_df %>%
  filter(
    !(CommonName == "Northern Gannet" &
        (
          (patch_name == "Patch1" & TotalObserved == 1136) |
            (patch_name == "Patch4" & TotalObserved == 3158) |
            (patch_name == "Patch3" & TotalObserved == 28)
        )
    )
  )

# ============================================================
# Parameters
# ============================================================

shape <- 2
mean_delay <- 20

start_date <- as.Date("2022-08-01")

# ============================================================
# Backfill 1: Funk Island (3158 deaths, reported Sept 15)
# ============================================================

total_deaths_1 <- 3158
report_date_1 <- as.Date("2022-09-15")

death_dates_1 <- seq(start_date, report_date_1, by = "day")
delay_days_1 <- as.numeric(report_date_1 - death_dates_1)

delay_prob_1 <- dnbinom(delay_days_1, size = shape, mu = mean_delay)
delay_prob_1 <- delay_prob_1 / sum(delay_prob_1)

set.seed(123)

allocated_deaths_1 <- as.vector(
  rmultinom(1, size = total_deaths_1, prob = delay_prob_1)
)

backfill_patch4 <- data.frame(
  DateObserved = death_dates_1,
  CommonName = "Northern Gannet",
  TotalObserved = allocated_deaths_1,
  MortBin = cut(
    allocated_deaths_1,
    breaks = c(1,10,50,100,500,1000,Inf),
    labels = c("1–10","10–50","50–100","100–500","500–1000","1000+"),
    include.lowest = TRUE
  ),
  patch_id = factor("4"),
  patch_name = "Patch4",
  Long = -53.90,
  Lat = 47.60,
  Source = "Survey_Corrected"
) %>%
  filter(TotalObserved > 0)

# ============================================================
# Backfill 2: Cape St. Mary's (1136 deaths, reported Sept 8)
# ============================================================

total_deaths_2 <- 1136
report_date_2 <- as.Date("2022-09-08")

death_dates_2 <- seq(start_date, report_date_2, by = "day")
delay_days_2 <- as.numeric(report_date_2 - death_dates_2)

delay_prob_2 <- dnbinom(delay_days_2, size = shape, mu = mean_delay)
delay_prob_2 <- delay_prob_2 / sum(delay_prob_2)

set.seed(456)

allocated_deaths_2 <- as.vector(
  rmultinom(1, size = total_deaths_2, prob = delay_prob_2)
)

backfill_patch1 <- data.frame(
  DateObserved = death_dates_2,
  CommonName = "Northern Gannet",
  TotalObserved = allocated_deaths_2,
  MortBin = cut(
    allocated_deaths_2,
    breaks = c(1,10,50,100,500,1000,Inf),
    labels = c("1–10","10–50","50–100","100–500","500–1000","1000+"),
    include.lowest = TRUE
  ),
  patch_id = factor("1"),
  patch_name = "Patch1",
  Long = -54.20,
  Lat = 46.82,
  Source = "Survey_Corrected"
) %>%
  filter(TotalObserved > 0)

# ============================================================
# Backfill 3: Baccalieu Island (28 deaths, reported Sept 14)
# ============================================================

total_deaths_3 <- 28
report_date_3 <- as.Date("2022-09-14")

death_dates_3 <- seq(start_date, report_date_3, by = "day")
delay_days_3 <- as.numeric(report_date_3 - death_dates_3)

delay_prob_3 <- dnbinom(delay_days_3, size = shape, mu = mean_delay)
delay_prob_3 <- delay_prob_3 / sum(delay_prob_3)

set.seed(789)

allocated_deaths_3 <- as.vector(
  rmultinom(1, size = total_deaths_3, prob = delay_prob_3)
)

backfill_patch3 <- data.frame(
  DateObserved = death_dates_3,
  CommonName = "Northern Gannet",
  TotalObserved = allocated_deaths_3,
  MortBin = cut(
    allocated_deaths_3,
    breaks = c(1,10,50,100,500,1000,Inf),
    labels = c("1–10","10–50","50–100","100–500","500–1000","1000+"),
    include.lowest = TRUE
  ),
  patch_id = factor("3"),
  patch_name = "Patch3",
  Long = -52.80,
  Lat = 48.10,
  Source = "Survey_Corrected"
) %>%
  filter(TotalObserved > 0)

# ============================================================
# Combine all data
# ============================================================

points_model_df <- bind_rows(
  points_model_df,
  backfill_patch4,
  backfill_patch1,
  backfill_patch3
) %>%
  mutate(
    DateObserved = as.Date(DateObserved),
    patch_id = factor(patch_id),
    patch_name = as.character(patch_name),
    CommonName = as.character(CommonName),
    Source = as.character(Source)
  ) %>%
  select(
    DateObserved,
    CommonName,
    TotalObserved,
    MortBin,
    patch_id,
    patch_name,
    Long,
    Lat,
    Source
  ) %>%
  arrange(DateObserved, patch_name, CommonName, Source)

# ============================================================
# Export
# ============================================================

write.csv(
  points_model_df,
  "NFLD_mortalities_7patch_with_backdistributed_gannets_NB.csv",
  row.names = FALSE
)

# ============================================================
# Plot NB distribution: Fig A2.2
# ============================================================

plot_df <- points_model_df %>%
  filter(
    CommonName == "Northern Gannet",
    patch_name %in% c("Patch1", "Patch3", "Patch4"),
    Source == "Survey_Corrected"
  ) %>%
  group_by(DateObserved, patch_name) %>%
  summarise(DailyDeaths = sum(TotalObserved), .groups = "drop")

ggplot(plot_df, aes(x = DateObserved, y = DailyDeaths, color = patch_name)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Back-distributed Northern Gannet Mortalities (Negative Binomial)",
    x = "Date",
    y = "Daily deaths",
    color = "Patch"
  ) +
  theme_minimal()

# ============================================================
# Plot NB distribution: Fig A2.2
# ============================================================

custom_patch_names <- c(
  Patch1 = "Cape St. Mary's",
  Patch2 = "Witless Bay",
  Patch3 = "Baccalieu Island",
  Patch4 = "Funk Island",
  Patch5 = "Hare Bay",
  Patch6 = "West Coast",
  Patch7 = "Lawn Islands"
)

plot_df <- points_model_df %>%
  filter(
    CommonName == "Northern Gannet",
    patch_name %in% c("Patch1", "Patch3", "Patch4"),
    Source == "Survey_Corrected"
  ) %>%
  group_by(DateObserved, patch_name) %>%
  summarise(DailyDeaths = sum(TotalObserved), .groups = "drop") %>%
  mutate(
    patch_name = factor(
      patch_name,
      levels = names(custom_patch_names),
      labels = custom_patch_names[names(custom_patch_names)]
    )
  )

p_nb <- ggplot(
  plot_df,
  aes(x = DateObserved, y = DailyDeaths, color = patch_name)
) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "Daily deaths",
    color = NULL
  ) +
  guides(
    color = guide_legend(nrow = 1)
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

ggsave(
  "Fig_A2.2.png",
  plot = p_nb,
  width = 7,
  height = 4.5,
  dpi = 600
)
