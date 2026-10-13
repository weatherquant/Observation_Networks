# ========================================== #
#
# NETWORK OPTIMISATION ANALYSIS
# JASON WEST
# BUREAU OF METEOROLOGY, AUSTRALIA
# 12 OCTOBER 2026
#
# ========================================== #

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(geosphere)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# ----------------------------------------------------------
# File paths
# ----------------------------------------------------------

f_sa1_shape <- paste0("SA1_2021_AUST_GDA2020.shp")

f_population <- paste0("SA1Pop.csv")

f_fsoi <- paste0("station_fsoi.txt")

out_dir <- paste0("outputs")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Read SA1 boundaries
# ----------------------------------------------------------

shp_data <- st_read(f_sa1_shape, quiet = TRUE)

# Retain SA1 identifiers as character.
# Do not convert 11-digit geographical identifiers to numeric.
shp_data <- shp_data |>
  mutate(SA1_CODE21 = as.character(SA1_CODE21))

# ----------------------------------------------------------
# Read population data
# ----------------------------------------------------------

pop <- fread(f_population, na.strings = c("", "NA", "N/A", "..", "-", "np"))

# Rename original SA1 field
setnames(pop, old = "SA1", new = "SA1_CODE21")

# Preserve identifier as character
pop[, SA1_CODE21 := as.character(SA1_CODE21)]

# Inspect values that cannot be converted directly
invalid_population_values <- pop[!is.na(Year_2024) &
                                   is.na(suppressWarnings(as.numeric(gsub(
                                     ",", "", trimws(as.character(Year_2024))
                                   )))), unique(Year_2024)]

print(invalid_population_values)

# Clean population values
pop[, Year_2024 := gsub(",", "", trimws(as.character(Year_2024)))]

pop[Year_2024 %in% c("", "NA", "N/A", "..", "-", "np"), Year_2024 := NA_character_]

pop[, Year_2024 := suppressWarnings(as.numeric(Year_2024))]

# Check missing values after cleaning
cat(
  "Population records:",
  nrow(pop),
  "\n",
  "Missing population values:",
  sum(is.na(pop$Year_2024)),
  "\n"
)

# ----------------------------------------------------------
# Join population to SA1 geography
# ----------------------------------------------------------

pop_dens <- shp_data |>
  left_join(pop[, .(SA1_CODE21, Year_2024)], by = "SA1_CODE21")

# Check join coverage
cat(
  "SA1 polygons:",
  nrow(pop_dens),
  "\n",
  "SA1 polygons with population:",
  sum(!is.na(pop_dens$Year_2024)),
  "\n",
  "SA1 polygons without population:",
  sum(is.na(pop_dens$Year_2024)),
  "\n"
)

# Do not automatically replace all missing values with zero until
# their meaning has been checked.

# ----------------------------------------------------------
# Retain populated SA1 regions
# ----------------------------------------------------------

populated_sa1 <- pop_dens |>
  filter(!is.na(Year_2024), Year_2024 > 0)

# Ensure valid geometries before point construction
populated_sa1 <- st_make_valid(populated_sa1)

# point_on_surface is preferable to centroid for irregular polygons
populated_points <- populated_sa1 |>
  st_point_on_surface()

# Transform to longitude/latitude for Haversine distances
populated_points_ll <- st_transform(populated_points, 4326)

pop_coordinates <- st_coordinates(populated_points_ll)

population_points <- populated_points_ll |>
  st_drop_geometry() |>
  mutate(pop_lon = pop_coordinates[, 1], pop_lat = pop_coordinates[, 2]) |>
  select(
    SA1_CODE21,
    SA2_CODE21,
    SA2_NAME21,
    SA3_CODE21,
    SA3_NAME21,
    STE_NAME21,
    AREASQKM21,
    Year_2024,
    pop_lon,
    pop_lat
  )

summary(population_points$Year_2024)

# ----------------------------------------------------------
# Read station FSOI file
# ----------------------------------------------------------

# setwd("/Data")

# f_fsoi <- paste0("fsoi_aus_network_synops.txt")

library(data.table)
library(dplyr)

data_dir <- paste0("Geography")

f_fsoi <- file.path(data_dir, "fsoi_aus_network_synops.txt")

# ----------------------------------------------------------
# Read a whitespace-delimited file
# comment.char removes any initial lines beginning with #
# sep = "" treats one or more spaces/tabs as the delimiter
# ----------------------------------------------------------

fsoi <- read.table(
  file = f_fsoi,
  header = TRUE,
  sep = "",
  comment.char = "#",
  quote = "",
  fill = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE,
  check.names = FALSE
)

fsoi <- as.data.table(fsoi)

print(names(fsoi))
print(head(fsoi))
str(fsoi)

required_fsoi_columns <- c("WMO_ID",
                           "Longitude",
                           "Latitude",
                           "SYNOP_Total_FSOI",
                           "SYNOP_obs",
                           "FSOI_per_obs")

missing_columns <- setdiff(required_fsoi_columns, names(fsoi))

if (length(missing_columns) > 0) {
  stop(paste("Missing columns:", paste(missing_columns, collapse = ", ")))
}

stations <- fsoi |>
  transmute(
    WMO_ID = as.character(WMO_ID),
    lon = as.numeric(Longitude),
    lat = as.numeric(Latitude),
    fsoi_total = as.numeric(SYNOP_Total_FSOI),
    synop_obs = as.numeric(SYNOP_obs),
    fsoi_per_obs_micro = as.numeric(FSOI_per_obs)
  ) |>
  filter(!is.na(WMO_ID), !is.na(lon), !is.na(lat), !is.na(fsoi_total)) |>
  distinct(WMO_ID, .keep_all = TRUE) |>
  mutate(
    # Negative FSOI represents forecast-error reduction
    fsoi_benefit = pmax(-fsoi_total, 0),
    
    # Positive FSOI represents forecast-error increase
    fsoi_degradation = pmax(fsoi_total, 0),
    
    # Beneficial impact per observation
    fsoi_per_obs_benefit =
      pmax(-fsoi_per_obs_micro, 0),
    
    fsoi_class = case_when(
      fsoi_total < 0 ~ "Beneficial",
      fsoi_total > 0 ~ "Forecast degrading",
      TRUE ~ "Neutral"
    )
  )

cat("Rows in raw FSOI file:", nrow(fsoi), "\n")
cat("Stations after cleaning:", nrow(stations), "\n")
cat("Unique WMO identifiers:", n_distinct(stations$WMO_ID), "\n")
cat("Beneficial stations:", sum(stations$fsoi_total < 0), "\n")
cat("Forecast-degrading stations:",
    sum(stations$fsoi_total > 0),
    "\n")

summary(stations$fsoi_total)
summary(stations$synop_obs)
summary(stations$fsoi_per_obs_micro)

duplicate_wmo <- fsoi |>
  count(WMO_ID, name = "records") |>
  filter(records > 1)

print(duplicate_wmo)

# ==========================================================
# PART 2
# Link each FSOI station to SA1 population information
# ==========================================================

# ----------------------------------------------------------
# Basic checks
# ----------------------------------------------------------

stopifnot(nrow(stations) > 0, nrow(population_points) > 0)

stopifnot(all(c("WMO_ID", "lon", "lat") %in% names(stations)), all(
  c("SA1_CODE21", "Year_2024", "pop_lon", "pop_lat") %in% names(population_points)
))

# Population values should already be numeric,
# but enforce this before calculations
population_points <- population_points |>
  mutate(
    Year_2024 = as.numeric(Year_2024),
    pop_lon = as.numeric(pop_lon),
    pop_lat = as.numeric(pop_lat)
  ) |>
  filter(!is.na(Year_2024),
         Year_2024 > 0,!is.na(pop_lon),!is.na(pop_lat))

cat("Populated SA1 points available:",
    nrow(population_points),
    "\n")

cat("Total population represented:",
    format(
      sum(population_points$Year_2024, na.rm = TRUE),
      big.mark = ",",
      scientific = FALSE
    ),
    "\n")

# ----------------------------------------------------------
# Coordinate matrices
# ----------------------------------------------------------

population_coordinates <- as.matrix(population_points |>
                                      select(pop_lon, pop_lat))

population_values <- population_points$Year_2024

# ----------------------------------------------------------
# Function to calculate station population measures
# ----------------------------------------------------------

calculate_station_population <- function(station_lon,
                                         station_lat,
                                         population_coordinates,
                                         population_values,
                                         population_points,
                                         decay_scale_km = 100) {
  station_coordinate <- c(station_lon, station_lat)
  
  distances_km <- geosphere::distHaversine(station_coordinate, population_coordinates) / 1000
  
  nearest_index <- which.min(distances_km)
  
  data.frame(
    dist_to_nearest_pop_km =
      distances_km[nearest_index],
    
    nearest_SA1_CODE21 =
      population_points$SA1_CODE21[nearest_index],
    
    nearest_SA2_NAME21 =
      population_points$SA2_NAME21[nearest_index],
    
    nearest_SA3_NAME21 =
      population_points$SA3_NAME21[nearest_index],
    
    nearest_state =
      population_points$STE_NAME21[nearest_index],
    
    nearest_SA1_population =
      population_values[nearest_index],
    
    population_25km = sum(population_values[distances_km <= 25], na.rm = TRUE),
    
    population_50km = sum(population_values[distances_km <= 50], na.rm = TRUE),
    
    population_100km = sum(population_values[distances_km <= 100], na.rm = TRUE),
    
    population_250km = sum(population_values[distances_km <= 250], na.rm = TRUE),
    
    population_exposure_100km = sum(
      population_values *
        exp(-distances_km /
              decay_scale_km),
      na.rm = TRUE
    )
  )
}

# ----------------------------------------------------------
# Run population calculations for all stations
# ----------------------------------------------------------

population_results <- map_dfr(seq_len(nrow(stations)), function(i) {
  if (i %% 50 == 0 ||
      i == 1 ||
      i == nrow(stations)) {
    cat("Processing station", i, "of", nrow(stations), "\n")
  }
  
  station_population <-
    calculate_station_population(
      station_lon =
        stations$lon[i],
      
      station_lat =
        stations$lat[i],
      
      population_coordinates =
        population_coordinates,
      
      population_values =
        population_values,
      
      population_points =
        population_points,
      
      decay_scale_km = 100
    )
  
  bind_cols(data.frame(WMO_ID =
                         stations$WMO_ID[i]), station_population)
})

# ----------------------------------------------------------
# Join population measures to the FSOI station data
# ----------------------------------------------------------

stations_analysis <- stations |>
  left_join(population_results, by = "WMO_ID")

# ----------------------------------------------------------
# Add observation-count reliability categories
# ----------------------------------------------------------

stations_analysis <- stations_analysis |>
  mutate(
    observation_sample = case_when(
      synop_obs < 1000 ~
        "Very low",
      
      synop_obs < 10000 ~
        "Low",
      
      synop_obs < 50000 ~
        "Moderate",
      
      TRUE ~
        "High"
    ),
    
    low_sample_flag =
      synop_obs < 10000,
    
    fsoi_review_class = case_when(
      fsoi_total < 0 &
        synop_obs >= 10000 ~
        "Beneficial",
      
      fsoi_total < 0 &
        synop_obs < 10000 ~
        "Beneficial but low sample",
      
      fsoi_total > 0 &
        synop_obs >= 10000 ~
        "Potentially forecast degrading",
      
      fsoi_total > 0 &
        synop_obs < 10000 ~
        "Positive FSOI but low sample",
      
      TRUE ~
        "Neutral"
    )
  )

# ----------------------------------------------------------
# Check the combined results
# ----------------------------------------------------------

cat("Stations in combined dataset:", nrow(stations_analysis), "\n")

cat("Stations missing population distance:", sum(is.na(
  stations_analysis$dist_to_nearest_pop_km
)), "\n")

summary(stations_analysis$dist_to_nearest_pop_km)

summary(stations_analysis$population_100km)

summary(stations_analysis$population_exposure_100km)

table(stations_analysis$observation_sample, useNA = "ifany")

table(stations_analysis$fsoi_review_class, useNA = "ifany")

# ----------------------------------------------------------
# Inspect the most geographically remote stations
# ----------------------------------------------------------

most_remote_stations <-
  stations_analysis |>
  arrange(desc(dist_to_nearest_pop_km)) |>
  select(
    WMO_ID,
    lon,
    lat,
    fsoi_total,
    fsoi_benefit,
    synop_obs,
    fsoi_per_obs_micro,
    dist_to_nearest_pop_km,
    nearest_SA2_NAME21,
    nearest_SA3_NAME21,
    nearest_state,
    nearest_SA1_population,
    population_100km,
    population_exposure_100km
  ) |>
  slice_head(n = 30)

print(most_remote_stations)

# ----------------------------------------------------------
# Inspect stations with the greatest beneficial total FSOI
# ----------------------------------------------------------

highest_fsoi_benefit <-
  stations_analysis |>
  arrange(desc(fsoi_benefit)) |>
  select(
    WMO_ID,
    lon,
    lat,
    fsoi_total,
    fsoi_benefit,
    synop_obs,
    fsoi_per_obs_micro,
    observation_sample,
    dist_to_nearest_pop_km,
    population_100km,
    population_exposure_100km
  ) |>
  slice_head(n = 30)

print(highest_fsoi_benefit)

# ----------------------------------------------------------
# Inspect stations with positive total FSOI
# ----------------------------------------------------------

forecast_degrading_stations <-
  stations_analysis |>
  filter(fsoi_total > 0) |>
  arrange(desc(fsoi_total)) |>
  select(
    WMO_ID,
    lon,
    lat,
    fsoi_total,
    fsoi_degradation,
    synop_obs,
    fsoi_per_obs_micro,
    observation_sample,
    dist_to_nearest_pop_km,
    population_100km,
    population_exposure_100km
  )

print(forecast_degrading_stations)

# ----------------------------------------------------------
# Examine associations before constructing a cost model
# ----------------------------------------------------------

association_data <-
  stations_analysis |>
  select(
    fsoi_total,
    fsoi_benefit,
    fsoi_per_obs_micro,
    synop_obs,
    dist_to_nearest_pop_km,
    population_50km,
    population_100km,
    population_250km,
    population_exposure_100km
  )

spearman_correlations <- cor(association_data, use = "pairwise.complete.obs", method = "spearman")

print(round(spearman_correlations, 3))

# ----------------------------------------------------------
# Export the combined station dataset
# ----------------------------------------------------------

fwrite(stations_analysis,
       file.path(out_dir, "fsoi_station_population_dataset.csv"))

fwrite(most_remote_stations,
       file.path(out_dir, "most_remote_fsoi_stations.csv"))

fwrite(highest_fsoi_benefit,
       file.path(out_dir, "highest_fsoi_benefit_stations.csv"))

fwrite(
  forecast_degrading_stations,
  file.path(out_dir, "positive_fsoi_stations_for_review.csv")
)

fwrite(
  as.data.table(spearman_correlations, keep.rownames = "variable"),
  file.path(out_dir, "fsoi_population_spearman_correlations.csv")
)

cat("Combined FSOI-population outputs written to:",
    normalizePath(out_dir),
    "\n")

############################################
# Check outliers

stations_analysis |>
  filter(dist_to_nearest_pop_km ==
           max(dist_to_nearest_pop_km, na.rm = TRUE)) |>
  select(
    WMO_ID,
    lon,
    lat,
    fsoi_total,
    synop_obs,
    dist_to_nearest_pop_km,
    nearest_SA1_CODE21,
    nearest_SA2_NAME21,
    nearest_state
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Check stations >1000km from a populated SA1
stations_analysis |>
  filter(dist_to_nearest_pop_km > 1000) |>
  arrange(desc(dist_to_nearest_pop_km)) |>
  select(
    WMO_ID,
    lon,
    lat,
    fsoi_total,
    synop_obs,
    dist_to_nearest_pop_km,
    nearest_SA2_NAME21,
    nearest_state
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ==========================================================
# INITIAL MAINTENANCE COST INDEX
# ==========================================================

# These are scenario parameters, not actual dollar costs
fixed_site_cost <- 1.0
distance_cost_rate <- 0.01
remote_threshold_km <- 250
remote_loading <- 1.0

stations_analysis <- stations_analysis |>
  mutate(
    maintenance_cost =
      fixed_site_cost +
      distance_cost_rate * dist_to_nearest_pop_km +
      if_else(
        dist_to_nearest_pop_km > remote_threshold_km,
        remote_loading,
        0
      )
  )

stations_analysis <- stations_analysis |>
  mutate(
    analysis_region = case_when(
      lon >= 112 & lon <= 154 &
        lat >= -44 & lat <= -10 ~
        "Mainland Australia and Tasmania",
      
      lat < -45 ~
        "Antarctic",
      
      TRUE ~
        "Offshore or external territory"
    )
  )

table(stations_analysis$analysis_region)

analysis_stations <- stations_analysis |>
  filter(analysis_region ==
           "Mainland Australia and Tasmania")

#  Normalise benefit and population measures
norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  
  if (!all(is.finite(r)) || diff(r) == 0) {
    return(rep(0, length(x)))
  }
  
  (x - r[1]) / diff(r)
}

analysis_stations <- analysis_stations |>
  mutate(
    forecast_score =
      norm01(fsoi_benefit),
    
    population_score =
      norm01(log1p(population_exposure_100km)),
    
    cost_score =
      norm01(maintenance_cost),
    
    composite_benefit =
      0.8 * forecast_score +
      0.2 * population_score,
    
    benefit_cost_ratio =
      composite_benefit /
      pmax(maintenance_cost, 1e-9)
  )

# First budget-selection algorithm
select_under_budget <- function(data,
                                budget_fraction,
                                weight_fsoi = 0.8,
                                weight_population = 0.2) {
  dat <- data |>
    mutate(
      decision_benefit =
        weight_fsoi * forecast_score +
        weight_population * population_score,
      
      decision_ratio =
        decision_benefit /
        pmax(maintenance_cost, 1e-9)
    ) |>
    arrange(desc(decision_ratio), desc(decision_benefit), WMO_ID)
  
  full_cost <- sum(dat$maintenance_cost, na.rm = TRUE)
  
  budget <- budget_fraction * full_cost
  
  dat <- dat |>
    mutate(
      cumulative_cost =
        cumsum(maintenance_cost),
      
      retained =
        cumulative_cost <= budget,
      
      decision =
        if_else(retained, "Retain", "Review for rationalisation")
    )
  
  list(
    budget_fraction = budget_fraction,
    budget = budget,
    full_cost = full_cost,
    station_results = dat,
    
    portfolio_summary = tibble(
      budget_fraction =
        budget_fraction,
      
      budget =
        budget,
      
      retained_cost =
        sum(dat$maintenance_cost[dat$retained]),
      
      stations_retained =
        sum(dat$retained),
      
      stations_reviewed =
        sum(!dat$retained),
      
      retained_fsoi_benefit =
        sum(dat$fsoi_benefit[dat$retained]),
      
      full_fsoi_benefit =
        sum(dat$fsoi_benefit),
      
      retained_population_exposure =
        sum(dat$population_exposure_100km[dat$retained])
    )
  )
}

# Run several budget scenarios
budget_fractions <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.70)

budget_results <- lapply(budget_fractions, function(b) {
  select_under_budget(
    data = analysis_stations,
    budget_fraction = b,
    weight_fsoi = 0.8,
    weight_population = 0.2
  )
})

head(budget_results, 1)

budget_summary <- map_dfr(budget_results, function(x) {
  x$portfolio_summary |>
    mutate(full_cost = x$full_cost)
}) |>
  mutate(
    fsoi_retained_percent =
      100 *
      retained_fsoi_benefit /
      full_fsoi_benefit,
    
    cost_retained_percent =
      100 *
      retained_cost /
      full_cost,
    
    population_exposure_retained_percent =
      100 *
      retained_population_exposure /
      budget_results[[1]]$portfolio_summary$retained_population_exposure
  )

print(as.data.frame(budget_summary), row.names = FALSE)

budget_station_decisions <- map_dfr(
  budget_results,
  function(x) {
    
    x$station_results |>
      transmute(
        budget_fraction =
          x$budget_fraction,
        
        WMO_ID,
        lon,
        lat,
        fsoi_total,
        fsoi_benefit,
        synop_obs,
        population_exposure_100km,
        dist_to_nearest_pop_km,
        maintenance_cost,
        decision_benefit,
        decision_ratio,
        retained,
        decision
      )
  }
)

fwrite(
  budget_summary,
  file.path(
    out_dir,
    "budget_scenario_summary.csv"
  )
)

fwrite(
  budget_station_decisions,
  file.path(
    out_dir,
    "station_decisions_by_budget.csv"
  )
)

# Figure 1: Cost–benefit retention curve
p_frontier <- ggplot(
  budget_summary,
  aes(
    x = cost_retained_percent,
    y = fsoi_retained_percent
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_line(
    colour = "#2166AC",
    linewidth = 1.1
  ) +
  geom_point(
    aes(
      size = stations_retained
    ),
    colour = "#2166AC"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title =
      "Forecast impact retained under alternative network budgets",
    x =
      "Network maintenance cost retained (%)",
    y =
      "Beneficial FSOI retained (%)",
    size =
      "Stations retained"
  )

ggsave(
  file.path(
    out_dir,
    "forecast_impact_budget_frontier.png"
  ),
  plot = p_frontier,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

# Figure 2: Station typology plot
p_typology <- ggplot(
  analysis_stations,
  aes(
    x = maintenance_cost,
    y = fsoi_benefit,
    size = population_exposure_100km,
    colour = fsoi_review_class
  )
) +
  geom_point(
    alpha = 0.7
  ) +
  scale_size_continuous(
    trans = "log1p",
    breaks = c(
      1e4,
      1e5,
      5e5,
      1e6,
      2e6
    ),
    labels = c(
      "10k",
      "100k",
      "500k",
      "1M",
      "2M"
    ),
    name = "Population exposure"
  ) +
  scale_colour_manual(
    values = c(
      "Beneficial" =
        "#2166AC",
      
      "Beneficial but low sample" =
        "#67A9CF",
      
      "Potentially forecast degrading" =
        "#B2182B",
      
      "Positive FSOI but low sample" =
        "#EF8A62",
      
      "Neutral" =
        "grey50"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title =
      "Forecast impact, population exposure and maintenance cost",
    x =
      "Maintenance cost index",
    y =
      "Beneficial total FSOI",
    colour =
      "Station classification"
  )

ggsave(
  file.path(
    out_dir,
    "station_value_typology.png"
  ),
  plot = p_typology,
  width = 8,
  height = 5.5,
  dpi = 300
)

# Figure 3: Decommissioning/review map at a selected budget
# For example, an 80% cost scenario:
selected_budget <- 0.80

map_data <- budget_station_decisions |>
  filter(
    budget_fraction ==
      selected_budget
  ) |>
  mutate(
    decision = factor(
      decision,
      levels = c(
        "Retain",
        "Review for rationalisation"
      )
    )
  )

p_map <- ggplot() +
  geom_sf(
    data = shp_data,
    fill = "grey97",
    colour = "grey75",
    linewidth = 0.15
  ) +
  geom_point(
    data = map_data,
    aes(
      x = lon,
      y = lat,
      colour = decision,
      size = fsoi_benefit
    ),
    alpha = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "Retain" = "#2166AC",
      "Review for rationalisation" =
        "#D6604D"
    )
  ) +
  coord_sf(
    xlim = c(112, 154),
    ylim = c(-44, -10),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title =
      "Indicative network decisions under an 80% cost budget",
    colour =
      "Portfolio decision",
    size =
      "Beneficial FSOI"
  )

ggsave(
  file.path(
    out_dir,
    "station_decision_map_80_percent_budget.png"
  ),
  plot = p_map,
  width = 8.5,
  height = 7,
  dpi = 300
)

# Figure 4: Regional station-to-population relationships
region_name <- "Greater Darwin"

regional_population <- populated_points_ll |>
  filter(
    SA3_NAME21 == region_name
  )

regional_stations <- analysis_stations |>
  filter(
    lon >= 129,
    lon <= 132,
    lat >= -14,
    lat <= -11
  )

regional_links <- regional_stations |>
  left_join(
    population_points |>
      select(
        SA1_CODE21,
        pop_lon,
        pop_lat
      ),
    by = c(
      "nearest_SA1_CODE21" =
        "SA1_CODE21"
    )
  )

# Build line geometry
link_geometry <- lapply(
  seq_len(nrow(regional_links)),
  function(i) {
    
    st_linestring(
      matrix(
        c(
          regional_links$lon[i],
          regional_links$lat[i],
          regional_links$pop_lon[i],
          regional_links$pop_lat[i]
        ),
        ncol = 2,
        byrow = TRUE
      )
    )
  }
)

regional_links_sf <- st_sf(
  regional_links,
  geometry = st_sfc(
    link_geometry,
    crs = 4326
  )
)

pmap2 <- ggplot() +
  geom_sf(
    data = regional_population,
    aes(
      size = Year_2024
    ),
    colour = "grey40",
    alpha = 0.5
  ) +
  geom_sf(
    data = regional_links_sf,
    colour = "grey60",
    linewidth = 0.4
  ) +
  geom_point(
    data = regional_stations,
    aes(
      x = lon,
      y = lat,
      colour = fsoi_benefit
    ),
    size = 3
  ) +
  scale_colour_viridis_c() +
  theme_minimal()

ggsave(
  file.path(
    out_dir,
    "link_build.png"
  ),
  plot = pmap2,
  width = 8.5,
  height = 7,
  dpi = 300
)

# ==========================================================
# IMPROVED FIGURE 4
# Greater Darwin population, station impact and service radii
# ==========================================================

library(sf)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------
# 1. Regional geographic extent
# ----------------------------------------------------------

darwin_extent <- st_bbox(c(
  xmin = 129.0,
  ymin = -14.0,
  xmax = 132.0,
  ymax = -11.0
), crs = st_crs(4326))

# Use the populated SA1 polygons rather than representative points
darwin_population <- populated_sa1 |>
  st_transform(4326) |>
  st_crop(darwin_extent) |>
  mutate(population_density =
           if_else(AREASQKM21 > 0, Year_2024 / AREASQKM21, NA_real_))

# ----------------------------------------------------------
# 2. Select stations in the regional extent
# ----------------------------------------------------------

darwin_stations <- analysis_stations |>
  filter(lon >= 129.0, lon <= 132.0, lat >= -14.0, lat <= -11.0)

darwin_stations_sf <- st_as_sf(
  darwin_stations,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# ----------------------------------------------------------
# 3. Construct service-radius buffers
#    Use Australian Albers so distances are in metres
# ----------------------------------------------------------

darwin_stations_projected <- st_transform(darwin_stations_sf, 3577)

station_buffer_50km <- st_buffer(darwin_stations_projected, dist = 50000) |>
  st_transform(4326)

station_buffer_100km <- st_buffer(darwin_stations_projected, dist = 100000) |>
  st_transform(4326)

# ----------------------------------------------------------
# 4. Station classification for plotting
# ----------------------------------------------------------

darwin_stations_sf <- darwin_stations_sf |>
  mutate(
    station_status = case_when(
      fsoi_total < 0 &
        synop_obs >= 10000 ~
        "Beneficial",
      
      fsoi_total < 0 &
        synop_obs < 10000 ~
        "Beneficial, low sample",
      
      fsoi_total >= 0 ~
        "Positive FSOI review",
      
      TRUE ~
        "Other"
    )
  )

# ----------------------------------------------------------
# 5. Improved regional figure
# ----------------------------------------------------------

pmap2 <- ggplot() +
  
  # SA1 population-density background
  geom_sf(
    data = darwin_population,
    aes(fill = population_density),
    colour = NA,
    alpha = 0.85
  ) +
  
  # 100-km station radii
  geom_sf(
    data = station_buffer_100km,
    fill = NA,
    colour = "#636363",
    linewidth = 0.35,
    linetype = "dashed",
    alpha = 0.45
  ) +
  
  # 50-km station radii
  geom_sf(
    data = station_buffer_50km,
    fill = NA,
    colour = "#252525",
    linewidth = 0.45,
    linetype = "solid",
    alpha = 0.55
  ) +
  
  # Station locations
  geom_sf(
    data = darwin_stations_sf,
    aes(
      colour = fsoi_benefit,
      size = synop_obs,
      shape = station_status
    ),
    alpha = 0.9,
    stroke = 0.7
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    trans = "log1p",
    name = paste0("Population density\n", "(people per km²)"),
    na.value = "grey95"
  ) +
  
  scale_colour_viridis_c(
    option = "D",
    direction = -1,
    name = paste0("Beneficial total FSOI\n", "(-FSOI, J kg\u207B\u00B9)")
  ) +
  
  scale_size_continuous(trans = "log1p",
                        range = c(2, 7),
                        name = "Assimilated\nobservations") +
  
  scale_shape_manual(
    values = c(
      "Beneficial" = 21,
      "Beneficial, low sample" = 24,
      "Positive FSOI review" = 4,
      "Other" = 3
    ),
    name = "Station status"
  ) +
  
  coord_sf(
    xlim = c(129.0, 132.0),
    ylim = c(-14.0, -11.0),
    expand = FALSE
  ) +
  
  guides(
    fill = guide_colourbar(order = 1),
    colour = guide_colourbar(order = 2),
    size = guide_legend(order = 3),
    shape = guide_legend(order = 4)
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    panel.grid.major =
      element_line(colour = "grey85", linewidth = 0.25),
    
    legend.position = "right",
    
    plot.title =
      element_text(face = "bold"),
    
    plot.subtitle =
      element_text(colour = "grey30")
  ) +
  
  labs(
    title =
      "Surface-observation value and population exposure in the Greater Darwin region",
    
    subtitle =
      "Solid and dashed boundaries denote 50-km and 100-km station radii",
    
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  filename = file.path(out_dir, "figure_4_greater_darwin_station_population.png"),
  plot = pmap2,
  width = 9.5,
  height = 7.2,
  dpi = 300
)

pmap2

#####################################
### Revised FIGURE 4
dt <- as.data.table(stations)
dt <- dt[fsoi_benefit > 0]
dt <- dt[order(-fsoi_benefit)]
dt[, station_rank := .I]
dt[, cum_stations := station_rank / .N]
dt[, cum_fsoi := cumsum(fsoi_benefit) / sum(fsoi_benefit)]

library(ggplot2)
library(scales)
library(viridis)

df <- stations_analysis[
  fsoi_benefit > 0 &
    maintenance_cost > 0 &
    population_exposure_100km > 0
]

df[, pop_exposure_m := population_exposure_100km / 1e6]
df[, pop_group :=
     cut(
       population_exposure_100km,
       breaks = quantile(
         population_exposure_100km,
         probs = c(0,0.25,0.5,0.75,1),
         na.rm = TRUE
       ),
       include.lowest = TRUE,
       labels = c(
         "Lowest 25%",
         "25-50%",
         "50-75%",
         "Highest 25%"
       )
     )]

# 75th percentile thresholds
cost_cut <- quantile(df$maintenance_cost, 0.75, na.rm = TRUE)
fsoi_cut <- quantile(df$fsoi_benefit, 0.75, na.rm = TRUE)

fig4 <- ggplot(
  df,
  aes(
    x = maintenance_cost,
    y = fsoi_benefit
  )
) +
  
  geom_point(
    aes(colour = pop_group),
    alpha = 0.75,
    size = 2.3
  ) +
  
  scale_colour_brewer(
    palette = "RdYlBu",
    direction = -1,
    name = "Population\nexposure"
  ) +

  scale_x_log10(
    labels = label_number()
  ) +
  
  scale_y_log10(
    labels = label_scientific(digits = 2)
  ) +

  geom_vline(
    xintercept = cost_cut,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = fsoi_cut,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  annotate(
    "label",
    x = 1.15,
    y = max(df$fsoi_benefit) * 0.75,
    label = "Priority\nretention",
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "label",
    x = max(df$maintenance_cost) * 0.75,
    y = min(df$fsoi_benefit[df$fsoi_benefit > 0]) * 15,
    label = "Review\ncandidates",
    fontface = "bold",
    size = 4
  ) +
  
  labs(
    title =
      "Forecast value versus cost",
    
    subtitle =
      "Colours show population exposure. Dashed lines show 75th-percentile thresholds.",
    
    x =
      "Maintenance cost index (log scale)",
    
    y =
      "Beneficial FSOI (log scale)"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

fig4

ggsave(
  filename = file.path(out_dir, "fig4.png"),
  plot = fig4,
  width = 9.5,
  height = 7.2,
  dpi = 300
)

# Robust classification across sensitivities
# A station should only be described as a potential rationalisation candidate if it is consistently excluded across:
# several budgets;
# several FSOI/population weights;
# linear and nonlinear cost functions;
# alternative population radii;
# inclusion and exclusion of low-sample stations;
# total versus per-observation FSOI;
# different remoteness loadings.

# 10. Generalisable relationships
# The paper should test relationships rather than only report a selected portfolio.
# Relationship 1: Remoteness and cost
# Ci=f(di)C_i=f(d_i)Ci.=f(di.)
# Test linear, quadratic and piecewise-linear specifications.
# Relationship 2: Population exposure and network density
# Test whether stations serving larger population exposure are systematically closer together and therefore more spatially redundant.
# Relationship 3: FSOI benefit and remoteness
# Test whether remote observations deliver disproportionately high forecast impact:
lm(
  log1p(fsoi_benefit) ~
    log1p(dist_to_nearest_pop_km) +
    log1p(synop_obs),
  data = analysis_stations
)

# Relationship 4: FSOI and observation volume
lm(
  fsoi_benefit ~
    log1p(synop_obs),
  data = analysis_stations
)

lm(
  log1p(fsoi_benefit) ~
    log1p(synop_obs) +
    log1p(dist_to_nearest_pop_km) +
    log1p(population_exposure_100km),
  data = analysis_stations
)

## Full regression outputs
# Refit and retain the three models
model_1 <- lm(
  log1p(fsoi_benefit) ~
    log1p(dist_to_nearest_pop_km) +
    log1p(synop_obs),
  data = analysis_stations
)

model_2 <- lm(
  fsoi_benefit ~
    log1p(synop_obs),
  data = analysis_stations
)

model_3 <- lm(
  log1p(fsoi_benefit) ~
    log1p(synop_obs) +
    log1p(dist_to_nearest_pop_km) +
    log1p(population_exposure_100km),
  data = analysis_stations
)

# Full coefficient tables, fit statistics and conventional diagnostics
summary(model_1)
summary(model_2)
summary(model_3)

confint(model_1)
confint(model_2)
confint(model_3)

AIC(model_1, model_2, model_3)
BIC(model_1, model_2, model_3)

# Now estimate robust standard errors (HC3)
robust_ci <- function(model) {
  b <- coef(model)
  se <- sqrt(diag(sandwich::vcovHC(model, type = "HC3")))
  
  data.frame(
    estimate = b,
    robust_se = se,
    lower_95 = b - qt(0.975, df.residual(model)) * se,
    upper_95 = b + qt(0.975, df.residual(model)) * se
  )
}

robust_ci(model_1)
robust_ci(model_2)
robust_ci(model_3)

# ==========================================================
# Budget portfolio results table
# ==========================================================

budget_summary <- purrr::map_dfr(budget_results, function(x) {
  # Station-level results for this budget
  dat <- x$station_results
  
  tibble(
    budget_fraction =
      x$budget_fraction,
    
    budget_limit =
      x$budget,
    
    full_network_cost =
      x$full_cost,
    
    retained_cost =
      sum(dat$maintenance_cost[dat$retained], na.rm = TRUE),
    
    unused_budget =
      x$budget -
      sum(dat$maintenance_cost[dat$retained], na.rm = TRUE),
    
    stations_retained =
      sum(dat$retained, na.rm = TRUE),
    
    stations_reviewed =
      sum(!dat$retained, na.rm = TRUE),
    
    retained_fsoi_benefit =
      sum(dat$fsoi_benefit[dat$retained], na.rm = TRUE),
    
    full_fsoi_benefit =
      sum(dat$fsoi_benefit, na.rm = TRUE),
    
    retained_population_exposure =
      sum(dat$population_exposure_100km[dat$retained], na.rm = TRUE),
    
    full_population_exposure =
      sum(dat$population_exposure_100km, na.rm = TRUE)
  )
}) |>
  mutate(
    budget_percent =
      100 * budget_fraction,
    
    cost_retained_percent =
      100 *
      retained_cost /
      full_network_cost,
    
    fsoi_retained_percent =
      100 *
      retained_fsoi_benefit /
      full_fsoi_benefit,
    
    population_exposure_retained_percent =
      100 *
      retained_population_exposure /
      full_population_exposure,
    
    unused_budget_percent =
      100 *
      unused_budget /
      full_network_cost,
    
    stations_retained_percent =
      100 *
      stations_retained /
      (stations_retained +
         stations_reviewed),
    
    fsoi_advantage_over_cost =
      fsoi_retained_percent -
      cost_retained_percent
  ) |>
  arrange(desc(budget_fraction))

print(
  as.data.frame(
    budget_summary
  ),
  row.names = FALSE
)

table_budget <- budget_summary |>
  transmute(
    `Budget limit (%)` =
      round(budget_percent, 0),
    
    `Cost index retained (%)` =
      round(cost_retained_percent, 1),
    
    `Unused budget (%)` =
      round(unused_budget_percent, 2),
    
    `Stations retained` =
      stations_retained,
    
    `Stations reviewed` =
      stations_reviewed,
    
    `Stations retained (%)` =
      round(stations_retained_percent, 1),
    
    `Beneficial FSOI retained (%)` =
      round(fsoi_retained_percent, 1),
    
    `Aggregate exposure score retained (%)` =
      round(population_exposure_retained_percent, 1),
    
    `FSOI retained minus cost retained (percentage points)` =
      round(fsoi_advantage_over_cost, 1)
  )

print(as.data.frame(table_budget), row.names = FALSE)

fwrite(
  as.data.table(
    table_budget
  ),
  file.path(
    out_dir,
    "budget_portfolio_results_table.csv"
  )
)

# Influence and residual diagnostics
par(mfrow = c(2, 2))
plot(model_1)

par(mfrow = c(2, 2))
plot(model_2)

par(mfrow = c(2, 2))
plot(model_3)

# Additional diagnostic statistics
install.packages(
  c("lmtest", "sandwich", "car"),
  dependencies = TRUE
)

library(lmtest)
library(sandwich)
library(car)

# Heteroskedasticity tests
bptest(model_1)
bptest(model_2)
bptest(model_3)

# Heteroskedasticity-robust coefficient inference
coeftest(
  model_1,
  vcov. = vcovHC(model_1, type = "HC3")
)

coeftest(
  model_2,
  vcov. = vcovHC(model_2, type = "HC3")
)

coeftest(
  model_3,
  vcov. = vcovHC(model_3, type = "HC3")
)

# Multicollinearity
vif(model_1)
vif(model_3)

# Influence thresholds
influence_summary <- function(model) {
  n <- nobs(model)
  p <- length(coef(model))
  
  data.frame(
    observations = n,
    parameters = p,
    residual_df = df.residual(model),
    influential_cooks =
      sum(cooks.distance(model) > 4 / n),
    high_leverage =
      sum(hatvalues(model) > 2 * p / n),
    large_studentised_residual =
      sum(abs(rstudent(model)) > 2)
  )
}

influence_summary(model_1)
influence_summary(model_2)
influence_summary(model_3)

# SPATIAL DIAGNOSTICS
install.packages(c("spdep", "sf"), dependencies = TRUE)

library(spdep)
library(sf)

station_sf <- st_as_sf(
  analysis_stations,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# Australian Albers projected coordinates
station_sf_projected <- st_transform(
  station_sf,
  3577
)

station_coordinates <- st_coordinates(
  station_sf_projected
)

# Eight-nearest-neighbour spatial structure
knn_8 <- knearneigh(
  station_coordinates,
  k = 8
)

neighbour_list <- knn2nb(knn_8)

spatial_weights <- nb2listw(
  neighbour_list,
  style = "W",
  zero.policy = TRUE
)

moran.test(
  residuals(model_1),
  spatial_weights,
  zero.policy = TRUE
)

moran.test(
  residuals(model_2),
  spatial_weights,
  zero.policy = TRUE
)

moran.test(
  residuals(model_3),
  spatial_weights,
  zero.policy = TRUE
)

