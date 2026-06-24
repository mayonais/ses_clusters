library(tidycensus)
library(tidyverse)

library(dplyr)
library(cluster)
library(factoextra)
library(tidyr)

library("PReMiuM")
library(sf)
library(tigris)

year <- 2018
survey <- "acs5"
state <- "CA"
county <- "Los Angeles"

#census_api_key("e3028f4138dfedd065a29bbe27a70a9", install = TRUE)

# --------------------------------------------------------------------

# vars <- load_variables(year, survey, cache = TRUE)
# View(vars)

acs <- get_acs(
  geography = "zcta",
  survey = survey,
  year = year,
  state = state,
  variables = c(
    total_pop = "B01003_001",
    
    race_total = "B03002_001",
    white_nh = "B03002_003",
    black_nh = "B03002_004",
    asian_nh = "B03002_006",
    hispanic = "B03002_012",
    
    median_age = "B01002_001",
    
    median_income = "B19013_001",
    
    poverty = "B17001_002",
    poverty_total = "B17001_001",
    
    unemployment = "B23025_005",
    labor_force = "B23025_003",
    
    total_education = "B15003_001",
    edu_001 = "B15003_002", edu_002 = "B15003_003", edu_003 = "B15003_004",
    edu_004 = "B15003_005", edu_005 = "B15003_006", edu_006 = "B15003_007",
    edu_007 = "B15003_008", edu_008 = "B15003_009", edu_009 = "B15003_010",
    edu_010 = "B15003_011", edu_011 = "B15003_012", edu_012 = "B15003_013",
    edu_013 = "B15003_014", edu_014 = "B15003_015", edu_015 = "B15003_016",
    
    total_insurance = "B27010_001",
    ins_001 = "B27010_017", ins_002 = "B27010_033",
    ins_003 = "B27010_050", ins_004 = "B27010_066",
    
    total_language = "B16004_001",
    eng_001 = "B16004_007",  eng_002 = "B16004_008",
    eng_003 = "B16004_012",  eng_004 = "B16004_013",
    eng_005 = "B16004_017",  eng_006 = "B16004_018",
    eng_007 = "B16004_022",  eng_008 = "B16004_023",
    eng_009 = "B16004_029",  eng_010 = "B16004_030",
    eng_011 = "B16004_034",  eng_012 = "B16004_035",
    eng_013 = "B16004_039",  eng_014 = "B16004_040",
    eng_015 = "B16004_044",  eng_016 = "B16004_045",
    eng_017 = "B16004_051",  eng_018 = "B16004_052",
    eng_019 = "B16004_056",  eng_020 = "B16004_057",
    eng_021 = "B16004_061",  eng_022 = "B16004_062",
    eng_023 = "B16004_066",  eng_024 = "B16004_067",
    
    total_disability = "B18101_001",
    dis_001 = "B18101_004", dis_002 = "B18101_007", dis_003 = "B18101_010",
    dis_004 = "B18101_013", dis_005 = "B18101_016", dis_006 = "B18101_019",
    dis_007 = "B18101_023", dis_008 = "B18101_026", dis_009 = "B18101_029",
    dis_010 = "B18101_032", dis_011 = "B18101_035", dis_012 = "B18101_038",
    
    owner_occ = "B25106_002",
    owner_30_1 = "B25106_006", owner_30_2 = "B25106_010",
    owner_30_3 = "B25106_014", owner_30_4 = "B25106_018",
    owner_30_5 = "B25106_022",
    
    renter_occ = "B25106_024",
    renter_30_1 = "B25106_028", renter_30_2 = "B25106_032",
    renter_30_3 = "B25106_036", renter_30_4 = "B25106_040",
    renter_30_5 = "B25106_044",
    
    housing_crowding_total = "B25014_001",
    crowd_001 = "B25014_007", crowd_002 = "B25014_013",
    
    outdoor_total = "C24050_001",
    agri_forestry = "C24050_002",
    construction_ind = "C24050_003"
  ),
  output = "wide",
  geometry = TRUE
)

acs_clean <- acs %>%
  filter(total_popE > 0) %>%
  mutate(
    total_popE,
    
    pct_white = white_nhE / race_totalE,
    pct_black = black_nhE / race_totalE,
    pct_asian = asian_nhE / race_totalE,
    pct_hispanic = hispanicE / race_totalE,
    
    poverty_rate = povertyE / poverty_totalE,
    unemployment_rate = unemploymentE / labor_forceE,
    
    education_low =
      (edu_001E + edu_002E + edu_003E + edu_004E + edu_005E + edu_006E +
         edu_007E + edu_008E + edu_009E + edu_010E + edu_011E +
         edu_012E + edu_013E + edu_014E + edu_015E) /
      total_educationE,
    
    uninsured_rate =
      (ins_001E + ins_002E + ins_003E + ins_004E) / total_insuranceE,
    
    language_isolation =
      (eng_001E + eng_002E + eng_003E + eng_004E + eng_005E + eng_006E +
         eng_007E + eng_008E + eng_009E + eng_010E + eng_011E + eng_012E +
         eng_013E + eng_014E + eng_015E + eng_016E + eng_017E + eng_018E +
         eng_019E + eng_020E + eng_021E + eng_022E + eng_023E + eng_024E) /
      total_languageE,
    
    disability_rate =
      (dis_001E + dis_002E + dis_003E + dis_004E + dis_005E + dis_006E +
         dis_007E + dis_008E + dis_009E + dis_010E + dis_011E + dis_012E) /
      total_disabilityE,
    
    housing_overcrowding =
      (crowd_001E + crowd_002E) / housing_crowding_totalE,
    
    owner_burden_rate =
      (owner_30_1E + owner_30_2E + owner_30_3E +
         owner_30_4E + owner_30_5E) / owner_occE,
    
    renter_burden_rate =
      (renter_30_1E + renter_30_2E + renter_30_3E +
         renter_30_4E + renter_30_5E) / renter_occE,
    
    outdoor_worker_rate = (agri_forestryE + construction_indE) / outdoor_totalE
   ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% #
  select(
    GEOID,
    total_popE,
    pct_white, pct_black, pct_asian, pct_hispanic,
    median_ageE,
    poverty_rate,
    median_incomeE,
    unemployment_rate,
    education_low,
    uninsured_rate,
    language_isolation,
    disability_rate,
    housing_overcrowding,
    owner_burden_rate,
    renter_burden_rate,
    outdoor_worker_rate
  ) %>%
  rename(
    median_age = median_ageE,
    median_income = median_incomeE,
    total_pop = total_popE
  )

# ----- COUNTY -----------------------------------------------

options(tigris_use_cache = TRUE)
#county_boundary <- tigris::counties(state = state, cb = TRUE) %>%
#  filter(NAME == county) %>%
#  st_make_valid() 

#zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
#  st_transform(st_crs(county_boundary)) %>% 
#  st_make_valid()

#county_zctas <- zcta_shapes %>%
#  st_filter(county_boundary) %>%
#  pull(ZCTA5CE20) 

#acs_clean <- acs_clean %>%
#  filter(GEOID %in% county_zctas) 

#nrow(acs_clean)
#stopifnot(all(acs_clean$GEOID %in% county_zctas))

# ----- STATE -------------------------------------------

state_boundary <- tigris::states(cb = TRUE) %>%
  filter(STUSPS == state) %>%
  st_make_valid()

zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
  st_transform(st_crs(state_boundary)) %>%
  st_make_valid()

state_zctas <- zcta_shapes %>%
  st_filter(state_boundary) %>%
  pull(ZCTA5CE20)

acs_clean <- acs_clean %>%
  filter(GEOID %in% state_zctas)

nrow(acs_clean)

# -------------------------------------------------------------

vars <- c(
  "pct_white", "pct_black", "pct_asian", "pct_hispanic",
  "median_age", "poverty_rate", "median_income",
  "unemployment_rate", "education_low",
  "uninsured_rate", "language_isolation",
  "disability_rate", "housing_overcrowding",
  "owner_burden_rate", "renter_burden_rate",
  "outdoor_worker_rate"
)

acs_zip_clean <- acs_clean %>%
  filter(if_all(all_of(vars), ~ !is.na(.) & is.finite(.)))

# ----- profile regression (PReMiuM) ---------------------------------

pr_data <- acs_zip_clean %>%
  st_drop_geometry() %>%
  select(all_of(vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ ntile(.x, 4))) %>%
  drop_na() %>%
  as.data.frame()

n_runs <- 3
diss_mats <- vector("list", n_runs)

set.seed(1)

for (i in seq_len(n_runs)) {
  set.seed(100 + i)
  
  prof <- profRegr(
    excludeY = TRUE,
    xModel = "Discrete",
    data = pr_data,
    covNames = vars,
    nBurn = 2000,
    nSweeps = 10000,
    nProgress = 1
  )
  diss_mats[[i]] <- calcDissimilarityMatrix(prof)
}
gc()

# ----- checking stability ---------------------------------------

v_list <- vector("list", n_runs)

for (i in 1:n_runs) {
  v_list[[i]] <- as.numeric(diss_mats[[i]]$disSimMat)
}

cor_vals <- c()
diff_vals <- c()

for (i in 1:(n_runs-1)) {
  for (j in (i+1):n_runs) {
    cor_vals <- c(cor_vals, cor(v_list[[i]], v_list[[j]]))
    diff_vals <- c(diff_vals, mean(abs(v_list[[i]] - v_list[[j]])))
  }
}

avg_dist <- rep(0, n_runs)

for (i in 1:n_runs) {
  avg_dist[i] <- mean(sapply(1:n_runs, function(j) mean(abs(v_list[[i]] - v_list[[j]]))))
}
best_run <- which.min(avg_dist)

rm(v_list)

# ----- labeling clusters -----------------------------------------------------

clusObj <- calcOptimalClustering(diss_mats[[best_run]])
pr_data$cluster <- clusObj$clustering

table(pr_data$cluster)

label_data <- acs_zip_clean %>%
  st_drop_geometry() %>%
  mutate(cluster = clusObj$clustering) %>%
  group_by(cluster) %>%
  summarise(
    across(all_of(vars), median),
    n = n()
  )

income_breaks <- quantile(acs_zip_clean$median_income, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
age_breaks <- quantile(acs_zip_clean$median_age, probs = c(0.25, 0.75), na.rm = TRUE)

income_label <- function(inc) {
  case_when(
    inc < income_breaks[1] ~ "Low income",
    inc < income_breaks[2] ~ "Mid-low income",
    inc < income_breaks[3] ~ "Mid-high income",
    TRUE                   ~ "High income"
  )
}

get_race_label <- function(pw, pb, pa, ph) {
  race_vals  <- c(White = pw, Black = pb, Asian = pa, Hispanic = ph)
  dominant_races <- names(race_vals[race_vals >= 0.40])
  if (length(dominant_races) == 0) {
    return("Mixed race")
  }
  paste(dominant_races, collapse = " & ")
}

cluster_profiles <- label_data %>%
  rowwise() %>%
  mutate(
    race_group   = get_race_label(pct_white, pct_black, pct_asian, pct_hispanic),
    age_group    = case_when(
      median_age <= age_breaks[1] ~ "Younger",
      median_age >= age_breaks[2] ~ "Older",
      TRUE                        ~ NA_character_
    ),
    income_group = income_label(median_income),
    label        = paste(na.omit(c(race_group, age_group, income_group)), collapse = " / ")
  ) %>%
  ungroup()

cluster_labels <- setNames(
  cluster_profiles$label,
  as.character(cluster_profiles$cluster)
)

# ---- SAVE ED-READY DATASET -----------------------------------

pr_data <- acs_zip_clean
pr_data$cluster <- clusObj$clustering
pr_data$cluster_label <- cluster_labels[as.character(pr_data$cluster)]

clustered_zctas <- pr_data %>%
  dplyr::select(GEOID, cluster, cluster_label, everything())

saveRDS(
  clustered_zctas,
  file = paste0("ed_analysis_ready/", gsub(" ", "_", state), "_zcta_clustered_for_ED.rds")
)

write.csv(
  clustered_zctas %>% st_drop_geometry(),
  file = paste0("ed_analysis_ready/", gsub(" ", "_", state), "_zcta_clustered_for_ED.csv"),
  row.names = FALSE
)

View(label_data)

# ----- PLOT MAP ---------------------------------------------------------------

zcta_cropped <- st_intersection(zcta_shapes, state_boundary)

zcta_cropped <- zcta_cropped %>%
  left_join(
    clustered_zctas %>%
      st_drop_geometry() %>%
      select(GEOID, cluster, cluster_label),
    by = c("ZCTA5CE20" = "GEOID")
  )

map <- ggplot(state_boundary) +
  geom_sf(fill = "grey90", color = "white", linewidth = 0.1) +
  geom_sf(data = zcta_cropped, aes(fill = cluster_label),
          color = "white", linewidth = 0.1) +
  coord_sf(
    xlim = c(st_bbox(state_boundary)["xmin"], st_bbox(state_boundary)["xmax"]),
    ylim = c(st_bbox(state_boundary)["ymin"], st_bbox(state_boundary)["ymax"]),
    expand = FALSE
  ) +
  scale_fill_viridis_d(
    option = "turbo",
    na.value = "transparent", 
    drop = FALSE
  ) +
  labs(title = paste("Socioeconomic Clusters of", state, "ZCTAs")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(
  filename = paste0("maps/", gsub(" ", "_", state), "_clusters.png"),
  plot = map,
  width = 8,
  height = 6,
  dpi = 300
)
rm(map)
gc()
# ----- adding to txt --------------------------------------------

sink(paste0("stability_outputs/", gsub(" ", "_", state), "_outputs.txt"))

cat("=== Correlations ===\n")
print(summary(cor_vals))

cat("\n=== Absolute Differences ===\n")
print(summary(diff_vals))

#cat("\n=== Cluster counts by run ===\n")
#run_cluster_k <- sapply(seq_along(diss_mats), function(i) {
#  cl <- PReMiuM::calcOptimalClustering(diss_mats[[i]])
#  length(unique(cl$clustering))
#})
#print(run_cluster_k)

cat("\nSelected run:", best_run, "out of", length(diss_mats), "\n")

cat("\n=== Cluster label counts ===\n")
cluster_label_table <- sort(table(pr_data$cluster_label), decreasing = TRUE)
print(cluster_label_table)

cat("\n=== ZCTA counts ===\n")
total_included <- nrow(pr_data)
excluded_ids <- setdiff(acs_clean$GEOID, pr_data$GEOID)
total_excluded <- length(excluded_ids)
cat("Included ZCTAs (clustered):", total_included, "\n")
cat("Excluded ZCTAs:", total_excluded, "\n")
cat("\nCluster sizes:\n")
print(table(pr_data$cluster))

cat("\n=== Excluded ZCTAs check ===\n")
excluded_ids <- setdiff(acs_clean$GEOID, pr_data$GEOID)
excluded_summary <- acs_clean %>%
  mutate(excluded_flag = ifelse(GEOID %in% excluded_ids, "excluded", "included")) %>%
  group_by(excluded_flag) %>%
  summarise(across(where(is.numeric), ~mean(is.na(.)))) %>%
  sf::st_drop_geometry()
print(excluded_summary, n = Inf, width = Inf)

sink()

saveRDS(diss_mats, paste0("rds/", gsub(" ", "_", state), "_diss_mats.rds"))
saveRDS(best_run,  paste0("rds/", gsub(" ", "_", state), "_best_run.rds"))

params <- list(
  seed_base = 100,
  n_runs = n_runs,
  nBurn = 2000,
  nSweeps = 10000,
  alpha_values = c(1)
)

saveRDS(params, paste0("rds/", gsub(" ", "_", state), "_model_params.rds"))
