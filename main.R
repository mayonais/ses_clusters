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

# ----- cleaning -------------------------------------------------------------

options(tigris_use_cache = TRUE)

# download county boundaries
county_boundary <- tigris::counties(state = state, cb = TRUE) %>%
  filter(NAME == county) %>%
  st_make_valid() # fix invalid geometries

zcta_shapes <- zctas(cb = TRUE, year = 2020) %>%
  st_transform(st_crs(county_boundary)) %>% # transform zctas to have same coordinate system
  st_make_valid()

county_zctas <- zcta_shapes %>%
  st_filter(county_boundary) %>%
  pull(ZCTA5CE20) # filter zcta by county

acs_clean <- acs_clean %>%
  filter(GEOID %in% county_zctas) # filter acs data to county zcta

nrow(acs_clean)

# confirm all selected zctas intersect county
stopifnot(all(acs_clean$GEOID %in% county_zctas))

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
  drop_na(all_of(vars))

X <- scale(st_drop_geometry(acs_zip_clean[, vars]))

colSums(is.na(X))

# ----- profile regression (PReMiuM) ---------------------------------

pr_data <- acs_zip_clean %>%
  st_drop_geometry() %>%
  select(all_of(vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ ntile(.x, 4))) %>%
  drop_na() %>%
  as.data.frame()

pr_data$y <- sample(0L:1L, nrow(pr_data), replace = TRUE)

n_runs <- 10
diss_mats <- vector("list", n_runs)

covNames <- setdiff(names(pr_data), "y")

pr_df <- data.frame(
  lapply(pr_data[, covNames], function(x) as.vector(drop(as.matrix(x)))),
  check.names = FALSE
)
pr_df$y <- as.vector(pr_data$y)
covNames <- setdiff(names(pr_df), "y")

set.seed(1)

for (i in seq_len(n_runs)) {
  set.seed(100 + i)
  
  prof <- profRegr(
    yModel = "Bernoulli",
    xModel = "Discrete",
    data = pr_df,
    covNames = covNames,
    outcome = "y",
    nBurn = 5000,
    nSweeps = 10000
  )
  diss_mats[[i]] <- calcDissimilarityMatrix(prof)
}

# ----- checking stability ---------------------------------------

v_list <- vector("list", 10)

for (i in 1:10) {
  v_list[[i]] <- as.numeric(diss_mats[[i]]$disSimMat)
}

cor_vals <- c()
diff_vals <- c()

for (i in 1:9) {
  for (j in (i+1):10) {
    cor_vals <- c(cor_vals, cor(v_list[[i]], v_list[[j]]))
    diff_vals <- c(diff_vals, mean(abs(v_list[[i]] - v_list[[j]])))
  }
}

avg_dist <- rep(0, 10)

for (i in 1:10) {
  avg_dist[i] <- mean(sapply(1:10, function(j) mean(abs(v_list[[i]] - v_list[[j]]))))
}
best_run <- which.min(avg_dist)

# ----- labeling and plotting the (most central) best run --------------------

clusObj <- calcOptimalClustering(diss_mats[[best_run]])
pr_data$cluster <- clusObj$clustering

table(pr_data$cluster)

cluster_profiles <- pr_data %>%
  group_by(cluster) %>%
  summarise(
    across(where(is.numeric), median),
    n = n()
  )

write.csv(
  cluster_profiles,
  paste0("cluster_profiles/", gsub(" ", "_", county), "_cluster_profiles.csv"),
  row.names = FALSE
)

cluster_labels <- c(
  "1"  = "Hispanic / Younger / Low income",
  "12" = "Hispanic / Mid-high income",
  
  "2"  = "White / Mid-high income",
  "5"  = "White / Older / Mid-low income",
  "8"  = "White / Older / High income",
  
  "3"  = "Mixed race / Mid-low income",
  "7"  = "Mixed race / High income",
  "11" = "Mixed race / Mid-high income",
  
  "4"  = "Black / Younger / Low income",
  "6"  = "Black / Low income",
  
  "9"  = "Asian / High income",
  "10" = "Asian / Mid-low income"
)

pr_data <- acs_zip_clean
pr_data$cluster <- clusObj$clustering
pr_data$cluster_label <- cluster_labels[as.character(pr_data$cluster)]

map <- ggplot(pr_data) +
  geom_sf(aes(fill = cluster_label),
          color = "white", linewidth = 0.1) +
  scale_fill_brewer(palette = "Set3", na.value = "grey80") +
  labs(title = paste("Socioeconomic Clusters of", county, "ZCTAs")) +
  theme_minimal()

ggsave(
  filename = paste0("maps/", gsub(" ", "_", county), "_clusters.png"),
  plot = map,
  width = 8,
  height = 6,
  dpi = 300
)

# ----- adding to txt --------------------------------------------

sink(paste0("stability_outputs/", gsub(" ", "_", county), "_outputs.txt"))

cat("=== Correlations ===\n")
print(summary(cor_vals))

cat("\n=== Absolute Differences ===\n")
print(summary(diff_vals))

cat("\n=== Cluster counts by run ===\n")
run_cluster_k <- sapply(seq_along(diss_mats), function(i) {
  tmp <- capture.output(
    cl <- calcOptimalClustering(diss_mats[[i]])
  )
  length(unique(cl$clustering))
})
print(run_cluster_k)

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
