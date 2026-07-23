library(tidycensus)
library(tidyverse)

library(dplyr)
library(cluster)
library(factoextra)
library(tidyr)

library("PReMiuM")
library(sf)
library(tigris)
library(car) # for vif()

year <- 2018
survey <- "acs5"
state <- "CA"
county <- "Los Angeles"

options(tigris_use_cache = TRUE)

#census_api_key("e3028f4138dfedd065a29bbe27a70a9", install = TRUE)

# ========================================================================
# VARIABLE SECTION
# ========================================================================

#vars <- load_variables(year, survey, cache = TRUE)
#View(vars)

acs <- get_acs(
  geography = "zcta",
  survey = survey,
  year = year,
  state = state, # always pull statewide, filter by county after
  variables = c(
    total_pop = "B01003_001",
    
    race_total = "B03002_001",
    white_nh = "B03002_003",
    black_nh = "B03002_004",
    asian_nh = "B03002_006",
    hispanic = "B03002_012",
    
    total_pop_sexage = "B01001_001",
    age17_male1 = "B01001_003", age17_male2 = "B01001_004",   
    age17_male3 = "B01001_005", age17_male4 = "B01001_006",  
    age17_female1 = "B01001_027", age17_female2 = "B01001_028", 
    age17_female3 = "B01001_029", age17_female4 = "B01001_030", 
    
    age65_male1 = "B01001_020", age65_male2 = "B01001_021",   
    age65_male3 = "B01001_022", age65_male4 = "B01001_023",   
    age65_male5 = "B01001_024", age65_male6 = "B01001_025", 
    age65_female1 = "B01001_044", age65_female2 = "B01001_045",
    age65_female3 = "B01001_046", age65_female4 = "B01001_047",
    age65_female5 = "B01001_048", age65_female6 = "B01001_049",
    
    #elderly_total = "B09021_022", elderly_alone = "B09021_023",
    
    median_income = "B19013_001", 
    
    poverty = "B17001_002",
    poverty_total = "B17001_001",
    
    #unemployment = "B23025_005",
    #labor_force = "B23025_003",
    
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
    
    #owner_occ = "B25106_002",
    #owner_30_1 = "B25106_006", owner_30_2 = "B25106_010",
    #owner_30_3 = "B25106_014", owner_30_4 = "B25106_018",
    #owner_30_5 = "B25106_022",
    
    renter_occ = "B25106_024",
    renter_30_1 = "B25106_028", renter_30_2 = "B25106_032",
    renter_30_3 = "B25106_036", renter_30_4 = "B25106_040",
    renter_30_5 = "B25106_044",
    
    housing_crowding_total = "B25014_001",
    crowd_001 = "B25014_007", crowd_002 = "B25014_013",
    
    #vehicles_total = "B08201_001",
    #no_vehicles = "B08201_002",
    
    housing_units_total = "B25024_001", 
    boat_rv_van = "B25024_011", mobile_home = "B25024_010",
    
    #phone_total = "B25043_001",
    #phone_owner_none = "B25043_007",
    #phone_renter_none = "B25043_016",
    
    #yrbuilt_total = "B25034_001",
    #built_1970_79 = "B25034_007",
    #built_1960_69 = "B25034_008",
    #built_1950_59 = "B25034_009",
    #built_1940_49 = "B25034_010",
    #built_pre1939 = "B25034_011",
    
    #tenure_total = "B25003_001",
    #renter_occ_tenure = "B25003_003",
    
    outdoor_total = "C24050_001",
    agri_forestry = "C24050_002",
    construction_ind = "C24050_003"),
  output = "wide",
  geometry = TRUE)

acs_clean <- acs %>% filter(total_popE > 0) %>% mutate(total_popE,
    
    pct_white = white_nhE / race_totalE,
    pct_black = black_nhE / race_totalE,
    pct_asian = asian_nhE / race_totalE,
    pct_hispanic = hispanicE / race_totalE,
    
    poverty_rate = povertyE / poverty_totalE,
    #unemployment_rate = unemploymentE / labor_forceE,
    #vehicle_rate = no_vehiclesE / vehicles_totalE,
    
    pct_age17 = (age17_male1E + age17_male2E + age17_male3E + age17_male4E +
                   age17_female1E + age17_female2E + age17_female3E + age17_female4E) 
    / total_pop_sexageE,
    
    pct_age65 = (age65_male1E + age65_male2E + age65_male3E + age65_male4E + 
                   age65_male5E + age65_male6E + age65_female1E + age65_female2E + 
                   age65_female3E + age65_female4E + age65_female5E + age65_female6E) 
    / total_pop_sexageE,
    
    #elderly_alone_rate = elderly_aloneE / elderly_totalE, 
    
    education_low =
      (edu_001E + edu_002E + edu_003E + edu_004E + edu_005E + edu_006E +
         edu_007E + edu_008E + edu_009E + edu_010E + edu_011E +
         edu_012E + edu_013E + edu_014E + edu_015E) /
      total_educationE,
    
    uninsured_rate = (ins_001E + ins_002E + ins_003E + ins_004E) / total_insuranceE,
    
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
    
    housing_overcrowding = (crowd_001E + crowd_002E) / housing_crowding_totalE,
    
    #owner_burden_rate = (owner_30_1E + owner_30_2E + owner_30_3E +
    #     owner_30_4E + owner_30_5E) / owner_occE,
    
    renter_burden_rate = (renter_30_1E + renter_30_2E + renter_30_3E +
         renter_30_4E + renter_30_5E) / renter_occE,
    
    alt_housing_rate = (mobile_homeE + boat_rv_vanE) / housing_units_totalE,

    #no_phone_rate = (phone_owner_noneE + phone_renter_noneE) / phone_totalE,   
    #renter_rate = renter_occ_tenureE / tenure_totalE,
    
    #old_housing_rate = (built_1970_79E + built_1960_69E + built_1950_59E +
    #                      built_1940_49E + built_pre1939E) / yrbuilt_totalE,
    
    outdoor_worker_rate = (agri_forestryE + construction_indE) / outdoor_totalE
   ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>% #
  select(GEOID, total_popE, pct_white, pct_black, pct_asian, pct_hispanic,
         poverty_rate, median_incomeE,
         education_low, uninsured_rate, language_isolation, disability_rate, 
         outdoor_worker_rate, pct_age17, pct_age65,
         renter_burden_rate, housing_overcrowding, alt_housing_rate) %>%
  rename(median_income = median_incomeE, total_pop = total_popE)

vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic",
  "poverty_rate", "renter_burden_rate", "housing_overcrowding",
  "median_income", "education_low", "pct_age17",
  "pct_age65", "uninsured_rate", "language_isolation", "disability_rate",
  "outdoor_worker_rate",
  "alt_housing_rate")

rm(acs)

# ===========================================================================
# COUNTY
# =======================================================================

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

# ===========================================================================
# STATE
# =======================================================================

ca_zcta_codes <- read_delim("datasets/2018_Gaz_zcta_national.txt", delim = "\t") %>%
  janitor::clean_names() %>%
  mutate(geoid = as.character(geoid))

ca_zcta_codes <- ca_zcta_codes %>%
  filter(as.numeric(geoid) >= 90001 & as.numeric(geoid) <= 96162) %>%
  pull(geoid)

acs_clean <- acs_clean %>% filter(GEOID %in% ca_zcta_codes)

nrow(acs_clean)
gc()

# =======================================================================
# PROFILE REGRESSION (PReMiuM)
# ========================================================================

pr_data <- acs_clean %>%
  st_drop_geometry() %>%
  select(GEOID, all_of(vars)) %>%
  mutate(across(all_of(vars), as.numeric)) %>%
  mutate(across(all_of(vars), ~ ntile(.x, 4))) %>%
  as.data.frame()

# LOAD EXISTING RUNS IF PRESENT
acs_diss_path <- paste0("rds/", gsub(" ", "_", state), "_ACS_diss_mats_MIN.rds")
acs_rho_path <- paste0("rds/", gsub(" ", "_", state), "_ACS_rho_list_MIN.rds")

acs_diss_mats <- if (file.exists(acs_diss_path)) readRDS(acs_diss_path) else list()
acs_rho_list <- if (file.exists(acs_rho_path))  readRDS(acs_rho_path)  else list()

n_existing <- length(acs_diss_mats)
n_new <- 1
n_runs <- n_existing + n_new

acs_diss_mats <- c(acs_diss_mats, vector("list", n_new))
acs_rho_list <- c(acs_rho_list,  vector("list", n_new))

for (i in (n_existing + 1):(n_existing + n_new)) {
  set.seed(100 + i)
  output_stem <- paste0("acs_prof_run", i)
  
  prof <- profRegr(
    excludeY = TRUE,
    xModel = "Discrete",
    data = pr_data %>% select(all_of(vars)),
    covNames = vars,
    varSelectType = "Continuous",
    output = output_stem,
    nBurn = 2000,
    nSweeps = 10000,
    nProgress = 1)
  
  acs_diss_mats[[i]] <- calcDissimilarityMatrix(prof)
  
  rho_samples <- read.table(paste0(output_stem, "_rho.txt"))
  colnames(rho_samples) <- vars
  acs_rho_list[[i]] <- colMeans(rho_samples)}

saveRDS(acs_diss_mats, acs_diss_path)
saveRDS(acs_rho_list,  acs_rho_path)

# ==========================================================================
# CHECK STABILITY
# ========================================================================

n_runs <- 3
v_list <- vector("list", n_runs)
for (i in 1:n_runs) {v_list[[i]] <- as.numeric(acs_diss_mats[[i]]$disSimMat)}

cor_vals <- c()
diff_vals <- c()

for (i in 1:(n_runs-1)) {
  for (j in (i+1):n_runs) {
    cor_vals <- c(cor_vals, cor(v_list[[i]], v_list[[j]]))
    diff_vals <- c(diff_vals, mean(abs(v_list[[i]] - v_list[[j]])))}}

avg_dist <- rep(0, n_runs)

for (i in 1:n_runs) {
  avg_dist[i] <- mean(sapply(1:n_runs, function(j) mean(abs(v_list[[i]] - v_list[[j]]))))}
best_run <- which.min(avg_dist)

rm(v_list)
rm(prof)
rm(rho_samples)

# =========================================================================
# LABEL CLUSTERS 
# ==========================================================================

clusObj <- calcOptimalClustering(acs_diss_mats[[best_run]])
pr_data$cluster <- clusObj$clustering

saveRDS(pr_data, "rds/pr_data_with_clusters_MIN.rds")
rm(clusObj)

cluster_profiles <- acs_clean %>% 
  st_drop_geometry() %>%
  left_join(pr_data %>% select(GEOID, cluster), by = "GEOID") %>%
  group_by(cluster) %>% 
  summarise(across(all_of(c(vars)), #"pct_white", "pct_black", "pct_asian", "pct_hispanic"
                   ~ median(.x, na.rm = TRUE)), n = n())

income_breaks <- quantile(acs_clean$median_income, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
age65_breaks <- quantile(acs_clean$pct_age65, probs = 0.75, na.rm = TRUE)
age17_breaks <- quantile(acs_clean$pct_age17, probs = 0.75, na.rm = TRUE)

get_race_label <- function(pw, pb, pa, ph) {
  race_vals <- c(White = pw, Black = pb, Asian = pa, Hispanic = ph)
  sorted_vals <- sort(race_vals, decreasing = TRUE)
  top_group <- names(sorted_vals)[1]
  second_group <- names(sorted_vals)[2]
  top_val <- sorted_vals[1]
  second_val <- sorted_vals[2]
  if (sum(race_vals > 0) == 2 && (top_val - second_val) < 0.15) {
    return(paste(sort(c(top_group, second_group)), collapse = " & "))
  } else if (top_val < 0.4 && sum(race_vals >= 0.25) >= 3)
  { return("Mixed racial composition") }
  else if (top_val >= 0.40 && second_val >= 0.40)
  { return(paste(sort(c(top_group, second_group)), collapse = " & ")) }
  else if (second_val >= 0.30 && (top_val - second_val) < 0.15)
  { return(paste(sort(c(top_group, second_group)), collapse = " & ")) }
  else if (top_val >= 0.40) { return(top_group) } else
    { return("Mixed racial composition") } }

income_label <- function(inc) {
  case_when(
    inc < income_breaks[1] ~ "Low income",
    inc < income_breaks[2] ~ "Mid-low income",
    inc < income_breaks[3] ~ "Mid-high income",
    TRUE ~ "High income")}

cluster_profiles <- cluster_profiles %>%
  rowwise() %>% mutate(
    race_group   = get_race_label(pct_white, pct_black, pct_asian, pct_hispanic),
    age_group = case_when(pct_age65 >= age65_breaks ~ "Older",
      pct_age17 >= age17_breaks ~ "Younger", TRUE ~ NA_character_),
    income_group = income_label(median_income),
    label = paste(na.omit(c(race_group, age_group, income_group)), collapse = " / ")
  ) %>% ungroup()

cluster_labels <- setNames(cluster_profiles$label,
  as.character(cluster_profiles$cluster))

# =======================================================================
# SAVE ED-READY DATASET
# ======================================================================

pr_data$cluster_label <- cluster_labels[as.character(pr_data$cluster)]
cluster_assignment <- pr_data %>% select(GEOID, cluster, cluster_label)

clustered_zctas <- acs_clean %>% st_drop_geometry() %>%
  left_join(cluster_assignment, by = "GEOID") %>%
  dplyr::select(GEOID, cluster, cluster_label, everything())

clustered_zctas <- clustered_zctas %>% left_join(
    cluster_profiles %>% select(cluster, all_of(vars)) %>% 
      rename_with(~paste0(.x, "_cluster_median"), all_of(vars)), by = "cluster") %>%
  mutate(across(all_of(vars), ~ if_else(is.na(.x), get(paste0(cur_column(), "_cluster_median")), .x))) %>%
  select(-ends_with("_cluster_median"))

saveRDS(clustered_zctas, file = paste0(
  "ed_analysis_ready/", gsub(" ", "_", state), "_ACS_zcta_clustered_for_ED_MIN.rds"))

write.csv(clustered_zctas, file = paste0(
  "ed_analysis_ready/", gsub(" ", "_", state), "_ACS_zcta_clustered_for_ED_MIN.csv"),
          row.names = FALSE)

# ===========================================================================
# MAP
# =========================================================================

zcta_cache_path <- "rds/zcta_shapes_2010_national.rds"

if (!file.exists(zcta_cache_path)) {
  cat("No cached 2010 ZCTA shapefile found — downloading once\n")
  zcta_shapes_acs <- zctas(cb = TRUE, year = 2010)
  saveRDS(zcta_shapes_acs, zcta_cache_path)} else {
    zcta_shapes_acs <- readRDS(zcta_cache_path)}

state_boundary <- tigris::states(cb = TRUE) %>%
  filter(STUSPS == state) %>% st_make_valid()

zcta_shapes_acs <- zcta_shapes_acs %>%
  filter(ZCTA5 %in% clustered_zctas$GEOID) %>%
  st_transform(st_crs(state_boundary)) %>% st_make_valid()

zcta_shapes_acs <- st_intersection(zcta_shapes_acs, state_boundary) %>%
  left_join(clustered_zctas %>% st_drop_geometry() %>% select(GEOID, cluster, cluster_label),
    by = c("ZCTA5" = "GEOID"))

map_acs <- ggplot(state_boundary) +
  geom_sf(fill = "grey90", color = "white", linewidth = 0.1) +
  geom_sf(data = zcta_shapes_acs, aes(fill = cluster_label),
          color = "white", linewidth = 0.1) + coord_sf(
            xlim = c(st_bbox(state_boundary)["xmin"], st_bbox(state_boundary)["xmax"]),
            ylim = c(st_bbox(state_boundary)["ymin"], st_bbox(state_boundary)["ymax"]),
            expand = FALSE  ) +
  scale_fill_viridis_d(option = "turbo", na.value = "transparent", drop = FALSE) +
  labs(title = paste("Socioeconomic Clusters of", state, "ZCTAs")) +
  theme_minimal() + theme(legend.position = "right", legend.direction = "vertical",
    legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7),
    legend.title = element_text(size = 9)) + guides(fill = guide_legend(ncol = 1))

ggsave(filename = paste0("maps/", state, "_ACS_clusters_MIN.png"),
       plot = map_acs, width = 8, height = 6, dpi = 300)
rm(map_acs)
rm(clustered_zctas)
gc()

# ============================================================
# CHECKS
# ============================================================

rho_summary_acs <- do.call(rbind, acs_rho_list) %>%
  as.data.frame() %>% setNames(vars) %>% mutate(run = 1:n_runs)

cor_long_acs <- pr_data %>% select(all_of(vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  as.table() %>% as.data.frame() %>% 
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  mutate(var1 = as.character(var1), var2 = as.character(var2)) %>%
  filter(var1 != var2) %>%
  mutate(pair_key = paste(pmin(var1, var2), pmax(var1, var2))) %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  select(var1, var2, correlation) %>% arrange(desc(abs(correlation))) %>%
  filter(abs(correlation) > 0.6, !is.nan(correlation)) %>% as.data.frame()

#vif_results <- data.frame(variable = vars, VIF = sapply(vars, function(v) {
#  model <- lm(as.formula(paste(v, "~ .")), data = pr_data %>% select(all_of(vars)))
#  max(vif(model), na.rm = TRUE)
#})) %>% arrange(desc(VIF))
pr_data$dummy_outcome <- rnorm(nrow(clustered_zctas))
model <- lm(dummy_outcome ~ ., data = pr_data[, c("dummy_outcome", vars)])
vif_results <- data.frame(variable = names(vif(model)), VIF = vif(model)) %>%
  arrange(desc(VIF))

rho_avg_acs <- colMeans(rho_summary_acs[, vars])

total_included <- nrow(pr_data)
total_excluded <- length(setdiff(acs_clean$GEOID, pr_data$GEOID))
missing_zcta_codes <- setdiff(ca_zcta_codes, acs_clean$GEOID)

missingness_summary <- lapply(vars, function(v) {
  x <- acs_clean[[v]]
  data.frame(variable = v, n_missing = sum(is.na(x)),
    pct_missing = round(mean(is.na(x)) * 100, 2),
    n_zero = sum(x == 0, na.rm = TRUE),
    pct_zero = round(mean(x == 0, na.rm = TRUE) * 100, 2))}) %>%
  bind_rows() %>% arrange(desc(pct_missing))

saveRDS(rho_summary_acs, paste0("rds/", state, "_ACS_rho_by_run_MIN.rds"))
saveRDS(acs_diss_mats, paste0("rds/", state, "_ACS_diss_mats_MIN.rds"))
saveRDS(best_run, paste0("rds/", state, "_ACS_best_run_MIN.rds"))

# ------ TXT FILE ------------------------------------------------------

sink(paste0("outputs/", state, "_ACS_outputs_MIN.txt"))
cat("=== Correlations ===\n")
print(summary(cor_vals))
cat("\n\n=== Absolute Differences ===\n")
print(summary(diff_vals))
cat("\n\nSelected run:", best_run, "out of", n_runs, "\n")

cat("\n\n=== ZCTA counts ===\n")
cat("Total CA ZCTAs:", length(ca_zcta_codes), "\n")
cat("ACS included:", nrow(acs_clean), " | excluded:", length(ca_zcta_codes) - nrow(acs_clean), "\n")
cat("\nMissing ZCTA codes (in reference list, not in acs_clean):\n")
print(missing_zcta_codes)

cat("\n\n=== Cluster label counts ===\n")
print(sort(table(pr_data$cluster_label), decreasing = TRUE))

cat("\n\n=== PReMiuM Posterior Inclusion Probabilities (mean rho, averaged across runs) ===\n")
print(sort(rho_avg_acs, decreasing = TRUE))
cat("\n\n=== Pairwise Correlations (|r| > 0.6) ===\n")
print(cor_long_acs)
cat("\n\n=== Variance Inflation Factor Check ===\n")
print(vif_results)
cat("\n\n=== Missingness summary ===\n")
print(missingness_summary)

cat("\n\n=== Completeness by ZCTA population ===\n")
print(acs_clean %>% st_drop_geometry() %>%
    mutate(complete = complete.cases(select(., all_of(vars)))) %>%
    group_by(complete) %>%
    summarise(mean_pop = mean(total_pop, na.rm = TRUE),
              median_pop = median(total_pop, na.rm = TRUE),
              n = n()))
sink()

rm(vif_check_acs)
rm(acs_diss_mats)
rm(acs_rho_list)
gc()
