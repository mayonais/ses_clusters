library(tidycensus)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(1)

# cluster assignments
clusters <- read.csv("ed_analysis_ready/CA_ACS_zcta_clustered_for_ED.csv",
                     colClasses = c(GEOID = "character"))

ca_pop <- get_acs(
  geography = "zcta",
  variables = "B01003_001",
  year = 2018,
  survey = "acs5",
  state = "CA"
) %>%
  filter(estimate > 0) %>%
  transmute(
    zip = GEOID,
    pop = estimate
  ) %>%
  semi_join(clusters, by = c("zip" = "GEOID"))

ses_clusters <- clusters %>%
  dplyr::select(GEOID, cluster, cluster_label)

ca_pop <- ca_pop %>%
  left_join(ses_clusters, by = c("zip" = "GEOID")) %>%
  filter(!is.na(cluster)) # drop any zips without a cluster assignment

n_zips <- nrow(ca_pop)
pop <- ca_pop$pop

# ----- study period and summer peak ------------------------------------------

dates <- do.call(c, lapply(2017:2018, function(y) {
  seq.Date(as.Date(paste0(y, "-05-01")), as.Date(paste0(y, "-10-31")), by = "day")
}))

year_of_date <- format(dates, "%Y")
doy <- as.numeric(format(dates, "%j"))

# real age-adjusted HRI ED visit rates, per 10,000 (input target)
# "2008" = 11.6, "2009" = 10.8, "2010" = 9.5,  "2011" = 10.1, "2012" = 12.9, "2013" = 13.1,
# "2014" = 12.9, "2015" = 13.5, "2016" = 13.3, 
real_ed_rates <- c("2017" = 19.6, "2018" = 14.9)

# probability of an ED visit per day, per person
peak_doy <- 200
season_sd <- 50
season_shape <- exp(-(doy - peak_doy)^2 / (2 * season_sd^2))

daily_rate_by_date <- numeric(length(dates))
for (yr in names(real_ed_rates)) {
  idx <- year_of_date == yr
  w <- season_shape[idx]
  w <- w / sum(w)
  daily_rate_by_date[idx] <- w * (real_ed_rates[yr] / 10000)}

# ----- vulnerability profile by cluster ---------------------------------------

source("analysis_functions.R")

cluster_profiles <- clusters %>%
  group_by(cluster, cluster_label) %>%
  summarise(
    across(
      c(
        pct_white, pct_black, pct_asian, pct_hispanic, pct_age17, pct_age65,
        poverty_rate, median_income, unemployment_rate, education_low,
        uninsured_rate, language_isolation,
        disability_rate, vehicle_rate, no_phone_rate,
        elderly_alone_rate, renter_rate,
        housing_overcrowding,
        owner_burden_rate, renter_burden_rate,
        alt_housing_rate, old_housing_rate,
        outdoor_worker_rate
      ),
      ~ median(.x, na.rm = TRUE)
    ),
    n = n(),
    .groups = "drop"
  )

cluster_profiles <- clusters %>%
  compute_vulnerability_score() %>%
  group_by(cluster, cluster_label) %>%
  summarise(vulnerability_score = mean(vulnerability_score, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  select(cluster, cluster_label, vulnerability_score)

# ----- build cluster-level risk multipliers -------------------------------

cluster_effect_scale <- 0.20

cluster_multiplier <- cluster_profiles %>%
  mutate(
    true_multiplier = exp(
      as.numeric(scale(vulnerability_score)) *
        cluster_effect_scale
    )
  ) %>%
  select(cluster, true_multiplier)

# explicit named vector for lookup — keyed by cluster id, not row position
mult_vec <- setNames(cluster_multiplier$true_multiplier,
                     as.character(cluster_multiplier$cluster))

zip_cluster_multiplier <- mult_vec[as.character(ca_pop$cluster)]

# zip-level random variation
zip_noise_raw <- pmax(rnorm(n_zips, mean = 1, sd = 0.3), 0.05)
zip_noise <- zip_noise_raw / mean(zip_noise_raw)
risk_multiplier <- zip_cluster_multiplier * zip_noise

# combines cluster-level risk differences with ZIP-level variation
cluster_truth <- cluster_multiplier %>%
  left_join(
    cluster_profiles %>% select(cluster, cluster_label, vulnerability_score),
    by = "cluster"
  )

# ----- simulate individual ED visits ----------------------------------

icd_codes <- c("276.51", "E86.0", # dehydration
               "584", "585", "586", "N17", "N18", "N19", # kidney failure
               "992", "T67") # heat illness/stroke
icd_probs <- c(rep(0.80 / 2, 2), rep(0.13 / 6, 6), rep(0.07 / 2, 2))

# zip-day-level random variation to introduce overdispersion
nb_theta <- 0.5
zip_visits <- map(seq_len(n_zips), function(z) { # a vector of expected visit counts per day
  base_lambda <- pop[z] * risk_multiplier[z] * daily_rate_by_date
  daily_variation <- rgamma(length(dates), shape = nb_theta, rate = nb_theta)
  lambda <- base_lambda * daily_variation
  daily_counts <- rpois(length(dates), lambda)
  total_visits <- sum(daily_counts)
  if (total_visits == 0) return(NULL)
    
  visit_dates <- rep(dates, daily_counts) # one entry per individual simulated visit, not per day 
  data.frame(
    zip = ca_pop$zip[z],
    date = visit_dates,
    icd_code = sample(icd_codes, length(visit_dates), replace = TRUE, prob = icd_probs))})

ed_visits <- bind_rows(zip_visits)
rm(zip_visits)
gc()

# checks
head(ed_visits)
nrow(ed_visits)
length(unique(ed_visits$zip)) # zips with >=1 visit

# ------ aggregate to ZCTA-day counts, fill zeros, merge pop --------------------------------------------------------

ed_combined <- expand.grid(
  GEOID = unique(ca_pop$zip),
  date  = dates,
  stringsAsFactors = FALSE) %>%
  left_join(
    ed_visits %>%
      mutate(date = as.Date(date)) %>%
      group_by(zip, date) %>%
      summarize(hri_ed = n(), .groups = "drop") %>%
      rename(GEOID = zip),
    by = c("GEOID", "date")) %>%
  mutate(hri_ed = replace_na(hri_ed, 0)) %>%
  left_join(ca_pop %>% dplyr::select(zip, pop) %>% rename(GEOID = zip), by = "GEOID")

gc()

write.csv(ed_combined, "CA_simulated_ed/daily_counts_2017_2018.csv", row.names = FALSE)
saveRDS(ed_combined, "CA_simulated_ed/daily_counts_2017_2018.rds")

# ----- CHECKS ------------------------------------------------------------------

# overall simulated rate vs. target real_ed_rates (9.5-19.6 per 10,000/yr)
zcta_summary <- ed_combined %>% group_by(GEOID, pop) %>%
  summarize(total_events = sum(hri_ed), .groups = "drop")

# rate
sum(zcta_summary$total_events) / sum(zcta_summary$pop) * 10000 / 2

# per-year calibration check
ed_combined <- ed_combined %>% mutate(year = format(date, "%Y")) 
year_summary <- ed_combined %>% group_by(year) %>% 
  summarize( total_events = sum(hri_ed), total_pop_days = sum(pop), 
             n_days = n_distinct(date), .groups = "drop" ) %>% 
  mutate( avg_pop = total_pop_days / n_days, 
          simulated_rate = total_events / avg_pop * 10000, 
          target_rate = real_ed_rates[year] )

print(year_summary %>% dplyr::select(year, simulated_rate, target_rate))

# verify cluster risk differences were incorporated
cluster_check <- ed_combined %>% left_join(ca_pop %>% select(zip, cluster) %>%
              mutate(cluster = as.character(cluster)),
              by = c("GEOID" = "zip")) %>%
  group_by(cluster) %>%  summarize(
    rate_per_10000_yr = sum(hri_ed) / sum(pop) * 10000 / 2, .groups = "drop") %>%
  left_join(cluster_truth %>% mutate(cluster = as.character(cluster)),
    by = "cluster")

write.csv(cluster_check, "CA_simulated_ed/cluster_multiplier.csv", row.names = FALSE)

cor(
  cluster_check$rate_per_10000_yr,
  cluster_check$true_multiplier,
  method = "spearman"
)
