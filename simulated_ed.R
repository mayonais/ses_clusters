library(tidycensus)
library(dplyr)
library(purrr)
library(MASS)
library(tidyr)

set.seed(1)

ca_pop <- get_acs(
  geography = "zcta",
  variables = "B01003_001",
  year = 2020,
  survey = "acs5"
) %>%
  filter(estimate > 0) %>%
  transmute(zip = GEOID, pop = estimate) %>%
  filter(as.numeric(zip) >= 90001 & as.numeric(zip) <= 96162)

n_zips <- nrow(ca_pop)
pop <- ca_pop$pop

dates <- do.call(c, lapply(2008:2018, function(y) {
  seq.Date(as.Date(paste0(y, "-05-01")), as.Date(paste0(y, "-10-31")), by = "day")
}))
n_days_per_year <- table(format(dates, "%Y"))
year_of_date <- format(dates, "%Y")

# real age-adjusted HRI ED visit rates, per 10,000 (input target)
real_ed_rates <- c(
  "2008" = 11.6, "2009" = 10.8, "2010" = 9.5,  "2011" = 10.1, "2012" = 12.9,
  "2013" = 13.1, "2014" = 12.9, "2015" = 13.5, "2016" = 13.3, "2017" = 19.6,
  "2018" = 14.9
)

# convert annual rate per 10,000 to daily rate per person
daily_rate_per_person <- (real_ed_rates / 10000) / as.numeric(n_days_per_year[names(real_ed_rates)])
daily_rate_by_date <- daily_rate_per_person[year_of_date]

# allows for zip level variation
sdlog_val <- 0.5
alpha <- rlnorm(n_zips, meanlog = -sdlog_val^2 / 2, sdlog = sdlog_val)

icd_codes <- c("276.51", "E86.0", # dehydration
               "584", "585", "586", "N17", "N18", "N19", # kidney failure
               "992", "T67") # heat illness/stroke

icd_probs <- c(
  rep(0.80 / 2, 2),
  rep(0.13 / 6, 6),  
  rep(0.07 / 2, 2)  
)

zip_visits <- map(seq_len(n_zips), function(z) {
  lambda <- pop[z] * alpha[z] * daily_rate_by_date # a vector of expected visit counts per day
  daily_counts <- rpois(length(dates), lambda)
  total_visits <- sum(daily_counts)
  if (total_visits == 0) return(NULL)
  
  visit_dates <- rep(dates, daily_counts) # one entry per individual simulated visit, not per day 
  data.frame(
    zip = ca_pop$zip[z],
    date = visit_dates,
    icd_code = sample(icd_codes, length(visit_dates), replace = TRUE)
  )
})

ed_visits <- bind_rows(zip_visits)
rm(zip_visits)
gc()

# checks
head(ed_visits)
nrow(ed_visits)
length(unique(ed_visits$zip)) # zips with >=1 visit
range(ed_visits$date)

#write.csv(ed_visits, "CA_ed_simulated_events.csv", row.names = FALSE)

# ---------------------------------------------------------------------------

ed_zcta_day <- ed_visits %>%
  mutate(date = as.Date(date)) %>%
  group_by(zip, date) %>%
  summarize(hri_ed = n(), .groups = "drop") %>%
  rename(GEOID = zip)

# aggregate ED events to ZCTA-day counts, combine full ZCTA × day,
# fill missing days with 0, and merge population
ed_combined <- expand.grid(
  GEOID = unique(ed_visits$zip),
  date  = seq.Date(as.Date("2008-05-01"), as.Date("2018-10-31"), by = "day"),
  stringsAsFactors = FALSE) %>%
  left_join(
    ed_visits %>%
      mutate(date = as.Date(date)) %>%
      group_by(zip, date) %>%
      summarize(hri_ed = n(), .groups = "drop") %>%
      rename(GEOID = zip),
    by = c("GEOID", "date")) %>%
  mutate(hri_ed = replace_na(hri_ed, 0)) %>%
  left_join(ca_pop %>% rename(GEOID = zip), by = "GEOID")

gc()
head(ed_combined)

write.csv(ed_combined, "ed_analysis_ready/CA_ed_simulated_counts.csv", row.names = FALSE)
saveRDS(ed_combined, "ed_analysis_ready/CA_ed_simulated_counts.rds")
