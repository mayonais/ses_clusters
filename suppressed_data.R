library(tidyverse)
library(lubridate)
library(dplyr)
library(tidyr)

state <- "CA"

clustered_zctas <- readRDS(paste0(
  "ed_analysis_ready/", gsub(" ", "_", state), "_ACS_zcta_clustered_for_ED_MIN.rds"))

daily_ed <- read_csv("datasets/Dummy_20241005_TotalERDaily.csv") %>%
  mutate(zip = as.character(zip), date = as.Date(date, format = "%m/%d/%Y"),
         month = sprintf("%02d", month(date)), doy = yday(date)) %>%
  filter(zip %in% clustered_zctas$GEOID, Population > 0)

monthly_ed <- read_csv("datasets/Dummy_20241005_TotalERMonthly.csv") %>%
  mutate(zip = as.character(zip), month = sprintf("%02d", as.numeric(month))) %>%
  filter(zip %in% clustered_zctas$GEOID, Population > 0)

missing_from_daily <- setdiff(unique(clustered_zctas$GEOID), unique(daily_ed$zip))
extra_in_daily <- setdiff(unique(daily_ed$zip), unique(clustered_zctas$GEOID))
length(missing_from_daily) # zctas missing from daily file
length(extra_in_daily) # exttra zips in daily file not in clusters
print(missing_from_daily)

# =========================================================================
# IMPUTE SUPPRESSED CELLS (<12), HEAT-SHAPED, BOUNDED 1-11, DAILY D2 only
# =========================================================================

set.seed(1)

impute_suppressed_no_heat <- function(raw_col) {
  suppressed <- str_detect(as.character(raw_col), "^<")
  raw_num <- suppressWarnings(as.numeric(raw_col))
  # assign suppressed cells uniformly between 1 and 11
  imputed <- sample(1:11, size = length(raw_col), replace = TRUE)
  if_else(suppressed, imputed, raw_num)}

daily_ed <- daily_ed %>%
  mutate(D2_suppressed = str_detect(as.character(D2Dx1), "^<"),
         D2 = impute_suppressed_no_heat(D2Dx1))

# heat-dependent mean: hotter days receive larger imputed counts
#expected_value <- 1 + 10 * (hi_scaled^6)
#noise <- rgamma(length(raw_col), shape = 1, rate = 0.5) - 0.5 # randomness
#imputed <- pmin(11, pmax(1, round(expected_value + noise))) # ensures 1-11 integer
#if_else(suppressed, imputed, raw_num)}

# rescales max heat index value to 0-1 (hi_scaled)
#hi_min_d <- min(daily_ed$Max_HI_Value, na.rm = TRUE)
#hi_max_d <- max(daily_ed$Max_HI_Value, na.rm = TRUE)
#daily_ed <- daily_ed %>% mutate(
#    HI_scaled = (Max_HI_Value - hi_min_d) / (hi_max_d - hi_min_d),
#    D2_suppressed = str_detect(as.character(D2Dx1), "^<"),
#    D2 = impute_suppressed(D2Dx1, HI_scaled))

# ============================================================================
# ENSURE DAILY COUNTS MATCH MONTHLY TOTALS 
# ============================================================================

monthly_ed <- monthly_ed %>% mutate(
  mo_suppressed = str_detect(as.character(D2Dx1), "^<"), # <12 is true; exact count false
  mo_upper = suppressWarnings(as.numeric(str_remove(as.character(D2Dx1), "^<"))), 
  mo_num = suppressWarnings(as.numeric(D2Dx1))) 

# final visits = baseline + rounded-down share + 1 only for the
# n_extra days with the biggest leftover decimals (others get +0)
allocate_remainder <- function(pool, share, baseline) {
  raw_alloc <- pool * share # each day's fair (decimal) share of the pool
  floor_alloc <- floor(raw_alloc) # round every day down first (under-allocates)
  remainder <- raw_alloc - floor_alloc 
  n_extra <- pool - sum(floor_alloc) # how many whole visits still need to be handed out
  rank_remainder <- rank(-remainder, ties.method = "first") # rank days by biggest leftover decimal
  baseline + floor_alloc + if_else(rank_remainder <= n_extra, 1, 0)}

daily_ed <- daily_ed %>% 
  left_join(monthly_ed %>% select(zip, month, mo_suppressed, mo_upper, mo_num),
            by = c("zip", "month")) %>%
  group_by(zip, month) %>% mutate(
    n_days = n(), n_suppressed_days = sum(D2_suppressed),
    daily_sum_raw = sum(D2, na.rm = TRUE),
    baseline = if_else(D2_suppressed, 1, 0), 
    share = if_else(daily_sum_raw > 0, D2 / daily_sum_raw, 1 / n_days),
                                  
# ===== exact monthly number given -> daily sum must match exactly =====
    target_exact = if_else(!is.na(mo_num), round(mo_num), NA_real_),
    pool_exact = if_else(!is.na(mo_num), pmax(target_exact - sum(baseline), 0), NA_real_),
    D2_exact_case = allocate_remainder(pool_exact, share, baseline),
                                  
# ===== suppressed monthly value (<12) -> daily sum stays under cap =====
    exceeds_cap = mo_suppressed & (daily_sum_raw >= mo_upper),
    pool_cap = if_else(exceeds_cap, pmax((mo_upper - 1) - sum(baseline), 0), NA_real_),
    D2_cap_case = allocate_remainder(pool_cap, share, baseline),
                                      
    D2 = case_when(!is.na(mo_num) ~ D2_exact_case,
    mo_suppressed & exceeds_cap ~ D2_cap_case,
    mo_suppressed & !exceeds_cap & sum(D2) == 0 ~ 
    if_else(row_number() == which.max(Max_HI_Value), 1, D2), 
    mo_suppressed & !exceeds_cap ~ pmax(D2, if_else(D2_suppressed, 1, D2)),
    TRUE ~ D2)) %>% ungroup() %>%
    select(-mo_suppressed, -mo_upper, -mo_num, -baseline, -share,
           -target_exact, -pool_exact, -D2_exact_case,
           -exceeds_cap, -pool_cap, -D2_cap_case)
gc()

# ============================================================
# CHECKS
# ============================================================

check_nonsuppressed <- daily_ed %>% group_by(zip, month) %>%
  summarise(daily_sum = sum(D2, na.rm = TRUE), .groups = "drop") %>%
  left_join(monthly_ed %>% select(zip, month, mo_num), by = c("zip", "month")) %>%
  filter(!is.na(mo_num)) %>% mutate(ok = daily_sum == round(mo_num))

# zip-months with a real monthly number that doesn't match (should be 0)
sum(!check_nonsuppressed$ok)

check_suppressed <- daily_ed %>% group_by(zip, month) %>%
  summarise(daily_sum = sum(D2, na.rm = TRUE), .groups = "drop") %>%
  inner_join(monthly_ed %>% filter(mo_suppressed) %>% select(zip, month, mo_upper),
             by = c("zip", "month")) %>%
  mutate(in_valid_range = daily_sum >= 1 & daily_sum < mo_upper)

print(sum(check_suppressed$daily_sum < 1)) # should be 0
print(sum(check_suppressed$daily_sum >= 12)) # should be 0

# D2 imputed value distribution
print(daily_ed %>% filter(D2_suppressed) %>%
        count(D2) %>% mutate(percent = n / sum(n)))

total_pop <- daily_ed %>% distinct(zip, Population) %>% summarise(sum(Population)) %>% pull()
rate <- sum(daily_ed$D2, na.rm = TRUE) / total_pop * 10000
round(rate, 2) # expecting ~14.9

daily_ed %>% mutate(heat_decile = ntile(Max_HI_Value, 10)) %>%
  group_by(heat_decile) %>%
  summarise(mean_HI = mean(Max_HI_Value), mean_D2 = mean(D2), n = n())

saveRDS(daily_ed, "rds/CA_daily_ed_no_heat.rds")
