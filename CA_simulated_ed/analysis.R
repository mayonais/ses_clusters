library(dplyr)
library(mgcv)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(gridExtra)

ses_clusters <- read.csv("../ed_analysis_ready/CA_zcta_clustered_for_ED.csv",
                         colClasses = c(GEOID = "character"))

ed_data <- read.csv("daily_counts_2008_2018.csv", colClasses = c(GEOID = "character"))

ses_clusters <- ses_clusters %>% 
  dplyr::select(GEOID, cluster, cluster_label)

ed_combined <- ed_data %>%
  inner_join(ses_clusters, by = "GEOID")

nrow(ed_combined)

ed_model <- ed_combined %>% 
  mutate(cluster = factor(cluster)) %>%
  filter(!is.na(cluster)) %>%
  mutate(
    date = as.Date(date),
    year = factor(year(date)),  # matches CHS's gamma_Year(d)
    day_of_week = factor(wday(date)), # matches CHS's delta_DOW(d)
    day_of_year = yday(date)
  ) %>%
  dplyr::select(hri_ed, cluster, pop, year, day_of_week, day_of_year)

nrow(ed_model)

rm(ed_data)
rm(ed_combined)
rm(ses_clusters)
gc()

# ----- choosing reference cluster --------------------------------------------------

cluster_truth <- read.csv("cluster_multiplier.csv")
reference_cluster <- cluster_truth %>%
  slice_min(vulnerability_score) %>%
  pull(cluster)

ed_model <- ed_model %>% mutate(cluster = factor(cluster),
         cluster = relevel(cluster, ref = as.character(reference_cluster)))

levels(ed_model$cluster)

# ------- models ---------------------------------------------------------

capture.output({
  cat("====================================================================\n")
  cat("POISSON GAM\n")
  cat("====================================================================\n")
  poisson <- bam( 
    hri_ed ~ cluster + year + day_of_week +
      # similar to CHS's second-order random walk f(DOY_d)
      # learns a smooth seasonal curve
      s(day_of_year, bs = "cr") + 
      offset(log(pop)), # tells the model to estimate rates instead of raw counts
    family = poisson(),
    data = ed_model,
    discrete = TRUE 
  )
  print(summary(poisson))
  cat("\nAIC:\n")
  print(AIC(poisson))
  rm(poisson)
  gc()
  
  cat("\n====================================================================\n")
  cat("NEGATIVE BINOMIAL GAM\n")
  cat("====================================================================\n")
  nb <- bam(
    hri_ed ~ cluster + year + day_of_week + s(day_of_year, bs = "cr") + 
      offset(log(pop)), family = nb(), data = ed_model, discrete = TRUE)
  print(summary(nb))
  cat("\nAIC:\n")
  print(AIC(nb))
  
  cat("\n\nZERO-INFLATION DIAGNOSTIC\n")
  cat("---------------------------------------------------------------------\n")
  observed_zero_prop <- mean(ed_model$hri_ed == 0)
  fitted_mu <- fitted(nb) # predicted mean ED count from NB model for every row
  theta_est <- nb$family$getTheta(TRUE)
  predicted_zero_prop <- mean(
    dnbinom(0, mu = fitted_mu, size = theta_est)
  ) # average probability of zero across all ZIP-days
  cat("Observed zeros:", observed_zero_prop, "\n") 
  cat("NB predicted zeros:", predicted_zero_prop, "\n")
  cat("Excess zeros:", observed_zero_prop - predicted_zero_prop, "\n") # if large, use zero-inflated NB
  
  cat("\n\nGAM CHECK\n")
  cat("---------------------------------------------------------------------\n")
  gam.check(nb, rep = 0)
  
  nb_p_table <- as.data.frame(summary(nb)$p.table)
  nb_p_table$term <- rownames(nb_p_table)
  saveRDS(nb_p_table, "nb_coef_table.rds")
  
  rm(nb)
  gc()
}, file = "model_summary.txt")
gc()
rm(fitted_mu)

# ----- extract incidence rate ratios ----------------------------------------------

colnames(nb_p_table)
nb_p_table <- readRDS("nb_coef_table.rds")

coef_table <- nb_p_table %>%
  rename(estimate  = 1, std.error = 2, p.value   = 4) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    ci_low_log = estimate - 1.96 * std.error,
    ci_high_log = estimate + 1.96 * std.error,
    IRR = exp(estimate),
    ci_low = exp(ci_low_log),
    ci_high = exp(ci_high_log),
    pct_change = (IRR - 1) * 100,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""),
    predictor = case_when(
      grepl("^cluster", term) ~ "cluster",
      grepl("^year", term) ~ "year",
      grepl("^day_of_week", term) ~ "day_of_week",
      TRUE ~ term),
    level = gsub("[^0-9]+", "", term)) %>%
  select(predictor, level, IRR, ci_low, ci_high, pct_change, p.value, sig) %>%
  mutate(across(c(IRR, ci_low, ci_high, pct_change), ~round(.x, 3)))

# attach cluster labels
clusters <- read.csv("../ed_analysis_ready/CA_zcta_clustered_for_ED.csv",
                     colClasses = c(GEOID = "character"))
cluster_labels <- clusters %>%
  distinct(cluster, cluster_label) %>%
  mutate(cluster = as.character(cluster))

coef_table <- coef_table %>%
  mutate(cluster_label = if_else(
    predictor == "cluster",
    cluster_labels$cluster_label[match(level, cluster_labels$cluster)],
    NA_character_))

# build reference row
reference_cluster <- as.character(levels(ed_model$cluster)[1])

reference_row <- data.frame(
  predictor = "cluster", level = reference_cluster,
  IRR = 1, ci_low = 1, ci_high = 1, pct_change = 0, 
  p.value = NA, sig = "ref",
  cluster_label = cluster_labels$cluster_label[match(reference_cluster, cluster_labels$cluster)])

# combine and sort
coef_table <- bind_rows(reference_row, coef_table) %>%
  mutate(
    predictor = factor(predictor, levels = c("cluster", "year", "day_of_week")),
    level_numeric = as.numeric(level)) %>%
  arrange(predictor, level_numeric) %>%
  select(-level_numeric)

print(coef_table, row.names = FALSE)
write.csv(coef_table, "incidence_rate_ratios", row.names = FALSE)

# ----- check -------------------------------------------------

reference_multiplier <- cluster_truth %>%
  filter(cluster == as.character(reference_cluster)) %>%
  pull(true_multiplier)

validation <- coef_table %>% filter(predictor == "cluster") %>%
  left_join(cluster_truth %>% mutate(cluster = as.character(cluster),
                     true_IRR = true_multiplier / reference_multiplier),
            by = c("level" = "cluster"))

validation_table <- validation %>%
  select(cluster = level, cluster_label.x, estimated_IRR = IRR, assigned_IRR = true_IRR)

irr_correlation <- cor(validation$IRR, validation$true_IRR)

irr_plot <- ggplot(validation, aes(x = true_IRR, y = IRR)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Assigned IRR from simulated ED data",
    y = "Estimated IRR from NB GAM",
    title = "Comparison of assigned and estimated IRRs (by cluster)",
    subtitle = paste0("Correlation = ", round(irr_correlation, 3))) +
  theme_minimal()

pdf("cluster_IRR_validation.pdf", width = 8, height = 8)
print(irr_plot)

grid::grid.newpage()
grid::grid.draw(
  gridExtra::tableGrob(validation_table))

dev.off()