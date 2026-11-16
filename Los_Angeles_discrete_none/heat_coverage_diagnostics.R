# Run from the project root: source("Los_Angeles_discrete_none/heat_coverage_diagnostics.R")
# Describes heat coverage and missing zone assignments without modifying data.
input_dir <- "Los_Angeles_discrete_none"
daily <- readRDS(file.path(input_dir, "Los Angeles_daily_ed_no_heat.rds"))
zones <- read.csv("datasets/cz_zip.csv", colClasses = "character")
stopifnot(!anyDuplicated(zones$zip))

study_dates <- seq(as.Date("2018-05-01"), as.Date("2018-10-31"), by = "day")
# Select the coverage records to display; this is not an analysis eligibility rule.
display_cutoff <- 0.80
daily$date <- as.Date(daily$date)
daily <- daily[daily$date %in% study_dates, ]
stopifnot(!anyNA(daily$zip), !anyDuplicated(daily[c("zip", "date")]))

coverage <- do.call(rbind, lapply(split(daily, as.character(daily$zip)), function(x) {
  # Zero is the source's unavailable-heat sentinel, as in prepare_heat_data().
  usable <- is.finite(x$Max_HI_Value) & x$Max_HI_Value != 0
  data.frame(zip = as.character(x$zip[1]),
             usable_days = length(unique(x$date[usable])))
}))
rownames(coverage) <- NULL
coverage$missing_days <- length(study_dates) - coverage$usable_days
coverage$coverage_pct <- 100 * coverage$usable_days / length(study_dates)
coverage$climate_zone <- zones$cz[match(coverage$zip, zones$zip)]
coverage <- coverage[order(coverage$usable_days, coverage$zip), ]
displayed <- coverage[coverage$coverage_pct < 100 * display_cutoff, ]
missing_zone <- coverage[is.na(coverage$climate_zone) |
                           !nzchar(trimws(coverage$climate_zone)), ]
missing_zone <- missing_zone[order(missing_zone$zip), ]

format_table <- function(x) {
  x$coverage_pct <- sprintf("%.1f%%", x$coverage_pct)
  capture.output(print(x, row.names = FALSE, right = FALSE))
}
columns <- c("zip", "usable_days", "missing_days", "coverage_pct")
report <- c(
  "Heat coverage and missing climate-zone diagnostics",
  "Input: Los_Angeles_discrete_none/Los Angeles_daily_ed_no_heat.rds",
  "Climate-zone lookup: datasets/cz_zip.csv",
  sprintf("Study period: %s through %s (%d expected days per ZIP)",
          min(study_dates), max(study_dates), length(study_dates)),
  "Usable heat: finite Max_HI_Value other than zero; zero denotes unavailable heat.",
  sprintf("Total ZIPs: %d", nrow(coverage)),
  "", sprintf("Selected ZIP heat-coverage records: %d", nrow(displayed)),
  format_table(displayed[columns]),
  "", sprintf("ZIPs missing climate-zone assignments: %d",
              nrow(missing_zone)),
  format_table(missing_zone[columns]),
  "", "Catalina (90704) heat coverage:",
  format_table(coverage[coverage$zip == "90704", columns]),
  "", "Descriptive report only; no analysis inputs were modified."
)
output_path <- file.path(input_dir, "heat_coverage_diagnostics.txt")
writeLines(report, output_path)
cat(paste(report, collapse = "\n"), "\n")
message("Report written to: ", output_path)
