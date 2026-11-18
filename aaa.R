
################################################################################
# 0) Setup
################################################################################
library(plm)
library(dplyr)
library(modelsummary)
library(sandwich)

# SET YOUR PATH HERE
base_path <- "C:\\Users\\nico_\\Dropbox\\Perez&Alvayay2025\\Data\\Output\\panel and spatial"
#base_path <- "C:/Users/Usuario/Dropbox/Research/Perez&Alvayay2025"

data_path <- file.path(base_path, "Data")
output_path <- file.path(data_path, "Output/spatial")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1) Load data
################################################################################
cat("=== Loading spatial panel ===\n")
#sf_data <- readRDS(file.path(data_path, "Output/spatial_panel.rds"))
sf_data <- st_read("C:\\Users\\nico_\\Dropbox\\Perez&Alvayay2025\\Data\\Output\\panel and spatial\\dtSf.geojson")
sf_data <- sf::st_drop_geometry(sf_data)

################################################################################
# 2) Helper function: 3-ring DID tables
################################################################################
run_did_table_3rings <- function(outcome, file_name) {
  
  cat("\n--- Running 3-ring spatial DID for:", outcome, "---\n")
  
  fml <- as.formula(paste0(
    outcome,
    " ~ D_T + WD_T + W2D_T + W3D_T + black + povline + female + median +
      white + owner_hper + vacant_rate + density + demo_rate + amenity"
  ))
  
  # Full sample
  m_all <- plm(fml, data = sf_data,
               index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  # Tercile 1
  m_t1 <- plm(fml, data = sf_data |> filter(demo_tercile == 1),
              index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  m_t2 <- plm(fml, data = sf_data |> filter(demo_tercile == 2),
              index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  # Tercile 3
  m_t3 <- plm(fml, data = sf_data |> filter(demo_tercile == 3),
              index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  # Q1 (comparison)
  m_q1 <- plm(fml, data = sf_data |> filter(demo_qbgr == 1),
              index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  # Q4 (comparison)
  m_q4 <- plm(fml, data = sf_data |> filter(demo_qbgr == 4),
              index = c("GEOID", "year"), model = "within", effect = "twoways")
  
  models <- list(
    "(1) Full"    = m_all,
    "(2) T1-Low"  = m_t1,
    "(3) T2-Med"  = m_t2,
    "(4) T3-High" = m_t3
  )
  
  coef_map <- c(
    "D_T"   = "ATT (Direct)",
    "WD_T"  = "AIT (1st order)",
    "W2D_T" = "AIT (2nd order)",
    "W3D_T" = "AIT (3rd order)"
  )
  
  vcov_list <- lapply(models, function(m) vcovHC(m, type = "HC1", cluster = "group"))
  
  modelsummary(
    models,
    coef_map = coef_map,
    vcov = vcov_list,
    fmt = 3,
    estimate  = "{estimate}{stars}",
    statistic = "({std.error})",
    stars = c("*" = .10, "**" = .05, "***" = .01),
    gof_map = c("nobs", "r.squared", "adj.r.squared"),
    output = file.path(output_path, file_name)
  )
  
  cat("  Saved:", file_name, "\n")
  
  # Print key results
  cat("\n  Key coefficients:\n")
  for (nm in names(models)) {
    b <- coef(models[[nm]])
    se <- sqrt(diag(vcov_list[[nm]]))
    cat(sprintf("  %s: ATT=%.3f (%.3f), AIT1=%.3f (%.3f), AIT2=%.3f (%.3f), AIT3=%.3f (%.3f)\n",
                nm, b["D_T"], se["D_T"], b["WD_T"], se["WD_T"],
                b["W2D_T"], se["W2D_T"], b["W3D_T"], se["W3D_T"]))
  }
}

################################################################################
# 3) Run tables
################################################################################

run_did_table_3rings(
  outcome   = "total_crime",
  file_name = "spatial_3rings_total_crime2.rtf"
)

run_did_table_3rings(
  outcome   = "property",
  file_name = "spatial_3rings_property2.rtf"
)

run_did_table_3rings(
  outcome   = "violent",
  file_name = "spatial_3rings_violent2.rtf"
)

run_did_table_3rings(
  outcome   = "log_light",
  file_name = "spatial_3rings_log_light2.rtf"
)

cat("\n=== 14_spatial_tables_3rings.R COMPLETE ===\n")
