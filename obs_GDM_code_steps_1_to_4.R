# ==============================================================================
# STEP 1. GDM FITTING
# ==============================================================================
# 1.1 Load packages, source functions, and input data
# ------------------------------------------------------------------------------
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(raster)
library(gdm)
library(nnls)
library(SDMTools)
library(devtools)
library(doSNOW)
library(foreach)
library(pbapply)
library(tibble)
library(doParallel)
library(Matrix)

source_dir <- "E:/BILBI - example dataset and code/ObsPair_GDM_v2-main/ObsPair_GDM_v2-main/src"
# Load the required R script files in sequence
source(file.path(source_dir, "utils.R"))
source(file.path(source_dir, "gdm_functions.R"))
source(file.path(source_dir, "site_aggregator.R"))
source(file.path(source_dir, "site_richness_extractor.R"))
source(file.path(source_dir, "obs_pair_sampler.R"))
source(file.path(source_dir, "gen_windows.R"))
source(file.path(source_dir, "run_chunked_env.R"))
source(file.path(source_dir, "plotting.R"))
source(file.path(source_dir, "predict_spatial.R"))
# Modified from obs_pair_sampler.R: restricts sampling to within 256 km and preserves the natural match:mismatch ratio in each distance bin, rather than forcing a 1:1 ratio.
source(file.path(source_dir, "obs_pair_sampler_stratified_odds.R"))
source(file.path(source_dir, "config.R")) # config$nMatch <- as.integer(env_or_default("RECA_NMATCH", "100000"))
load_packages()



# Observation records with environmental covariates for all sampled cells/sites.
# 0 = forest-associated species; 1 = forest-dependent species.
combined_data <- readRDS("E:/Weicheng_Code_5/observation_data/NT_TSMF_combined.rds")%>% 
  filter(Forest_Dependency == 0)

#The main purpose of this section (through line 132) is to organize and clean citizen-science observation records into a structure suitable for subsequent obs-GDM analysis.
# 1.2 Standardize observation data for Andrew's spatiotemporal workflow
# Temporal fields are added only to satisfy the required input structure.
dat <- combined_data %>%
  rename(
    decimalLongitude = x,
    decimalLatitude = y,
    scientificName = species
  ) %>%
  mutate(
    eventDate = "2000-01-01",
    coordinateUncertaintyInMeters = 1
  )
# Spatial resolution (~1 km)
res    <- 0.009
# Analysis extent
box <- extent(min(combined_data$x), max(combined_data$x), min(combined_data$y), max(combined_data$y))
# Aggregate records to grid cells
datRED <- siteAggregator(dat, res, box)


# w is the mismatch-to-match ratio from a random sample of observation-record pairs.
# It corrects the artificial 1:1 match:mismatch fitting design before predicted
# mismatch probabilities are converted to community-level compositional dissimilarity.
#source(file.path(source_dir, "utils.R"))
cat("Estimating mismatch/match ratio (w)...\n")
w <- estimate_w(datRED, nSamples = config$w_estimation_samples)


# 1.3 Prepare site-by-species data and sampling objects
cat("\n--- Step 1.3: Prepare data for sampling ---\n")
# Prepare site × species data
data <- data.frame(
  ID        = datRED$ID,
  Latitude  = datRED$latID,
  Longitude = datRED$lonID,
  species   = datRED$gen_spec,
  nRecords  = datRED$nRecords,
  nRecords.exDateLocDups = datRED$nRecords.exDateLocDups,
  nSiteVisits = datRED$nSiteVisits,
  richness  = datRED$richness,
  stringsAsFactors = FALSE
)

# # Count records per unique site × species for m1 expansion
# LocDups <- as.factor(paste(datRED$ID, datRED$gen_spec, sep = ":"))
# ones    <- rep(1, nrow(datRED))
# count   <- bySum(ones, LocDups)

# 1.4 Build the site-by-species matrix (observed richness)
cat("\n--- Step 1.4: Build site × species matrix ---\n")
cl <- makeCluster(config$cores_to_use)
registerDoSNOW(cl)
frog.auGrid <- site.richness.extractor.bigData(frog.auGrid = data)
stopCluster(cl)

# # The proxy code below requires these objects. Some source versions create them
# # internally; construct them here only when they are absent.
# if (!exists("site_sp_matrix", inherits = FALSE)) {
#   site_levels <- unique(as.character(data$ID))
#   species_levels <- unique(as.character(data$species))
#   site_sp_matrix <- Matrix::sparseMatrix(i = match(as.character(data$ID), site_levels), j = match(as.character(data$species), species_levels), x = 1, dims = c(length(site_levels), length(species_levels)), dimnames = list(site_levels, species_levels))
#   site_sp_matrix@x[] <- 1
# }
# if (!exists("site_levels", inherits = FALSE)) site_levels <- rownames(site_sp_matrix)

# Restore full dataset
frog.auGrid <- data.frame(
  ID        = datRED$ID,
  Latitude  = datRED$latID,
  Longitude = datRED$lonID,
  species   = datRED$gen_spec,
  eventDate = as.character(datRED$eventDate),
  nRecords  = datRED$nRecords,
  nRecords.exDateLocDups = datRED$nRecords.exDateLocDups,
  nSiteVisits = datRED$nSiteVisits,
  richness  = datRED$richness,
  stringsAsFactors = FALSE,
  Site.Richness = datRED$richness
)

# Build a full-record proxy mapping each frog.auGrid row to its site.
# This replaces m2 = m1[row.vect, ], which caused out-of-memory errors.
full_site_map <- match(as.character(frog.auGrid$ID), site_levels)
m2 <- site_species_proxy(site_sp_matrix, full_site_map)
if (exists("m1", inherits = FALSE)) rm(m1)
gc()


# 1.5 Sample observation pairs in six distance bins within 256 km
# Aim for a similar total number of pairs per bin, but do not force a 1:1
# match:mismatch ratio. Mismatches naturally become more common with distance;
# preserving this pattern lets the model estimate distance decay.
# Keep the relative sampling odds of matches and mismatches the same across bins.
# nMatch sets candidates per iteration; target_nMatch sets the total-pair target
# per bin; max_iterations limits sampling attempts.

result <- obsPairSampler.distanceStratified_2(
  frog.auGrid      = frog.auGrid,
  m1               = m2,
  nMatch           = config$nMatch,
  target_nMatch    = 20000,
  richness         = TRUE,
  speciesThreshold = 500,
  max_iterations   = 50,
  seed             = 123)



obsPairs_out <- result$pairs


# 1.6 Attach coordinates and environmental covariates to sampled pairs
coord_map <- unique(data.frame(
  lonID = datRED$lonID,
  latID = datRED$latID,
  RAW_lon = datRED$RAW_longdec,
  RAW_lat = datRED$RAW_latdec))

# Match original coordinates for Site 1
match1 <- match(paste(obsPairs_out$Lon1, obsPairs_out$Lat1), 
                paste(coord_map$lonID, coord_map$latID))
obsPairs_out$RAW_Lon1 <- coord_map$RAW_lon[match1]
obsPairs_out$RAW_Lat1 <- coord_map$RAW_lat[match1]

# Match original coordinates for Site 2
match2 <- match(paste(obsPairs_out$Lon2, obsPairs_out$Lat2), 
                paste(coord_map$lonID, coord_map$latID))
obsPairs_out$RAW_Lon2 <- coord_map$RAW_lon[match2]
obsPairs_out$RAW_Lat2 <- coord_map$RAW_lat[match2]

# Match environmental covariates directly from dat
env_vars <- c("BD","CLAY","PH","SILT","TWI","PTA","EAAS","WDI","WDX","EPI","TXX","TNI","RUG","TRX","ELEV","FLII")

# Match Site 1 environmental covariates
match_env1 <- match(paste(round(obsPairs_out$RAW_Lon1, 8), round(obsPairs_out$RAW_Lat1, 8)),
                    paste(round(dat$decimalLongitude, 8), round(dat$decimalLatitude, 8)))

for(var in env_vars) {
  obsPairs_out[[var]] <- dat[[var]][match_env1]
}

# Match Site 2 environmental covariates (_2 suffix)
match_env2 <- match(paste(round(obsPairs_out$RAW_Lon2, 8), round(obsPairs_out$RAW_Lat2, 8)),
                    paste(round(dat$decimalLongitude, 8), round(dat$decimalLatitude, 8)))

for(var in env_vars) {
  obsPairs_out[[paste0(var, "_2")]] <- dat[[var]][match_env2]
}

# 1.7 Construct the GDM training table
mod_data <- data.frame(
  Match = obsPairs_out$Match,
  Lon1 = obsPairs_out$RAW_Lon1,
  Lat1 = obsPairs_out$RAW_Lat1,
  Lon2 = obsPairs_out$RAW_Lon2,
  Lat2 = obsPairs_out$RAW_Lat2
)

# Add environmental covariates
for(var in env_vars) {
  mod_data[[var]] <- obsPairs_out[[var]]
  mod_data[[paste0(var, "_2")]] <- obsPairs_out[[paste0(var, "_2")]]
}


# Add remaining response and richness variables
mod_data$Richness.S1 <- obsPairs_out$Richness.S1
mod_data$Richness.S2 <- obsPairs_out$Richness.S2
# mod_data$SORENSON <- obsPairs_out$Sorenson
# mod_data$TARGET <- obsPairs_out$TARGET



target_order <- c(
  "Match", "Lon1", "Lat1", "Lon2", "Lat2",
  "BD", "CLAY", "PH", "SILT", "TWI", "PTA", "EAAS", "WDI", "WDX", 
  "EPI", "TXX", "TNI", "RUG", "TRX", "ELEV", "FLII",
  "BD_2", "CLAY_2", "PH_2", "SILT_2", "TWI_2", "PTA_2", "EAAS_2", 
  "WDI_2", "WDX_2", "EPI_2", "TXX_2", "TNI_2", "RUG_2", "TRX_2", "ELEV_2", "FLII_2",
  "Richness.S1", "Richness.S2"
)

existing_cols <- target_order[target_order %in% names(mod_data)]


# Reorder columns for GDM fitting
mod_data <- mod_data[, existing_cols]
print(head(mod_data))

# ==============================================================================
# 1.8 Fit the two-stage obs-pair GDM
# Stage 1: environmental predictors, including ELEV and FLII
# Stage 2: geographic distance with Stage-1 linear predictor as an offset
# ==============================================================================

gdm_output_path <- "E:/Weicheng_Code_5/GDM_fitting"
dir.create(gdm_output_path, showWarnings = FALSE, recursive = TRUE)

env_vars <- c("BD", "CLAY", "PH", "SILT", "TWI", "PTA", "EAAS", "WDI", "WDX", "EPI", "TXX", "TNI", "RUG", "TRX", "ELEV", "FLII")

# Calculate geographic distance in metres, following Andrew's original code
mod_data$geoDist <- geosphere::distGeo(
  mod_data[, c("Lon1", "Lat1")],
  mod_data[, c("Lon2", "Lat2")]
)

required_columns <- c("Match", env_vars, paste0(env_vars, "_2"), "geoDist")
mod_data_env <- mod_data[complete.cases(mod_data[, required_columns]), required_columns, drop = FALSE]

cat("Number of pairs used for fitting:", nrow(mod_data_env), "\n")
cat("Match pairs:", sum(mod_data_env$Match == 0), "\n")
cat("Mismatch pairs:", sum(mod_data_env$Match == 1), "\n")
cat("Geographic-distance range (km):", paste(round(range(mod_data_env$geoDist) / 1000, 3), collapse = " to "), "\n")


# ------------------------------------------------------------------------------
# 1.8.1 Transform environmental predictors using I-splines
# ------------------------------------------------------------------------------

toSpline <- mod_data_env[, c(env_vars, paste0(env_vars, "_2")), drop = FALSE]
splined_env <- splineData(toSpline)

mod_ready_env <- data.frame(Match = mod_data_env$Match, as.data.frame(splined_env), check.names = FALSE)

cat("Environmental spline columns:", ncol(mod_ready_env) - 1L, "\n")


# ------------------------------------------------------------------------------
# 1.8.2 Transform geographic distance using I-splines
# ------------------------------------------------------------------------------

geo_quantiles <- as.numeric(quantile(mod_data_env$geoDist, probs = c(0, 0.5, 1), na.rm = TRUE, names = FALSE))

if (length(unique(geo_quantiles)) < 3L) {
  stop("Geographic-distance quantiles are not unique; geo I-splines cannot be fitted.")
}

geo_splines <- data.frame(
  geoDist_spl1 = I_spline(mod_data_env$geoDist, geo_quantiles[1], geo_quantiles[1], geo_quantiles[2]),
  geoDist_spl2 = I_spline(mod_data_env$geoDist, geo_quantiles[1], geo_quantiles[2], geo_quantiles[3]),
  geoDist_spl3 = I_spline(mod_data_env$geoDist, geo_quantiles[2], geo_quantiles[3], geo_quantiles[3])
)

mod_ready_geo <- data.frame(Match = mod_data_env$Match, geo_splines, check.names = FALSE)

cat("Geographic-distance knots (km):", paste(round(geo_quantiles / 1000, 3), collapse = ", "), "\n")


# ------------------------------------------------------------------------------
# 1.8.3 Fit Stage 1: environmental predictors
# ------------------------------------------------------------------------------

formula_env <- reformulate(names(mod_ready_env)[-1], response = "Match")

cat("\nStage 1: fitting environmental predictors...\n")
# Do not use case weights in fitGDM() for this binary match/mismatch analysis.
# Case weights change how much each pair contributes to fitting and may alter
# the intercept and I-spline response curves. Their effect in this formulation
# has not yet been validated, so retain the explicit pair-sampling design.
#
# The w used later in ObsTrans(p0, w, p) is different: it is a post-fit
# correction used to convert predicted mismatch probabilities into calibrated
# dissimilarity. It does not change the fitted model. Estimate w from the
# sampling design, not from the overall ratio in the final distance-balanced sample
obsGDM_env <- fitGDM(
  formula = formula_env,
  data = mod_ready_env,
  weights = NULL,
  lambda = config$ridge_lambda)

# summary(obsGDM_env)
stage1_eta <- as.numeric(obsGDM_env$linear.predictors)
stage1_probability <- as.numeric(fitted(obsGDM_env))
stage1_D2 <- (obsGDM_env$null.deviance - obsGDM_env$deviance) / obsGDM_env$null.deviance
stage1_R2 <- RsqGLM(model = obsGDM_env)

cat("Stage-1 intercept:", round(coef(obsGDM_env)[1], 6), "\n")
cat("Stage-1 deviance explained:", round(stage1_D2, 4), "\n")
cat("Stage-1 Nagelkerke R2:", round(stage1_R2$Nagelkerke, 4), "\n")


# ------------------------------------------------------------------------------
# 1.8.4 Fit Stage 2: geographic distance with Stage-1 eta as an offset
# ------------------------------------------------------------------------------

formula_geo <- reformulate(names(mod_ready_geo)[-1], response = "Match")

cat("\nStage 2: fitting geographic distance with Stage-1 eta as offset...\n")
obsGDM_geo <- fitGDM(
  formula = formula_geo,
  data = mod_ready_geo,
  weights = NULL,
  lambda = config$ridge_lambda,
  offset = stage1_eta
)
# summary(obsGDM_geo)
full_eta <- as.numeric(obsGDM_geo$linear.predictors)
full_probability <- as.numeric(fitted(obsGDM_geo))

cat("Stage-2 intercept:", round(coef(obsGDM_geo)[1], 6), "\n")


# ------------------------------------------------------------------------------
# 1.8.5 Combine Stage-1 and Stage-2 parameters
# ------------------------------------------------------------------------------

env_coefficients <- coef(obsGDM_env)[-1]
geo_coefficients <- coef(obsGDM_geo)[-1]

env_coefficients[is.na(env_coefficients)] <- 0
geo_coefficients[is.na(geo_coefficients)] <- 0

combined_coefficients <- c(env_coefficients, geo_coefficients)
combined_intercept <- as.numeric(coef(obsGDM_env)[1] + coef(obsGDM_geo)[1])

cat("Combined intercept:", round(combined_intercept, 6), "\n")
cat("Environmental coefficient sum:", round(sum(env_coefficients), 6), "\n")
cat("Geographic coefficient sum:", round(sum(geo_coefficients), 6), "\n")


# ------------------------------------------------------------------------------
# 1.8.6 Store environmental and geographic knots
# ------------------------------------------------------------------------------

env_quantiles <- unlist(lapply(env_vars, function(var) {
  quantile(c(mod_data_env[[var]], mod_data_env[[paste0(var, "_2")]]), probs = c(0, 0.5, 1), na.rm = TRUE, names = FALSE)
}), use.names = FALSE)

combined_quantiles <- c(env_quantiles, geo_quantiles)


# ------------------------------------------------------------------------------
# 1.8.7 Build the combined prediction-ready fit object
# ------------------------------------------------------------------------------

fit <- list(
  intercept = combined_intercept,
  predictors = c(env_vars, "geoDist"),
  splines = rep(3L, length(env_vars) + 1L),
  coefficients = as.numeric(combined_coefficients),
  quantiles = as.numeric(combined_quantiles),
  predicted = full_probability,
  ecological = full_eta - combined_intercept,
  sample = nrow(mod_ready_env),
  species_group = "NT_TSMF",
  n_pairs = nrow(mod_data_env),
  w_ratio = w,
  biAverage = FALSE,
  fit_mode = "two_stage_offset_geo",
  add_geo_distance = TRUE,
  offset_geo_after_env = TRUE,
  geo_distance_units = "metres",
  stage1_intercept = as.numeric(coef(obsGDM_env)[1]),
  stage2_intercept = as.numeric(coef(obsGDM_geo)[1]),
  stage1_predictors = env_vars,
  stage1_coefficients = as.numeric(env_coefficients),
  stage1_quantiles = as.numeric(env_quantiles),
  stage1_splines = rep(3L, length(env_vars)),
  stage1_predicted = stage1_probability,
  stage1_D2 = stage1_D2,
  stage1_nagelkerke_R2 = stage1_R2$Nagelkerke,
  geo_quantiles = geo_quantiles,
  ridge_lambda = config$ridge_lambda,
  case_weights_applied = FALSE,
  run_timestamp = Sys.time()
)


# ------------------------------------------------------------------------------
# 1.8.8 Validate the fit object and combined linear predictor
# ------------------------------------------------------------------------------

stopifnot(
  length(fit$predictors) == length(fit$splines),
  length(fit$coefficients) == sum(fit$splines),
  length(fit$quantiles) == sum(fit$splines),
  length(fit$stage1_predictors) == length(fit$stage1_splines),
  length(fit$stage1_coefficients) == sum(fit$stage1_splines),
  length(fit$stage1_quantiles) == sum(fit$stage1_splines)
)

reconstructed_env <- as.numeric(as.matrix(splined_env) %*% env_coefficients)
reconstructed_geo <- as.numeric(as.matrix(geo_splines) %*% geo_coefficients)
reconstructed_eta <- combined_intercept + reconstructed_env + reconstructed_geo
max_eta_difference <- max(abs(reconstructed_eta - full_eta), na.rm = TRUE)

cat("Maximum reconstructed-eta difference:",
    format(max_eta_difference, scientific = TRUE), "\n")

if (max_eta_difference > 1e-6) {
  warning("The combined fit object does not exactly reproduce the Stage-2 linear predictor.")
}


# ------------------------------------------------------------------------------
# 1.8.9 Create parameter and model-summary tables
# ------------------------------------------------------------------------------

fit_parameters <- data.frame(
  predictor = rep(fit$predictors, fit$splines),
  spline = unlist(lapply(fit$splines, seq_len)),
  quantile = fit$quantiles,
  coefficient = fit$coefficients
)

fit_summary <- data.frame(
  fit_mode = fit$fit_mode,
  stage1_intercept = fit$stage1_intercept,
  stage2_intercept = fit$stage2_intercept,
  combined_intercept = fit$intercept,
  stage1_D2 = fit$stage1_D2,
  stage1_nagelkerke_R2 = fit$stage1_nagelkerke_R2,
  environmental_coefficient_sum = sum(env_coefficients),
  geographic_coefficient_sum = sum(geo_coefficients),
  w_ratio = fit$w_ratio,
  n_pairs = fit$n_pairs,
  maximum_eta_reconstruction_error = max_eta_difference,
  run_timestamp = as.character(fit$run_timestamp)
)

print(fit_parameters)
print(fit_summary)


# ------------------------------------------------------------------------------
# 1.8.10 Save fitting objects and summaries
# ------------------------------------------------------------------------------

saveRDS(fit, file.path(gdm_output_path, "GDM_fit_NT_TSMF_FA_two_stage.rds"))
saveRDS(obsGDM_env, file.path(gdm_output_path, "obsGDM_env_NT_TSMF_FA.rds"))
saveRDS(obsGDM_geo, file.path(gdm_output_path, "obsGDM_geo_NT_TSMF_FA.rds"))
saveRDS(mod_data_env, file.path(gdm_output_path, "mod_data_env_NT_TSMF_FA.rds"))
saveRDS(mod_ready_env, file.path(gdm_output_path, "mod_ready_env_NT_TSMF_FA.rds"))
saveRDS(mod_ready_geo, file.path(gdm_output_path, "mod_ready_geo_NT_TSMF_FA.rds"))
saveRDS(result, file.path(gdm_output_path, "sampling_result_NT_TSMF_FA.rds"))

write.csv(fit_parameters, file.path(gdm_output_path, "GDM_fit_parameters_NT_TSMF_FA_two_stage.csv"), row.names = FALSE)
write.csv(fit_summary, file.path(gdm_output_path, "GDM_fit_summary_NT_TSMF_FA_two_stage.csv"), row.names = FALSE)



# ==============================================================================
# STEP 2. ENVIRONMENTAL MATCHING OF INTACT AND DEGRADED FOREST PIXELS
# One degraded pixel is selected per distance window for each intact pixel.
# ==============================================================================

library(data.table)
library(geosphere)
library(MatchIt)

input_file <- "E:/Jiahua_ForWeicheng1/data_forMATCHING/NT_TSMF_combined_data.rds"
output_dir <- "E:/Weicheng_Code_5/GDM_fitting"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)
max_intact_per_band <- 1000L
max_complete_intact_per_band <- 500L
correlation_cutoff <- 0.70
chunk_size <- 50000L
prescreen_radius_km <- 256
max_elevation_difference_m <- 200

# These are the narrow windows used in the previous script.
distance_windows <- list("0.5-8km" = list(c(0.5, 5), c(7, 9)), "15-17km" = list(c(15, 17)), "31-33km" = list(c(31, 33)), "63-65km" = list(c(63, 65)), "127-129km" = list(c(127, 129)), "255-256km" = list(c(255, 256)))
distance_labels <- c("0.5-8km" = 8, "15-17km" = 16, "31-33km" = 32, "63-65km" = 64, "127-129km" = 128, "255-256km" = 256)
ecological_priority <- c("ELEV", "PTA", "TXX", "TNI", "WDI", "WDX", "EPI", "EAAS", "TWI", "RUG", "TRX", "PH", "CLAY", "SILT", "BD")

# ------------------------------------------------------------------------------
# 2.1 Load and prepare forest pixels
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(input_file))
required <- c("cell", "lon", "lat", "FLII", ecological_priority)
missing_columns <- setdiff(required, names(dt))
if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

dt <- unique(dt, by = "cell")
dt <- dt[is.finite(FLII) & FLII >= 0 & FLII <= 10 & is.finite(lon) & is.finite(lat) & complete.cases(dt[, ..ecological_priority])]
if (nrow(dt) < 2L) stop("At least two complete forest pixels are required.")

# Divide this region's observed absolute-elevation range into three equal-width bands.
elevation_range <- range(dt$ELEV)
if (!all(is.finite(elevation_range)) || diff(elevation_range) <= 0) stop("ELEV must contain a non-zero finite range to define three elevation bands.")
elevation_breaks <- seq(elevation_range[1], elevation_range[2], length.out = 4L)
elevation_band_definition <- data.table(elev_band = c("low", "middle", "high"), lower_m = elevation_breaks[1:3], upper_m = elevation_breaks[2:4])
dt[, elev_band := cut(ELEV, breaks = elevation_breaks, labels = c("low", "middle", "high"), include.lowest = TRUE, right = TRUE)]
cat(sprintf("Regional elevation range: %.2f to %.2f m; equal-width breaks: %.2f and %.2f m\n", elevation_range[1], elevation_range[2], elevation_breaks[2], elevation_breaks[3]))

intact <- dt[FLII >= 9.6]
degraded <- dt[FLII < 9.6]

if (!nrow(intact) || !nrow(degraded)) stop("At least one integrity class is empty.")

# ------------------------------------------------------------------------------
# 2.2 Screen correlated matching variables
# ------------------------------------------------------------------------------

p <- length(ecological_priority)
n_cov <- 0
sum_x <- numeric(p)
sum_xx <- matrix(0, p, p)

for (first in seq.int(1L, nrow(dt), by = chunk_size)) {
  last <- min(first + chunk_size - 1L, nrow(dt))
  x <- as.matrix(dt[first:last, ..ecological_priority])
  storage.mode(x) <- "double"
  
  n_cov <- n_cov + nrow(x)
  sum_x <- sum_x + colSums(x)
  sum_xx <- sum_xx + crossprod(x)
}

cov_all <- (sum_xx - tcrossprod(sum_x) / n_cov) / (n_cov - 1)
dimnames(cov_all) <- list(ecological_priority, ecological_priority)

sd_all <- sqrt(diag(cov_all))
valid <- ecological_priority[is.finite(sd_all) & sd_all > 0]

if (!length(valid)) {
  stop("All matching covariates are constant.")
}

cor_all <- cov_all[valid, valid, drop = FALSE] /
  outer(sd_all[valid], sd_all[valid])
diag(cor_all) <- 1

retained <- character()
selection_log <- vector("list", length(ecological_priority))

for (i in seq_along(ecological_priority)) {
  v <- ecological_priority[i]
  
  conflicts <- if (v %in% valid && length(retained)) {
    retained[abs(cor_all[v, retained]) > correlation_cutoff]
  } else {
    character()
  }
  
  accepted <- v %in% valid && !length(conflicts)
  
  reason <- if (!(v %in% valid)) {
    "zero variance"
  } else if (length(conflicts)) {
    paste("correlated with", paste(conflicts, collapse = ", "))
  } else {
    "retained"
  }
  
  selection_log[[i]] <- data.table(
    variable = v,
    retained = accepted,
    reason = reason
  )
  
  if (accepted) retained <- c(retained, v)
}

selection_log <- rbindlist(selection_log)

# Check global covariance for linear dependence. MatchIt calculates its own
# Mahalanobis distances later from each anchor's eligible candidate pool.
repeat {
  sigma <- cov_all[retained, retained, drop = FALSE]
  chol_sigma <- tryCatch(chol(sigma), error = function(e) NULL)
  
  if (!is.null(chol_sigma)) break
  if (length(retained) == 1L) {
    stop("The matching covariance is singular.")
  }
  
  dropped <- tail(retained, 1L)
  retained <- head(retained, -1L)
  
  selection_log[
    variable == dropped,
    `:=`(
      retained = FALSE,
      reason = "singular covariance after correlation screening"
    )
  ]
}

sigma_inverse <- chol2inv(chol_sigma)

if (!("ELEV" %in% retained)) {
  stop("ELEV must be retained among the matching variables.")
}

# For example, if retained contains ELEV, PTA and WDI, this becomes:
# treat ~ ELEV + PTA + WDI
matching_formula <- reformulate(retained, response = "treat")

cat("MatchIt formula:\n")
print(matching_formula)

# ------------------------------------------------------------------------------
# 2.3 Randomly select up to 1,000 intact pixels per elevation band
# ------------------------------------------------------------------------------

anchor_list <- lapply(levels(dt$elev_band), function(band) {
  pool <- intact[elev_band == band]
  if (!nrow(pool)) return(pool)
  pool[sample.int(.N, min(.N, max_intact_per_band))]
})

anchors <- rbindlist(anchor_list, use.names = TRUE)

cat("Randomly selected intact pixels by elevation band:\n")
print(anchors[, .(selected_intact = .N), by = elev_band])

degraded_by_band <- lapply(levels(dt$elev_band), function(band) {
  pool <- copy(degraded[elev_band == band])
  setorder(pool, lat)
  pool
})

names(degraded_by_band) <- levels(dt$elev_band)

window_names <- names(distance_windows)
result_list <- vector("list", nrow(anchors) * length(window_names))

no_candidate <- matrix(
  0L,
  nrow = length(levels(dt$elev_band)),
  ncol = length(window_names),
  dimnames = list(levels(dt$elev_band), window_names)
)

mark_missing <- function(band, labels = window_names) {
  no_candidate[band, labels] <<- no_candidate[band, labels] + 1L
}

# ------------------------------------------------------------------------------
# 2.4 Match one degraded pixel in each of the six distance windows
# ------------------------------------------------------------------------------

for (i in seq_len(nrow(anchors))) {
  a <- anchors[i]
  band <- as.character(a$elev_band)
  pool <- degraded_by_band[[band]]
  
  if (!nrow(pool)) {
    mark_missing(band)
    next
  }
  
  # Initial spatial screening.
  lat_delta <- prescreen_radius_km / 110.574
  lo <- findInterval(max(-90, a$lat - lat_delta), pool$lat) + 1L
  hi <- findInterval(min(90, a$lat + lat_delta), pool$lat)
  
  if (lo > hi || lo > nrow(pool)) {
    mark_missing(band)
    next
  }
  
  candidates <- pool[lo:hi]
  
  max_abs_lat <- min(
    89.999,
    max(abs(c(a$lat - lat_delta, a$lat + lat_delta)))
  )
  
  lon_km_per_degree <- 111.320 * cos(max_abs_lat * pi / 180)
  
  lon_delta <- if (lon_km_per_degree <= 0.001) {
    180
  } else {
    min(180, prescreen_radius_km / lon_km_per_degree)
  }
  
  candidates <- candidates[
    abs(((lon - a$lon + 180) %% 360) - 180) <= lon_delta
  ]
  
  if (!nrow(candidates)) {
    mark_missing(band)
    next
  }
  
  # Calculate exact geographic distance in metres.
  geo_m <- geosphere::distGeo(
    cbind(
      rep(a$lon, nrow(candidates)),
      rep(a$lat, nrow(candidates))
    ),
    as.matrix(candidates[, .(lon, lat)])
  )
  
  groups <- rep(NA_character_, length(geo_m))
  
  for (label in window_names) {
    for (interval in distance_windows[[label]]) {
      inside <- geo_m >= interval[1] * 1000 &
        geo_m <= interval[2] * 1000
      groups[inside] <- label
    }
  }
  
  # This is the binding 200-metre elevation restriction.
  # MatchIt never receives candidates that exceed it.
  keep <- which(
    !is.na(groups) &
      geo_m <= prescreen_radius_km * 1000 &
      abs(candidates$ELEV - a$ELEV) <=
      max_elevation_difference_m
  )
  
  if (!length(keep)) {
    mark_missing(band)
    next
  }
  
  candidates <- candidates[keep]
  geo_m <- geo_m[keep]
  groups <- groups[keep]
  
  missing_windows <- setdiff(window_names, unique(groups))
  
  if (length(missing_windows)) {
    mark_missing(band, missing_windows)
    next
  }
  
  # Create six copies of the intact pixel, one for each distance window.
  treated_names <- paste0("intact_", seq_along(window_names))
  control_names <- paste0("degraded_", seq_len(nrow(candidates)))
  
  match_input <- as.data.frame(rbindlist(list(
    a[rep(1L, length(window_names)), ..retained],
    candidates[, ..retained]
  )))
  
  match_input$treat <- c(
    rep(1L, length(window_names)),
    rep(0L, nrow(candidates))
  )
  
  match_input$window <- factor(
    c(window_names, groups),
    levels = window_names
  )
  
  rownames(match_input) <- c(treated_names, control_names)
  
  # ========================================================================
  # MATCHIT: MatchIt calculates the Mahalanobis distance itself.
  #
  # matching_formula uses all retained environmental variables.
  # exact = ~ window restricts matches to the same distance window.
  #
  # No caliper is passed here: all candidates already satisfy the
  # <= 200-metre elevation restriction above.
  # ========================================================================
  matched <- MatchIt::matchit(
    formula = matching_formula,
    data = match_input,
    method = "nearest",
    distance = "mahalanobis",
    exact = ~ window,
    replace = TRUE,
    ratio = 1,
    m.order = "data"
  )
  
  chosen_controls <- matched$match.matrix[treated_names, 1L]
  unmatched_windows <- window_names[is.na(chosen_controls)]
  
  if (length(unmatched_windows)) {
    mark_missing(band, unmatched_windows)
    
    message(sprintf(
      "Intact cell %s: MatchIt could not match %s; skipping.",
      as.character(a$cell),
      paste(unmatched_windows, collapse = ", ")
    ))
    
    next
  }
  
  best_by_window <- match(
    as.character(chosen_controls),
    control_names
  )
  
  if (anyNA(best_by_window)) {
    stop(
      "MatchIt returned an unknown control for intact cell ",
      a$cell
    )
  }
  
  # Global-covariance distance is a diagnostic only.
  # It is not the distance MatchIt used to select these controls.
  selected_difference <- sweep(
    as.matrix(candidates[best_by_window, ..retained]),
    2,
    as.numeric(unlist(a[, ..retained])),
    FUN = "-"
  )
  
  mahalanobis_d2_global <- rowSums(
    (selected_difference %*% sigma_inverse) *
      selected_difference
  )
  
  for (j in seq_along(window_names)) {
    label <- window_names[j]
    best <- best_by_window[j]
    
    anchor_row <- copy(a)
    matched_row <- copy(candidates[best])
    
    setnames(anchor_row, paste0(names(anchor_row), "_1"))
    setnames(matched_row, paste0(names(matched_row), "_2"))
    
    pair <- cbind(anchor_row, matched_row)
    
    pair[, `:=`(
      geoDist = geo_m[best],
      distance_km = geo_m[best] / 1000,
      dist_range = label,
      dist_label = unname(distance_labels[label]),
      mahalanobis_d2_global = mahalanobis_d2_global[j],
      elevation_difference_m = abs(ELEV_1 - ELEV_2),
      flii_diff = abs(FLII_1 - FLII_2)
    )]
    
    result_list[[
      (i - 1L) * length(window_names) + j
    ]] <- pair
  }
  
  if (i %% 100L == 0L) {
    message(sprintf(
      "Matched/processed: %d/%d",
      i,
      nrow(anchors)
    ))
  }
}

candidate_match <- rbindlist(result_list, fill = TRUE)

if (!nrow(candidate_match)) {
  stop("No complete six-window matches were found.")
}

# Check the actual output, including the 200-metre elevation limit.
stopifnot(
  uniqueN(
    candidate_match,
    by = c("cell_1", "dist_range")
  ) == nrow(candidate_match),
  
  all(candidate_match$FLII_1 >= 9.6),
  all(candidate_match$FLII_2 < 9.6),
  
  all(
    candidate_match$elev_band_1 ==
      candidate_match$elev_band_2
  ),
  
  all(candidate_match$geoDist <= 256000),
  
  all(
    candidate_match$elevation_difference_m <=
      max_elevation_difference_m
  ),
  
  all(mapply(
    function(d, label) {
      any(vapply(
        distance_windows[[label]],
        function(interval) {
          d >= interval[1] & d <= interval[2]
        },
        logical(1)
      ))
    },
    candidate_match$distance_km,
    candidate_match$dist_range
  ))
)

# ------------------------------------------------------------------------------
# 2.5 Retain up to 500 intact pixels with matches in all six windows
# ------------------------------------------------------------------------------

complete_anchors <- candidate_match[
  ,
  .(windows_matched = uniqueN(dist_range)),
  by = .(
    cell_1,
    elev_band = as.character(elev_band_1)
  )
][windows_matched == length(window_names)]

set.seed(124)
chosen_anchors <- complete_anchors[
  , .SD[sample.int(.N, min(.N, max_complete_intact_per_band))],
  by = elev_band
]

result_match <- candidate_match[
  cell_1 %in% chosen_anchors$cell_1
]

stopifnot(
  nrow(result_match) ==
    nrow(chosen_anchors) * length(window_names),
  
  result_match[
    ,
    all(.N == length(window_names)),
    by = cell_1
  ][
    ,
    all(V1)
  ]
)

# ------------------------------------------------------------------------------
# 2.6 Summarize matching results and balance
# ------------------------------------------------------------------------------

band_levels <- levels(dt$elev_band)

matching_summary <- data.table(
  elev_band = band_levels,
  elevation_lower_m = elevation_band_definition$lower_m,
  elevation_upper_m = elevation_band_definition$upper_m,
  
  initial_intact = as.integer(table(
    factor(anchors$elev_band, levels = band_levels)
  )),
  
  complete_six_windows = as.integer(table(
    factor(complete_anchors$elev_band, levels = band_levels)
  )),
  
  retained_intact = as.integer(table(
    factor(chosen_anchors$elev_band, levels = band_levels)
  )),
  
  initial_missing_windows = as.integer(
    rowSums(no_candidate)
  )
)

matching_summary[
  ,
  `:=`(
    excluded_incomplete =
      initial_intact - complete_six_windows,
    
    excluded_after_random_selection =
      complete_six_windows - retained_intact,
    
    final_pairs =
      retained_intact * length(window_names)
  )
]

distance_summary <- CJ(
  elev_band = band_levels,
  dist_range = window_names,
  sorted = FALSE
)

distance_summary[
  ,
  `:=`(
    dist_label = unname(
      distance_labels[dist_range]
    ),
    
    target_pairs =
      matching_summary$retained_intact[
        match(elev_band, matching_summary$elev_band)
      ],
    
    initial_no_candidate = as.integer(
      no_candidate[cbind(
        match(elev_band, rownames(no_candidate)),
        match(dist_range, colnames(no_candidate))
      )]
    )
  )
]

actual <- result_match[
  ,
  .(
    pairs = .N,
    distinct_intact = uniqueN(cell_1),
    distinct_degraded = uniqueN(cell_2),
    
    median_mahalanobis_d2_global =
      median(mahalanobis_d2_global),
    
    median_elevation_difference_m =
      median(elevation_difference_m),
    
    p90_elevation_difference_m =
      as.numeric(
        quantile(elevation_difference_m, 0.9)
      )
  ),
  by = .(
    elev_band = as.character(elev_band_1),
    dist_range
  )
]

distance_summary <- merge(
  distance_summary,
  actual,
  by = c("elev_band", "dist_range"),
  all.x = TRUE,
  sort = FALSE
)

distance_summary[
  is.na(pairs),
  `:=`(
    pairs = 0L,
    distinct_intact = 0L,
    distinct_degraded = 0L
  )
]

setorder(distance_summary, elev_band, dist_label)

# Global SMD for each environmental variable.
balance <- if (nrow(result_match) > 1L) {
  rbindlist(lapply(ecological_priority, function(v) {
    x <- result_match[[paste0(v, "_1")]]
    y <- result_match[[paste0(v, "_2")]]
    
    pooled_sd <- sqrt(
      (stats::var(x) + stats::var(y)) / 2
    )
    
    data.table(
      variable = v,
      retained = v %in% retained,
      
      smd_degraded_minus_intact =
        if (is.finite(pooled_sd) && pooled_sd > 0) {
          (mean(y) - mean(x)) / pooled_sd
        } else {
          NA_real_
        }
    )
  }))
} else {
  data.table()
}

# ------------------------------------------------------------------------------
# 2.7 Save matching outputs
# ------------------------------------------------------------------------------

matching_output_file <- file.path(
  output_dir,
  "result_match_NT_TSMF_FA_single_degraded.rds"
)
saveRDS(result_match, matching_output_file)

fwrite(selection_log,file.path(output_dir, "matching_variable_selection_NT_TSMF_FA.csv"))

fwrite(elevation_band_definition, file.path(output_dir, "matching_elevation_band_definition_NT_TSMF_FA.csv"))

fwrite(matching_summary,file.path(output_dir, "matching_elevation_summary_NT_TSMF_FA.csv"))

fwrite(distance_summary,file.path(output_dir, "matching_distance_summary_NT_TSMF_FA.csv"))

fwrite(balance,file.path(output_dir, "matching_balance_NT_TSMF_FA.csv"))

cat("Final matched pairs:", nrow(result_match), "\n")
print(matching_summary)
print(distance_summary)



# ==============================================================================
# STEP 3. PREDICT DISSIMILARITY FOR THE MATCHED PAIRS
# ==============================================================================

library(data.table)
library(geosphere)

# ------------------------------------------------------------------------------
# 3.1 Load the fitted GDM and matched pairs
# ------------------------------------------------------------------------------

gdm_output_path <- "E:/Weicheng_Code_5/GDM_fitting"
source_dir <- "E:/BILBI - example dataset and code/ObsPair_GDM_v2-main/ObsPair_GDM_v2-main/src"
matching_file <- file.path(
  "E:/Weicheng_Code_5/GDM_fitting",
  "result_match_NT_TSMF_FA_single_degraded.rds"
)

source(file.path(source_dir, "utils.R"))
source(file.path(source_dir, "gdm_functions.R"))

fit <- readRDS(file.path(gdm_output_path, "GDM_fit_NT_TSMF_FA_two_stage.rds"))
result_match <- as.data.table(readRDS(matching_file))


# ------------------------------------------------------------------------------
# 3.2 Define the two-stage pairwise GDM prediction function
# ------------------------------------------------------------------------------

gdm_pairwise_dissimilarity <- function(fit, pair_df, lon1_col = NULL, lat1_col = NULL, lon2_col = NULL, lat2_col = NULL, verbose = TRUE) {
  
  if (verbose) cat("=== Two-stage pairwise GDM prediction ===\n")
  
  pair_df <- as.data.frame(pair_df)
  
  if (!identical(fit$fit_mode, "two_stage_offset_geo")) {
    stop("The supplied fit object is not a two-stage offset-geo model.")
  }
  
  if (!identical(fit$geo_distance_units, "metres")) {
    stop("This prediction function expects geo-distance knots in metres.")
  }
  
  if (nrow(pair_df) == 0L) stop("pair_df contains no rows.")
  
  find_coordinate_column <- function(data, supplied, candidates) {
    if (!is.null(supplied)) {
      if (!supplied %in% names(data)) stop("Coordinate column not found: ", supplied)
      return(supplied)
    }
    
    available <- candidates[candidates %in% names(data)]
    if (length(available) == 0L) stop("None of these coordinate columns were found: ", paste(candidates, collapse = ", "))
    available[1]
  }
  
  lon1_col <- find_coordinate_column(pair_df, lon1_col, c("lon_1", "Lon1", "RAW_Lon1", "longitude_1"))
  lat1_col <- find_coordinate_column(pair_df, lat1_col, c("lat_1", "Lat1", "RAW_Lat1", "latitude_1"))
  lon2_col <- find_coordinate_column(pair_df, lon2_col, c("lon_2", "Lon2", "RAW_Lon2", "longitude_2"))
  lat2_col <- find_coordinate_column(pair_df, lat2_col, c("lat_2", "Lat2", "RAW_Lat2", "latitude_2"))
  
  env_predictors <- fit$stage1_predictors
  env_splines <- fit$stage1_splines
  env_csp <- c(0L, cumsum(env_splines))
  n_pairs <- nrow(pair_df)
  
  environmental_distance <- rep(0, n_pairs)
  
  for (i in seq_along(env_predictors)) {
    predictor <- env_predictors[i]
    variable_1 <- paste0(predictor, "_1")
    variable_2 <- paste0(predictor, "_2")
    
    if (!variable_1 %in% names(pair_df)) stop("Missing Site-1 predictor column: ", variable_1)
    if (!variable_2 %in% names(pair_df)) stop("Missing Site-2 predictor column: ", variable_2)
    
    ns <- env_splines[i]
    parameter_index <- (env_csp[i] + 1L):env_csp[i + 1L]
    coefficient_i <- fit$stage1_coefficients[parameter_index]
    quantile_i <- fit$stage1_quantiles[parameter_index]
    x1 <- pair_df[[variable_1]]
    x2 <- pair_df[[variable_2]]
    predictor_contribution <- rep(0, n_pairs)
    
    for (s in seq_len(ns)) {
      if (s == 1L) {
        q1 <- quantile_i[1]; q2 <- quantile_i[1]; q3 <- quantile_i[2]
      } else if (s == ns) {
        q1 <- quantile_i[ns - 1L]; q2 <- quantile_i[ns]; q3 <- quantile_i[ns]
      } else {
        q1 <- quantile_i[s - 1L]; q2 <- quantile_i[s]; q3 <- quantile_i[s + 1L]
      }
      
      spline_1 <- I_spline(x1, q1, q2, q3)
      spline_2 <- I_spline(x2, q1, q2, q3)
      predictor_contribution <- predictor_contribution + coefficient_i[s] * abs(spline_1 - spline_2)
    }
    
    environmental_distance <- environmental_distance + predictor_contribution
  }
  
  
  # Geographic distance is a pair-level predictor, not a Site-1/Site-2 variable
  geographic_distance_m <- geosphere::distGeo(
    cbind(pair_df[[lon1_col]], pair_df[[lat1_col]]),
    cbind(pair_df[[lon2_col]], pair_df[[lat2_col]])
  )
  
  geo_index <- match("geoDist", fit$predictors)
  if (is.na(geo_index)) stop("The fit object does not contain the geoDist predictor.")
  
  full_csp <- c(0L, cumsum(fit$splines))
  geo_parameter_index <- (full_csp[geo_index] + 1L):full_csp[geo_index + 1L]
  geo_coefficients <- fit$coefficients[geo_parameter_index]
  geo_quantiles <- fit$quantiles[geo_parameter_index]
  geo_ns <- fit$splines[geo_index]
  geographic_contribution <- rep(0, n_pairs)
  
  for (s in seq_len(geo_ns)) {
    if (s == 1L) {
      q1 <- geo_quantiles[1]; q2 <- geo_quantiles[1]; q3 <- geo_quantiles[2]
    } else if (s == geo_ns) {
      q1 <- geo_quantiles[geo_ns - 1L]; q2 <- geo_quantiles[geo_ns]; q3 <- geo_quantiles[geo_ns]
    } else {
      q1 <- geo_quantiles[s - 1L]; q2 <- geo_quantiles[s]; q3 <- geo_quantiles[s + 1L]
    }
    
    geographic_contribution <- geographic_contribution + geo_coefficients[s] * I_spline(geographic_distance_m, q1, q2, q3)
  }
  
  ecological_distance <- environmental_distance + geographic_contribution
  linear_predictor <- fit$intercept + ecological_distance
  predicted_mismatch_probability <- inv.logit(linear_predictor)
  p0 <- inv.logit(fit$intercept)
  dissimilarity <- ObsTrans(p0, obs_trans_weight(fit), predicted_mismatch_probability)$out
  
  output <- data.frame(
    cell_1 = pair_df$cell_1,
    cell_2 = pair_df$cell_2,
    gdm_geographic_distance_km = geographic_distance_m / 1000,
    gdm_environmental_contribution = environmental_distance,
    gdm_geographic_contribution = geographic_contribution,
    gdm_ecological_distance = ecological_distance,
    gdm_linear_predictor = linear_predictor,
    gdm_predicted_mismatch_probability = predicted_mismatch_probability,
    dissimilarity = dissimilarity
  )
  
  if (verbose) {
    cat("Environmental distance:", paste(round(range(environmental_distance, na.rm = TRUE), 6), collapse = " to "), "\n")
    cat("Geographic contribution:", paste(round(range(geographic_contribution, na.rm = TRUE), 6), collapse = " to "), "\n")
    cat("Total ecological distance:", paste(round(range(ecological_distance, na.rm = TRUE), 6), collapse = " to "), "\n")
    cat("Dissimilarity:", paste(round(range(dissimilarity, na.rm = TRUE), 6), collapse = " to "), "\n")
  }
  
  output
}


# ------------------------------------------------------------------------------
# 3.3 Predict compositional dissimilarity
# ------------------------------------------------------------------------------

diss_result <- gdm_pairwise_dissimilarity(
  fit = fit,
  pair_df = result_match,
  lon1_col = "lon_1",
  lat1_col = "lat_1",
  lon2_col = "lon_2",
  lat2_col = "lat_2",
  verbose = TRUE
)


# ------------------------------------------------------------------------------
# 3.4 Combine predictions with the matched-pair data
# ------------------------------------------------------------------------------

prediction_columns <- c(
  "gdm_geographic_distance_km",
  "gdm_environmental_contribution",
  "gdm_geographic_contribution",
  "gdm_ecological_distance",
  "gdm_linear_predictor",
  "gdm_predicted_mismatch_probability",
  "dissimilarity"
)

diss_result_match <- cbind(
  as.data.frame(result_match),
  diss_result[, prediction_columns, drop = FALSE]
)

diss_result_match <- as.data.table(diss_result_match)
if (anyDuplicated(names(diss_result_match)) > 0L) stop("Duplicated columns after prediction: ", paste(unique(names(diss_result_match)[duplicated(names(diss_result_match))]), collapse = ", "))


# ------------------------------------------------------------------------------
# 3.5 Validate the prediction output
# ------------------------------------------------------------------------------

stopifnot(nrow(diss_result_match) == nrow(result_match))

cat("\nPrediction summary:\n")
print(summary(diss_result_match[, ..prediction_columns]))

# ------------------------------------------------------------------------------
# 3.6 Save prediction results
# ------------------------------------------------------------------------------

saveRDS(
  diss_result_match,
  file.path(gdm_output_path, "diss_result_match_NT_TSMF_FA_two_stage.rds")
)
write.csv(
  diss_result_match,
  file.path(gdm_output_path, "diss_result_match_NT_TSMF_FA_two_stage.csv"),
  row.names = FALSE
)



# ==============================================================================
# APPENDIX A. SUPERSEDED GAM/GAMM ANALYSIS (DISABLED)
# ==============================================================================
if (FALSE) {
# The superseded two-figure GAM code is retained below only for provenance.
# It is not executed. The corrected three-model analysis follows this block.
#为每一个距离梯度构建一个gamm的计算和绘图代码，所有距离梯度组合在一张图中，还有许多gam质量控制的图
#-------------------------------------------------------------------------------
library(data.table)
library(mgcv)
library(ggplot2)

## Ensure that the beta-regression response is strictly within (0, 1)
diss_result_match <- readRDS("E:/Weicheng_Code_5/GDM_fitting/diss_result_match_NT_TSMF_FA_two_stage.rds")

analysis_dt <- copy(diss_result_match)
setDT(analysis_dt)


analysis_dt <- as.data.table(copy(diss_result_match))

# ==============================================================================
# A.1.1 Construct pairwise environmental differences
# ==============================================================================

# FLII is excluded because flii_diff is the focal predictor.
candidate_env <- c("BD", "CLAY", "PH", "SILT", "TWI", "PTA", "EAAS", "WDI", "WDX", "EPI", "TXX", "TNI", "RUG", "TRX", "ELEV")
required_pair_columns <- c(paste0(candidate_env, "_1"), paste0(candidate_env, "_2"))
missing_columns <- setdiff(required_pair_columns, names(analysis_dt))

if (length(missing_columns) > 0L) stop("Missing environmental columns: ", paste(missing_columns, collapse = ", "))

# Calculate the absolute environmental difference between the two paired sites.
for (v in candidate_env) analysis_dt[, (paste0(v, "_diff")) := abs(get(paste0(v, "_1")) - get(paste0(v, "_2")))]

candidate_adjustment_vars <- paste0(candidate_env, "_diff")


# ==============================================================================
# A.1.2 Define ecological retention priority
# ==============================================================================

# Variables listed earlier have higher priority. For a pair with |r| > 0.7,
# the higher-priority variable is retained and the other variable is removed.
ecological_priority <- c(
  "ELEV_diff", "PTA_diff", "TXX_diff", "TNI_diff", "WDI_diff",
  "WDX_diff", "EPI_diff", "EAAS_diff", "TWI_diff", "RUG_diff",
  "TRX_diff", "PH_diff", "CLAY_diff", "SILT_diff", "BD_diff"
)

if (!setequal(ecological_priority, candidate_adjustment_vars)) {
  stop("ecological_priority must contain every candidate adjustment variable exactly once.")
}


# ==============================================================================
# A.1.3 Remove variables with insufficient or zero variation
# ==============================================================================

variation_table <- rbindlist(lapply(candidate_adjustment_vars, function(v) {
  x <- analysis_dt[[v]]
  data.table(variable = v, valid_n = sum(is.finite(x)), unique_n = uniqueN(x[is.finite(x)]), standard_deviation = sd(x, na.rm = TRUE))
}))

variation_table[, usable := valid_n >= 20L & unique_n >= 5L & is.finite(standard_deviation) & standard_deviation > 0]

excluded_low_variation <- variation_table[usable == FALSE]
candidate_adjustment_vars <- variation_table[usable == TRUE, variable]
ecological_priority <- ecological_priority[ecological_priority %in% candidate_adjustment_vars]

cat("\nVariables excluded because of insufficient variation:\n")
print(excluded_low_variation)


# ==============================================================================
# A.1.4 Calculate correlations
# ==============================================================================

correlation_method <- "pearson"
correlation_cutoff <- 0.70

correlation_matrix <- cor(analysis_dt[, ..candidate_adjustment_vars], use = "pairwise.complete.obs", method = correlation_method)

cat("\nCorrelation matrix:\n")
print(round(correlation_matrix, 2))


# ==============================================================================
# A.1.5 Extract variable pairs with |r| > 0.7
# ==============================================================================

correlation_pairs <- as.data.table(as.table(correlation_matrix))
setnames(correlation_pairs, c("variable_1", "variable_2", "correlation"))

correlation_pairs <- correlation_pairs[variable_1 != variable_2]
correlation_pairs[, `:=`(
  variable_1 = as.character(variable_1),
  variable_2 = as.character(variable_2),
  absolute_correlation = abs(correlation)
)]

correlation_pairs[, pair_id := vapply(seq_len(.N), function(i) paste(sort(c(variable_1[i], variable_2[i])), collapse = "~"), character(1))]
correlation_pairs <- correlation_pairs[!duplicated(pair_id) & absolute_correlation > correlation_cutoff]
correlation_pairs[, pair_id := NULL]
setorder(correlation_pairs, -absolute_correlation)

cat("\nVariable pairs with |r| > ", correlation_cutoff, ":\n", sep = "")
print(correlation_pairs)


# ==============================================================================
# A.1.6 Filter correlated variables according to ecological priority
# ==============================================================================

selected_covariates <- character(0)

filter_log <- data.table(
  variable = character(),
  priority = integer(),
  decision = character(),
  conflicting_variable = character(),
  correlation = numeric(),
  reason = character()
)

for (i in seq_along(ecological_priority)) {
  current_variable <- ecological_priority[i]
  
  if (length(selected_covariates) == 0L) {
    selected_covariates <- current_variable
    filter_log <- rbind(filter_log, data.table(
      variable = current_variable, priority = i, decision = "Retained",
      conflicting_variable = NA_character_, correlation = NA_real_,
      reason = "Highest-priority available variable"
    ))
    next
  }
  
  current_correlations <- correlation_matrix[current_variable, selected_covariates, drop = TRUE]
  current_correlations <- current_correlations[is.finite(current_correlations)]
  conflicting_variables <- names(current_correlations)[abs(current_correlations) > correlation_cutoff]
  
  if (length(conflicting_variables) == 0L) {
    selected_covariates <- c(selected_covariates, current_variable)
    filter_log <- rbind(filter_log, data.table(
      variable = current_variable, priority = i, decision = "Retained",
      conflicting_variable = NA_character_, correlation = NA_real_,
      reason = paste0("No correlation above |r| = ", correlation_cutoff, " with retained variables")
    ))
  } else {
    strongest_conflict <- conflicting_variables[which.max(abs(current_correlations[conflicting_variables]))]
    strongest_correlation <- current_correlations[strongest_conflict]
    
    filter_log <- rbind(filter_log, data.table(
      variable = current_variable, priority = i, decision = "Removed",
      conflicting_variable = strongest_conflict, correlation = as.numeric(strongest_correlation),
      reason = paste0("Removed because |r| > ", correlation_cutoff, " with the higher-priority variable ", strongest_conflict)
    ))
  }
}

removed_covariates <- filter_log[decision == "Removed", variable]

cat("\nRetained adjustment covariates:\n")
print(selected_covariates)

cat("\nRemoved adjustment covariates:\n")
print(removed_covariates)

cat("\nComplete filtering log:\n")
print(filter_log)


# ==============================================================================
# A.1.7 Prepare the final analysis dataset
# ==============================================================================

keep_vars <- c("dissimilarity", "flii_diff", "dist_label", selected_covariates)
analysis_dt <- analysis_dt[complete.cases(analysis_dt[, ..keep_vars])]
analysis_dt[, dissimilarity := pmin(pmax(dissimilarity, 1e-6), 1 - 1e-6)]

dist_levels <- c(8, 16, 32, 64, 128, 256)
analysis_dt[, dist_label := factor(dist_label, levels = dist_levels)]
analysis_dt <- analysis_dt[!is.na(dist_label)]

cat("\nFinal sample size:", nrow(analysis_dt), "\n")
cat("Final adjustment covariates:", paste(selected_covariates, collapse = ", "), "\n")



#============================================================================================================================================================
#============================================================================================================================================================
##==============GAMM============================================================
# ==============================================================================
# DISTANCE-SPECIFIC GAM ANALYSIS AND PLOTTING
# ==============================================================================

library(data.table)
library(mgcv)
library(ggplot2)

output_path <- "E:/Weicheng_Code_5/GDM_fitting"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# A.1.8 Check input data
# ==============================================================================

required_vars <- c("dissimilarity", "flii_diff", "dist_label", selected_covariates)
missing_vars <- setdiff(required_vars, names(analysis_dt))
if (length(missing_vars) > 0L) stop("Missing variables: ", paste(missing_vars, collapse = ", "))

analysis_dt <- analysis_dt[complete.cases(analysis_dt[, ..required_vars])]
analysis_dt[, dissimilarity := pmin(pmax(dissimilarity, 1e-6), 1 - 1e-6)]

dist_levels <- c("8", "16", "32", "64", "128", "256")
analysis_dt[, dist_label := factor(as.character(dist_label), levels = dist_levels)]
analysis_dt <- analysis_dt[!is.na(dist_label)]

cat("Final sample size:", nrow(analysis_dt), "\n")
cat("Selected adjustment covariates:", paste(selected_covariates, collapse = ", "), "\n")
print(analysis_dt[, .N, by = dist_label])


# ==============================================================================
# A.1.9 Construct the formula automatically
# ==============================================================================

# Separate nonlinear FLII smooths are fitted for each distance class, while all
# automatically selected environmental differences are included as adjustments.
adjustment_terms <- sprintf("s(%s, bs = 'cr', k = 5)", selected_covariates)
model_terms <- c("dist_label", "s(flii_diff, by = dist_label, bs = 'cr', k = 5)", adjustment_terms)
gamm_formula <- as.formula(paste("dissimilarity ~", paste(model_terms, collapse = " + ")))

cat("\nFinal GAM formula:\n")
print(gamm_formula)


# ==============================================================================
# A.1.10 Fit the full GAM
# ==============================================================================

mod_full <- gam(formula = gamm_formula, data = analysis_dt, family = betar(link = "logit"), method = "fREML", select = TRUE)

cat("\nModel summary:\n")
print(summary(mod_full))


# ==============================================================================
# A.1.11 GAM quality control
# ==============================================================================

cat("\nGAM basis-dimension and residual checks:\n")
gam.check(mod_full)

cat("\nConcurvity diagnostics:\n")
print(concurvity(mod_full, full = TRUE))

cat("\nModel convergence:", mod_full$converged, "\n")

pdf(file.path(output_path, "GAM_quality_control.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
gam.check(mod_full)
dev.off()


# ==============================================================================
# A.1.12 Extract FLII smooth statistics for each distance class
# ==============================================================================

model_summary <- summary(mod_full)
smooth_table <- as.data.frame(model_summary$s.table)
smooth_table$smooth_term <- rownames(model_summary$s.table)
rownames(smooth_table) <- NULL

flii_smooth_table <- smooth_table[grepl("^s\\(flii_diff\\):dist_label", smooth_table$smooth_term), , drop = FALSE]
if (nrow(flii_smooth_table) != length(dist_levels)) warning("The number of FLII smooths does not match the number of distance classes.")

p_column <- grep("^p", names(flii_smooth_table), value = TRUE, ignore.case = TRUE)[1]
statistic_column <- intersect(c("F", "Chi.sq"), names(flii_smooth_table))[1]

if (is.na(p_column)) stop("No p-value column was found in the GAM smooth table.")
if (is.na(statistic_column)) stop("No F or Chi.sq column was found in the GAM smooth table.")

flii_statistics <- data.table(
  dist_label = factor(dist_levels, levels = dist_levels),
  smooth_term = flii_smooth_table$smooth_term,
  edf = as.numeric(flii_smooth_table[["edf"]]),
  reference_df = as.numeric(flii_smooth_table[["Ref.df"]]),
  test_statistic = as.numeric(flii_smooth_table[[statistic_column]]),
  p_value = as.numeric(flii_smooth_table[[p_column]])
)

flii_statistics[, significance := fcase(
  is.na(p_value), "p = NA",
  p_value < 0.001, "p < 0.001",
  p_value < 0.01, "p < 0.01",
  p_value < 0.05, "p < 0.05",
  default = paste0("p = ", sprintf("%.3f", p_value))
)]

flii_statistics[, label := paste0("edf = ", sprintf("%.2f", edf), "\n", significance)]

cat("\nFLII smooth statistics by distance class:\n")
print(flii_statistics)


# ==============================================================================
# A.1.13 Construct prediction data automatically
# ==============================================================================

# Adjustment covariates are held at their overall medians.
adjustment_values <- vapply(selected_covariates, function(v) median(analysis_dt[[v]], na.rm = TRUE), numeric(1))

prediction_list <- lapply(dist_levels, function(current_dist) {
  current_data <- analysis_dt[dist_label == current_dist]
  pred_dt <- data.table(flii_diff = seq(min(current_data$flii_diff, na.rm = TRUE), max(current_data$flii_diff, na.rm = TRUE), length.out = 300), dist_label = factor(current_dist, levels = dist_levels))
  for (v in selected_covariates) pred_dt[, (v) := adjustment_values[v]]
  pred_link <- predict(mod_full, newdata = pred_dt, type = "link", se.fit = TRUE)
  link_inverse <- mod_full$family$linkinv
  pred_dt[, `:=`(
    predicted_dissimilarity = link_inverse(pred_link$fit),
    lower = link_inverse(pred_link$fit - 1.96 * pred_link$se.fit),
    upper = link_inverse(pred_link$fit + 1.96 * pred_link$se.fit)
  )]
  pred_dt
})

prediction_dt <- rbindlist(prediction_list)
prediction_dt[, dist_label := factor(as.character(dist_label), levels = dist_levels)]


# ==============================================================================
# A.1.14 Sample original points for plotting only
# ==============================================================================

# The GAM uses all observations; sampling only reduces rendering time.
set.seed(123)
point_dt <- analysis_dt[, .SD[sample.int(.N, min(.N, 5000L))], by = dist_label]


# ==============================================================================
# A.1.15 Prepare panel annotations
# ==============================================================================

panel_statistics <- analysis_dt[, .(
  sample_size = .N,
  flii_min = min(flii_diff, na.rm = TRUE),
  flii_max = max(flii_diff, na.rm = TRUE)
), by = dist_label]

panel_statistics <- merge(panel_statistics, flii_statistics, by = "dist_label", all.x = TRUE, sort = FALSE)
panel_statistics[, dist_label := factor(as.character(dist_label), levels = dist_levels)]

x_limits <- c(0, max(analysis_dt$flii_diff, na.rm = TRUE))
y_limits <- c(0, 1)

panel_statistics[, `:=`(
  label_x = x_limits[1] + 0.04 * diff(x_limits),
  label_y = 0.97,
  panel_label = paste0("n = ", format(sample_size, big.mark = ",", scientific = FALSE), "\n", "edf = ", sprintf("%.2f", edf), "\n", significance)
)]


# ==============================================================================
# A.1.16 Define colours and facet labels
# ==============================================================================

dist_cols <- c("8" = "#9ECAE1", "16" = "#6BAED6", "32" = "#4292C6", "64" = "#2171B5", "128" = "#08519C", "256" = "#08306B")
facet_labels <- setNames(paste0(dist_levels, " km"), dist_levels)


# ==============================================================================
# A.1.17 Plot the six distance-specific nonlinear relationships
# ==============================================================================

p_combined <- ggplot() +
  geom_point(data = point_dt, aes(flii_diff, dissimilarity, colour = dist_label), alpha = 0.14, size = 0.65, show.legend = FALSE) +
  geom_ribbon(data = prediction_dt, aes(flii_diff, ymin = lower, ymax = upper, fill = dist_label), alpha = 0.22, colour = NA, show.legend = FALSE) +
  geom_line(data = prediction_dt, aes(flii_diff, predicted_dissimilarity, colour = dist_label), linewidth = 1.15, show.legend = FALSE) +
  geom_text(data = panel_statistics, aes(label_x, label_y, label = panel_label), hjust = 0, vjust = 1, size = 3.4, lineheight = 1.05) +
  facet_wrap(~dist_label, ncol = 3, nrow = 2, scales = "fixed", labeller = as_labeller(facet_labels)) +
  scale_colour_manual(values = dist_cols, drop = FALSE) +
  scale_fill_manual(values = dist_cols, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, ceiling(x_limits[2]), by = 2), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
  labs(
    title = "FLII difference and predicted compositional dissimilarity across distance classes",
    subtitle = "Distance-specific nonlinear smooths adjusted for automatically selected environmental differences",
    x = "FLII difference",
    y = "Predicted compositional dissimilarity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "grey93", colour = "grey35"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(8, 12, 8, 8)
  )

p_combined


# # ==============================================================================
# # A.1.18 Save outputs
# # ==============================================================================
# 
# ggsave(file.path(output_path, "FLII_dissimilarity_by_distance_adjusted_GAM.png"), p_combined, width = 12, height = 8, dpi = 400, bg = "white")
# ggsave(file.path(output_path, "FLII_dissimilarity_by_distance_adjusted_GAM.pdf"), p_combined, width = 12, height = 8)
# 
# fwrite(prediction_dt, file.path(output_path, "FLII_dissimilarity_GAM_predictions.csv"))
# fwrite(flii_statistics, file.path(output_path, "FLII_smooth_statistics_by_distance.csv"))
# fwrite(data.table(selected_covariate = selected_covariates, fixed_value = adjustment_values), file.path(output_path, "GAM_selected_adjustment_covariates.csv"))

# saveRDS(mod_full, file.path(output_path, "FLII_dissimilarity_adjusted_GAM.rds"))





# ==============================================================================
# FLII DIFFERENCE × ELEVATION INTERACTION
# ==============================================================================

# ==============================================================================
# GAM: FLII difference × elevation class × distance class
# ==============================================================================

library(data.table)
library(mgcv)
library(ggplot2)
library(ggpubr)

analysis_dt <- as.data.table(copy(diss_result_match))

# ------------------------------------------------------------------------------
# A.2.1 Define the background elevation of each pair
# ------------------------------------------------------------------------------

analysis_dt[, ELEV_mean := (ELEV_1 + ELEV_2) / 2]

dist_levels <- c(8, 16, 32, 64, 128, 256)
elev_levels <- c("low", "middle", "high")

analysis_dt[, dist_label := factor(dist_label, levels = dist_levels)]
# Reuse the region-specific equal-width elevation bands created during matching.
analysis_dt[, ELEV_group := factor(as.character(elev_band_1), levels = elev_levels)]

# ------------------------------------------------------------------------------
# A.2.2 Prepare adjustment covariates
# ------------------------------------------------------------------------------

# selected_covariates comes from the previous correlation-filtering procedure.
# ELEV_diff may remain because it controls elevation differences within pairs,
# whereas ELEV_group describes the overall elevational setting of each pair.

selected_covariates <- selected_covariates[selected_covariates %in% names(analysis_dt)]

required_vars <- c("dissimilarity", "flii_diff", "dist_label", "ELEV_mean", "ELEV_group", selected_covariates)
analysis_dt <- analysis_dt[complete.cases(analysis_dt[, ..required_vars])]
analysis_dt[, dissimilarity := pmin(pmax(as.numeric(dissimilarity), 1e-6), 1 - 1e-6)]

# Remove unused distance and elevation levels after complete-case filtering
analysis_dt[, dist_label := droplevels(dist_label)]
analysis_dt[, ELEV_group := droplevels(ELEV_group)]

# Each level represents one distance-class × elevation-class combination
analysis_dt[, dist_elev_group := interaction(dist_label, ELEV_group, sep = " km | ", drop = TRUE)]
analysis_dt[, dist_elev_group := factor(dist_elev_group)]

# ------------------------------------------------------------------------------
# A.2.3 Examine sample size and actual elevation range
# ------------------------------------------------------------------------------

elevation_summary <- analysis_dt[, .(
  n_pairs = .N,
  elevation_min = min(ELEV_mean),
  elevation_median = median(ELEV_mean),
  elevation_max = max(ELEV_mean),
  flii_min = min(flii_diff),
  flii_max = max(flii_diff),
  dissimilarity_mean = mean(dissimilarity)
), by = .(dist_label, ELEV_group)]

setorder(elevation_summary, dist_label, ELEV_group)
print(elevation_summary)

# Warn about strata with limited sample sizes
small_strata <- elevation_summary[n_pairs < 100]
if (nrow(small_strata) > 0L) {
  warning("Some distance × elevation strata contain fewer than 100 pairs.")
  print(small_strata)
}

# ------------------------------------------------------------------------------
# A.2.4 Construct the GAM formula automatically
# ------------------------------------------------------------------------------

# dist_elev_group provides a separate intercept for every combination.
# The factor-by smooth provides a separate nonlinear FLII curve for every
# distance × elevation combination.

main_terms <- c(
  "dist_elev_group",
  "s(flii_diff, by = dist_elev_group, bs = 'cr', k = 5)"
)

adjustment_terms <- if (length(selected_covariates) > 0L) {
  paste0("s(", selected_covariates, ", bs = 'cr', k = 5)")
} else {
  character(0)
}

gamm_formula <- as.formula(paste("dissimilarity ~", paste(c(main_terms, adjustment_terms), collapse = " + ")))
print(gamm_formula)

# ------------------------------------------------------------------------------
# A.2.5 Fit the model
# ------------------------------------------------------------------------------

mod_full <- gam(
  formula = gamm_formula,
  data = analysis_dt,
  family = betar(link = "logit"),
  method = "REML",
  select = TRUE
)

print(summary(mod_full))
gam.check(mod_full)

# Additional concurvity diagnostic
concurvity_results <- concurvity(mod_full, full = TRUE)
print(concurvity_results)

# ------------------------------------------------------------------------------
# A.2.6 Extract model smooth statistics
# ------------------------------------------------------------------------------

smooth_table <- as.data.table(summary(mod_full)$s.table, keep.rownames = "smooth")

rename_map <- c("Ref.df" = "reference_df", "F" = "F_statistic", "Chi.sq" = "chi_square", "p-value" = "p_value")
old_names <- intersect(names(rename_map), names(smooth_table))
setnames(smooth_table, old = old_names, new = unname(rename_map[old_names]))

if (!"p_value" %in% names(smooth_table)) stop("No p-value column was found. Available columns: ", paste(names(smooth_table), collapse = ", "))

smooth_table[, significance := fifelse(
  p_value < 0.001, "p < 0.001",
  fifelse(p_value < 0.01, "p < 0.01", fifelse(p_value < 0.05, "p < 0.05", sprintf("p = %.3f", as.numeric(p_value))))
)]

print(smooth_table)

# ------------------------------------------------------------------------------
# A.2.7 Construct prediction data
# ------------------------------------------------------------------------------

# Adjustment covariates are fixed at their medians during prediction.
adjustment_medians <- setNames(lapply(selected_covariates, function(v) median(analysis_dt[[v]], na.rm = TRUE)), selected_covariates)

prediction_list <- lapply(seq_len(nrow(elevation_summary)), function(i) {
  current_dist <- elevation_summary$dist_label[i]
  current_elev <- elevation_summary$ELEV_group[i]
  current_data <- analysis_dt[dist_label == current_dist & ELEV_group == current_elev]
  
  prediction_data <- data.table(
    flii_diff = seq(min(current_data$flii_diff), max(current_data$flii_diff), length.out = 200),
    dist_label = factor(current_dist, levels = levels(analysis_dt$dist_label)),
    ELEV_group = factor(current_elev, levels = levels(analysis_dt$ELEV_group))
  )
  
  prediction_data[, dist_elev_group := interaction(dist_label, ELEV_group, sep = " km | ", drop = TRUE)]
  prediction_data[, dist_elev_group := factor(dist_elev_group, levels = levels(analysis_dt$dist_elev_group))]
  
  for (v in selected_covariates) prediction_data[, (v) := adjustment_medians[[v]]]
  
  prediction_link <- predict(mod_full, newdata = prediction_data, type = "link", se.fit = TRUE)
  prediction_data[, `:=`(
    predicted_dissimilarity = plogis(as.numeric(prediction_link$fit)),
    lower = plogis(as.numeric(prediction_link$fit) - 1.96 * as.numeric(prediction_link$se.fit)),
    upper = plogis(as.numeric(prediction_link$fit) + 1.96 * as.numeric(prediction_link$se.fit))
  )]
  
  prediction_data
})

pred_dt <- rbindlist(prediction_list, fill = TRUE)
pred_dt[, dist_label := factor(dist_label, levels = dist_levels)]
pred_dt[, ELEV_group := factor(ELEV_group, levels = elev_levels)]

# ------------------------------------------------------------------------------
# A.2.8 Prepare panel labels
# ------------------------------------------------------------------------------

panel_statistics <- analysis_dt[, .(
  n_pairs = .N,
  elevation_min = min(ELEV_mean),
  elevation_max = max(ELEV_mean)
), by = .(dist_label, ELEV_group)]

panel_statistics[, sample_label := paste0(
  ELEV_group, ": n = ", format(n_pairs, big.mark = ","),
  " (", round(elevation_min), "–", round(elevation_max), " m)"
)]

distance_labels <- panel_statistics[, .(
  panel_label = paste(sample_label, collapse = "\n")
), by = dist_label]

# ------------------------------------------------------------------------------
# A.2.9 Plot six distance panels with three elevation curves
# ------------------------------------------------------------------------------

elevation_colours <- c(
  "low" = "#2CA25F",
  "middle" = "#F39C34",
  "high" = "#8C2D04"
)

plot_list <- lapply(levels(droplevels(analysis_dt$dist_label)), function(d) {
  plot_data <- analysis_dt[dist_label == d]
  prediction_data <- pred_dt[dist_label == d]
  annotation_text <- distance_labels[dist_label == d, panel_label]
  
  ggplot() +
    geom_point(
      data = plot_data,
      aes(x = flii_diff, y = dissimilarity, colour = ELEV_group),
      alpha = 0.10,
      size = 0.55
    ) +
    geom_ribbon(
      data = prediction_data,
      aes(x = flii_diff, ymin = lower, ymax = upper, fill = ELEV_group, group = ELEV_group),
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(
      data = prediction_data,
      aes(x = flii_diff, y = predicted_dissimilarity, colour = ELEV_group, group = ELEV_group),
      linewidth = 1.1
    ) +
    annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.05, vjust = 1.15, size = 3, lineheight = 1.1) +
    scale_colour_manual(values = elevation_colours, drop = FALSE) +
    scale_fill_manual(values = elevation_colours, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.6),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(7, 10, 7, 7)
    ) +
    labs(title = paste0(d, " km"), x = "FLII difference", y = "Predicted compositional dissimilarity")
})

p_combined <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, align = "hv", common.legend = TRUE, legend = "bottom")

p_final <- annotate_figure(
  p_combined,
  top = text_grob("FLII difference and compositional dissimilarity across elevation and distance classes", face = "bold", size = 15),
  bottom = text_grob("Elevation class is based on the mean elevation of the two paired sites", size = 10)
)

print(p_final)
}


# ==============================================================================
# STEP 4. DESCRIPTIVE GAM SMOOTHING AND VISUALISATION (UNDER DEVELOPMENT)
# GAMs are used only as descriptive plotting smoothers, not as an independent
# inferential test of the FLII effect. This section is still being debugged.
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Load packages and prediction data
# ------------------------------------------------------------------------------

library(data.table)
library(mgcv)
library(ggplot2)

gdm_output_path <- "E:/Weicheng_Code_5/GDM_fitting"
gam_output_path <- file.path(gdm_output_path, "GAM_results_NT_TSMF_FA")
dir.create(gam_output_path, recursive = TRUE, showWarnings = FALSE)

prediction_file <- file.path(
  gdm_output_path,
  "diss_result_match_NT_TSMF_FA_two_stage.rds"
)

if (!file.exists(prediction_file)) {
  stop("The two-stage GDM prediction file does not exist: ", prediction_file)
}

gam_data <- as.data.table(readRDS(prediction_file))

# ------------------------------------------------------------------------------
# 4.2 Prepare the GAM dataset and environmental-difference covariates
# ------------------------------------------------------------------------------

required_pair_columns <- c(
  "cell_1", "dissimilarity", "flii_diff", "distance_km",
  "dist_label", "elev_band_1"
)

missing_pair_columns <- setdiff(required_pair_columns, names(gam_data))
if (length(missing_pair_columns)) {
  stop("Missing GAM input columns: ", paste(missing_pair_columns, collapse = ", "))
}

# The same intact forest pixel occurs once in each distance window. This factor
# identifies those repeated observations and is fitted as a random intercept.
gam_data[, intact_cell := factor(cell_1)]

distance_levels <- c("8", "16", "32", "64", "128", "256")
elevation_levels <- c("low", "middle", "high")

gam_data[, dist_label := factor(as.character(dist_label), levels = distance_levels)]
# Elevation classes were defined during matching and are reused unchanged here.
# The GAM does not calculate a new elevation classification.
gam_data[, elev_band := factor(as.character(elev_band_1), levels = elevation_levels)]

# Beta regression requires a response strictly between zero and one.
gam_data[, dissimilarity := pmin(
  pmax(as.numeric(dissimilarity), 1e-6),
  1 - 1e-6
)]

# Construct absolute environmental differences. FLII is excluded because it is
# the focal predictor, while distance and mean elevation enter models explicitly.
environmental_variables <- c(
  "BD", "CLAY", "PH", "SILT", "TWI", "PTA", "EAAS", "WDI",
  "WDX", "EPI", "TXX", "TNI", "RUG", "TRX", "ELEV"
)

environmental_pair_columns <- c(
  paste0(environmental_variables, "_1"),
  paste0(environmental_variables, "_2")
)

missing_environmental_columns <- setdiff(
  environmental_pair_columns,
  names(gam_data)
)

if (length(missing_environmental_columns)) {
  stop(
    "Missing paired environmental columns: ",
    paste(missing_environmental_columns, collapse = ", ")
  )
}

for (variable in environmental_variables) {
  difference_name <- paste0(variable, "_diff")
  gam_data[, (difference_name) := abs(
    get(paste0(variable, "_1")) - get(paste0(variable, "_2"))
  )]
}

environmental_priority <- c(
  "ELEV_diff", "PTA_diff", "TXX_diff", "TNI_diff", "WDI_diff",
  "WDX_diff", "EPI_diff", "EAAS_diff", "TWI_diff", "RUG_diff",
  "TRX_diff", "PH_diff", "CLAY_diff", "SILT_diff", "BD_diff"
)

# ------------------------------------------------------------------------------
# 4.3 Screen adjustment covariates and prepare model inputs
# ------------------------------------------------------------------------------

# Remove unusable covariates and then apply the same deterministic ecological
# priority rule used elsewhere: a lower-priority variable is excluded whenever
# its absolute correlation with any retained variable exceeds 0.70.
variation_table_gam <- rbindlist(lapply(environmental_priority, function(v) {
  values <- gam_data[[v]]
  data.table(
    variable = v,
    valid_n = sum(is.finite(values)),
    unique_n = uniqueN(values[is.finite(values)]),
    standard_deviation = sd(values, na.rm = TRUE)
  )
}))

usable_adjustments <- variation_table_gam[
  valid_n >= 20L & unique_n >= 5L &
    is.finite(standard_deviation) & standard_deviation > 0,
  variable
]

adjustment_correlation <- cor(
  gam_data[, ..usable_adjustments],
  use = "pairwise.complete.obs",
  method = "pearson"
)

selected_adjustments <- character()
adjustment_selection_log <- vector("list", length(environmental_priority))

for (i in seq_along(environmental_priority)) {
  variable <- environmental_priority[i]
  usable <- variable %in% usable_adjustments
  conflicts <- character()

  if (usable && length(selected_adjustments)) {
    correlations <- adjustment_correlation[
      variable,
      selected_adjustments,
      drop = TRUE
    ]
    conflicts <- names(correlations)[
      is.finite(correlations) & abs(correlations) > 0.70
    ]
  }

  retained_variable <- usable && !length(conflicts)
  if (retained_variable) {
    selected_adjustments <- c(selected_adjustments, variable)
  }

  adjustment_selection_log[[i]] <- data.table(
    variable = variable,
    retained = retained_variable,
    reason = if (!usable) {
      "Insufficient variation"
    } else if (length(conflicts)) {
      paste("Correlated above 0.70 with", paste(conflicts, collapse = ", "))
    } else {
      "Retained"
    }
  )
}

adjustment_selection_log <- rbindlist(adjustment_selection_log)
fwrite(
  adjustment_selection_log,
  file.path(gam_output_path, "GAM_adjustment_variable_selection.csv")
)

cat("GAM adjustment covariates:",
    paste(selected_adjustments, collapse = ", "), "\n")

complete_variables <- c(
  "dissimilarity", "flii_diff", "distance_km", "dist_label",
  "elev_band", "intact_cell", selected_adjustments
)

gam_data <- gam_data[
  complete.cases(gam_data[, ..complete_variables]) &
    !is.na(dist_label) & !is.na(elev_band)
]

gam_data[, intact_cell := droplevels(intact_cell)]
gam_data[, dist_label := droplevels(dist_label)]
gam_data[, elev_band := droplevels(elev_band)]

if (!nrow(gam_data)) stop("No complete observations remain for GAM fitting.")

adjustment_terms <- if (length(selected_adjustments)) {
  sprintf("s(%s, bs = 'cr', k = 5)", selected_adjustments)
} else {
  character()
}

random_effect_term <- "s(intact_cell, bs = 're')"
adjustment_medians <- vapply(
  selected_adjustments,
  function(v) median(gam_data[[v]], na.rm = TRUE),
  numeric(1)
)

add_prediction_adjustments <- function(newdata) {
  for (v in selected_adjustments) {
    newdata[, (v) := adjustment_medians[[v]]]
  }
  newdata
}

# ------------------------------------------------------------------------------
# 4.4 Define conditional-average prediction with intact-pixel random effects
# ------------------------------------------------------------------------------

# Produce a conditional-average curve that explicitly includes the estimated
# intact-pixel random intercepts. Each curve-grid row is crossed with the
# observed intact cells supplied in profiles, predictions are calculated with
# s(intact_cell) included, and response-scale predictions are then averaged.
# The lower and upper bounds average the cell-specific model intervals; they
# should be interpreted as approximate conditional-average uncertainty bands.
predict_random_effect_average <- function(model, curve_grid, profiles) {
  curve_grid <- copy(as.data.table(curve_grid))
  profiles <- unique(copy(as.data.table(profiles)))

  if (!nrow(profiles)) stop("No intact-cell profiles were supplied for prediction.")

  curve_grid[, curve_row := .I]
  duplicated_profile_names <- intersect(names(curve_grid), names(profiles))
  if (length(duplicated_profile_names)) {
    stop(
      "Prediction curve and profile tables contain duplicate columns: ",
      paste(duplicated_profile_names, collapse = ", ")
    )
  }

  inverse_link <- model$family$linkinv
  curve_chunks <- split(
    seq_len(nrow(curve_grid)),
    ceiling(seq_len(nrow(curve_grid)) / 25L)
  )

  averaged <- rbindlist(lapply(curve_chunks, function(rows) {
    grid_chunk <- curve_grid[rows]
    expanded_curve <- grid_chunk[
      rep(seq_len(.N), each = nrow(profiles))
    ]
    expanded_profiles <- profiles[
      rep(seq_len(.N), times = nrow(grid_chunk))
    ]
    expanded <- cbind(expanded_curve, expanded_profiles)
    expanded[, intact_cell := factor(
      as.character(intact_cell),
      levels = levels(gam_data$intact_cell)
    )]
    expanded <- add_prediction_adjustments(expanded)

    prediction <- predict(
      model,
      newdata = expanded,
      type = "link",
      se.fit = TRUE
    )

    expanded[, `:=`(
      conditional_prediction = inverse_link(as.numeric(prediction$fit)),
      conditional_lower = inverse_link(as.numeric(prediction$fit) -
                                         1.96 * as.numeric(prediction$se.fit)),
      conditional_upper = inverse_link(as.numeric(prediction$fit) +
                                         1.96 * as.numeric(prediction$se.fit))
    )]

    expanded[, .(
      predicted = mean(conditional_prediction),
      lower = mean(conditional_lower),
      upper = mean(conditional_upper),
      random_effect_cells = uniqueN(intact_cell)
    ), by = curve_row]
  }))

  output <- merge(curve_grid, averaged, by = "curve_row", sort = FALSE)
  output[, curve_row := NULL]
  output
}

# ------------------------------------------------------------------------------
# 4.5 GAM 1 and Figure 1: FLII curves by distance class
# ------------------------------------------------------------------------------

formula_gam1 <- as.formula(paste(
  "dissimilarity ~",
  paste(c(
    "dist_label",
    "s(flii_diff, by = dist_label, bs = 'cr', k = 5)",
    adjustment_terms,
    random_effect_term
  ), collapse = " + ")
))

gam1_distance <- bam(
  formula = formula_gam1,
  data = gam_data,
  family = betar(link = "logit"),
  method = "fREML",
  select = TRUE,
  discrete = TRUE
)

prediction_gam1 <- rbindlist(lapply(levels(gam_data$dist_label), function(d) {
  subset_data <- gam_data[dist_label == d]
  newdata <- data.table(
    flii_diff = seq(
      min(subset_data$flii_diff),
      max(subset_data$flii_diff),
      length.out = 250L
    ),
    dist_label = factor(d, levels = levels(gam_data$dist_label))
  )
  profiles <- unique(subset_data[, .(intact_cell)])
  predict_random_effect_average(gam1_distance, newdata, profiles)
}))

set.seed(123)
plot_points_gam1 <- gam_data[
  , .SD[sample.int(.N, min(.N, 3000L))],
  by = dist_label
]

figure1 <- ggplot() +
  geom_point(
    data = plot_points_gam1,
    aes(flii_diff, dissimilarity),
    alpha = 0.12,
    size = 0.55,
    colour = "grey35"
  ) +
  geom_ribbon(
    data = prediction_gam1,
    aes(flii_diff, ymin = lower, ymax = upper),
    fill = "#3182BD",
    alpha = 0.20
  ) +
  geom_line(
    data = prediction_gam1,
    aes(flii_diff, predicted),
    colour = "#08519C",
    linewidth = 1.05
  ) +
  facet_wrap(~dist_label, ncol = 3, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  labs(
    title = "FLII difference and compositional dissimilarity by distance",
    subtitle = "Curves average conditional predictions across intact-pixel random intercepts",
    x = "FLII difference",
    y = "Predicted compositional dissimilarity"
  )

# ------------------------------------------------------------------------------
# 4.6 GAM 2 and Figure 2: FLII curves by distance and elevation class
# ------------------------------------------------------------------------------

gam_data[, distance_elevation := interaction(
  dist_label,
  elev_band,
  sep = " | ",
  drop = TRUE
)]

formula_gam2 <- as.formula(paste(
  "dissimilarity ~",
  paste(c(
    "distance_elevation",
    "s(flii_diff, by = distance_elevation, bs = 'cr', k = 5)",
    adjustment_terms,
    random_effect_term
  ), collapse = " + ")
))

gam2_distance_elevation <- bam(
  formula = formula_gam2,
  data = gam_data,
  family = betar(link = "logit"),
  method = "fREML",
  select = TRUE,
  discrete = TRUE
)

observed_strata <- unique(gam_data[, .(dist_label, elev_band)])

prediction_gam2 <- rbindlist(lapply(seq_len(nrow(observed_strata)), function(i) {
  d <- observed_strata$dist_label[i]
  e <- observed_strata$elev_band[i]
  subset_data <- gam_data[dist_label == d & elev_band == e]

  newdata <- data.table(
    flii_diff = seq(
      min(subset_data$flii_diff),
      max(subset_data$flii_diff),
      length.out = 200L
    ),
    dist_label = factor(d, levels = levels(gam_data$dist_label)),
    elev_band = factor(e, levels = levels(gam_data$elev_band))
  )

  newdata[, distance_elevation := interaction(
    dist_label,
    elev_band,
    sep = " | ",
    drop = TRUE
  )]
  newdata[, distance_elevation := factor(
    distance_elevation,
    levels = levels(gam_data$distance_elevation)
  )]

  profiles <- unique(subset_data[, .(intact_cell)])
  predict_random_effect_average(gam2_distance_elevation, newdata, profiles)
}), fill = TRUE)

figure2 <- ggplot() +
  geom_point(
    data = gam_data,
    aes(flii_diff, dissimilarity, colour = elev_band),
    alpha = 0.08,
    size = 0.5
  ) +
  geom_ribbon(
    data = prediction_gam2,
    aes(
      flii_diff,
      ymin = lower,
      ymax = upper,
      fill = elev_band,
      group = elev_band
    ),
    alpha = 0.14,
    colour = NA
  ) +
  geom_line(
    data = prediction_gam2,
    aes(flii_diff, predicted, colour = elev_band, group = elev_band),
    linewidth = 1.0
  ) +
  facet_wrap(~dist_label, ncol = 3, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "FLII difference and compositional dissimilarity by elevation and distance",
    subtitle = "Elevation bands come from matching; curves average intact-pixel random effects",
    x = "FLII difference",
    y = "Predicted compositional dissimilarity",
    colour = "Elevation",
    fill = "Elevation"
  )

# ------------------------------------------------------------------------------
# 4.7 GAM 3 and Figure 3: global FLII curve with elevation-distance interaction
# ------------------------------------------------------------------------------

# dist_label and elev_band are fixed factors, and their interaction allows the
# background dissimilarity to differ among all distance-by-elevation strata.
# The single s(flii_diff) term estimates the requested overall FLII relationship.
# The random intercept accounts for repeated pairs from each intact pixel.
formula_gam3 <- as.formula(paste(
  "dissimilarity ~",
  paste(c(
    "s(flii_diff, bs = 'cr', k = 5)",
    "dist_label * elev_band",
    adjustment_terms,
    random_effect_term
  ), collapse = " + ")
))

gam3_global_interaction <- bam(
  formula = formula_gam3,
  data = gam_data,
  family = betar(link = "logit"),
  method = "fREML",
  select = TRUE,
  discrete = TRUE
)

prediction_gam3 <- data.table(
  flii_diff = seq(
    min(gam_data$flii_diff),
    max(gam_data$flii_diff),
    length.out = 300L
  )
)

# Average across the observed distance, elevation and intact-cell profiles.
# Thus the plotted curve includes both fixed stratum effects and the fitted
# intact-cell random effects instead of selecting one arbitrary reference cell.
global_profiles <- unique(gam_data[, .(
  intact_cell,
  dist_label,
  elev_band
)])

prediction_gam3 <- predict_random_effect_average(
  gam3_global_interaction,
  prediction_gam3,
  global_profiles
)

set.seed(124)
plot_points_gam3 <- gam_data[
  sample.int(.N, min(.N, 15000L))
]

figure3 <- ggplot() +
  geom_point(
    data = plot_points_gam3,
    aes(flii_diff, dissimilarity),
    alpha = 0.08,
    size = 0.5,
    colour = "grey40"
  ) +
  geom_ribbon(
    data = prediction_gam3,
    aes(flii_diff, ymin = lower, ymax = upper),
    fill = "#31A354",
    alpha = 0.20
  ) +
  geom_line(
    data = prediction_gam3,
    aes(flii_diff, predicted),
    colour = "#006D2C",
    linewidth = 1.15
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Global FLII difference and compositional dissimilarity",
    subtitle = "Averaged across observed elevation-distance strata and intact-pixel random effects",
    x = "FLII difference",
    y = "Predicted compositional dissimilarity"
  )

# ------------------------------------------------------------------------------
# 4.8 Print diagnostics and save GAM models, predictions, and figures
# ------------------------------------------------------------------------------

cat("\nGAM 1 summary:\n")
print(summary(gam1_distance))
cat("\nGAM 2 summary:\n")
print(summary(gam2_distance_elevation))
cat("\nGAM 3 summary:\n")
print(summary(gam3_global_interaction))

saveRDS(gam1_distance, file.path(
  gam_output_path,
  "GAM1_FLII_by_distance_random_intact_cell.rds"
))
saveRDS(gam2_distance_elevation, file.path(
  gam_output_path,
  "GAM2_FLII_by_distance_elevation_random_intact_cell.rds"
))
saveRDS(gam3_global_interaction, file.path(
  gam_output_path,
  "GAM3_global_FLII_elevation_distance_interaction.rds"
))

fwrite(prediction_gam1, file.path(
  gam_output_path,
  "GAM1_predictions_by_distance.csv"
))
fwrite(prediction_gam2, file.path(
  gam_output_path,
  "GAM2_predictions_by_distance_elevation.csv"
))
fwrite(prediction_gam3, file.path(
  gam_output_path,
  "GAM3_global_predictions.csv"
))

ggsave(
  file.path(gam_output_path, "Figure1_FLII_by_distance.png"),
  figure1,
  width = 12,
  height = 8,
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(gam_output_path, "Figure2_FLII_by_distance_elevation.png"),
  figure2,
  width = 12,
  height = 8,
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(gam_output_path, "Figure3_global_FLII_interaction.png"),
  figure3,
  width = 8,
  height = 6,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(gam_output_path, "Figure1_FLII_by_distance.pdf"),
  figure1,
  width = 12,
  height = 8
)
ggsave(
  file.path(gam_output_path, "Figure2_FLII_by_distance_elevation.pdf"),
  figure2,
  width = 12,
  height = 8
)
ggsave(
  file.path(gam_output_path, "Figure3_global_FLII_interaction.pdf"),
  figure3,
  width = 8,
  height = 6
)

print(figure1)
print(figure2)
print(figure3)
