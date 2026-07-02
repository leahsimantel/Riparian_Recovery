################################################################################
###    CH01 Script D: Partial Dependence Plots 
###    40-meter spacing 
###    03/17/2026

# ==============================================================
# Libraries
# ==============================================================
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(randomForest)
  library(cluster)      
  library(RColorBrewer)
  library(corrplot)
  library(car)
  library(forcats)
  library(lme4)
  library(purrr)
  library(patchwork)
  library(readr)
  library(stringr)
})
# --------

#          -------------------------------------------                         #
#########     LOAD ENVIRONMENTAL BLOCKING RESULTS (SKIP RECOMPUTE)   ###########
####          -------------------------------------------                   ####
# Load latest file. 
rds_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD5CWD30_twi20260307_095447.rds"

envblock_results <- readRDS(rds_path)

# Re-name by YSF (critical for downstream code)
envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}

envblock_results_named <- envblock_by_ysf(envblock_results)

cat("Loaded environmental blocking results from:\n", rds_path, "\n")
# ------------------------------------------------------------------------------
cwd_zscore_var <- "cwd_5yr_zscore_08"

# ==============================================================
# 3. RANDOM FOREST (RF) — SEV_NUM + LEAVE-ONE-OUT MODEL COMPARISON
# =====================================================================
# LOCO (Leave-One-Covariate-Out) predictor sets: ---------------------
rf_vars_full <- c(
  "sev_num",
  "twi",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  "cwd_30YrAvg_TT",
  cwd_zscore_var
)

rf_var_sets <- c(
  list(full = rf_vars_full),
  lapply(rf_vars_full, function(v) setdiff(rf_vars_full, v))
)

names(rf_var_sets) <- c(
  "full",
  paste0("drop_", rf_vars_full)
)

# Create RF function -----------------------------------------------
rf_run_one_ysf_loco <- function(envblock_result,
                                rf_vars,
                                trees   = ntree,
                                seed_rf = 23) {
  
  stopifnot(is.list(envblock_result),
            "data_with_folds_env" %in% names(envblock_result))
  
  target_ysf <- envblock_result$YSF
  df <- envblock_result$data_with_folds_env
  k_folds <- length(levels(df$fold))
  
  message("  → YSF = ", target_ysf,
          " | predictors = ", paste(rf_vars, collapse = ", "))
  
  fold_results <- vector("list", length = k_folds)
  names(fold_results) <- paste0("Fold_", levels(df$fold))
  
  for (kf in levels(df$fold)) {
    
    message("    • Fold ", kf, "/", k_folds)
    
    train_data <- df %>% dplyr::filter(fold != kf)
    test_data  <- df %>% dplyr::filter(fold == kf)
    
    # ---- numeric median imputation (train median) ----
    for (col in rf_vars) {
      med <- suppressWarnings(stats::median(train_data[[col]], na.rm = TRUE))
      if (is.na(med)) med <- 0
      train_data[[col]][is.na(train_data[[col]])] <- med
      test_data[[col]][is.na(test_data[[col]])]   <- med
    }
    
    set.seed(seed_rf + as.integer(as.character(kf)))
    
    rf_formula <- as.formula(
      paste("delta_ndvi_min ~", paste(rf_vars, collapse = " + "))
    )
    
    rf_model <- randomForest::randomForest(
      rf_formula,
      data  = train_data,
      ntree = trees
    )
    
    preds <- predict(rf_model, newdata = test_data)
    obs   <- test_data$delta_ndvi_min
    N     <- length(obs)
    
    mse  <- if (N) mean((preds - obs)^2) else NA_real_
    rmse <- if (is.na(mse)) NA_real_ else sqrt(mse)
    
    denom <- sum((obs - mean(obs))^2)
    r2    <- if (denom == 0) NA_real_ else 1 - sum((obs - preds)^2) / denom
    
    fold_results[[paste0("Fold_", kf)]] <- list(
      model = rf_model,
      preds = preds,
      obs   = obs,
      N     = N,
      rmse  = rmse,
      r2    = r2
    )
  }
  
  # ---- Aggregate performance ----------------------------------------
  rmse_vals <- sapply(fold_results, `[[`, "rmse")
  r2_vals   <- sapply(fold_results, `[[`, "r2")
  Ns        <- sapply(fold_results, `[[`, "N")
  
  perf <- tibble::tibble(
    YSF = target_ysf,
    Mean_RMSE_unweighted     = mean(rmse_vals, na.rm = TRUE),
    Mean_PseudoR2_unweighted = mean(r2_vals,   na.rm = TRUE),
    Mean_RMSE_weighted       = stats::weighted.mean(rmse_vals, Ns, na.rm = TRUE),
    Mean_PseudoR2_weighted   = stats::weighted.mean(r2_vals,   Ns, na.rm = TRUE)
  )
  
  list(
    YSF          = target_ysf,
    rf_vars     = rf_vars,
    fold_models = fold_results,
    perf        = perf
  )
}

# ---- RUN LOCO RF MODELS ---------------------------------------------
envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}
envblock_results_named <- envblock_by_ysf(envblock_results)

#rf_loco_results <- purrr::imap(
rf_var_sets,
function(vars, model_name) {
  
  message("\n==============================")
  message("LOCO RF model: ", model_name)
  message("==============================")
  
  res <- lapply(
    as.character(ysf_set),
    function(ysf_key) {
      rf_run_one_ysf_loco(
        envblock_result = envblock_results_named[[ysf_key]],
        rf_vars = vars,
        trees   = ntree,
        seed_rf = 168
      )
    }
  )
  
  names(res) <- as.character(ysf_set)
  
  list(
    model_name = model_name,
    rf_vars    = vars,
    results    = res
  )
}
)

# ---- PERFORMANCE COMPARISON (ΔRMSE, Δpseudo-R²) ----------------------

perf_loco_all <- purrr::map_dfr(
  rf_loco_results,
  function(m) {
    dplyr::bind_rows(lapply(m$results, `[[`, "perf")) %>%
      dplyr::mutate(
        model  = m$model_name,
        n_vars = length(m$rf_vars)
      )
  }
)

perf_loco_summary <- perf_loco_all %>%
  
  # Join full-model weighted metrics (baseline per YSF)
  dplyr::left_join(
    perf_loco_all %>%
      dplyr::filter(model == "full") %>%
      dplyr::select(
        YSF,
        RMSE_full = Mean_RMSE_weighted,
        R2_full   = Mean_PseudoR2_weighted
      ),
    by = "YSF"
  ) %>%
  
  # Compute deltas relative to full model
  dplyr::mutate(
    delta_RMSE = Mean_RMSE_weighted     - RMSE_full,
    delta_R2   = Mean_PseudoR2_weighted - R2_full
  ) %>%
  
  # Rank models WITHIN each YSF
  dplyr::group_by(YSF) %>%
  dplyr::arrange(
    Mean_RMSE_weighted,                 # primary: lower RMSE is better
    dplyr::desc(Mean_PseudoR2_weighted),# secondary: higher R² is better
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    rank_within_YSF = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

print(perf_loco_summary, n = Inf)
# -------------------------------------------------------
# SAVE LOCO RF RESULTS (timestamped + rolling "latest") ------------------------
out_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch01_RF_LOCO"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

loco_save <- list(
  meta = list(
    created      = Sys.time(),
    ysf_set      = ysf_set,
    predictors   = rf_vars_full,
    ntree        = ntree,
    seed_rf      = 168,
    severity_mode = "sev_num",
    description  = "LOCO Random Forest results with env-block CV"
  ),
  rf_loco_results   = rf_loco_results,
  perf_loco_all     = perf_loco_all,
  perf_loco_summary = perf_loco_summary
)

# ---- Timestamped archive 
#saveRDS(
#  loco_save,
#  file = file.path(out_dir, paste0("rf_loco_results_", stamp, ".rds"))
#)

# ---- Rolling latest version 
saveRDS(
  loco_save,
  file = file.path(out_dir, "rf_loco_results_latest.rds")
)

cat("✅ LOCO RF results saved to:\n", out_dir, "\n")
# ====================================================

# LOAD LOCO RF RESULTS (skip re-running LOCO models)
# ==============================================================

rds_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch01_RF_LOCO/rf_loco_results_cwd5_cwd30_twi.rds"
loco_save <- readRDS(rds_path)

rf_loco_results   <- loco_save$rf_loco_results
perf_loco_all     <- loco_save$perf_loco_all
perf_loco_summary <- loco_save$perf_loco_summary

str(perf_loco_summary)

# Double check - need envbloc_results_named
stopifnot(exists("envblock_results_named"))
# Reconstruct YSF set from loaded env-blocking results
ysf_set <- as.numeric(names(envblock_results_named))
stopifnot(length(ysf_set) > 0)




# ========================================================
# Single-variable PDPs for LOCO models -----------------------------------------
# ========================================================
# ---- controls ----
pdp_grid_n   <- 60
pdp_q_lo     <- 0.01
pdp_q_hi     <- 0.99
pdp_row_cap  <- 2500
pdp_seed     <- 999

cwd_col  <- "cwd_5yr_zscore_08"
sev_mode <- "sev_num"   # LOCO RFs ALWAYS use numeric severity

# ---- predictors to plot PDP for ----
base_vars_effect <- c(
  "tmeanz_JJA",
  "twi",
  "veg_climate_index_08",
  "pptz_JJA",
  "swez_Apr",
  cwd_col
)

# Include severity PDP only when numeric
vars_effect <- c("sev_num", base_vars_effect)

# ---- FUNCTION: robust grid helper (numeric) ----
.master_grid_cont <- function(x, n, q_lo, q_hi) {
  xnum <- suppressWarnings(as.numeric(x))
  if (!length(xnum) || length(unique(xnum[is.finite(xnum)])) < 2) return(numeric(0))
  ql <- suppressWarnings(stats::quantile(xnum, probs = q_lo, na.rm = TRUE))
  qh <- suppressWarnings(stats::quantile(xnum, probs = q_hi, na.rm = TRUE))
  if (!is.finite(ql) || !is.finite(qh) || ql == qh) return(numeric(0))
  seq(ql, qh, length.out = n)
}
# ---- FUNCTION: get fold training data (mirror RF prep, dynamic factor vs numeric) ----
.get_train_data_for_fold <- function(ysf_key, fold_name,
                                     rf_numeric_cols,
                                     rf_factor_cols = character(0),
                                     row_cap = pdp_row_cap,
                                     seed = pdp_seed) {
  kf_val <- sub("^Fold_", "", fold_name)
  
  df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env
  df2    <- df_all %>% dplyr::mutate(.fold_chr = as.character(fold))
  
  train_raw <- df2 %>% dplyr::filter(.fold_chr != kf_val)
  if (!nrow(train_raw)) train_raw <- df2
  
  train <- train_raw %>%
    dplyr::select(-.fold_chr) %>%
    sf::st_drop_geometry()
  
  # factor handling
  for (col in intersect(rf_factor_cols, names(train))) {
    train[[col]] <- factor(as.character(train[[col]]))
  }
  
  # numeric imputation (median from full training fold)
  for (col in intersect(rf_numeric_cols, names(train))) {
    train[[col]] <- suppressWarnings(as.numeric(train[[col]]))
    med <- suppressWarnings(stats::median(train[[col]], na.rm = TRUE))
    if (is.na(med)) med <- 0
    train[[col]][is.na(train[[col]])] <- med
  }
  
  # speed cap (after imputation)
  if (is.finite(row_cap) && nrow(train) > row_cap) {
    set.seed(seed)
    train <- train[sample.int(nrow(train), row_cap), , drop = FALSE]
  }
  
  train
}
# ---- FUNCTION: brute-force PD at a grid value ----
.predict_grid_cont <- function(model, train_df, var, grid_x) {
  if (is.null(model) || !is.data.frame(train_df) || !length(grid_x)) return(tibble())
  purrr::map_dfr(grid_x, function(xx) {
    newdat <- train_df
    newdat[[var]] <- xx
    preds <- try(predict(model, newdata = newdat), silent = TRUE)
    if (inherits(preds, "try-error")) return(tibble())
    pnum <- as.numeric(preds)
    tibble::tibble(
      x         = xx,
      fold_mean = mean(pnum, na.rm = TRUE),
      fold_sd   = stats::sd(pnum, na.rm = TRUE)
    )
  })
}
# ---- FUNCTION: build CV PDPs across YSFs ----
build_cv_pdp <- function(ysf_set,
                         rf_results_envblock,
                         envblock_results_named,
                         cwd_col,
                         vars_effect,
                         pdp_grid_n,
                         sev_mode = "sev_num") {
  
  rf_factor_cols  <- character(0)
  rf_numeric_cols <- c(
    "sev_num",
    "twi",
    "pptz_JJA",
    "tmeanz_JJA",
    "swez_Apr",
    "veg_climate_index_08",
    cwd_col
  )
  
  set.seed(pdp_seed)
  
  purrr::map_dfr(as.character(ysf_set), function(ysf_key) {
    
    rf_obj <- rf_results_envblock[[as.character(ysf_key)]]
    if (is.null(rf_obj) || is.null(rf_obj$fold_models)) return(tibble())
    
    df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
      sf::st_drop_geometry()
    
    # grids per variable
    master_grids <- lapply(vars_effect, function(vv) {
      if (!vv %in% names(df_all)) return(NULL)
      .master_grid_cont(df_all[[vv]], pdp_grid_n, pdp_q_lo, pdp_q_hi)
    })
    names(master_grids) <- vars_effect
    
    fold_names <- names(rf_obj$fold_models)
    
    fold_curves <- purrr::map_dfr(fold_names, function(fnm) {
      
      rf_model <- rf_obj$fold_models[[fnm]]$model
      if (is.null(rf_model)) return(tibble())
      
      tr_df <- .get_train_data_for_fold(
        ysf_key,
        fnm,
        rf_numeric_cols,
        rf_factor_cols,
        pdp_row_cap,
        pdp_seed
      )
      
      vars_here <- intersect(vars_effect, names(tr_df))
      
      purrr::map_dfr(vars_here, function(vv) {
        grid_x <- master_grids[[vv]]
        if (is.null(grid_x) || !length(grid_x)) return(tibble())
        out <- .predict_grid_cont(rf_model, tr_df, vv, grid_x)
        if (!nrow(out)) return(tibble())
        out %>%
          dplyr::mutate(
            YSF      = as.numeric(ysf_key),
            variable = vv,
            Fold     = fnm
          )
      })
    })
    
    if (!nrow(fold_curves)) return(tibble())
    
    fold_curves %>%
      dplyr::group_by(YSF, variable, x) %>%
      dplyr::summarise(
        y_mean = mean(fold_mean, na.rm = TRUE),
        # total uncertainty = between-fold + within-fold variance
        sd_between = if (dplyr::n() > 1) stats::sd(fold_mean, na.rm = TRUE) else 0,
        mean_within_var = {
          v <- fold_sd^2
          v <- v[is.finite(v)]
          if (length(v)) mean(v) else 0
        },
        y_sd = sqrt(pmax(0, sd_between^2 + mean_within_var)),
        n_folds = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(variable, x)
  })
}
# ========================================================

## Build PDPs; create pdp_loco_cv ========================
message("\n==============================")
message("Building LOCO single-variable PDPs")
message("==============================")

pdp_loco_results <- purrr::imap(
  rf_loco_results,
  function(m, model_name) {
    
    message("\n--- PDPs for LOCO model: ", model_name, " ---")
    
    build_cv_pdp(
      ysf_set                = ysf_set,
      rf_results_envblock    = m$results,
      envblock_results_named = envblock_results_named,
      cwd_col                = cwd_col,
      vars_effect            = intersect(vars_effect, m$rf_vars),
      pdp_grid_n             = pdp_grid_n,
      sev_mode               = "sev_num"
    ) %>%
      dplyr::mutate(model = model_name)
  }
)
pdp_loco_cv <- dplyr::bind_rows(pdp_loco_results)
stopifnot(nrow(pdp_loco_cv) > 0)

# ---- factor prep 
model_levels <- names(rf_loco_results)

pdp_loco_cv <- pdp_loco_cv %>%
  dplyr::mutate(
    YSF      = factor(YSF, levels = sort(unique(YSF), decreasing = TRUE)),
    model    = factor(model, levels = model_levels),
    variable = factor(variable, levels = vars_effect)
  )

# ---- QA
pdp_loco_cv %>%
  dplyr::group_by(model, YSF, variable) %>%
  dplyr::summarise(
    n_x = dplyr::n(),
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE),
    min_sd = min(y_sd, na.rm = TRUE),
    max_sd = max(y_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# PLOTTING =========
vars_to_plot <- c(
  #"tmeanz_JJA",
  "swez_Apr"
  #"pptz_JJA"
  #"cwd_5yr_zscore_08",
  #"veg_climate_index_08",
  #"sev_num",
  #"twi"
)
# ========================================================
vars_to_plot <- intersect(vars_to_plot, unique(pdp_loco_cv$variable))

pdp_plot_df <- pdp_loco_cv %>%
  dplyr::filter(variable %in% vars_to_plot) %>%
  dplyr::mutate(variable = factor(variable, levels = vars_to_plot))

var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "twi"                  = "Topographic Wetness Index",
  "veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

ggplot(
  pdp_plot_df,
  aes(x = x, y = y_mean, color = YSF, group = YSF)
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = y_mean - sd_between,
        ymax = y_mean + sd_between,
        fill = YSF),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  facet_grid(model ~ variable,
             scales = "free_x",
             labeller = as_labeller(var_labels)) +
  labs(
    title = "LOCO Partial Dependence of Climate Predictors on ΔNDVI",
    subtitle = "Comparison across Leave-One-Covariate-Out models",
    x = NULL,
    y = "Model-Predicted ΔNDVI"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 13),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )
################################################################################





###############################################################################
# 2-VARIABLE PDPs FOR LOCO RF RESULTS
# - Fire severity interactions across multiple YSF groups
# - Uses LOCO RF object structure
# - Uses the top-ranked LOCO model separately for each YSF by default
###############################################################################

### Load libraries ============
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(scales)
})
## =================================

# =============================================================================
# REQUIRED OBJECTS
# =============================================================================
stopifnot(exists("rf_loco_results"))
stopifnot(exists("perf_loco_summary"))
stopifnot(exists("envblock_results_named"))
stopifnot(exists("ysf_set"))

# =============================================================================
# USER CONTROLS
# =============================================================================
pdp2_target_ysf_vec <- c(5, 10, 15, 20)
pdp2_grid_res   <- 35
pdp2_q_lo       <- 0.01
pdp2_q_hi       <- 0.99
pdp2_row_cap    <- 3500
pdp2_seed       <- 999

# Use the top-ranked LOCO model separately for each YSF
use_top_ranked_model <- TRUE

# If FALSE above, provide one named model per YSF
# manual_model_by_ysf <- c(
#   "5"  = "drop_twi",
#   "10" = "drop_twi",
#   "15" = "drop_",
#   "20" = "drop_twi"
# )

# Variables to plot as rows (in this order)
severity_partner_vars_panel <- c(
  "tmeanz_JJA",
  "pptz_JJA",
  "swez_Apr",
  "cwd_30yrAvg_TT"
)

# Variable Labels
var_labels <- c(
  "sev_num"           = "Fire Severity",
  "tmeanz_JJA"        = "Temperature",
  "pptz_JJA"          = "Precipitation",
  "swez_Apr"          = "SWE",
  "cwd_5yr_zscore_08" = "CWD5",
  "cwd_30yrAvg_TT"    = "CWD30",
  "twi"               = "TWI"
)

model_labels <- c(
  full                    = "Full model",
  drop_sev_num            = "Drop: Fire severity",
  drop_twi                = "Drop: TWI",
  drop_pptz_JJA           = "Drop: Summer precip",
  drop_tmeanz_JJA         = "Drop: Summer temp",
  drop_swez_Apr           = "Drop: April SWE",
  drop_cwd_30yrAvg_TT     = "Drop: CWD 30Yr Normal",
  drop_cwd_5yr_zscore_08  = "Drop: Postfire CWD"
)

# =============================================================================
# CHOOSE MODEL FOR EACH YSF
# =============================================================================
get_model_for_ysf <- function(ysf_val) {
  if (use_top_ranked_model) {
    model_name_use <- perf_loco_summary %>%
      dplyr::filter(YSF == ysf_val, rank_within_YSF == 1) %>%
      dplyr::pull(model)
    
    if (length(model_name_use) == 0) {
      stop("No top-ranked LOCO model found for YSF = ", ysf_val)
    }
    
    model_name_use[1]
  } else {
    manual_model_by_ysf[[as.character(ysf_val)]]
  }
}

model_by_ysf <- setNames(
  vapply(pdp2_target_ysf_vec, get_model_for_ysf, character(1)),
  pdp2_target_ysf_vec
)

cat("Models by YSF:\n")
print(model_by_ysf)
# =============================================================================
# FUNCTIONS
# =============================================================================

.master_grid_cont2 <- function(x, n = 35, q_lo = 0.01, q_hi = 0.99) {
  xnum <- suppressWarnings(as.numeric(x))
  xnum <- xnum[is.finite(xnum)]
  
  if (length(xnum) < 2 || length(unique(xnum)) < 2) return(numeric(0))
  
  ql <- suppressWarnings(stats::quantile(xnum, probs = q_lo, na.rm = TRUE))
  qh <- suppressWarnings(stats::quantile(xnum, probs = q_hi, na.rm = TRUE))
  
  if (!is.finite(ql) || !is.finite(qh) || ql == qh) return(numeric(0))
  
  seq(ql, qh, length.out = n)
}

get_train_df_for_fold_loco <- function(ysf_key,
                                       fold_name,
                                       rf_vars_model,
                                       row_cap = Inf,
                                       seed = 999) {
  
  kf_val <- sub("^Fold_", "", fold_name)
  
  df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
    dplyr::mutate(.fold_chr = as.character(fold)) %>%
    dplyr::filter(!is.na(delta_ndvi_min)) %>%
    sf::st_drop_geometry()
  
  train <- df_all %>%
    dplyr::filter(.fold_chr != kf_val) %>%
    dplyr::select(-.fold_chr)
  
  if (!nrow(train)) {
    train <- df_all %>% dplyr::select(-.fold_chr)
  }
  
  keep_cols <- intersect(c("delta_ndvi_min", rf_vars_model), names(train))
  train <- train[, keep_cols, drop = FALSE]
  
  # numeric coercion + median imputation to mirror RF prep
  for (cc in intersect(rf_vars_model, names(train))) {
    train[[cc]] <- suppressWarnings(as.numeric(train[[cc]]))
    med <- suppressWarnings(stats::median(train[[cc]], na.rm = TRUE))
    if (is.na(med)) med <- 0
    train[[cc]][is.na(train[[cc]])] <- med
  }
  
  if (is.finite(row_cap) && nrow(train) > row_cap) {
    set.seed(seed)
    train <- train[sample.int(nrow(train), row_cap), , drop = FALSE]
  }
  
  train
}

build_cv_pdp2_one_pair_loco <- function(ysf_val,
                                        rf_results_this_model,
                                        envblock_results_named,
                                        rf_vars_model,
                                        var1,
                                        var2,
                                        x_grid,
                                        y_grid,
                                        row_cap = Inf,
                                        seed = 999) {
  
  ysf_key <- as.character(ysf_val)
  rf_obj  <- rf_results_this_model[[ysf_key]]
  
  if (is.null(rf_obj) || is.null(rf_obj$fold_models)) {
    return(tibble())
  }
  
  if (!length(x_grid) || !length(y_grid)) {
    return(tibble())
  }
  
  grid_df <- tidyr::expand_grid(x = x_grid, y = y_grid)
  fold_names <- names(rf_obj$fold_models)
  
  fold_surfaces <- purrr::map_dfr(fold_names, function(fnm) {
    
    mdl <- rf_obj$fold_models[[fnm]]$model
    if (is.null(mdl)) return(tibble())
    
    train_df <- get_train_df_for_fold_loco(
      ysf_key       = ysf_key,
      fold_name     = fnm,
      rf_vars_model = rf_vars_model,
      row_cap       = row_cap,
      seed          = seed
    )
    
    if (!(var1 %in% names(train_df)) || !(var2 %in% names(train_df))) {
      return(tibble())
    }
    
    pd_fold <- purrr::pmap_dfr(grid_df, function(x, y) {
      newdat <- train_df
      newdat[[var1]] <- x
      newdat[[var2]] <- y
      
      preds <- try(predict(mdl, newdata = newdat), silent = TRUE)
      if (inherits(preds, "try-error")) return(tibble())
      
      tibble(
        x = x,
        y = y,
        yhat = mean(as.numeric(preds), na.rm = TRUE)
      )
    })
    
    if (!nrow(pd_fold)) return(tibble())
    
    pd_fold %>%
      dplyr::mutate(
        YSF  = as.numeric(ysf_key),
        Fold = fnm
      )
  })
  
  if (!nrow(fold_surfaces)) {
    return(tibble())
  }
  
  fold_surfaces %>%
    dplyr::group_by(YSF, x, y) %>%
    dplyr::summarise(
      yhat_mean = mean(yhat, na.rm = TRUE),
      yhat_sd   = if (dplyr::n() > 1) stats::sd(yhat, na.rm = TRUE) else 0,
      n_folds   = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      var1 = var1,
      var2 = var2,
      .before = 1
    )
}
# =============================================================================
# BUILD SHARED MASTER GRIDS 
# (dynamic, based on User Controls)
# =============================================================================
all_panel_data <- purrr::map_dfr(
  as.character(pdp2_target_ysf_vec),
  function(ysf_key) {
    envblock_results_named[[ysf_key]]$data_with_folds_env %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(!is.na(delta_ndvi_min)) %>%
      dplyr::mutate(YSF = as.numeric(ysf_key))
  }
)

# shared x-grid for sev_num across all selected YSFs
x_grid_shared <- .master_grid_cont2(
  all_panel_data$sev_num,
  n    = pdp2_grid_res,
  q_lo = pdp2_q_lo,
  q_hi = pdp2_q_hi
)

# shared y-grid for each row variable across all selected YSFs
y_grid_lookup <- setNames(
  lapply(severity_partner_vars_panel, function(v2) {
    if (!v2 %in% names(all_panel_data)) {
      stop("Variable ", v2, " not found in pooled panel data.")
    }
    
    gy <- .master_grid_cont2(
      all_panel_data[[v2]],
      n    = pdp2_grid_res,
      q_lo = pdp2_q_lo,
      q_hi = pdp2_q_hi
    )
    
    if (!length(gy)) {
      stop("Could not build shared y-grid for variable: ", v2)
    }
    
    gy
  }),
  severity_partner_vars_panel
)

cat("Shared x grid range:", min(x_grid_shared), "to", max(x_grid_shared), "\n")
cat("Shared x grid length:", length(x_grid_shared), "\n")

print(
  tibble::tibble(
    var2  = names(y_grid_lookup),
    y_min = purrr::map_dbl(y_grid_lookup, min),
    y_max = purrr::map_dbl(y_grid_lookup, max),
    n_y   = purrr::map_int(y_grid_lookup, length)
  )
)
# =============================================================================
# BUILD ALL PDPs FOR ALL YSF x ROW VARIABLES
# =============================================================================
pdp2_all <- purrr::map_dfr(
  pdp2_target_ysf_vec,
  function(ysf_now) {
    
    model_name_now <- model_by_ysf[[as.character(ysf_now)]]
    model_obj_now  <- rf_loco_results[[model_name_now]]
    
    if (is.null(model_obj_now)) {
      stop("Model ", model_name_now, " not found in rf_loco_results for YSF ", ysf_now)
    }
    
    rf_vars_model_now         <- model_obj_now$rf_vars
    rf_results_this_model_now <- model_obj_now$results
    
    if (!"sev_num" %in% rf_vars_model_now) {
      stop("Model ", model_name_now, " for YSF ", ysf_now, " does not include sev_num.")
    }
    
    missing_vars <- setdiff(severity_partner_vars_panel, rf_vars_model_now)
    if (length(missing_vars) > 0) {
      stop(
        "Model ", model_name_now, " for YSF ", ysf_now,
        " is missing required variables: ",
        paste(missing_vars, collapse = ", ")
      )
    }
    
    purrr::map_dfr(
      severity_partner_vars_panel,
      function(v2) {
        message("Building 2D PDP: sev_num x ", v2, " | YSF = ", ysf_now,
                " | model = ", model_name_now)
        
        build_cv_pdp2_one_pair_loco(
          ysf_val                = ysf_now,
          rf_results_this_model  = rf_results_this_model_now,
          envblock_results_named = envblock_results_named,
          rf_vars_model          = rf_vars_model_now,
          var1                   = "sev_num",
          var2                   = v2,
          x_grid                 = x_grid_shared,
          y_grid                 = y_grid_lookup[[v2]],
          row_cap                = pdp2_row_cap,
          seed                   = pdp2_seed
        ) %>%
          dplyr::mutate(model_name = model_name_now)
      }
    )
  }
)

if (!nrow(pdp2_all)) {
  stop("No multi-YSF 2-variable PDP output was generated.")
}

pdp2_all <- pdp2_all %>%
  dplyr::mutate(
    yhat_mean = signif(yhat_mean, 2),
    yhat_sd   = signif(yhat_sd, 2)
  )

## Check: outputs should match for each YSF
pdp2_all %>%
  dplyr::filter(var1 == "sev_num", var2 %in% severity_partner_vars_panel) %>%
  dplyr::group_by(var2, YSF) %>%
  dplyr::summarise(
    x_min = min(x, na.rm = TRUE),
    x_max = max(x, na.rm = TRUE),
    n_x   = dplyr::n_distinct(x),
    y_min = min(y, na.rm = TRUE),
    y_max = max(y, na.rm = TRUE),
    n_y   = dplyr::n_distinct(y),
    .groups = "drop"
  ) %>%
  print(n = Inf)
# =============================================================================
# UNIFORM AXIS LIMITS FOR PANEL COMPARISON
# - x limits shared across all plots (sev_num)
# - y limits shared within each row variable
# =============================================================================

# shared x-axis limits across all panels
x_limits_shared <- range(pdp2_all$x, na.rm = TRUE)

# shared y-axis limits for each row variable
y_limits_by_var <- pdp2_all %>%
  dplyr::group_by(var2) %>%
  dplyr::summarise(
    y_min = min(y, na.rm = TRUE),
    y_max = max(y, na.rm = TRUE),
    .groups = "drop"
  )

# turn into a named list for easy lookup inside the plotting loop
y_limits_lookup <- setNames(
  lapply(seq_len(nrow(y_limits_by_var)), function(i) {
    c(y_limits_by_var$y_min[i], y_limits_by_var$y_max[i])
  }),
  y_limits_by_var$var2
)

cat("Shared x limits:", paste(x_limits_shared, collapse = " to "), "\n")
print(y_limits_by_var)
# =============================================================================
# SHARED LEGEND / COLOR SCALE
# =============================================================================
step_fill  <- 0.025
step_label <- 0.05

global_min <- min(pdp2_all$yhat_mean, na.rm = TRUE)
global_max <- max(pdp2_all$yhat_mean, na.rm = TRUE)

global_limits <- c(
  floor(global_min / step_fill) * step_fill,
  ceiling(global_max / step_fill) * step_fill
)

global_limits <- c(
  global_limits[1] - step_fill,
  global_limits[2] + step_fill
)

global_limits <- round(global_limits, 10)

global_breaks <- seq(global_limits[1], global_limits[2], by = step_fill)
global_breaks <- unique(round(global_breaks, 10))

legend_breaks <- seq(
  ceiling(global_limits[1] / step_label) * step_label,
  floor(global_limits[2] / step_label) * step_label,
  by = step_label
)
legend_breaks <- unique(round(legend_breaks, 10))

fill_scale_pdp_panel <- scale_fill_steps2(
  low = "#4B3621",
  mid = "beige",
  high = "darkgreen",
  midpoint = 0,
  limits = global_limits,
  breaks = legend_breaks,
  labels = format(legend_breaks, trim = TRUE, scientific = FALSE),
  oob = scales::squish,
  name = "Predicted ΔNDVI",
  guide = guide_coloursteps(
    title.position = "top",
    title.hjust = 0.5,
    even.steps = TRUE,
    show.limits = FALSE,
    title.theme = element_text(
      family = "Times New Roman",
      face   = "bold",
      size   = 18,
      margin = ggplot2::margin(b = 10)
    ),
    label.theme = element_text(
      family = "Times New Roman",
      size   = 12
    ),
    barheight = unit(5.0, "cm"),
    barwidth  = unit(0.75, "cm"),
    ticks = TRUE,
    label.position = "right"
  )
)

cat("Data range:", global_min, "to", global_max, "\n")
cat("Legend limits:", global_limits[1], "to", global_limits[2], "\n")
# =============================================================================
# PLOT FUNCTION
# =============================================================================
.make_pdp2_plot_single <- function(df_pair,
                                   fill_scale = fill_scale_pdp_panel,
                                   var_labels = var_labels,
                                   show_legend = FALSE,
                                   show_xlab   = TRUE,
                                   show_ylab   = TRUE,
                                   col_title   = NULL,
                                   x_limits    = NULL,
                                   y_limits    = NULL) {
  
  stopifnot(nrow(df_pair) > 0)
  
  v1 <- df_pair$var1[1]
  v2 <- df_pair$var2[1]
  
  v1_lab <- if (!is.null(var_labels[[v1]])) var_labels[[v1]] else v1
  v2_lab <- if (!is.null(var_labels[[v2]])) var_labels[[v2]] else v2
  
  # reserve equal title / label space across plots
  plot_title_use <- if (!is.null(col_title)) col_title else " "
  x_lab_use      <- if (show_xlab) v1_lab else " "
  y_lab_use      <- if (show_ylab) v2_lab else " "
  
  p <- ggplot(df_pair, aes(x = x, y = y, fill = yhat_mean)) +
    geom_raster() +
    geom_contour(
      aes(z = yhat_mean),
      color = "black",
      alpha = 0.35
    ) +
    fill_scale +
    labs(
      title = plot_title_use,
      x = x_lab_use,
      y = y_lab_use
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman"),
      
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      
      axis.title.x = element_text(
        size = 13,
        face = "bold",
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 13,
        face = "bold",
        margin = ggplot2::margin(r = 8)
      ),
      
      axis.text = element_text(size = 11),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = ggplot2::margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"),
      panel.grid = element_blank(),
      aspect.ratio = 1
    )
  
  if (!is.null(x_limits) || !is.null(y_limits)) {
    p <- p + coord_cartesian(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE
    )
  }
  
  p
}
# =============================================================================
# BUILD 16 INDIVIDUAL PANELS (4 rows x 4 columns)
# rows = variables, columns = YSF
# =============================================================================
row_vars <- severity_partner_vars_panel
col_ysf  <- pdp2_target_ysf_vec

plots_16 <- vector("list", length = length(row_vars) * length(col_ysf))
k <- 1

for (r in seq_along(row_vars)) {
  for (c in seq_along(col_ysf)) {
    
    v2      <- row_vars[r]
    ysf_now <- col_ysf[c]
    
    df_pair <- pdp2_all %>%
      dplyr::filter(YSF == ysf_now, var1 == "sev_num", var2 == v2)
    
    if (!nrow(df_pair)) {
      stop("Missing PDP surface for YSF = ", ysf_now, " and variable = ", v2)
    }
    
    show_xlab <- r == length(row_vars)
    show_ylab <- c == 1
    col_title <- if (r == 1) paste0(ysf_now, " YSF") else NULL
    
    plots_16[[k]] <- .make_pdp2_plot_single(
      df_pair,
      fill_scale  = fill_scale_pdp_panel,
      var_labels  = var_labels,
      show_legend = (k == 1),
      show_xlab   = show_xlab,
      show_ylab   = show_ylab,
      col_title   = col_title,
      x_limits    = x_limits_shared,
      y_limits    = y_limits_lookup[[v2]]
    )
    
    k <- k + 1
  }
}

shared_legend <- cowplot::get_legend(plots_16[[1]])
plots_16_noleg <- lapply(plots_16, function(p) p + theme(legend.position = "none"))

# ----- build the panel grid
panel_grid <- cowplot::plot_grid(
  plotlist = plots_16_noleg,
  ncol = length(col_ysf),
  align = "hv",
  axis = "tblr"
)

# =============================================================================
# TITLE + SUBTITLE
# =============================================================================
title_grob <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Partial Dependence Plots",
    fontfamily = "Times New Roman",
    fontface   = "bold",
    size       = 24,
    hjust      = 0.5
  )

subtitle_grob <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Fire severity interactions across post-fire intervals",
    fontfamily = "Times New Roman",
    fontface   = "plain",
    size       = 18,
    hjust      = 0.5,
    y          = 0.9
  )
# =============================================================================
# FINAL PANEL
# =============================================================================
pdp2_panel_all_ysf <- cowplot::plot_grid(
  cowplot::plot_grid(
    title_grob,
    subtitle_grob,
    cowplot::plot_grid(
      panel_grid,
      shared_legend,
      ncol = 2,
      rel_widths = c(1, 0.12)
    ),
    ncol = 1,
    rel_heights = c(0.08, 0.06, 1)
  ),
  ncol = 1
)

print(pdp2_panel_all_ysf)
# ==============================================================================





