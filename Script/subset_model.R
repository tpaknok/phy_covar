source("Script/random_matching.R")

matched <- random_matching(all_SES_pair_long_clim, n_iter = 100)

library(glmmTMB)
library(dplyr)
library(broom.mixed)

# --- define models as a named list of formulas + family ---
model_specs <- list(
  cover_diff = list(
    formula     = cover_diff_adj ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor)) + scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) + (1|plotID_full),
    family      = beta_family(link = "logit")
  ),
  sync = list(
    formula     = async ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor))*scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) + (1|plotID_full),
    family      = gaussian()
  ),
  pop = list(
    formula     = pop_stab_diff ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor))*scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor))*scale(native_sim) + scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) + scale(year) + scale(log(gamma)) + (1|plotID_full),
    family      = gaussian()
  )
)

# --- run all models on all null iterations ---
run_null_models <- function(matched_df, model_specs, alpha = 0.05) {
  
  iters <- split(matched_df, matched_df$iter)
  
  # --- run models across all iterations ---
  results <- lapply(names(iters), function(it) {
    message("iter ", it)
    iter_df <- iters[[it]]
    
    lapply(names(model_specs), function(mod_name) {
      spec <- model_specs[[mod_name]]
      
      fit <- tryCatch(
        glmmTMB(formula     = spec$formula,
                dispformula = spec$dispformula,
                data        = iter_df,
                family      = spec$family,
                control     = glmmTMBControl(parallel = 6),
                REML        = TRUE),
        error = function(e) {
          message("  model ", mod_name, " iter ", it, " failed: ", e$message)
          return(NULL)
        }
      )
      
      if (is.null(fit)) return(NULL)
      
      converged <- fit$fit$convergence == 0
      
      tidy(fit, effects = "fixed") %>%
        mutate(iter      = as.integer(it),
               model     = mod_name,
               sig       = p.value < alpha,
               converged = converged)
    }) %>% bind_rows()
    
  }) %>% bind_rows()
  
  # --- convergence diagnostics ---
  conv_summary <- results %>%
    distinct(model, iter, converged) %>%
    group_by(model) %>%
    summarise(
      total       = n(),
      n_converged = sum(converged),
      n_dropped   = sum(!converged),
      pct_kept    = round(100 * mean(converged), 1),
      .groups     = "drop"
    )
  
  message("\n--- Convergence summary ---")
  print(conv_summary)
  
  # --- significance summary (converged iterations only) ---
  sig_summary <- results %>%
    filter(!is.na(p.value),
           converged) %>%
    group_by(model, term) %>%
    summarise(
      n_iter        = n(),
      prop_sig      = mean(sig, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate   = sd(estimate, na.rm = TRUE),
      mean_se       = mean(std.error, na.rm = TRUE),
      .groups       = "drop"
    ) %>%
    arrange(model, desc(prop_sig))
  
  list(
    raw_coefs    = results,
    sig_summary  = sig_summary,
    conv_summary = conv_summary
  )
}

# --- usage ---
null_out     <- run_null_models(matched, model_specs)

raw_coefs    <- null_out$raw_coefs      # full coefficient table (all iters)
sig_summary  <- null_out$sig_summary    # proportion significant per term (converged only)
conv_summary <- null_out$conv_summary   # convergence diagnostics per model
