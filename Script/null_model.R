# ============================================================
#  NULL MODEL PIPELINE — PARALLELIZED
# ============================================================

library(furrr)
library(future)

extract_glmmtmb_coefs <- function(fit, mod_name) {
  
  cond_tidy <- broom.mixed::tidy(fit, effects = "fixed", component = "cond") %>%
    mutate(component = "cond")
  
  disp_coefs <- summary(fit)$coefficients$disp
  disp_tidy  <- tibble(
    term      = rownames(disp_coefs),
    estimate  = disp_coefs[, "Estimate"],
    std.error = disp_coefs[, "Std. Error"],
    statistic = disp_coefs[, "z value"],
    p.value   = disp_coefs[, "Pr(>|z|)"],
    component = "disp"
  )
  
  bind_rows(cond_tidy, disp_tidy) %>%
    mutate(model = mod_name)
}

compare_observed_to_null <- function(observed_model_results, null_raw_coefs) {
  # observed_model_results: tidy output from your real model (broom::tidy or similar)
  # null_raw_coefs: the raw_coefs tibble from run_null_models()
  
  null_raw_coefs %>%
    # Join on whatever identifies a coefficient uniquely in your setup
    inner_join(
      observed_model_results %>%
        select(model, component, term, obs_estimate = estimate),
      by = c("model", "component", "term")
    ) %>%
    group_by(model, component, term, obs_estimate) %>%
    summarise(
      n_null_total    = n(),
      n_gt_obs        = sum(estimate > obs_estimate,  na.rm = TRUE),
      n_lt_obs        = sum(estimate < obs_estimate,  na.rm = TRUE),
      n_eq_obs        = sum(estimate == obs_estimate, na.rm = TRUE),
      prop_gt_obs     = n_gt_obs / n_null_total,
      prop_lt_obs     = n_lt_obs / n_null_total,
      # Two-sided empirical p-value
      p_empirical     = 2 * min(prop_gt_obs, prop_lt_obs),
      mean_null       = mean(estimate, na.rm = TRUE),
      sd_null         = sd(estimate,   na.rm = TRUE),
      # How many SDs is the observed from the null mean?
      z_vs_null       = (obs_estimate - mean_null) / sd_null,
      .groups         = "drop"
    ) %>%
    mutate(
      direction = case_when(
        prop_gt_obs > prop_lt_obs ~ "obs_below_null",  # most null > observed
        prop_lt_obs > prop_gt_obs ~ "obs_above_null",  # most null < observed
        TRUE                      ~ "symmetric"
      )
    ) %>%
    arrange(model, component, p_empirical)
}
# ------------------------------------------------------------
# 1. DATA GENERATION FUNCTION (one randomized iteration)
# ------------------------------------------------------------

get_null_df <- function() {
  
  null_result <- list()
  
  for (i in 1:length(unique(Cover_subset$ID))) {
    
    Cover_df_subset <- Cover_subset %>%
      ungroup() %>%
      filter(ID == unique(Cover_subset$ID)[[i]]) %>%
      mutate(Epithet = word(Accepted_SPNAME, 2, 2, "_")) %>%
      filter(Epithet != "sp.")
    
    Cover_df_subset$Genus <- word(Cover_df_subset$Accepted_SPNAME, 1, sep = "_")
    
    foo_df <- Cover_df_subset %>%
      group_by(Accepted_SPNAME) %>%
      dplyr::count() %>%
      as.data.frame()
    
    Cover_df_subset <- Cover_df_subset %>%
      left_join(foo_df, join_by(Accepted_SPNAME))
    
    comp_df <- Cover_df_subset %>%
      select(year, Accepted_SPNAME, max_cover) %>%
      group_by(year, Accepted_SPNAME) %>%
      dplyr::summarize(max_cover = sum(max_cover), .groups = "drop") %>%
      pivot_wider(id_cols     = year,
                  names_from  = Accepted_SPNAME,
                  values_from = max_cover,
                  values_fill = 0)
    
    no_var_sp <- names(which(apply(comp_df[, -1], 2, FUN = sd) == 0))
    comp_df   <- comp_df %>% select(-all_of(no_var_sp))
    
    plot_info           <- unique(Cover_df_subset[, c("site_code", "block", "plot", "subplot", "trt")])
    plot_info$year      <- length(unique(comp_df$year))
    plot_info$gamma     <- ncol(comp_df[, -1]) + length(no_var_sp)
    plot_info$alpha     <- mean(vegan::specnumber(comp_df[, -1]))
    plot_info$no_var_sp <- length(no_var_sp)
    
    sp_info_base <- unique(Cover_df_subset[, c("Accepted_SPNAME", "Family", "Genus", "n")])
    sp_info_base <- sp_info_base %>% filter(!Accepted_SPNAME %in% no_var_sp)
    
    # --- shuffle species columns (null randomization) ---
    sp_cols                <- comp_df[, -1]
    shuf_idx               <- sample(ncol(sp_cols))
    sp_cols_null           <- sp_cols[, shuf_idx]
    colnames(sp_cols_null) <- colnames(sp_cols)
    
    # --- drop zero-variance columns after shuffling ---
    no_var_sp_null <- names(which(apply(sp_cols_null, 2, FUN = sd) == 0))
    if (length(no_var_sp_null) > 0) {
      sp_cols_null <- sp_cols_null[, !colnames(sp_cols_null) %in% no_var_sp_null]
    }
    
    sp_info_null <- sp_info_base %>%
      filter(!Accepted_SPNAME %in% no_var_sp_null)
    
    sp_info_null$avg_cover <- colMeans(sp_cols_null)
    sp_info_null$var_cover <- apply(sp_cols_null, 2, var)
    sp_info_null$pop_stab  <- ifelse(
      sp_info_null$var_cover > 0,
      sp_info_null$avg_cover / sqrt(sp_info_null$var_cover),
      NA_real_
    )
    
    null_result[[i]] <- calc_temporal_metrics(
      sp_cols_null, comp_df$year,
      plot_info = plot_info,
      sp_info   = sp_info_null
    )
  }
  
  # --- bind pair_df from each plot ---
  null_result_df <- do.call(rbind, lapply(null_result, `[[`, "pair.df"))
  
  null_result_df$label           <- paste0(null_result_df$sp1, "_", null_result_df$sp2)
  null_result_df$async           <- null_result_df$obs.cor * -1
  null_result_df$pop_stab_diff   <- abs(log(null_result_df$pop_stab.x) - log(null_result_df$pop_stab.y))
  null_result_df$cover_diff      <- abs(null_result_df$avg_cover.x / 100 - null_result_df$avg_cover.y / 100)
  null_result_df$plotID_full     <- paste0(null_result_df$site_code, "_",
                                           null_result_df$block,     "_",
                                           null_result_df$plot,      "_",
                                           null_result_df$subplot)
  null_result_df$label_for_merge <- paste0(null_result_df$sp1, "_",
                                           null_result_df$sp2, "_",
                                           null_result_df$plotID_full)
  
  all_SES_pair_long_clim$label_for_merge <- paste0(all_SES_pair_long_clim$sp1, "_",
                                                   all_SES_pair_long_clim$sp2, "_",
                                                   all_SES_pair_long_clim$plotID_full)
  
  null_result_df <- null_result_df %>%
    left_join(
      all_SES_pair_long_clim[, c("native_sim", "phy.cor",
                                 "label_for_merge",
                                 "MAT_v2", "MAP_v2", "TEMP_VAR_v2", "MAP_VAR_v2")],
      by = join_by(label_for_merge)
    )
  
  null_result_df <- subset(null_result_df, phy.cor <= 400)
  
  n_rows <- nrow(null_result_df)
  null_result_df$cover_diff_adj <- (null_result_df$cover_diff * (n_rows - 1) + 0.5) / n_rows
  
  return(null_result_df)
}


# ------------------------------------------------------------
# 2. MODEL SPECIFICATIONS
# ------------------------------------------------------------

model_specs <- list(
  
  cover_diff = list(
    formula     = cover_diff_adj ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor)) + scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) + (1 | plotID_full),
    family      = beta_family(link = "logit")
  ),
  
  sync = list(
    formula     = async ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor)) * scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) + (1 | plotID_full),
    family      = gaussian()
  ),
  
  pop = list(
    formula     = pop_stab_diff ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) +
      (scale(sqrt(phy.cor)) * scale(native_sim) || plotID_full),
    dispformula = ~ scale(sqrt(phy.cor)) * scale(native_sim) +
      scale(MAP_VAR_v2) + scale(MAT_v2) + scale(TEMP_VAR_v2) + scale(MAP_v2) +
      scale(year) + scale(log(gamma)) + (1 | plotID_full),
    family      = gaussian()
  )
)


# ------------------------------------------------------------
# 3. SINGLE ITERATION WORKER (called inside furrr)
# ------------------------------------------------------------

run_one_iteration <- function(attempt_id, get_null_df, model_specs, alpha) {
  
  # --- generate data ---
  iter_df <- tryCatch(
    get_null_df(),
    error = function(e) {
      message("  [attempt ", attempt_id, "] Data generation failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(iter_df)) return(NULL)
  
  # --- fit all three models ---
  iter_coefs <- lapply(names(model_specs), function(mod_name) {
    
    spec <- model_specs[[mod_name]]
    
    fit <- tryCatch(
      glmmTMB(formula     = spec$formula,
              dispformula = spec$dispformula,
              data        = iter_df,
              family      = spec$family,
              # no internal parallelism — we are already parallelizing across iterations
              control     = glmmTMBControl(parallel = 1),
              REML        = TRUE),
      error = function(e) {
        message("  [attempt ", attempt_id, "] model [", mod_name, "] failed: ", e$message)
        NULL
      }
    )
    
    if (is.null(fit)) {
      message("  [attempt ", attempt_id, "] model [", mod_name, "] returned NULL")
      return(NULL)
    }
    
    if (fit$fit$convergence != 0) {
      message("  [attempt ", attempt_id, "] model [", mod_name, "] non-zero convergence code")
      return(NULL)
    }
    
    if (isTRUE(fit$sdr$pdHess == FALSE)) {
      message("  [attempt ", attempt_id, "] model [", mod_name, "] non-positive-definite Hessian")
      return(NULL)
    }
    
    cond_tidy  <- broom.mixed::tidy(fit, effects = "fixed", component = "cond") %>%
      mutate(component = "cond")
    
    disp_coefs <- summary(fit)$coefficients$disp
    disp_tidy  <- tibble(
      term      = rownames(disp_coefs),
      estimate  = disp_coefs[, "Estimate"],
      std.error = disp_coefs[, "Std. Error"],
      statistic = disp_coefs[, "z value"],
      p.value   = disp_coefs[, "Pr(>|z|)"],
      component = "disp"
    )
    
    bind_rows(cond_tidy, disp_tidy) %>%
      mutate(model = mod_name)
  })
  
  # --- all three must converge ---
  if (any(sapply(iter_coefs, is.null))) return(NULL)
  
  bind_rows(iter_coefs) %>%
    mutate(attempt = attempt_id,
           sig     = p.value < alpha)
}


# ------------------------------------------------------------
# 4. NULL MODEL RUNNER (parallelized with furrr)
# ------------------------------------------------------------

run_null_models <- function(
    get_null_df,
    model_specs,
    n_null   = 100,
    alpha    = 0.05,
    max_try  = n_null * 10,
    n_cores  = 1
) {
  
  converged_results <- list()
  n_converged       <- 0
  n_tried           <- 0
  
  if (n_cores == 1) {
    
    # -------------------------------------------------------
    # SEQUENTIAL — simple while loop, stops the moment
    # n_null converged iterations are collected
    # -------------------------------------------------------
    message("Running sequentially (n_cores = 1)")
    
    while (n_converged < n_null && n_tried < max_try) {
      
      n_tried <- n_tried + 1
      message(sprintf("Attempt %d  |  Converged so far: %d / %d",
                      n_tried, n_converged, n_null))
      
      res <- run_one_iteration(n_tried, get_null_df, model_specs, alpha)
      
      if (!is.null(res)) {
        n_converged <- n_converged + 1
        converged_results[[n_converged]] <- res %>% mutate(iter = n_converged)
        message(sprintf("  >>> Accepted iteration %d / %d", n_converged, n_null))
      }
    }
    
  } else {
    
    # -------------------------------------------------------
    # PARALLEL — batched furrr, sends n_cores * 2 attempts
    # at a time across workers
    # -------------------------------------------------------
    plan(multisession, workers = n_cores)
    on.exit(plan(sequential), add = TRUE)
    message(sprintf("Running with %d parallel workers", n_cores))
    
    batch_size <- n_cores * 2
    
    while (n_converged < n_null && n_tried < max_try) {
      
      this_batch  <- min(batch_size, max_try - n_tried)
      attempt_ids <- (n_tried + 1):(n_tried + this_batch)
      
      message(sprintf("Launching batch of %d attempts (tried: %d, converged: %d / %d)",
                      this_batch, n_tried, n_converged, n_null))
      
      batch_results <- furrr::future_map(
        attempt_ids,
        ~ run_one_iteration(.x, get_null_df, model_specs, alpha),
        .options = furrr_options(seed = TRUE)
      )
      
      n_tried <- n_tried + this_batch
      
      for (res in batch_results) {
        if (!is.null(res) && n_converged < n_null) {
          n_converged <- n_converged + 1
          converged_results[[n_converged]] <- res %>% mutate(iter = n_converged)
          message(sprintf("  >>> Accepted iteration %d / %d", n_converged, n_null))
        }
      }
    }
  }
  
  # --- warn if cap hit ---
  if (n_converged < n_null) {
    warning(sprintf(
      "Only %d / %d converged iterations after %d attempts. Consider increasing max_try.",
      n_converged, n_null, max_try
    ))
  }
  
  if (n_converged == 0) {
    warning("No iterations converged. Returning empty results.")
    return(list(
      raw_coefs       = tibble(),
      sig_summary     = tibble(),
      attempt_summary = tibble(n_null = n_null, n_converged = 0L,
                               n_attempted = n_tried, pct_kept = 0)
    ))
  }
  
  results <- bind_rows(converged_results)
  
  sig_summary <- results %>%
    filter(!is.na(p.value)) %>%
    group_by(model, component, term) %>%
    summarise(
      n_iter        = n(),
      prop_sig      = mean(sig,       na.rm = TRUE),
      mean_estimate = mean(estimate,  na.rm = TRUE),
      sd_estimate   = sd(estimate,    na.rm = TRUE),
      mean_se       = mean(std.error, na.rm = TRUE),
      .groups       = "drop"
    ) %>%
    arrange(model, component, desc(prop_sig))
  
  attempt_summary <- tibble(
    n_null      = n_null,
    n_converged = n_converged,
    n_attempted = n_tried,
    pct_kept    = round(100 * n_converged / n_tried, 1)
  )
  
  message("\n--- Attempt summary ---")
  print(attempt_summary)
  
  list(
    raw_coefs       = results,
    sig_summary     = sig_summary,
    attempt_summary = attempt_summary
  )
}

# ------------------------------------------------------------
# 5. RUN
# ------------------------------------------------------------
set.seed(999)

start <- Sys.time()
null_out <- run_null_models(
  get_null_df = get_null_df,
  model_specs = model_specs,
  n_null      = 100,
  alpha       = 0.05,
  max_try     = 150,
  n_cores     = 6
)
end <- Sys.time()
end-start

saveRDS(null_out, "null_out_100iter.rds")

null_out$attempt_summary
null_out$raw_coefs

observed_tidy <- bind_rows(
  
  # --- model 1: if m1_cover_diff is a plain glm/glmmTMB ---
  broom.mixed::tidy(m1_cover_diff, effects = "fixed", component = "cond") %>%   # drops component = "cond" override
    mutate(model = "cover_diff"),                            # keep component as-is from tidy()
  
  # --- model 2 ---
  broom.mixed::tidy(m1_sync, effects = "fixed", component = "cond") %>%
    mutate(model = "sync"),                                   # same — do NOT overwrite component
  broom.mixed::tidy(m1_pop, effects = "fixed", component = "cond") %>%
    mutate(model = "pop")                                   # same — do NOT overwrite component
  
)

observed_tidy <- bind_rows(
  extract_glmmtmb_coefs(m1_cover_diff, "cover_diff"),
  extract_glmmtmb_coefs(m1_sync,       "sync"),
  extract_glmmtmb_coefs(m1_pop,       "pop")
)

comparison <- unique(compare_observed_to_null(observed_tidy, null_out$raw_coefs))
