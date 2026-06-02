null_result_df <- do.call(rbind, lapply(seq_along(null_result), function(s) {
  do.call(rbind, lapply(seq_along(null_result[[s]]), function(i) {
    cbind(null_result[[s]][[i]]$pair.df, site = s, iter = i)
  }))
}))


null_result_df$label <- paste0(null_result_df$sp1,"_",null_result_df$sp2) 
null_result_df$async <- null_result_df$obs.cor*-1
null_result_df$pop_stab_diff <- abs(log(null_result_df$pop_stab.x)-log(null_result_df$pop_stab.y))
null_result_df$cover_diff <- abs(null_result_df$avg_cover.x/100-null_result_df$avg_cover.y/100)
null_result_df$plotID_full <- paste0(null_result_df$site_code,"_",null_result_df$block,"_",null_result_df$plot,"_",null_result_df$subplot)
null_result_df$label_for_merge <- paste0(null_result_df$sp1,"_",null_result_df$sp2,"_",null_result_df$plotID_full) 
all_SES_pair_long_clim$label_for_merge <- paste0(all_SES_pair_long_clim$sp1,"_",all_SES_pair_long_clim$sp2,"_",all_SES_pair_long_clim$plotID_full) 

null_result_df <- null_result_df %>%
  left_join(all_SES_pair_long_clim[,c("native_sim","phy.cor",
                                      "label_for_merge",
                                      "MAT_v2","MAP_v2","TEMP_VAR_v2","MAP_VAR_v2")],
            by=join_by(label_for_merge))


null_result_df <- subset(null_result_df,phy.cor <= 400)

null_result_df$cover_diff_adj <- (null_result_df$cover_diff * (unique(table(null_result_df$iter)) - 1) + 0.5) / unique(table(null_result_df$iter))


### model ----

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
run_null_models <- function(null_df, model_specs, alpha = 0.05) {
  
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
      
      # --- extract conditional fixed effects ---
      cond_tidy <- tidy(fit, effects = "fixed", component = "cond") %>%
        mutate(component = "cond")
      
      # --- extract dispersion fixed effects manually ---
      disp_coefs <- summary(fit)$coefficients$disp
      disp_tidy <- tibble(
        term      = rownames(disp_coefs),
        estimate  = disp_coefs[, "Estimate"],
        std.error = disp_coefs[, "Std. Error"],
        statistic = disp_coefs[, "z value"],
        p.value   = disp_coefs[, "Pr(>|z|)"],
        component = "disp"
      )
      
      bind_rows(cond_tidy, disp_tidy) %>%
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
    group_by(model, component, term) %>%
    summarise(
      n_iter        = n(),
      prop_sig      = mean(sig, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate   = sd(estimate, na.rm = TRUE),
      mean_se       = mean(std.error, na.rm = TRUE),
      .groups       = "drop"
    ) %>%
    arrange(model, component, desc(prop_sig))
  
  list(
    raw_coefs    = results,
    sig_summary  = sig_summary,
    conv_summary = conv_summary
  )
}

null_out     <- run_null_models(null_result_df, model_specs)



###

for (iter in 1:n_null) {
  message(iter)
  
  null_result<-list()
  for (i in 1:length(unique(Cover_subset$ID))) {
    message(i)
    
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
      pivot_wider(id_cols = year,
                  names_from  = Accepted_SPNAME,
                  values_from = max_cover,
                  values_fill = 0)
    
    no_var_sp <- names(which(apply(comp_df[,-1],2,FUN=sd) == 0))
    comp_df   <- comp_df %>% select(-all_of(no_var_sp))
    
    plot_info <- unique(Cover_df_subset[, c("site_code", "block", "plot", "subplot", "trt")])
    plot_info$year  <- length(unique(comp_df$year))
    plot_info$gamma <- ncol(comp_df[, -1])+length(no_var_sp)
    plot_info$alpha <- mean(vegan::specnumber(comp_df[, -1]))
    plot_info$no_var_sp <- length(no_var_sp)
    
    sp_info_base <- unique(Cover_df_subset[, c("Accepted_SPNAME", "Family", "Genus", "n")])
    sp_info_base <- sp_info_base %>% filter(!Accepted_SPNAME %in% no_var_sp)
    
    # --- null expectations ---
    sp_cols      <- comp_df[, -1]
    shuf_idx     <- sample(ncol(sp_cols))
    sp_cols_null <- sp_cols[, shuf_idx]
    colnames(sp_cols_null) <- colnames(sp_cols)
      
    sp_info_null <- sp_info_base
    sp_info_null$avg_cover <- colMeans(sp_cols_null)
    sp_info_null$var_cover <- apply(sp_cols_null, 2, var)
    sp_info_null$pop_stab  <- sp_info_null$avg_cover / sqrt(sp_info_null$var_cover)
      
    null_result[[i]] <- calc_temporal_metrics(
        sp_cols_null, comp_df$year,
        plot_info = plot_info, sp_info = sp_info_null
      )
  }
  
  null_result_df <- do.call(rbind,null_result)

  null_result_df$label <- paste0(null_result_df$sp1,"_",null_result_df$sp2) 
  null_result_df$async <- null_result_df$obs.cor*-1
  null_result_df$pop_stab_diff <- abs(log(null_result_df$pop_stab.x)-log(null_result_df$pop_stab.y))
  null_result_df$cover_diff <- abs(null_result_df$avg_cover.x/100-null_result_df$avg_cover.y/100)
  null_result_df$plotID_full <- paste0(null_result_df$site_code,"_",null_result_df$block,"_",null_result_df$plot,"_",null_result_df$subplot)
  null_result_df$label_for_merge <- paste0(null_result_df$sp1,"_",null_result_df$sp2,"_",null_result_df$plotID_full) 
  all_SES_pair_long_clim$label_for_merge <- paste0(all_SES_pair_long_clim$sp1,"_",all_SES_pair_long_clim$sp2,"_",all_SES_pair_long_clim$plotID_full) 
  
  null_result_df <- null_result_df %>%
    left_join(all_SES_pair_long_clim[,c("native_sim","phy.cor",
                                        "label_for_merge",
                                        "MAT_v2","MAP_v2","TEMP_VAR_v2","MAP_VAR_v2")],
              by=join_by(label_for_merge))
  
  
  null_result_df <- subset(null_result_df,phy.cor <= 400)
  
  null_result_df$cover_diff_adj <- (null_result_df$cover_diff * (unique(nrow(null_result_df)) - 1) + 0.5) / unique(nrow(null_result_df))
  
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
  
  
}
