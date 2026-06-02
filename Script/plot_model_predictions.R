library(tidyverse)
library(ggpubr)
library(glmmTMB)

plot_model_predictions <- function(
    model,
    data,
    y_label           = "Response",
    x_label           = "Phylogenetic distance",
    family            = c("gaussian", "gamma", "beta", "tweedie"),
    x_var             = "phy.cor",
    facet_var         = NULL,
    facet_label       = "Historical sympatry",
    n_x               = 100,
    n_facet           = 3,
    linetype_mean     = "dashed",
    linetype_var      = "solid",
    ci_level          = 0.95,
    ncol              = 1,
    nrow              = 2,
    dispersion_metric = c("disp", "sd", "cv"),
    log_transform = F
) {
  
  family            <- match.arg(family)
  dispersion_metric <- match.arg(dispersion_metric)
  z_mult            <- qnorm(1 - (1 - ci_level) / 2)
  
  # ---- 1. Extract variables from glmmTMB model -------------------------
  fixed_vars <- all.vars(nobars(model$modelInfo$allForm$formula))
  disp_vars  <- all.vars(nobars(model$modelInfo$allForm$dispformula))
  all_vars   <- unique(c(fixed_vars, disp_vars))
  
  response_var   <- fixed_vars[1]
  predictor_vars <- setdiff(all_vars, response_var)
  
  re_vars <- tryCatch({
    unique(unlist(lapply(findbars(model$modelInfo$allForm$formula), function(b) {
      all.vars(b[[3]])
    })))
  }, error = function(e) character(0))
  
  message("Response:   ", response_var)
  message("Predictors: ", paste(predictor_vars, collapse = ", "))
  message("RE vars:    ", paste(re_vars, collapse = ", "))
  
  # ---- 2. Build prediction grid ----------------------------------------
  grid_list <- lapply(setNames(predictor_vars, predictor_vars), function(v) {
    if (v == x_var) {
      seq(round(min(data[[v]], na.rm = TRUE),2), round(max(data[[v]], na.rm = TRUE),2), length.out = n_x)
    } else if (!is.null(facet_var) && v == facet_var) {
      seq(round(min(data[[v]], na.rm = TRUE),2), round(max(data[[v]], na.rm = TRUE),2), length.out = n_facet)
    } else if (v %in% re_vars || !v %in% names(data)) {
      NA
    } else if (is.numeric(data[[v]])) {
      mean(data[[v]], na.rm = TRUE)
    } else {
      lvl <- names(sort(table(data[[v]]), decreasing = TRUE))[1]
      if (is.null(lvl) || length(lvl) == 0) NA else lvl
    }
  })
  
  grid_list <- Filter(Negate(is.null), grid_list)
  newdat    <- expand.grid(grid_list, stringsAsFactors = FALSE)
  
  # ---- 3. Mean predictions ---------------------------------------------
  mean_pred <- cbind(
    newdat,
    predict(model, newdat, type = "response", se.fit = TRUE, re.form = NA)
  ) %>%
    dplyr::rename(mean = fit)
  
  # ---- 4. Dispersion predictions ---------------------------------------
  disp_pred <- predict(model, newdata = newdat, type = "disp", re.form = NA)
  
  # Tweedie power parameter
  tweedie_p <- if (family == "tweedie") {
    tryCatch(
      model$modelInfo$family$power,
      error = function(e) stop("Could not extract Tweedie power parameter.")
    )
  } else NULL
  
  # Residual SD — confounded with mean for non-Gaussian
  mean_pred$residual_sd <- switch(family,
                                  gaussian = disp_pred,
                                  gamma    = mean_pred$mean / sqrt(disp_pred),
                                  beta     = sqrt(mean_pred$mean * (1 - mean_pred$mean) / (1 + disp_pred)),
                                  tweedie  = sqrt(disp_pred * mean_pred$mean ^ tweedie_p)
  )
  
  # CV — partially removes mean contribution
  mean_pred$residual_cv <- switch(family,
                                  gaussian = mean_pred$residual_sd / mean_pred$mean,
                                  gamma    = 1 / sqrt(disp_pred),
                                  beta     = 1 / sqrt(1 + disp_pred),
                                  tweedie  = sqrt(disp_pred * mean_pred$mean ^ (tweedie_p - 2))
  )
  
  # Dispersion — pure dispformula effect, mean-independent
  mean_pred$disp_plot_var <- switch(dispersion_metric,
                                    disp = switch(family,
                                                  gaussian = disp_pred,
                                                  gamma    = 1 / sqrt(disp_pred),
                                                  beta     = 1 / disp_pred,
                                                  tweedie  = disp_pred
                                    ),
                                    sd = mean_pred$residual_sd,
                                    cv = mean_pred$residual_cv
  )
  
  disp_y_label <- switch(dispersion_metric,
                         disp = switch(family,
                                       gaussian = "Dispersion",
                                       gamma    = "Dispersion",
                                       beta     = "Dispersion",
                                       tweedie  = "Dispersion"
                         ),
                         sd = "Model residual SD",
                         cv = "Coefficient of variation"
  )
  
  # ---- 5. Ribbon bounds — mean uncertainty only ------------------------
  ribbon_lo_expr <- switch(family,
                           gaussian = if (log_transform & family == "gaussian") expression(exp(mean - z_mult * se.fit)) else expression(mean - z_mult * se.fit),
                           gamma    = expression(exp(log(mean) - z_mult * se.fit / mean)),
                           beta     = expression(plogis(qlogis(mean) - z_mult * se.fit / (mean * (1 - mean)))),
                           tweedie  = expression(exp(log(mean) - z_mult * se.fit / mean))
  )
  
  ribbon_hi_expr <- switch(family,
                           gaussian = if (log_transform & family == "gaussian") expression(exp(mean + z_mult * se.fit)) else expression(mean + z_mult * se.fit),
                           gamma    = expression(exp(log(mean) + z_mult * se.fit / mean)),
                           beta     = expression(plogis(qlogis(mean) + z_mult * se.fit / (mean * (1 - mean)))),
                           tweedie  = expression(exp(log(mean) + z_mult * se.fit / mean))
  )
  
  mean_pred <- mean_pred %>%
    mutate(
      ribbon_lo = eval(ribbon_lo_expr),
      ribbon_hi = eval(ribbon_hi_expr),
      mean = if (log_transform & family == "gaussian") exp(mean) else mean
    )
  
  
  # ---- 6. Optional facet -----------------------------------------------
  if (!is.null(facet_var)) {
    mean_pred$facet_label <- paste0(facet_label," = ", mean_pred[[facet_var]])
  }
  
  facet_layer <- if (!is.null(facet_var)) facet_wrap(~facet_label) else NULL
  group_var   <- if (is.null(facet_var)) 1 else mean_pred[[facet_var]]
  
  # ---- 7. Mean / CI plot -----------------------------------------------
  p_mean <- ggplot(mean_pred, aes(x = .data[[x_var]], y = mean, group = group_var)) +
    geom_ribbon(aes(ymin = ribbon_lo, ymax = ribbon_hi), alpha = 0.2) +
    geom_line(linetype = linetype_mean) +
    facet_layer +
    labs(y = y_label, x = x_label) +
    theme_bw()
  
  # ---- 8. Dispersion plot ----------------------------------------------
  p_var <- ggplot(mean_pred, aes(x = .data[[x_var]], y = disp_plot_var, group = group_var)) +
    geom_line(linetype = linetype_var) +
    facet_layer +
    labs(y = disp_y_label, x = x_label) +
    theme_bw()
  
  # ---- 9. Combine & return ---------------------------------------------
  p_combined <- ggarrange(p_mean, p_var, ncol = ncol, nrow = nrow)
  
  return(list(
    mean_plot = p_mean,
    var_plot  = p_var,
    combined  = p_combined,
    data      = mean_pred,
    grid_vars = grid_list
  ))
}