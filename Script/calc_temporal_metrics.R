library(picante)
library(dplyr)

calc_temporal_metrics <- function(X, year, plot_info, sp_info) {
  X <- as.data.frame(X)
  
  # strip any inherited row names from both info objects
  rownames(plot_info) <- NULL
  rownames(sp_info)   <- NULL
  
  n_sp <- ncol(X)
  
  # --- Observed covariances ---
  cor.ob <- cor(X)
  var    <- diag(cov(X))
  diag(cor.ob) <- NA
  
  # Lower-triangle pairs only
  idx <- which(lower.tri(cor.ob), arr.ind = TRUE)   # faster than expand.grid + filter
  pos_df <- data.frame(Var1 = idx[, 1], Var2 = idx[, 2])
  
  pos_df$sp1_var <- var[pos_df$Var1]
  pos_df$sp2_var <- var[pos_df$Var2]
  
  # Vectorised pairwise total variance (no loop)
  pair_sums  <- vapply(
    seq_len(nrow(pos_df)),
    function(i) var(X[, pos_df$Var1[i]] + X[, pos_df$Var2[i]]),
    numeric(1)
  )
  pos_df$tvar_pair  <- pair_sums
  pos_df$log_Vratio <- log(pair_sums / (pos_df$sp1_var + pos_df$sp2_var))
  pos_df$phi        <- pair_sums / (sqrt(pos_df$sp1_var) + sqrt(pos_df$sp2_var))^2
  
  # --- Observed abundance trends (all species in one lm call) ---
  Y        <- as.numeric(year - min(year))
  fit      <- lm(as.matrix(X) ~ Y)          # single call, ncol(X) regressions
  coef.obs <- coef(fit)[2, ]                 # slope row
  names(coef.obs) <- colnames(X)
  
  # --- Pair-level output ---
  sp_name <- data.frame(
    sp1 = colnames(X)[pos_df$Var1],
    sp2 = colnames(X)[pos_df$Var2],
    stringsAsFactors = FALSE
  )
  
  SES.pair.df <- data.frame(
    plot_info,
    sp_name,
    obs.cor    = cor.ob[lower.tri(cor.ob)],
    log_Vratio = pos_df$log_Vratio,
    phi        = pos_df$phi,
    Var1       = pos_df$Var1,
    Var2       = pos_df$Var2
  ) |>
    left_join(sp_info, by = join_by(sp1 == Accepted_SPNAME)) |>
    left_join(sp_info, by = join_by(sp2 == Accepted_SPNAME))
  
  # Average observed correlation per species (column means of lower-tri matrix)
  cor.ob.full           <- cor.ob
  cor.ob.full[upper.tri(cor.ob.full)] <- t(cor.ob.full)[upper.tri(cor.ob.full)]
  obs.cor.avg <- colMeans(cor.ob.full, na.rm = TRUE)
  
  # --- Species-level output ---
  species.df <- data.frame(
    sp             = colnames(X),
    var            = unname(var),
    mean_abundance = unname(colMeans(X)),
    obs.cor.avg    = unname(obs.cor.avg),
    coef.obs       = unname(coef.obs),
    plot_info
  ) |>
    left_join(sp_info, by = join_by(sp == Accepted_SPNAME))
  
  list(
    pair.df    = SES.pair.df,
    species.df = species.df
  )
}