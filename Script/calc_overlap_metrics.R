calc_overlap_metrics <- function(df, sp_col = "sp") {
  sp_names <- df[[sp_col]]
  mat <- as.matrix(df[, !names(df) %in% sp_col])
  rownames(mat) <- sp_names
  
  # Core matrix operations (no loop)
  range_sizes <- rowSums(mat)
  shared_mat  <- tcrossprod(mat)  # mat %*% t(mat) = pairwise intersections
  
  # Get upper triangle indices for all pairs
  idx <- which(upper.tri(shared_mat), arr.ind = TRUE)
  i   <- idx[, 1]
  j   <- idx[, 2]
  
  A      <- range_sizes[i]
  B      <- range_sizes[j]
  shared <- shared_mat[idx]
  union  <- A + B - shared
  
  data.frame(
    sp1            = sp_names[i],
    sp2            = sp_names[j],
    shared         = shared,
    range_sp1      = A,
    range_sp2      = B,
    sympatry_min   = ifelse(pmin(A, B) == 0, NA_real_, shared / pmin(A, B)),
    sympatry_union = ifelse(union == 0,       NA_real_, shared / union),
    range_symmetry = ifelse((A + B) == 0,     NA_real_, pmin(A, B) / (A + B))
  )
}