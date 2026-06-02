random_matching <- function(df, n_iter = 100) {
  
  plots <- split(df, df$plotID_full)
  
  results <- lapply(seq_len(n_iter), function(iter) {
    
    iter_result <- lapply(plots, function(plot_df) {
      
      # shuffle for randomness
      plot_df <- plot_df[sample(nrow(plot_df)), ]
      
      used_sp  <- character(0)
      selected <- logical(nrow(plot_df))
      
      for (k in seq_len(nrow(plot_df))) {
        s1 <- plot_df$sp1[k]
        s2 <- plot_df$sp2[k]
        
        if (!s1 %in% used_sp & !s2 %in% used_sp) {
          selected[k] <- TRUE
          used_sp <- c(used_sp, s1, s2)
        }
      }
      
      plot_df[selected, ]
    })
    
    bind_rows(iter_result)
  })
  
  bind_rows(results, .id = "iter")
}

