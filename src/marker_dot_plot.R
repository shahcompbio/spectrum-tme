
## marker dot plots showing logFC and pct.expr per cluster -------------------------

## input: seurat object and two-column cluster/marker table

calc_fc <- function(data_tbl, coi, cluster_var) {
  cluster_var <- enquo(cluster_var)
  data_tbl %>%
    mutate(!!cluster_var := ifelse(!!cluster_var == coi, coi, "Other")) %>%
    gather(gene, value, -!!cluster_var, -sample) %>%
    group_by(!!cluster_var, gene) %>%
    summarise(mean_value = mean(value),
              .groups = "drop") %>%
    spread(!!cluster_var, mean_value) %>%
    mutate(logFC = pull(., 2) - pull(., 3)) %>% 
    mutate(logFC = case_when(
      logFC < -2 ~ -2,
      logFC > 2 ~ 2,
      T ~ logFC
    )) %>% 
    mutate(!!cluster_var := coi) %>%
    select(gene, !!cluster_var, logFC)
}

marker_dot_plot_preprocess <- function(seu_obj, marker_tbl, cluster_var = cluster_label, 
                                       sample_ncells = 10000, sample_seed = 123, top_n = 5) {
  
  marker_tbl <- marker_tbl %>% 
    group_by(cluster_marker) %>% 
    slice(1:top_n) %>% 
    ungroup()
  
  cluster_var <- enquo(cluster_var)
  cois <- unique(pull(marker_tbl, cluster_marker))
    
  expr_data <- FetchData(seu_obj, c(as_label(cluster_var), "sample", unique(marker_tbl$gene)), 
                         slot = "data") %>% 
    as_tibble()
  
  if (is.numeric(sample_ncells)) {
    set.seed(sample_seed)
    expr_data <- sample_n(expr_data, sample_ncells)
  }
  
  logfc_data <- bind_rows(lapply(cois, function(x) calc_fc(expr_data, x, !!cluster_var)))
  
  plot_data <- expr_data %>% 
    gather(gene, value, -sample, -!!cluster_var) %>% 
    group_by(!!cluster_var, gene) %>% 
    summarise(mean_value = mean(value, na.rm = T),
              n_cells = n(),
              pct_expr = sum(value > 0)/n(),
              .groups = "drop") %>% 
    left_join(select(marker_tbl, gene, cluster_marker), by = "gene") %>% 
    left_join(logfc_data, by = c(as_label(cluster_var), "gene")) %>% 
    mutate(!!cluster_var := ordered(!!cluster_var, levels = cois)) %>% 
    mutate(cluster_marker = ordered(cluster_marker, levels = cois))
  
  return(plot_data)
  
}

plot_marker_dots <- function(dot_plot_data, super_set, cluster_var) {
  
  cluster_var <- enquo(cluster_var)
  
  rdbu_breaks <- c(seq(min(dot_plot_data$logFC, na.rm = T), 0, length.out = 5)[-5], 
                   0, seq(0, max(dot_plot_data$logFC, na.rm = T), length.out = 5)[-1])
  
  rdbu <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
  
  dot_plot_data <- dot_plot_data %>% 
    mutate(cluster_numeric_var = as.numeric(!!cluster_var),
           cluster_numeric_marker = as.numeric(cluster_marker))
  
  dot_plot <- ggplot(dot_plot_data) + 
    geom_point(aes(gene, !!cluster_var, color = logFC, size = pct_expr)) +
    facet_grid(vars(cluster_numeric_var), vars(cluster_numeric_marker), 
               scales = "free", space = "free", switch = "y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          strip.text.y.left = element_text(angle = 0),
          strip.background = element_rect()) +
    scale_size_continuous(range = c(0, 5)) + 
    scale_color_gradientn(colors = rdbu, values = scales::rescale(rdbu_breaks)) +
    labs(color = "logFC", size = "Expr. [%]", x = "Gene", y = "Cluster")
  
  g <- ggplot_gtable(ggplot_build(dot_plot))
  
  strip_elements <- which(grepl('strip-', g$layout$name))
  fills <- rep(clrs[[as_label(cluster_var)]][[super_set]], times = 2)
  k <- 1
  for (i in strip_elements) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  return(ggdraw() + draw_grob(g))
  
}

