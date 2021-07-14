

kde2d_tidy <- function(data_tbl, x, y, kernel_group, kernel_subset,  
                       n = 200, h = 0.4) {
  
  x <- enquo(x)
  y <- enquo(y)
  kernel_group <- enquo(kernel_group)
  
  data_tbl <- filter(data_tbl, !!kernel_group %in% kernel_subset)
  
  kde_data <- MASS::kde2d(data_tbl[[as_label(x)]], data_tbl[[as_label(y)]], n = n, h = h) 
  
  result_tbl <- kde_data$z %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "x") %>% 
    gather(y, value, -x) %>% 
    mutate(y = str_remove_all(y, "V") %>% as.numeric,
           x = as.numeric(x),
           !!kernel_group := paste0(kernel_subset, collapse = "_")) %>% 
    mutate(x = scales::rescale(x),
           y = scales::rescale(y))
  
  return(result_tbl)
  
}

kde2d_contrast <- function(data_tbl, x, y, kernel_group, kernel_subsets_a, kernel_subsets_b, 
                           n = 200, h = 0.2, quench = 0.015) {
  
  x <- enquo(x)
  y <- enquo(y)
  kernel_group <- enquo(kernel_group)
  kernel_a <- kde2d_tidy(data_tbl, !!x, !!y, !!kernel_group, kernel_subsets_a, n = n, h = h)
  kernel_b <- kde2d_tidy(data_tbl, !!x, !!y, !!kernel_group, kernel_subsets_b, n = n, h = h)
  
  kernel_tbl <- kernel_a %>% 
    left_join(kernel_b, by = c("x", "y")) %>% 
    mutate(value = value.x - value.y) %>% 
    mutate(delta_group = value > 0) %>% 
    group_by(delta_group) %>% 
    mutate(value_quenched = ifelse(value > quench, quench, ifelse(value < -quench, -quench, value))) %>% 
    mutate(value_scaled = ifelse(delta_group, scales::rescale(value_quenched, c(0, 1)), scales::rescale(value_quenched, c(-1, 0)))) %>% 
    ungroup()
  
  return(kernel_tbl)
  
}

# kernel_tbl <- kde2d_contrast(
#   data_tbl = plot_data_sub %>% filter(cell_type_myeloid == "Macrophage"), 
#   x = umapharmony_1, y = umapharmony_2, kernel_group = tumor_megasite, 
#   kernel_subsets_a = c("Adnexa"), kernel_subsets_b = c("Other", "Ascites")) %>% 
#   filter(value.x > 0.01 & value.y > 0.01
# )


