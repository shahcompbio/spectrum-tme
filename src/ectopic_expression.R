## compute ectopic expression across clusters

find_markers_wrapper <- function(seu_sub, group_var = tumor_megasite, 
                                 split_var = cluster_label) {
  
  group_var <- enquo(group_var)
  split_var <- enquo(split_var)
  
  seu_sub[[as_label(group_var)]] <- deframe(select(scrna_meta_tbl, sample, !!group_var))[seu_sub$sample]
  Idents(seu_sub) <- seu_sub[[as_label(group_var)]]
  
  all_var_names <- names(clrs[[as_label(group_var)]])
  sub_var_names <- all_var_names[all_var_names %in% unique(Idents(seu_sub))]
  
  comparisons_mat <- combn(rev(sub_var_names), 2)
  
  marker_list <- list()
  split_value <- unique(unlist(seu_sub[[as_label(split_var)]]))
  
  for (i in 1:ncol(comparisons_mat)) {
    
    id1 <- comparisons_mat[1,i]
    id2 <- comparisons_mat[2,i]
    
    marker_list[[paste0(id1, "_vs_", id2)]] <- FindMarkers(
      seu_sub, ident.1 = id1, ident.2 = id2) %>%
      as_tibble(rownames = "gene") %>% 
      mutate(comparison = paste0(id1, "_vs_", id2),
             group1 = id1, group2 = id2,
             !!split_var := split_value)
    
    print(paste0(id1 " vs " id2, " in ", split_value, ": DONE\n"))
    
  }
  
  return(bind_rows(marker_list))
  
}

load_and_compute_wrapper <- function(super_set, group_var, split_var) {
  
  seu_obj <- read_rds(paste0("/work/shah/users/uhlitzf/data/SPECTRUM/freeze/v7/", 
                             super_set, "_processed_filtered_annotated.rds"))
  
  group_var <- enquo(group_var)
  split_var <- enquo(split_var)
  
  seu_split <- SplitObject(seu_obj, as_label(split_var))
  
  find_markers_wrapper_wrapper <- function(seu_obj) find_markers_wrapper(seu_obj, group_var = !!group_var, split_var = !!split_var)
  
  marker_list <- parallel::mclapply(seu_split, find_markers_wrapper_wrapper, mc.cores = 10)
  marker_tbl <- bind_rows(marker_list[unlist(lapply(marker_list, is.data.frame))])
  
  write_tsv(marker_tbl, paste0("/work/shah/users/uhlitzf/data/SPECTRUM/freeze/v7/", 
                               super_set, "_ectopic_expression_markers.tsv"))
  
  rm(seu_obj)
  rm(seu_split)
  
  print(paste0(super_set, " COMPLETE"))
  
  
}

load_and_compute_wrapper("Ovarian.cancer.super", tumor_megasite, cluster_label)
load_and_compute_wrapper("Myeloid.super", tumor_megasite, cluster_label)
load_and_compute_wrapper("T.super", tumor_megasite, cluster_label_sub)
load_and_compute_wrapper("Fibroblast.super", tumor_megasite, cluster_label)
load_and_compute_wrapper("Endothelial.super", tumor_megasite, cluster_label)

