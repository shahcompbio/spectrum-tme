---
title: "MSK SPECTRUM data freeze: major cell type embeddings"
author: "Florian Uhlitz"
date: "2021-05-06"
output:
  bookdown::html_document2:
    code_folding: hide
editor_options: 
  chunk_output_type: console
---




```r
library(magrittr)
library(tidyverse)
library(Seurat)
library(readxl)
library(cowplot)
library(colorblindr)
library(viridis)
library(magick, lib.loc = "/home/uhlitzf/miniconda3/lib/R/library")
library(ggpubr)
library(knitr)
```

# Figure 2 upper


```r
## load global vars: 
source("src/global_vars.R")

# scrna_meta_tbl
# clrs
# markers_v7
# markers_v7_super
# cell_type_super_lookup

names(clrs$cell_type) <- str_replace_all(names(clrs$cell_type), "\\.", " ")
names(clrs$cell_type) <- str_replace_all(names(clrs$cell_type), "Ovarian", "Ov")

## load data --------------------------------------

## load cohort embeddings
seu_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/outs_pre/cells.tsv") %>% 
  mutate(cell_type = ifelse(cell_type == "Monocyte", "Myeloid.cell", cell_type)) %>% 
  mutate(cell_type = ifelse(cell_type == "Ovarian.cancer.cell", "Ov.cancer.cell", cell_type))

# seu_tbl <- read_tsv("/work/shah/isabl_data_lake/analyses/68/74/6874/cells.tsv") %>%
#   mutate(cell_type = ifelse(cell_type == "Monocyte", "Myeloid.cell", cell_type)) %>%
#   mutate(cell_type = ifelse(cell_type == "Ovarian.cancer.cell", "Ov.cancer.cell", cell_type))

## load consensus data
consOV_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/consensusOV/SPECTRUM_freeze_v7_consensusOV.tsv") %>%
  mutate(consensusOV = ordered(consensusOV, levels = names(clrs$consensusOV))) %>% 
  select(cell_id, consensusOV)

## join data
seu_tbl_full <- seu_tbl %>% 
  left_join(scrna_meta_tbl, by = "sample") %>% 
  filter(therapy == "pre-Rx", cell_type != "Other", tumor_supersite != "Unknown") %>% 
  rename(UMAP_1 = umap50_1, UMAP_2 = umap50_2) %>% 
  ## flip UMAP coordinates?
  # mutate(UMAP_1 = -UMAP_1, UMAP_2 = -UMAP_2) %>% 
  mutate(cell_type_super = cell_type_super_lookup[cell_type]) %>% 
  mutate(cell_type = ordered(str_replace_all(cell_type, "\\.", " "), 
                             levels = names(clrs$cell_type))) %>% 
  mutate(sort_short_x = ifelse(sort_short == "U" & cell_type_super == "Immune", 
                               "CD45+", ifelse(sort_short == "U" & cell_type_super == "Stromal", 
                                               "CD45-", sort_short))) %>% 
  left_join(consOV_tbl, by = "cell_id")
```


```r
# ## extract cell cycle tbl
# cell_cycle_list <- as.list(rep(NA, times = length(names(clrs$cluster_label)))) %>%
#   setNames(names(clrs$cluster_label))
# 
# for (i in 1:length(cell_cycle_list)) {
#   seu_obj <- read_rds(
#     paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/outs_pre/",
#            names(cell_cycle_list)[i], "_seurat_0.1.rds"))
#   cell_cycle_list[[i]] <- as_tibble(FetchData(seu_obj, c("cell_id", "S.Score", "G2M.Score", "CC.Diff", "Phase")))
#   rm(seu_obj)
# }
# 
# cell_cycle_tbl <- bind_rows(cell_cycle_list)
# 
# write_tsv(cell_cycle_tbl, "/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/cell_cycle_tbl.tsv")

cell_cycle_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/cell_cycle_tbl.tsv")

seu_tbl_full <- seu_tbl_full %>%
  left_join(cell_cycle_tbl, by = "cell_id")

# seu_tbl_full <- seu_tbl_full %>%
#   sample_n(10000)
```

## UMAPs


```r
dpi <- 72
rw <- 10
rh <- 10

base_umap <- ggplot(seu_tbl_full) +
  coord_fixed() +
  NoAxes() +
  theme(legend.position = c(0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = ggplot2::margin(1, 1, 1, 1),
        #panel.border = element_rect(linetype = 1, color = "black", size = 1),
        legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 22))

pt.size <- 0.1
pt.size.mini <- 0.01
pt.alpha <- 0.05
pt.alpha.mini <- 0.02

# median_tbl <- seu_tbl_full %>% 
#   filter(cell_type != "Other") %>% 
#   group_by(cell_type, reduction) %>% 
#   summarise(UMAP_1 = median(UMAP_1),
#             UMAP_2 = median(UMAP_2))

ncol1_legend <- list(
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1, 
                              label.position = "right")),
  theme(legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 9.5, 0), color = "white"),
        legend.spacing.x = unit(0, "npc"))
)

umap_site <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = tumor_supersite), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$tumor_supersite) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 2, 
                              label.position = "right"))


umap_site_mini <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = tumor_supersite), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$tumor_supersite) +
  guides(color = F)

umap_site_legend_ncol1 <- cowplot::get_legend(umap_site + ncol1_legend)
umap_site_legend <- cowplot::get_legend(umap_site)

umap_cell_type <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = fct_rev(cell_type)), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$cell_type) +
  ncol1_legend

umap_cell_type_mini <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = cell_type), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$cell_type) +
  guides(color = F)

umap_cell_type_legend <- cowplot::get_legend(umap_cell_type)

umap_patient <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = patient_id_short), size = pt.size, alpha = pt.alpha, raster.dpi = 200, raster.width = 10, raster.height = 10) +
  scale_color_manual(values = clrs$patient_id_short) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              nrow = 17, 
                              label.position = "right")) +
  theme(legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc"))

umap_patient_mini <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = patient_id_short), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$patient_id_short) +
  guides(color = F)

# umap_patient_legend <- cowplot::get_legend(umap_patient)

umap_consOV <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = consensusOV), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$consensusOV, labels = c("Immuno.", "Mesench.", "Prolif.", "Differen.")) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 2, 
                              label.position = "right")) +
  theme(legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc"))

umap_consOV_mini <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = consensusOV), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$consensusOV) +
  guides(color = F)

umap_consOV_legend <- cowplot::get_legend(umap_consOV)

umap_mutsig <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = consensus_signature), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$consensus_signature) + 
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 2, 
                              label.position = "right")) +
  theme(legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc"))

umap_mutsig_mini <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = consensus_signature), size = pt.size.mini, alpha = pt.alpha.mini, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$consensus_signature) +
  guides(color = F)

umap_mutsig_legend <- cowplot::get_legend(umap_mutsig)

umap_phase <- base_umap +
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = Phase), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  scale_color_manual(values = clrs$Phase) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1,
                              label.position = "right")) +
  theme(legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc")) +
  labs(title = "Phase", color = "")

continuous_umap_layers <- list(
  scale_color_viridis_c(),
  guides(color = guide_colorbar(label.position = "right",
                                title.position = "top",
                                title.hjust = 0,
                                title.vjust = 1,
                                direction = "vertical")),
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.03, "npc"),
        legend.position = c(-0.07, 1),
        legend.justification = c("left", "top"),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)))

umap_umis <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = log2(nCount_RNA)), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  labs(title = "log2(#UMIs)", color = "") +
  continuous_umap_layers

umap_genes <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = log2(nFeature_RNA)), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  labs(title = "log2(#Genes)", color = "") +
  continuous_umap_layers

umap_mitos <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = percent.mt), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  labs(title = "MT reads [%]", color = "") +
  continuous_umap_layers

umap_doublet_score <- base_umap + 
  ggrastr::geom_point_rast(aes(UMAP_1, UMAP_2, color = doublet_score), size = pt.size, alpha = pt.alpha, raster.dpi = dpi, raster.width = rw, raster.height = rh) +
  labs(title = "Doublet score", color = "") +
  continuous_umap_layers
```

## total numbers 


```r
total_numbers_plot <- function(column_var) {
  column_var <- enquo(column_var)
  plot_data <- seu_tbl_full %>% 
    filter(cell_type != "Other", therapy == "pre-Rx") %>%
    group_by(!!column_var) %>%
    tally %>%
    arrange(n) %>%
    mutate(column_label = paste0(!!column_var, " (", format(n, trim=T, big.mark=","), ")")) %>% 
    mutate(!!column_var := ordered(!!column_var, levels = unique(!!column_var)))
  
  common_layers <- list(
    scale_y_continuous(expand = c(0, 0)),
    guides(fill = F, color = F),
    coord_flip(clip = "off"),
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(face = "plain"))
  )
  
  p1 <- ggplot(plot_data) +
    geom_point(aes(!!column_var, 0, color = !!column_var), size = 4) +
    scale_color_manual(values = clrs[[as_label(column_var)]]) +
    common_layers 
    
  p2 <- ggplot(plot_data) +
    geom_bar(aes(!!column_var, n, fill = !!column_var), stat = "identity", width = 0.3) +
    geom_text(aes(!!column_var, 1, label = column_label), hjust = 0, vjust = 0, nudge_x = 0.2) +
    scale_fill_manual(values = clrs[[as_label(column_var)]]) +
    common_layers 
  
  plot_grid(p1, p2, align = "v", ncol = 2, rel_widths = c(0.03, 0.97))
  
}

numbers_plot_cell_type <- total_numbers_plot(cell_type)
numbers_plot_tumor_supersite <- total_numbers_plot(tumor_supersite)
```

## Mixing per cell type


```r
pm_score_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/patient_mixing_scores.tsv") %>%
  mutate(cell_type = str_replace_all(cell_type, "Monocyte", "Myeloid cell")) %>% 
  mutate(cell_type = str_replace_all(cell_type, "Ovarian.cancer.cell", "Ov cancer cell")) %>% 
  mutate(patient_id_short = str_sub(cell_id, 13, 15)) %>% 
  left_join(distinct(scrna_meta_tbl, patient_id_short, .keep_all = T), by = "patient_id_short") %>% 
  filter(cell_type != "Other")

pm_score_tbl_ordered <- pm_score_tbl %>%
  filter(cell_type != "Other") %>%
  group_by(cell_type) %>%
  mutate(median = median(score)) %>%
  ungroup() %>%
  arrange(-median) %>%
  mutate(cell_type = ordered(cell_type, levels = unique(cell_type)))

mixing_plot_cell_type <- ggplot(pm_score_tbl_ordered) +
  geom_violin(aes(cell_type, score, fill = cell_type),
              adjust = 3, color = "white", width = 1) +
  geom_boxplot(aes(cell_type, score, color = cell_type),
               width = 0.2, size = 1, fill = NA, outlier.shape = NA) +
  geom_boxplot(aes(cell_type, score, fill = cell_type),
               width = 0.2, color = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "black") +
  scale_fill_manual(values = clrs$cell_type) +
  scale_color_manual(values = clrs$cell_type) +
  guides(fill = F, color = F) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-5, 7),
                     breaks = c(-4, -2, 0, 2, 4, 6)) +
  labs(y = "Patient\nspecificity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_blank())


ggplot(pm_score_tbl) +
  geom_violin(aes(consensus_signature, score, fill = consensus_signature),
              adjust = 3, color = "white", width = 1) +
  geom_boxplot(aes(consensus_signature, score, color = consensus_signature),
               width = 0.2, size = 1, fill = NA, outlier.size = 0.01) +
  geom_boxplot(aes(consensus_signature, score, fill = consensus_signature),
               width = 0.2, color = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "black") +
  scale_fill_manual(values = clrs$consensus_signature) +
  scale_color_manual(values = clrs$consensus_signature) +
  facet_wrap(~cell_type, scales = "free") + 
  guides(fill = F, color = F) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-5, 7),
                     breaks = c(-4, -2, 0, 2, 4, 6)) +
  labs(y = "Patient\nspecificity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_110-1.png" width="576" />

```r
ggsave_nodb("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.pdf", width = 6, height = 9)
ggsave("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.png", width = 6, height = 9)
```

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.pdf]("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.pdf")

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.png]("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig.png")


```r
ggplot(pm_score_tbl) +
  geom_violin(aes(patient_id_short, score, fill = consensus_signature),
              adjust = 3, color = "white", width = 1) +
  geom_boxplot(aes(patient_id_short, score, color = consensus_signature),
               width = 0.2, size = 1, fill = NA, outlier.size = 0.01) +
  geom_boxplot(aes(patient_id_short, score, fill = consensus_signature),
               width = 0.2, color = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "black") +
  scale_fill_manual(values = clrs$consensus_signature) +
  scale_color_manual(values = clrs$consensus_signature) +
  facet_wrap(~cell_type, scales = "free") + 
  guides(fill = F, color = F) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(-5, 7),
                     breaks = c(-4, -2, 0, 2, 4, 6)) +
  labs(y = "Patient\nspecificity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_111-1.png" width="1728" />

```r
ggsave_nodb("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.pdf", width = 18, height = 9)
ggsave("figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.png", width = 18, height = 9)
```

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.pdf](figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.pdf)

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.png](figures/200_cohort_plotting_umap_v7/002_patient_mixing_mutsig_patient.png)


## pairwise comparisons


```r
### pairwise comparisons

pm_score_summary <- pm_score_tbl_ordered %>%
  separate(cell_id, into = c("patient_id", "sample_suffix"), remove = F, extra = "drop", sep = "_") %>%
  group_by(patient_id, cell_type) %>%
  summarise(mean_score = mean(score))

ptt_mat <- pairwise.t.test(pm_score_summary$mean_score, pm_score_summary$cell_type, p.adjust.method = "fdr") %$% 
  p.value

ptt_mat_full <- matrix(nrow = 9, ncol = 9)
rownames(ptt_mat_full) <- unique(c(colnames(ptt_mat), rownames(ptt_mat)))
colnames(ptt_mat_full) <- unique(c(colnames(ptt_mat), rownames(ptt_mat)))
ptt_mat_full[upper.tri(ptt_mat_full)] <- 1
ptt_mat_full[lower.tri(ptt_mat_full)] <- ptt_mat[lower.tri(ptt_mat, diag = T)]
diag(ptt_mat_full) <- 1

ptt_tbl <- ptt_mat_full %>% 
  as_tibble(rownames = "cell_type_1") %>%
  gather(cell_type_2, pvalue, -cell_type_1) %>%
  mutate(cell_type_2 = ordered(cell_type_2, levels = levels(pm_score_tbl_ordered$cell_type))) %>%
  mutate(cell_type_1 = ordered(cell_type_1, levels = levels(pm_score_tbl_ordered$cell_type))) %>%
  mutate(log10p = -log10(pvalue),
         sig = pvalue < 0.1,
         label = ifelse(sig, "*", "")) %>%
  replace_na(list(pvalue = 1, log10p = 0, sig = F, label = ""))

pm_comparisons <- ggplot(ptt_tbl) +
  # geom_tile(aes(fct_rev(cell_type_1), fct_rev(cell_type_2), fill = log10p)) +
  geom_point(aes(fct_rev(cell_type_1), fct_rev(cell_type_2), size = log10p), color = "black", shape = 21) +
  geom_point(aes(fct_rev(cell_type_1), fct_rev(cell_type_2), fill = log10p, size = log10p), shape = 21) +
  # geom_point(aes(fct_rev(cell_type_1), fct_rev(cell_type_2), shape = sig),
  #           color = "white", size = 5) +
  # geom_point(aes(fct_rev(cell_type_1), fct_rev(cell_type_2), shape = sig),
  #           color = "black", size = 3) +
  scale_fill_gradientn(colours = viridis(9)) +
  scale_size_continuous(guide = F, range = c(-1, 5)) +
  # scale_color_gradientn(colours = viridis(9)) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        plot.title = element_text(face = "plain", size = 12)) +
  scale_shape_manual(values = c("", "*"), guide = F) +
  # scale_x_discrete(expand = c(0, 0)) +
  # scale_y_discrete(expand = c(0, 0)) +
  labs(fill = "-log10(q-value)", title = "Pairwise comparisons\nof patient specificity")

pm_comparisons
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_115-1.png" width="480" />

```r
print(ptt_tbl, n = 100)
```

```
## # A tibble: 81 x 6
##    cell_type_1      cell_type_2        pvalue  log10p sig   label
##    <ord>            <ord>               <dbl>   <dbl> <lgl> <chr>
##  1 Ov cancer cell   Ov cancer cell   1   e+ 0  0      FALSE ""   
##  2 Fibroblast       Ov cancer cell   9.02e-46 45.0    TRUE  "*"  
##  3 Myeloid cell     Ov cancer cell   2.03e-62 61.7    TRUE  "*"  
##  4 T cell           Ov cancer cell   3.89e-79 78.4    TRUE  "*"  
##  5 Dendritic cell   Ov cancer cell   1.12e-86 86.0    TRUE  "*"  
##  6 Plasma cell      Ov cancer cell   5.60e-85 84.3    TRUE  "*"  
##  7 Endothelial cell Ov cancer cell   2.90e-84 83.5    TRUE  "*"  
##  8 B cell           Ov cancer cell   2.35e-89 88.6    TRUE  "*"  
##  9 Mast cell        Ov cancer cell   9.17e-98 97.0    TRUE  "*"  
## 10 Ov cancer cell   Fibroblast       1   e+ 0  0      FALSE ""   
## 11 Fibroblast       Fibroblast       1   e+ 0  0      FALSE ""   
## 12 Myeloid cell     Fibroblast       1.02e- 4  3.99   TRUE  "*"  
## 13 T cell           Fibroblast       1.11e-14 14.0    TRUE  "*"  
## 14 Dendritic cell   Fibroblast       4.81e-21 20.3    TRUE  "*"  
## 15 Plasma cell      Fibroblast       1.75e-19 18.8    TRUE  "*"  
## 16 Endothelial cell Fibroblast       8.01e-19 18.1    TRUE  "*"  
## 17 B cell           Fibroblast       1.33e-23 22.9    TRUE  "*"  
## 18 Mast cell        Fibroblast       3.67e-32 31.4    TRUE  "*"  
## 19 Ov cancer cell   Myeloid cell     1   e+ 0  0      FALSE ""   
## 20 Fibroblast       Myeloid cell     1   e+ 0  0      FALSE ""   
## 21 Myeloid cell     Myeloid cell     1   e+ 0  0      FALSE ""   
## 22 T cell           Myeloid cell     7.51e- 5  4.12   TRUE  "*"  
## 23 Dendritic cell   Myeloid cell     4.93e- 9  8.31   TRUE  "*"  
## 24 Plasma cell      Myeloid cell     6.16e- 8  7.21   TRUE  "*"  
## 25 Endothelial cell Myeloid cell     1.74e- 7  6.76   TRUE  "*"  
## 26 B cell           Myeloid cell     6.57e-11 10.2    TRUE  "*"  
## 27 Mast cell        Myeloid cell     1.03e-17 17.0    TRUE  "*"  
## 28 Ov cancer cell   T cell           1   e+ 0  0      FALSE ""   
## 29 Fibroblast       T cell           1   e+ 0  0      FALSE ""   
## 30 Myeloid cell     T cell           1   e+ 0  0      FALSE ""   
## 31 T cell           T cell           1   e+ 0  0      FALSE ""   
## 32 Dendritic cell   T cell           6.07e- 2  1.22   TRUE  "*"  
## 33 Plasma cell      T cell           1.59e- 1  0.797  FALSE ""   
## 34 Endothelial cell T cell           2.19e- 1  0.659  FALSE ""   
## 35 B cell           T cell           9.26e- 3  2.03   TRUE  "*"  
## 36 Mast cell        T cell           1.45e- 6  5.84   TRUE  "*"  
## 37 Ov cancer cell   Dendritic cell   1   e+ 0  0      FALSE ""   
## 38 Fibroblast       Dendritic cell   1   e+ 0  0      FALSE ""   
## 39 Myeloid cell     Dendritic cell   1   e+ 0  0      FALSE ""   
## 40 T cell           Dendritic cell   1   e+ 0  0      FALSE ""   
## 41 Dendritic cell   Dendritic cell   1   e+ 0  0      FALSE ""   
## 42 Plasma cell      Dendritic cell   6.59e- 1  0.181  FALSE ""   
## 43 Endothelial cell Dendritic cell   5.30e- 1  0.276  FALSE ""   
## 44 B cell           Dendritic cell   5.04e- 1  0.298  FALSE ""   
## 45 Mast cell        Dendritic cell   3.76e- 3  2.42   TRUE  "*"  
## 46 Ov cancer cell   Plasma cell      1   e+ 0  0      FALSE ""   
## 47 Fibroblast       Plasma cell      1   e+ 0  0      FALSE ""   
## 48 Myeloid cell     Plasma cell      1   e+ 0  0      FALSE ""   
## 49 T cell           Plasma cell      1   e+ 0  0      FALSE ""   
## 50 Dendritic cell   Plasma cell      1   e+ 0  0      FALSE ""   
## 51 Plasma cell      Plasma cell      1   e+ 0  0      FALSE ""   
## 52 Endothelial cell Plasma cell      8.36e- 1  0.0777 FALSE ""   
## 53 B cell           Plasma cell      2.58e- 1  0.589  FALSE ""   
## 54 Mast cell        Plasma cell      7.89e- 4  3.10   TRUE  "*"  
## 55 Ov cancer cell   Endothelial cell 1   e+ 0  0      FALSE ""   
## 56 Fibroblast       Endothelial cell 1   e+ 0  0      FALSE ""   
## 57 Myeloid cell     Endothelial cell 1   e+ 0  0      FALSE ""   
## 58 T cell           Endothelial cell 1   e+ 0  0      FALSE ""   
## 59 Dendritic cell   Endothelial cell 1   e+ 0  0      FALSE ""   
## 60 Plasma cell      Endothelial cell 1   e+ 0  0      FALSE ""   
## 61 Endothelial cell Endothelial cell 1   e+ 0  0      FALSE ""   
## 62 B cell           Endothelial cell 1.91e- 1  0.719  FALSE ""   
## 63 Mast cell        Endothelial cell 3.82e- 4  3.42   TRUE  "*"  
## 64 Ov cancer cell   B cell           1   e+ 0  0      FALSE ""   
## 65 Fibroblast       B cell           1   e+ 0  0      FALSE ""   
## 66 Myeloid cell     B cell           1   e+ 0  0      FALSE ""   
## 67 T cell           B cell           1   e+ 0  0      FALSE ""   
## 68 Dendritic cell   B cell           1   e+ 0  0      FALSE ""   
## 69 Plasma cell      B cell           1   e+ 0  0      FALSE ""   
## 70 Endothelial cell B cell           1   e+ 0  0      FALSE ""   
## 71 B cell           B cell           1   e+ 0  0      FALSE ""   
## 72 Mast cell        B cell           2.98e- 2  1.53   TRUE  "*"  
## 73 Ov cancer cell   Mast cell        1   e+ 0  0      FALSE ""   
## 74 Fibroblast       Mast cell        1   e+ 0  0      FALSE ""   
## 75 Myeloid cell     Mast cell        1   e+ 0  0      FALSE ""   
## 76 T cell           Mast cell        1   e+ 0  0      FALSE ""   
## 77 Dendritic cell   Mast cell        1   e+ 0  0      FALSE ""   
## 78 Plasma cell      Mast cell        1   e+ 0  0      FALSE ""   
## 79 Endothelial cell Mast cell        1   e+ 0  0      FALSE ""   
## 80 B cell           Mast cell        1   e+ 0  0      FALSE ""   
## 81 Mast cell        Mast cell        1   e+ 0  0      FALSE ""
```

```r
ggsave_nodb("figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.pdf", pm_comparisons, width = 5, height = 4)
ggsave("figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.png", pm_comparisons, width = 5, height = 4)
```

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.pdf](figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.pdf)

[figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.png](figures/200_cohort_plotting_umap_v7/002_patient_mixing_tests.png)

## consensusOV summary


```r
consOV_tbl_summary <- seu_tbl_full %>%
  select(cell_type, consensusOV, patient_id_short, tumor_supersite) %>%
  mutate(consOV_cell_type = cell_type) %>%
  # mutate(consOV_cell_type = ifelse(S.Score > 0.5, "Cycling cell", as.character(cell_type))) %>%
  group_by(consOV_cell_type, consensusOV) %>%
  tally %>%
  group_by(consOV_cell_type) %>%
  mutate(nrel = n/(sum(n))*100) %>%
  ungroup %>%
  mutate(consensusOV = ordered(consensusOV, levels = c(names(clrs$consensusOV)))) %>%
  arrange(consensusOV, nrel)

heat_layers <- list(
  geom_tile(),
  scale_fill_gradientn(colors = magma(9)[c(-1)], limits = c(0, 100), breaks = c(0, 50, 100)),
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "bottom",
    aspect.ratio = 4/9
    ),
  guides(fill = guide_colorbar(title.position = "top"))
)

tcga_heat <- consOV_tbl_summary %>%
  filter(consOV_cell_type != "Other") %>%
  mutate(cell_type = ordered(consOV_cell_type, levels = unique(c("Myeloid cell", "Fibroblast", "Cycling cell", names(clrs$cell_type))))) %>%
  mutate(consensusOV = ordered(consensusOV, levels = names(clrs$consensusOV))) %>%
  ggplot(aes(cell_type, consensusOV, fill = nrel)) +
  labs(fill = "Fraction\nof cells [%]") +
  heat_layers

tcga_heat_wide <- consOV_tbl_summary %>%
  filter(consOV_cell_type != "Other") %>%
  mutate(cell_type = ordered(consOV_cell_type, levels = unique(c("Myeloid cell", "Fibroblast", "Cycling cell", names(clrs$cell_type))))) %>%
  mutate(consensusOV = ordered(consensusOV, levels = names(clrs$consensusOV))) %>%
  ggplot(aes(consensusOV, cell_type, fill = nrel)) +
  labs(fill = "Fraction\nof cells [%]") +
  heat_layers +
  theme(aspect.ratio = 9/4, legend.position = "right")

# tcga_heat + guides(fill = F)
tcga_heat_legend <- get_legend(tcga_heat)

print(consOV_tbl_summary, n = 36)
```

```
## # A tibble: 36 x 4
##    consOV_cell_type consensusOV         n   nrel
##    <ord>            <ord>           <int>  <dbl>
##  1 Fibroblast       Immunoreactive    809  0.499
##  2 Ov cancer cell   Immunoreactive   3125  1.24 
##  3 Endothelial cell Immunoreactive    260  1.40 
##  4 Plasma cell      Immunoreactive   1129  5.39 
##  5 Mast cell        Immunoreactive     99  7.88 
##  6 B cell           Immunoreactive   1857  9.94 
##  7 T cell           Immunoreactive  43985 17.6  
##  8 Dendritic cell   Immunoreactive   1038 21.6  
##  9 Myeloid cell     Immunoreactive 154098 76.6  
## 10 Ov cancer cell   Mesenchymal      1366  0.542
## 11 Mast cell        Mesenchymal        13  1.03 
## 12 Dendritic cell   Mesenchymal        63  1.31 
## 13 B cell           Mesenchymal       724  3.88 
## 14 Endothelial cell Mesenchymal       950  5.13 
## 15 T cell           Mesenchymal     13311  5.32 
## 16 Myeloid cell     Mesenchymal     12984  6.45 
## 17 Plasma cell      Mesenchymal      1550  7.40 
## 18 Fibroblast       Mesenchymal    101439 62.6  
## 19 Myeloid cell     Proliferative    1329  0.660
## 20 Mast cell        Proliferative      50  3.98 
## 21 Dendritic cell   Proliferative     375  7.79 
## 22 T cell           Proliferative   21149  8.45 
## 23 Fibroblast       Proliferative   19280 11.9  
## 24 Ov cancer cell   Proliferative   49008 19.5  
## 25 Plasma cell      Proliferative    4117 19.7  
## 26 B cell           Proliferative    4339 23.2  
## 27 Endothelial cell Proliferative    4388 23.7  
## 28 Myeloid cell     Differentiated  32806 16.3  
## 29 Fibroblast       Differentiated  40550 25.0  
## 30 B cell           Differentiated  11753 62.9  
## 31 Plasma cell      Differentiated  14148 67.6  
## 32 T cell           Differentiated 171890 68.7  
## 33 Dendritic cell   Differentiated   3338 69.3  
## 34 Endothelial cell Differentiated  12933 69.8  
## 35 Ov cancer cell   Differentiated 198338 78.8  
## 36 Mast cell        Differentiated   1095 87.1
```

## Figure 2 upper panel grid

### v7.1


```r
arrow <- arrow(angle = 20, type = "closed", length = unit(0.1, "npc"))
umap_coord_anno <- ggplot(tibble(group = c("UMAP1", "UMAP2"),
                                 x = c(0, 0), xend = c(1, 0),
                                 y = c(0, 0), yend = c(0, 1),
                                 lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                                 angle = c(0, 90))) +
  geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
               arrow = arrow, size = 1, lineend = "round") +
  geom_text(aes(lx, ly, label = group, angle = angle), size = 4) +
  theme_void() +
  coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))

add_umap_coord <- function(gg_obj) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno, x = -0.015, y = -0.02, width = 0.2, height = 0.2)
  return(p)
}

full_grid1 <- ggdraw() +
  ## left
  draw_label("Patient", x = 0.25, y = 0.98, size = 20) +
  draw_plot(add_umap_coord(umap_patient), x = 0, y = 0, width = 0.5, height = 1) +
  draw_label(paste0("n = ", nrow(seu_tbl_full)), 
             x = 0.25, y = 0.03, size = 16, hjust = 0) + 
  ## right upper
  draw_plot(mixing_plot_cell_type, x = 0.5+0.5*0, y = 0.66, width = 0.5*0.3, height = 0.33) +
  draw_label("Cell type", x = 0.71, y = 0.98, size = 20, hjust = 0) +
  draw_plot(umap_cell_type_mini, x = 0.83, y = 0.68, width = 0.5*1/3, height = 0.31) +
  draw_plot(numbers_plot_cell_type, x = 0.7, y = 0.651, width = 0.13, height = 0.315) +
  ## right middle
  draw_label("TCGA", x = 0.53+0.5*0, y = 0.6, size = 20, hjust = 0) +
  draw_plot(umap_consOV_mini, x = 0.5+0.5*0, y = 0.3, width = 0.5*1/3, height = 0.31) +
  draw_label("Signature", x = 0.53+0.5*0.33, y = 0.6, size = 20, hjust = 0) +
  draw_plot(umap_mutsig_mini, x = 0.5+0.5*0.33, y = 0.3, width = 0.5*1/3, height = 0.31) +
  draw_label("Site", x = 0.53+0.5*0.66, y = 0.6, size = 20, hjust = 0) + 
  draw_plot(umap_site_mini, x = 0.5+0.5*0.66, y = 0.3, width = 0.5*1/3, height = 0.31) +
  ## right lower
  draw_grob(umap_consOV_legend, x = 0.52+0.5*0, y = -0.67, hjust = 0, vjust = 0) +
  draw_plot(tcga_heat + guides(fill = F), x = 0.49+0.5*0, y = -0.03, width = 0.18, height = 0.3) +
  draw_grob(tcga_heat_legend, x = 0.68, y = -0.35) +
  draw_grob(umap_mutsig_legend, x = 0.52+0.5*0.33, y = -0.67, hjust = 0, vjust = 0) +
  draw_grob(umap_site_legend, x = 0.52+0.5*0.66, y = -0.67, hjust = 0, vjust = 0)

full_grid1
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_060-1.png" width="1920" />

```r
full_grid1_rds <- list(rds = full_grid1, width = 20, height = 10)
# write_rds(full_grid1, "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.rds")
ggsave(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.png", full_grid1, width = 20, height = 10)
ggsave_nodb(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.pdf", full_grid1, width = 20, height = 10)
```

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.pdf](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.pdf)

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.png](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.1.png)

### v7.2


```r
full_grid2 <- ggdraw() +
  ## left
  draw_label("Patient", x = 0.25, y = 0.98, size = 20) +
  draw_plot(add_umap_coord(umap_patient), x = 0, y = 0, width = 0.5, height = 1) +
  draw_label(paste0("n = ", nrow(seu_tbl_full)), 
             x = 0.25, y = 0.03, size = 16, hjust = 0) + 
  ## right upper
  draw_plot(mixing_plot_cell_type, x = 0.5+0.5*0, y = 0.66, width = 0.5*0.3, height = 0.33) +
  draw_label("Cell type", x = 0.725, y = 0.98, size = 20, hjust = 0) +
  draw_plot(umap_cell_type_mini, x = 0.83, y = 0.68, width = 0.5*1/3, height = 0.31) +
  draw_plot(numbers_plot_cell_type, x = 0.715, y = 0.651, width = 0.13, height = 0.315) +
  ## right middle
  draw_label("TCGA", x = 0.53+0.5*0, y = 0.6, size = 20, hjust = 0) +
  draw_plot(umap_consOV_mini, x = 0.5+0.5*0, y = 0.3, width = 0.5*1/3, height = 0.31) +
  draw_label("Site", x = 0.725, y = 0.6, size = 20, hjust = 0) +
  draw_plot(numbers_plot_tumor_supersite, x = 0.715, y = 0.34, width = 0.13, height = 0.25) +
  draw_plot(umap_site_mini, x = 0.5+0.5*0.66, y = 0.3, width = 0.5*1/3, height = 0.31) +
  ## right lower
  draw_plot(tcga_heat + guides(fill = F), x = 0.49+0.5*0, y = -0.03, width = 0.18, height = 0.3) +
  draw_grob(tcga_heat_legend, x = 0.68, y = -0.35) +
  draw_grob(umap_consOV_legend, x = 0.52+0.5*0, y = -0.67, hjust = 0, vjust = 0)

full_grid2
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_065-1.png" width="1920" />

```r
full_grid2_rds <- list(rds = full_grid2, width = 20, height = 10)
# write_rds(full_grid2, "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.rds")
ggsave(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.png", full_grid2, width = 20, height = 10)
ggsave_nodb(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.pdf", full_grid2, width = 20, height = 10)
```

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.pdf](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.pdf)

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.png](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.2.png)

### v7.3


```r
full_grid3 <- ggdraw() +
  ## left
  draw_label("Patient", x = 0.25, y = 0.98, size = 20) +
  draw_plot(add_umap_coord(umap_patient), x = 0, y = 0, width = 0.5, height = 1) +
  draw_label(paste0("n = ", nrow(seu_tbl_full)), 
             x = 0.25, y = 0.03, size = 16, hjust = 0) + 
  ## right upper
  draw_label("Cell type", x = 0.53+0.5*0, y = 0.98, size = 20, hjust = 0) +
  draw_plot(numbers_plot_cell_type, x = 0.52+0.5*0 , y = 0.651, width = 0.13, height = 0.315) +
  draw_plot(umap_cell_type_mini, x = 0.64, y = 0.68, width = 0.5*1/3, height = 0.31) +
  draw_plot(mixing_plot_cell_type, x = 0.83, y = 0.66, width = 0.5*0.3, height = 0.33) +
  ## right middle
  draw_label("TCGA", x = 0.53+0.5*0, y = 0.6, size = 20, hjust = 0) +
  draw_plot(umap_consOV_mini, x = 0.5+0.5*0, y = 0.3, width = 0.5*1/3, height = 0.31) +
  draw_label("Site", x = 0.725, y = 0.6, size = 20, hjust = 0) +
  draw_plot(numbers_plot_tumor_supersite, x = 0.715, y = 0.34, width = 0.13, height = 0.25) +
  draw_plot(umap_site_mini, x = 0.5+0.5*0.66, y = 0.3, width = 0.5*1/3, height = 0.31) +
  ## right lower
  draw_plot(tcga_heat + guides(fill = F), x = 0.49+0.5*0, y = -0.03, width = 0.18, height = 0.3) +
  draw_grob(tcga_heat_legend, x = 0.68, y = -0.35) +
  draw_grob(umap_consOV_legend, x = 0.52+0.5*0, y = -0.67, hjust = 0, vjust = 0)

full_grid3
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_070-1.png" width="1920" />

```r
full_grid3_rds <- list(rds = full_grid3, width = 20, height = 10)
# write_rds(full_grid3, "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.rds")
ggsave(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.png", full_grid3, width = 20, height = 10)
ggsave_nodb(filename = "figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.pdf", full_grid3, width = 20, height = 10)
```

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.pdf](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.pdf)

[figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.png](figures/200_cohort_plotting_umap_v7/002_cohort_umaps_v7.3.png)

## S1


```r
## QC histograms
qc_tbl <- seu_tbl_full %>%
  mutate(`log2(#UMIs)` = log2(nCount_RNA),
         `log2(#Genes)` = log2(nFeature_RNA),
         `log2(#UMIs/#Genes)` = `log2(#UMIs)`-`log2(#Genes)`,
         `MT reads [%]` = percent.mt) %>%
  select(`log2(#UMIs)`,
         `log2(#Genes)`,
         # `log2(#UMIs/#Genes)`,
         `MT reads [%]`, cell_type) %>%
  gather(key, value, -cell_type) %>%
  group_by(key) %>%
  mutate(median_value = median(value))

qc_tbl <- bind_rows(qc_tbl, mutate(qc_tbl, cell_type = "All cell types")) %>%
  mutate(cell_type = ordered(cell_type,
                             levels = rev(c("All cell types", names(clrs$cell_type)))))

qc_plot <- ggplot(qc_tbl) +
  geom_violin(aes(cell_type, value, fill = cell_type),
              color = "white", width = 1) +
  geom_boxplot(aes(cell_type, value, color = cell_type),
               width = 0.15, size = 1, fill = NA, outlier.shape = NA) +
  geom_boxplot(aes(cell_type, value), fill = NA,
               width = 0.15, color = "white", outlier.shape = NA) +
  geom_hline(aes(yintercept = median_value),
             data = distinct(qc_tbl, key, median_value)) +
  facet_wrap(~key, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c(clrs$cell_type, `All cell types` = "grey10"),
                    guide = F) +
  scale_color_manual(values = c(clrs$cell_type, `All cell types` = "grey10"),
                     guide = F) +
  coord_flip() +
  theme(axis.title = element_blank(),
        strip.placement = "outside")
```


```r
custom_guide <- guide_legend(override.aes = list(size = 4, alpha = 1), 
                             ncol = 1, label.position = "right")

mutsig_legend <- cowplot::get_legend(umap_mutsig + theme(legend.title=element_text()) + labs(color = "Mutational\nsignature") + guides(color = custom_guide))
consOV_legend <- cowplot::get_legend(umap_consOV + theme(legend.title=element_text()) + labs(color = "TCGA subtype") + guides(color = custom_guide))
cell_type_legend <- cowplot::get_legend(umap_cell_type + theme(legend.title=element_text(), legend.text = element_text(size = 14, margin = ggplot2::margin(0, 10, 0, 0), color = "black")) + labs(color = "Cell type") + guides(color = custom_guide))
site_legend <- cowplot::get_legend(umap_site + theme(legend.title=element_text()) + labs(color = "Site") + guides(color = custom_guide))

figS1_panel <- ggdraw() +
  # draw_plot(add_umap_coord(umap_mutsig + guides(color = F) + labs(title = "Mutational signature")), width = 0.24, height = 0.48, x = 0.01, y = 0.5) +
  draw_plot(add_umap_coord(umap_umis), width = 0.24, height = 0.48, x = 0.01, y = 0.5) +
  draw_plot(add_umap_coord(umap_genes), width = 0.24, height = 0.48, x = 0.25, y = 0.5) +
  draw_plot(add_umap_coord(umap_mitos), width = 0.24, height = 0.48, x = 0.01, y = 0) +
  draw_plot(add_umap_coord(umap_phase), width = 0.24, height = 0.48, x = 0.25, y = 0) +
  draw_plot(qc_plot, width = 0.5, height = 0.48, x = 0.5, y = 0.5) +
  draw_grob(mutsig_legend, x = 0.55, y = 0.48, hjust = 0, vjust = 1) +
  draw_grob(consOV_legend, x = 0.55, y = 0.22, hjust = 0, vjust = 1) +
  draw_grob(cell_type_legend, x = 0.7, y = 0.48, hjust = 0, vjust = 1) +
  draw_grob(site_legend, x = 0.85, y = 0.48, hjust = 0, vjust = 1) +
  draw_label("A", x = 0.01, y = 0.98, fontface = "bold", size = 22) +
  draw_label("B", x = 0.52, y = 0.98, fontface = "bold", size = 22) +
  draw_label("C", x = 0.52, y = 0.48, fontface = "bold", size = 22)

figS1_panel
```

<img src="200_cohort_plotting_umap_v7_files/figure-html/chunk_105-1.png" width="1536" />

```r
figS1_panel_rds <- list(rds = figS1_panel, width = 16, height = 9)
ggsave("figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.png", figS1_panel, width = 16, height = 9)
ggsave_nodb("figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.pdf", figS1_panel, width = 16, height = 9)
# write_rds(figS1_panel_rds, "figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.rds")
```

[figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.pdf](figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.pdf)

[figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.png](figures/200_cohort_plotting_umap_v7/002_supplementary_umap_v7.png)

## session info 


```r
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 3.6.2 (2019-12-12)
##  os       Debian GNU/Linux 10 (buster)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       Etc/UTC                     
##  date     2021-05-06                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package      * version    date       lib
##    abind          1.4-5      2016-07-21 [1]
##    ape            5.3        2019-03-17 [2]
##    assertthat     0.2.1      2019-03-21 [2]
##    backports      1.1.10     2020-09-15 [1]
##    beeswarm       0.2.3      2016-04-25 [2]
##    bibtex         0.4.2.2    2020-01-02 [2]
##    Biobase        2.46.0     2019-10-29 [2]
##    BiocGenerics   0.32.0     2019-10-29 [2]
##    bookdown       0.22       2021-04-22 [1]
##    broom          0.7.6      2021-04-05 [1]
##    bslib          0.2.4      2021-01-25 [1]
##    Cairo          1.5-12.2   2020-07-07 [1]
##    callr          3.6.0      2021-03-28 [1]
##    car            3.0-8      2020-05-21 [1]
##    carData        3.0-4      2020-05-22 [1]
##    cellranger     1.1.0      2016-07-27 [1]
##    cli            2.5.0      2021-04-26 [1]
##    cluster        2.1.1      2021-02-14 [1]
##    codetools      0.2-18     2020-11-04 [1]
##    colorblindr  * 0.1.0      2020-01-13 [2]
##    colorspace   * 2.0-0      2020-11-11 [1]
##    cowplot      * 1.1.1      2020-12-30 [1]
##    crayon         1.4.1      2021-02-08 [1]
##    curl           4.3        2019-12-02 [2]
##    data.table     1.14.0     2021-02-21 [1]
##    DBI            1.1.0      2019-12-15 [2]
##    dbplyr         2.0.0      2020-11-03 [1]
##    desc           1.3.0      2021-03-05 [1]
##    devtools       2.2.1      2019-09-24 [2]
##    digest         0.6.25     2020-02-23 [1]
##    dplyr        * 1.0.2      2020-08-18 [1]
##    ellipsis       0.3.1      2020-05-15 [1]
##    evaluate       0.14       2019-05-28 [1]
##    fansi          0.4.1      2020-01-08 [2]
##    farver         2.0.3      2020-01-16 [1]
##    fitdistrplus   1.0-14     2019-01-23 [2]
##    forcats      * 0.5.1      2021-01-27 [1]
##    foreign        0.8-74     2019-12-26 [3]
##    fs             1.5.0      2020-07-31 [1]
##    future         1.15.1     2019-11-25 [2]
##    future.apply   1.4.0      2020-01-07 [2]
##    gbRd           0.4-11     2012-10-01 [2]
##    generics       0.1.0      2020-10-31 [1]
##    ggbeeswarm     0.6.0      2017-08-07 [2]
##    ggplot2      * 3.3.3      2020-12-30 [1]
##    ggpubr       * 0.4.0      2020-06-27 [1]
##    ggrastr        0.1.9      2020-06-20 [1]
##    ggrepel        0.9.1      2021-01-15 [1]
##    ggridges       0.5.2      2020-01-12 [2]
##    ggsignif       0.6.0      2019-08-08 [1]
##    globals        0.12.5     2019-12-07 [2]
##    glue           1.3.2      2020-03-12 [1]
##    gridExtra      2.3        2017-09-09 [1]
##    gtable         0.3.0      2019-03-25 [1]
##    haven          2.3.1      2020-06-01 [1]
##    highr          0.8        2019-03-20 [1]
##    hms            1.0.0      2021-01-13 [1]
##    htmltools      0.5.1.1    2021-01-22 [1]
##    htmlwidgets    1.5.3      2020-12-10 [1]
##    httr           1.4.2      2020-07-20 [1]
##    ica            1.0-2      2018-05-24 [2]
##    igraph         1.2.6      2020-10-06 [1]
##    irlba          2.3.3      2019-02-05 [2]
##    jquerylib      0.1.3      2020-12-17 [1]
##    jsonlite       1.7.2      2020-12-09 [1]
##    KernSmooth     2.23-18    2020-10-29 [1]
##    knitr        * 1.31       2021-01-27 [1]
##    labeling       0.4.2      2020-10-20 [1]
##    lattice        0.20-41    2020-04-02 [1]
##    lazyeval       0.2.2      2019-03-15 [1]
##    leiden         0.3.1      2019-07-23 [2]
##    lifecycle      0.2.0      2020-03-06 [1]
##    listenv        0.8.0      2019-12-05 [2]
##    lmtest         0.9-37     2019-04-30 [2]
##    lsei           1.2-0      2017-10-23 [2]
##    lubridate      1.7.10     2021-02-26 [1]
##  P magick       * 2.4.0      2020-06-23 [?]
##    magrittr     * 2.0.1      2020-11-17 [1]
##    MASS           7.3-53.1   2021-02-12 [1]
##    Matrix         1.3-2      2021-01-06 [1]
##    memoise        1.1.0      2017-04-21 [2]
##    metap          1.2        2019-12-08 [2]
##    mnormt         1.5-5      2016-10-15 [2]
##    modelr         0.1.8      2020-05-19 [1]
##    multcomp       1.4-12     2020-01-10 [2]
##    multtest       2.42.0     2019-10-29 [2]
##    munsell        0.5.0      2018-06-12 [1]
##    mutoss         0.1-12     2017-12-04 [2]
##    mvtnorm        1.0-12     2020-01-09 [2]
##    nlme           3.1-152    2021-02-04 [1]
##    npsurv         0.4-0      2017-10-14 [2]
##    numDeriv       2016.8-1.1 2019-06-06 [1]
##    openxlsx       4.1.5      2020-05-06 [1]
##    pbapply        1.4-2      2019-08-31 [2]
##    pillar         1.6.0      2021-04-13 [1]
##    pkgbuild       1.0.6      2019-10-09 [2]
##    pkgconfig      2.0.3      2019-09-22 [1]
##    pkgload        1.0.2      2018-10-29 [2]
##    plotly         4.9.3      2021-01-10 [1]
##    plotrix        3.7-7      2019-12-05 [2]
##    plyr           1.8.6      2020-03-03 [1]
##    png            0.1-7      2013-12-03 [1]
##    prettyunits    1.1.1      2020-01-24 [1]
##    processx       3.5.0      2021-03-23 [1]
##    ps             1.3.2      2020-02-13 [1]
##    purrr        * 0.3.4      2020-04-17 [1]
##    R.methodsS3    1.7.1      2016-02-16 [2]
##    R.oo           1.23.0     2019-11-03 [2]
##    R.utils        2.9.2      2019-12-08 [2]
##    R6             2.4.1      2019-11-12 [1]
##    RANN           2.6.1      2019-01-08 [2]
##    rappdirs       0.3.3      2021-01-31 [1]
##    RColorBrewer   1.1-2      2014-12-07 [1]
##    Rcpp           1.0.4      2020-03-17 [1]
##    RcppAnnoy      0.0.16     2020-03-08 [1]
##    RcppParallel   4.4.4      2019-09-27 [2]
##    Rdpack         0.11-1     2019-12-14 [2]
##    readr        * 1.4.0      2020-10-05 [1]
##    readxl       * 1.3.1      2019-03-13 [1]
##    remotes        2.3.0      2021-04-01 [1]
##    reprex         2.0.0      2021-04-02 [1]
##    reshape2       1.4.4      2020-04-09 [1]
##    reticulate     1.14       2019-12-17 [2]
##    rio            0.5.16     2018-11-26 [1]
##    rlang          0.4.8      2020-10-08 [1]
##    rmarkdown      2.7        2021-02-19 [1]
##    ROCR           1.0-11     2020-05-02 [1]
##    rprojroot      2.0.2      2020-11-15 [1]
##    rstatix        0.6.0      2020-06-18 [1]
##    rstudioapi     0.13       2020-11-12 [1]
##    rsvd           1.0.3      2020-02-17 [1]
##    Rtsne          0.15       2018-11-10 [2]
##    rvest          0.3.6      2020-07-25 [1]
##    sandwich       2.5-1      2019-04-06 [2]
##    sass           0.3.1      2021-01-24 [1]
##    scales         1.1.1      2020-05-11 [1]
##    sctransform    0.2.1      2019-12-17 [2]
##    SDMTools       1.1-221.2  2019-11-30 [2]
##    sessioninfo    1.1.1      2018-11-05 [2]
##    Seurat       * 3.1.2      2019-12-12 [2]
##    sn             1.5-4      2019-05-14 [2]
##    stringi        1.5.3      2020-09-09 [1]
##    stringr      * 1.4.0      2019-02-10 [1]
##    survival       3.2-10     2021-03-16 [1]
##    testthat       2.3.2      2020-03-02 [1]
##    TFisher        0.2.0      2018-03-21 [2]
##    TH.data        1.0-10     2019-01-21 [2]
##    tibble       * 3.0.4      2020-10-12 [1]
##    tidyr        * 1.1.2      2020-08-27 [1]
##    tidyselect     1.1.0      2020-05-11 [1]
##    tidyverse    * 1.3.0      2019-11-21 [1]
##    tsne           0.1-3      2016-07-15 [2]
##    usethis        1.5.1      2019-07-04 [2]
##    utf8           1.1.4      2018-05-24 [2]
##    uwot           0.1.5      2019-12-04 [2]
##    vctrs          0.3.5      2020-11-17 [1]
##    vipor          0.4.5      2017-03-22 [2]
##    viridis      * 0.6.0      2021-04-15 [1]
##    viridisLite  * 0.4.0      2021-04-13 [1]
##    withr          2.3.0      2020-09-22 [1]
##    xfun           0.22       2021-03-11 [1]
##    xml2           1.3.2      2020-04-23 [1]
##    yaml           2.2.1      2020-02-01 [1]
##    zip            2.1.1      2020-08-27 [1]
##    zoo            1.8-7      2020-01-10 [2]
##  source                                 
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  Bioconductor                           
##  Bioconductor                           
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  Github (clauswilke/colorblindr@1ac3d4d)
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  Bioconductor                           
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
## 
## [1] /home/uhlitzf/R/lib
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
## 
##  P ── Loaded and on-disk path mismatch.
```

