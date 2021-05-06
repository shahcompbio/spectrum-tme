---
title: "MSK SPECTRUM freeze Macrophages"
author: "Florian Uhlitz"
date: "2021-05-06"
output: 
  html_document: 
    df_print: kable
    number_sections: yes
    toc: yes
    toc_float: true
    code_folding: hide
editor_options: 
  chunk_output_type: console
params: 
  cell_type_super: "Myeloid.super"
  cell_type_major: ["Myeloid.cell", "Dendritic.cell", "Mast.cell"]
  cell_sort: "CD45+"
  louvain_resolution: 0.3
  louvain_cluster: "RNA_snn_res.0.3"
  pcut: 0.05
---



# Figure 6 lower


```r
library(magrittr)
library(tidyverse)
library(Seurat)
library(readxl)
library(cowplot)
library(colorblindr)
library(viridis)
library(ComplexHeatmap)
library(grid)
library(ggpubr)

coi <- params$cell_type_super
cell_sort <- params$cell_sort
cell_type_major <- params$cell_type_major
louvain_resolution <- params$louvain_resolution
louvain_cluster <- params$louvain_cluster
pcut <- params$pcut
```


```r
## load global vars: 
source("src/global_vars.R")

# scrna_meta_tbl
# clrs
# markers_v7
# markers_v7_super
# cell_type_super_lookup
```


```r
seu_obj_ml <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Myeloid.super_processed_filtered_annotated.rds")

seu_obj_mp <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Macrophages_processed.rds")

seu_obj_dc <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/DCs_processed.rds")

seu_obj_cc <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Ovarian.cancer.super_processed_filtered_annotated.rds")

seu_obj_tc <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/T.super_processed_filtered_annotated.rds")


gois_mp <- c("S100A8", "S100A9", "CXCL10", "CD274", "C1QC", "FN1", "MARCO", "SELENOP")
gois_dc <- grep("^IRF", rownames(seu_obj_dc), value = T)
gois <- c(gois_mp, gois_dc)

seu_obj_dc <- AddModuleScore(seu_obj_dc, list(IRF.genes.module = c("IRF1", "IRF7", "IRF8", "IRF9")), name = "IRF.genes.module")
seu_obj_dc$IRF.genes.module <- seu_obj_dc$IRF.genes.module1
seu_obj_dc$IRF.genes.module1 <- NULL

module_names <- grep("pathway|CD8|module", colnames(seu_obj_dc@meta.data), value = T)

module_names_sub <- c("MARCO", "CXCL10", "CD274", "JAK.STAT.pathway")

module_names_dc <- c("Interferon.regulators.module", "IRF.genes.module", "IRF1", "IRF7", "IRF8", "IRF9")

myfeatures <- c("cell_id", "umapharmony_1", "umapharmony_2", "sample", "cell_type", "cluster_label", "Phase")

my_subtypes <- names(clrs$cluster_label[[coi]])

plot_data <- as_tibble(FetchData(seu_obj_mp, c(myfeatures, module_names, gois))) %>%
  gather(key, value, -c(1:length(myfeatures))) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>%
  mutate(cell_type_super = cell_type_super_lookup[cell_type]) %>%
  mutate(sort_short_x = ifelse(sort_short == "U" & cell_type_super == "Immune",
                               "CD45+", ifelse(sort_short == "U" & cell_type_super == "Stromal", "CD45-", sort_short))) %>% 
  mutate(cluster_label = ordered(cluster_label, levels = my_subtypes),
         goi_expressed = value > 0) %>% 
  group_by(sample, key) %>% 
  mutate(n = length(sample)) %>% 
  mutate(pct.expr = sum(goi_expressed)/n*100) %>% 
  ungroup

plot_data_uniq <- distinct(plot_data, cell_id, .keep_all = T)

plot_data_dc <- as_tibble(FetchData(seu_obj_dc, c(myfeatures, module_names, gois))) %>%
  gather(key, value, -c(1:length(myfeatures))) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>%
  mutate(cell_type_super = cell_type_super_lookup[cell_type]) %>%
  mutate(sort_short_x = ifelse(sort_short == "U" & cell_type_super == "Immune",
                               "CD45+", ifelse(sort_short == "U" & cell_type_super == "Stromal", "CD45-", sort_short))) %>% 
  mutate(cluster_label = ordered(cluster_label, levels = my_subtypes),
         goi_expressed = value > 0) %>% 
  group_by(sample, key) %>% 
  mutate(n = length(sample)) %>% 
  mutate(pct.expr = sum(goi_expressed)/n*100) %>% 
  ungroup

plot_data_uniq_dc <- distinct(plot_data_dc, cell_id, .keep_all = T)
```


## Macrophage umap 


```r
dc_plot_cluster <- ggplot(plot_data_uniq) +
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = cluster_label), 
                           size = 0.01, alpha = 1, raster.dpi = 150) +
  # theme_void() +
  scale_color_manual(values = clrs$cluster_label[[coi]]) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  labs(color = "Cluster")

plot_data_gois <- plot_data %>% 
  filter(key %in% gois_mp) %>% 
  mutate(key = ordered(key, levels = gois_mp)) %>% 
  mutate(value = ifelse(value > 5, 5, value))

dc_plot_genes <- ggplot() +
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2), 
                           size = 0.1, alpha = 1, color = "grey80", raster.dpi = 75,
                           data = filter(plot_data_gois, value == 0)) +
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = value), 
                           size = 0.1, alpha = 1, raster.dpi = 75,
                           data = filter(plot_data_gois, value > 0)) +
  facet_wrap(~key, nrow = 2) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +  
  scale_color_gradientn(colors = viridis(9), breaks = c(0, 5), 
                        labels = c("low", "high")) +
  theme(aspect.ratio = 1) +
  labs(color = "Expression")

dc_plot_genes
```

<img src="620_Macrophage_umaps_and_boxplots_files/figure-html/chunk_120-1.png" width="768" />

```r
ggsave("figures/620_Macrophage_umaps_and_boxplots/005_macrophage_gois_umap.pdf", dc_plot_genes, width = 8, height = 4)
ggsave("figures/620_Macrophage_umaps_and_boxplots/005_macrophage_gois_umap.png", dc_plot_genes, width = 8, height = 4)
```

## boxplots

### macrophages


```r
mutsig_boxplot_data <- plot_data %>%
  filter(key %in% module_names_sub) %>%
  filter(!(consensus_signature %in% c("Undetermined", "HRD"))) %>%
  filter(tumor_supersite != "Ascites") %>%
  mutate(key = str_replace_all(key, "JAK.STAT.pathway", "JAK/STAT\npathway"))

mutsig_boxplot_data_stats <- mutsig_boxplot_data %>%
  group_by(sample, key) %>%
  mutate(value = mean(value),
         pct.expr = mean(pct.expr)) %>%
  distinct(sample, key, .keep_all = T)

mutsig_boxplot <- mutsig_boxplot_data %>%
  filter(key != "JAK/STAT\npathway") %>%
  ggplot(aes(consensus_signature, pct.expr)) +
  geom_violin(aes(consensus_signature, pct.expr, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
  geom_boxplot(aes(consensus_signature, pct.expr, color = consensus_signature),
               width = 0.5, size = 0.75, outlier.shape = NA) +
  geom_boxplot(aes(consensus_signature, pct.expr, fill = consensus_signature),
               color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
  # geom_boxplot(aes(consensus_signature, pct.expr, color = consensus_signature),
  #              fill = "white") +
  stat_compare_means(ref.group = ".all.",
                     data = filter(mutsig_boxplot_data_stats,
                                   key != "JAK/STAT\npathway"),
                     label = "p.signif", label.y = 85, hide.ns = T) +
  stat_compare_means(label.y = 95, label.x = 1.5, label.sep = "\n\n",
                     data = filter(mutsig_boxplot_data_stats,
                                   key != "JAK/STAT\npathway")) +
  facet_wrap(~key, ncol = 5) +
  scale_color_manual(values = clrs$consensus_signature) +
  scale_fill_manual(values = clrs$consensus_signature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) +
  # coord_cartesian(ylim = c(-1, 2)) +
  labs(y = "% Pos. Macrophages", color = "Mutational\nsignature",
       fill = "Mutational\nsignature")

mutsig_boxplot_jakstat <- mutsig_boxplot_data %>%
  filter(key == "JAK/STAT\npathway") %>%
  ggplot(aes(consensus_signature, value)) +
  geom_violin(aes(consensus_signature, value, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
  geom_boxplot(aes(consensus_signature, value, color = consensus_signature),
               width = 0.5, size = 0.75, outlier.shape = NA) +
  geom_boxplot(aes(consensus_signature, value, fill = consensus_signature),
               color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
  # geom_boxplot(aes(consensus_signature, value, color = consensus_signature),
  #              fill = "white") +
  stat_compare_means(ref.group = ".all.",
                     data = filter(mutsig_boxplot_data_stats,
                                   key == "JAK/STAT\npathway"),
                     label = "p.signif", label.y = 4.2, hide.ns = T) +
  stat_compare_means(label.y = 4.8, label.x = 1.5, label.sep = "\n\n",
                     data = filter(mutsig_boxplot_data_stats,
                                   key == "JAK/STAT\npathway")) +
  facet_wrap(~key, ncol = 5) +
  scale_color_manual(values = clrs$consensus_signature) +
  scale_fill_manual(values = clrs$consensus_signature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) +
  coord_cartesian(ylim = c(-1, 5)) +
  labs(y = "PROGENy score", color = "Mutational\nsignature",
       fill = "Mutational\nsignature")

mutsig_boxplot_grid <- plot_grid(
  mutsig_boxplot + guides(color = F, fill = F),
  mutsig_boxplot_jakstat + guides(color = F, fill = F),
  nrow = 1, rel_widths = c(0.72, 0.28))

mutsig_boxplot_grid
```

<img src="620_Macrophage_umaps_and_boxplots_files/figure-html/chunk_140-1.png" width="624" />

```r
ggsave("figures/620_Macrophage_umaps_and_boxplots/005_macrophage_boxplots.pdf", mutsig_boxplot_grid, width = 6.5, height = 3.5)
```

### DCs


```r
mutsig_boxplot_data_dc <- plot_data_dc %>%
  filter(key %in% module_names_dc) %>%
  filter(!(consensus_signature %in% c("Undetermined", "HRD"))) %>%
  filter(tumor_supersite != "Ascites") %>%
  mutate(key = str_replace_all(key, "Interferon.regulators.module", "IFN\nregulators"),
         key = str_replace_all(key, "IRF.genes.module", "All IRFs"))

mutsig_boxplot_data_stats_dc <- mutsig_boxplot_data_dc %>%
  group_by(sample, key) %>%
  mutate(value = mean(value),
         pct.expr = mean(pct.expr)) %>%
  distinct(sample, key, .keep_all = T)

mutsig_boxplot_dc <- mutsig_boxplot_data_dc %>%
  filter(key %in% gois_dc, key %in% module_names_dc) %>%
  ggplot(aes(consensus_signature, pct.expr)) +
  geom_violin(aes(consensus_signature, pct.expr, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
  geom_boxplot(aes(consensus_signature, pct.expr, color = consensus_signature),
               width = 0.5, size = 0.75, outlier.shape = NA) +
  geom_boxplot(aes(consensus_signature, pct.expr, fill = consensus_signature),
               color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
  # geom_boxplot(aes(consensus_signature, pct.expr, color = consensus_signature),
  #              fill = "white") +
  stat_compare_means(ref.group = ".all.",
                     data = filter(mutsig_boxplot_data_stats_dc,
                                   key %in% gois_dc, key %in% module_names_dc),
                     label = "p.signif", label.y = 85, hide.ns = T) +
  stat_compare_means(label.y = 95, label.x = 1.5, label.sep = "\n\n",
                     data = filter(mutsig_boxplot_data_stats_dc,
                                   key %in% gois_dc, key %in% module_names_dc)) +
  facet_wrap(~key, ncol = 5) +
  scale_color_manual(values = clrs$consensus_signature) +
  scale_fill_manual(values = clrs$consensus_signature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) +
  # coord_cartesian(ylim = c(-1, 2)) +
  labs(y = "% Pos. DCs", color = "Mutational\nsignature",
       fill = "Mutational\nsignature")

mutsig_boxplot_ifnreg <- mutsig_boxplot_data_dc %>%
  filter(key == "IFN\nregulators") %>%
  ggplot(aes(consensus_signature, value)) +
  geom_violin(aes(consensus_signature, value, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
  geom_boxplot(aes(consensus_signature, value, color = consensus_signature),
               width = 0.5, size = 0.75, outlier.shape = NA) +
  geom_boxplot(aes(consensus_signature, value, fill = consensus_signature),
               color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
  # geom_boxplot(aes(consensus_signature, value, color = consensus_signature),
  #              fill = "white") +
  stat_compare_means(ref.group = ".all.",
                     data = filter(mutsig_boxplot_data_stats_dc,
                                   key == "IFN\nregulators"),
                     label = "p.signif", label.y = 1.4, hide.ns = T) +
  stat_compare_means(label.y = 1.9, label.x = 1.5, label.sep = "\n\n",
                     data = filter(mutsig_boxplot_data_stats_dc,
                                   key == "IFN\nregulators")) +
  facet_wrap(~key, ncol = 5) +
  scale_color_manual(values = clrs$consensus_signature) +
  scale_fill_manual(values = clrs$consensus_signature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) +
  coord_cartesian(ylim = c(-1, 2)) +
  labs(y = "Module score", color = "Mutational\nsignature",
       fill = "Mutational\nsignature")

mutsig_boxplot_grid_dc <- plot_grid(
  mutsig_boxplot_dc + NoLegend(),
  mutsig_boxplot_ifnreg + NoLegend(),
  nrow = 1, rel_widths = c(0.75, 0.25))

mutsig_boxplot_grid_dc
```

<img src="620_Macrophage_umaps_and_boxplots_files/figure-html/chunk_145-1.png" width="720" />

```r
ggsave("figures/620_Macrophage_umaps_and_boxplots/005_DC_boxplots.pdf", mutsig_boxplot_grid_dc, width = 7.5, height = 3.5)
```


```r
patient_boxplot_wrapper <- function(xkey) {

  xkey <- enquo(xkey)

  patient_lvls <- plot_data %>%
    mutate(key = !!xkey) %>%
    group_by(patient_id_short) %>%
    summarise(median = median(value)) %>%
    ungroup() %>%
    arrange(median) %>%
    pull(patient_id_short)

  plot_data %>%
    mutate(patient_id_short = ordered(patient_id_short, levels = patient_lvls)) %>%
    ggplot() +
    geom_violin(aes(patient_id_short, !!xkey, fill = consensus_signature),
                color = "white") +
    # geom_boxplot(aes(consensus_signature, value, color = consensus_signature),
    #              width = 0.2, size = 2) +
    geom_boxplot(aes(patient_id_short, !!xkey, fill = consensus_signature),
                 color = "white", width = 0.2, outlier.shape = NA) +
    scale_color_manual(values = clrs$consensus_signature) +
    scale_fill_manual(values = clrs$consensus_signature) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank()) +
    labs(y = "Module score", color = "Mutational\nsignature", title = xkey)
}

# patient_boxplot_list <- lapply(module_names_sub, patient_boxplot_wrapper) %>%
#   setNames(module_names_sub)
```

## correlation


```r
plot_data_y <- FetchData(seu_obj_cc, c("sample", "cell_type", "JAK.STAT.pathway")) %>%
  bind_rows(FetchData(seu_obj_tc, c("sample", "cell_type", "JAK.STAT.pathway"))) %>% 
  bind_rows(FetchData(seu_obj_mp, c("sample", "cell_type", "JAK.STAT.pathway"))) %>% 
  as_tibble() %>% 
  gather(key_y, value, -sample, -cell_type) %>%
  left_join(scrna_meta_tbl, by = "sample") %>%
  filter((sort_short == "CD45+" & cell_type %in% c("T.cell", "Monocyte")) | (sort_short == "CD45-" & cell_type == "Ovarian.cancer.cell")) %>% 
  group_by(sample, cell_type, key_y, patient_id_short, tumor_supersite, consensus_signature) %>%
  summarise(mean_value_y = mean(value)) %>%
  ungroup %>%
  mutate(sample = str_remove_all(sample, "_CD45N|_CD45P"))

plot_data_cor <- mutsig_boxplot_data_stats_dc %>% 
  select(sample, key_x = key, value) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>%
  filter(sort_short == "CD45+") %>%
  group_by(sample, key_x) %>%
  summarise(mean_value_x = mean(value)) %>%
  ungroup %>%
  mutate(sample = str_remove_all(sample, "_CD45N|_CD45P")) %>%
  inner_join(plot_data_y, by = "sample")
 
cor_plot1 <- plot_data_cor %>%
  filter(key_x == "IFN\nregulators", key_y == "JAK.STAT.pathway") %>%
  mutate(cell_type = case_when(cell_type == "Ovarian.cancer.cell" ~ "Ov cancer cells",
                               cell_type == "T.cell" ~ "T cells",
                               cell_type == "Monocyte" ~ "Macrophages",
                               T ~ "NA")) %>% 
  mutate(cell_type = ordered(cell_type, levels = c("Ov cancer cells", "T cells", "Macrophages"))) %>% 
  ggplot() +
  geom_point(aes(mean_value_x, mean_value_y, color = consensus_signature)) +
  geom_smooth(aes(mean_value_x, mean_value_y), method = "lm", se = T,
              linetype = 2, size = 0.5, color = "black") +
  stat_cor(aes(mean_value_x, mean_value_y),
           method = "spearman", color = "black") +
  scale_color_manual(values = clrs$consensus_signature) +
  facet_grid(~cell_type) + 
  labs(y = "JAK/STAT pathway",
       x = "IFN regulators module score\n(Dendritic cells)",
       color = "Mutational\nsignature")
# 
# cor_plot2 <- {cor_plot1 + 
#     labs(y = "JAK-STAT score\n(CD8+ T cells)", 
#          x = "JAK-STAT score\n(cancer cells)",
#          color = "Mutational\nsignature")
# } %+%
#   filter(plot_data_cor, key_tc == "JAK.STAT.pathway", 
#          key_cc == "JAK.STAT.pathway")
# 
# cor_plot_grid <- plot_grid(cor_plot1 + guides(color = F), cor_plot2,
#                            align = "hv", axis = "x", 
#                            rel_widths = c(0.4, 0.6))
# 
# cor_plot_grid
# 

cor_plot1
```

<img src="620_Macrophage_umaps_and_boxplots_files/figure-html/chunk_170-1.png" width="768" />

```r
ggsave("figures/620_Macrophage_umaps_and_boxplots/005_DC_cancer_cell_correlation.pdf", cor_plot1, width = 8, height = 3.25)
```


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
##  package        * version    date       lib
##  abind            1.4-5      2016-07-21 [1]
##  ape              5.3        2019-03-17 [2]
##  assertthat       0.2.1      2019-03-21 [2]
##  backports        1.1.10     2020-09-15 [1]
##  beeswarm         0.2.3      2016-04-25 [2]
##  bibtex           0.4.2.2    2020-01-02 [2]
##  Biobase          2.46.0     2019-10-29 [2]
##  BiocGenerics     0.32.0     2019-10-29 [2]
##  bookdown         0.22       2021-04-22 [1]
##  broom            0.7.6      2021-04-05 [1]
##  bslib            0.2.4      2021-01-25 [1]
##  Cairo            1.5-12.2   2020-07-07 [1]
##  callr            3.6.0      2021-03-28 [1]
##  car              3.0-8      2020-05-21 [1]
##  carData          3.0-4      2020-05-22 [1]
##  cellranger       1.1.0      2016-07-27 [1]
##  circlize         0.4.10     2020-06-15 [1]
##  cli              2.5.0      2021-04-26 [1]
##  clue             0.3-57     2019-02-25 [1]
##  cluster          2.1.1      2021-02-14 [1]
##  codetools        0.2-18     2020-11-04 [1]
##  colorblindr    * 0.1.0      2020-01-13 [2]
##  colorspace     * 2.0-0      2020-11-11 [1]
##  ComplexHeatmap * 2.2.0      2019-10-29 [1]
##  cowplot        * 1.1.1      2020-12-30 [1]
##  crayon           1.4.1      2021-02-08 [1]
##  curl             4.3        2019-12-02 [2]
##  data.table       1.14.0     2021-02-21 [1]
##  DBI              1.1.0      2019-12-15 [2]
##  dbplyr           2.0.0      2020-11-03 [1]
##  desc             1.3.0      2021-03-05 [1]
##  devtools         2.2.1      2019-09-24 [2]
##  digest           0.6.25     2020-02-23 [1]
##  dplyr          * 1.0.2      2020-08-18 [1]
##  ellipsis         0.3.1      2020-05-15 [1]
##  evaluate         0.14       2019-05-28 [1]
##  fansi            0.4.1      2020-01-08 [2]
##  farver           2.0.3      2020-01-16 [1]
##  fitdistrplus     1.0-14     2019-01-23 [2]
##  forcats        * 0.5.1      2021-01-27 [1]
##  foreign          0.8-74     2019-12-26 [3]
##  fs               1.5.0      2020-07-31 [1]
##  future           1.15.1     2019-11-25 [2]
##  future.apply     1.4.0      2020-01-07 [2]
##  gbRd             0.4-11     2012-10-01 [2]
##  generics         0.1.0      2020-10-31 [1]
##  GetoptLong       1.0.2      2020-07-06 [1]
##  ggbeeswarm       0.6.0      2017-08-07 [2]
##  ggplot2        * 3.3.3      2020-12-30 [1]
##  ggpubr         * 0.4.0      2020-06-27 [1]
##  ggrastr          0.1.9      2020-06-20 [1]
##  ggrepel          0.9.1      2021-01-15 [1]
##  ggridges         0.5.2      2020-01-12 [2]
##  ggsignif         0.6.0      2019-08-08 [1]
##  GlobalOptions    0.1.2      2020-06-10 [1]
##  globals          0.12.5     2019-12-07 [2]
##  glue             1.3.2      2020-03-12 [1]
##  gridExtra        2.3        2017-09-09 [1]
##  gtable           0.3.0      2019-03-25 [1]
##  haven            2.3.1      2020-06-01 [1]
##  highr            0.8        2019-03-20 [1]
##  hms              1.0.0      2021-01-13 [1]
##  htmltools        0.5.1.1    2021-01-22 [1]
##  htmlwidgets      1.5.3      2020-12-10 [1]
##  httr             1.4.2      2020-07-20 [1]
##  ica              1.0-2      2018-05-24 [2]
##  igraph           1.2.6      2020-10-06 [1]
##  irlba            2.3.3      2019-02-05 [2]
##  jquerylib        0.1.3      2020-12-17 [1]
##  jsonlite         1.7.2      2020-12-09 [1]
##  KernSmooth       2.23-18    2020-10-29 [1]
##  knitr            1.31       2021-01-27 [1]
##  labeling         0.4.2      2020-10-20 [1]
##  lattice          0.20-41    2020-04-02 [1]
##  lazyeval         0.2.2      2019-03-15 [1]
##  leiden           0.3.1      2019-07-23 [2]
##  lifecycle        0.2.0      2020-03-06 [1]
##  listenv          0.8.0      2019-12-05 [2]
##  lmtest           0.9-37     2019-04-30 [2]
##  lsei             1.2-0      2017-10-23 [2]
##  lubridate        1.7.10     2021-02-26 [1]
##  magrittr       * 2.0.1      2020-11-17 [1]
##  MASS             7.3-53.1   2021-02-12 [1]
##  Matrix           1.3-2      2021-01-06 [1]
##  memoise          1.1.0      2017-04-21 [2]
##  metap            1.2        2019-12-08 [2]
##  mgcv             1.8-34     2021-02-16 [1]
##  mnormt           1.5-5      2016-10-15 [2]
##  modelr           0.1.8      2020-05-19 [1]
##  multcomp         1.4-12     2020-01-10 [2]
##  multtest         2.42.0     2019-10-29 [2]
##  munsell          0.5.0      2018-06-12 [1]
##  mutoss           0.1-12     2017-12-04 [2]
##  mvtnorm          1.0-12     2020-01-09 [2]
##  nlme             3.1-152    2021-02-04 [1]
##  npsurv           0.4-0      2017-10-14 [2]
##  numDeriv         2016.8-1.1 2019-06-06 [1]
##  openxlsx         4.1.5      2020-05-06 [1]
##  pbapply          1.4-2      2019-08-31 [2]
##  pillar           1.6.0      2021-04-13 [1]
##  pkgbuild         1.0.6      2019-10-09 [2]
##  pkgconfig        2.0.3      2019-09-22 [1]
##  pkgload          1.0.2      2018-10-29 [2]
##  plotly           4.9.3      2021-01-10 [1]
##  plotrix          3.7-7      2019-12-05 [2]
##  plyr             1.8.6      2020-03-03 [1]
##  png              0.1-7      2013-12-03 [1]
##  prettyunits      1.1.1      2020-01-24 [1]
##  processx         3.5.0      2021-03-23 [1]
##  ps               1.3.2      2020-02-13 [1]
##  purrr          * 0.3.4      2020-04-17 [1]
##  R.methodsS3      1.7.1      2016-02-16 [2]
##  R.oo             1.23.0     2019-11-03 [2]
##  R.utils          2.9.2      2019-12-08 [2]
##  R6               2.4.1      2019-11-12 [1]
##  RANN             2.6.1      2019-01-08 [2]
##  rappdirs         0.3.3      2021-01-31 [1]
##  RColorBrewer     1.1-2      2014-12-07 [1]
##  Rcpp             1.0.4      2020-03-17 [1]
##  RcppAnnoy        0.0.16     2020-03-08 [1]
##  RcppParallel     4.4.4      2019-09-27 [2]
##  Rdpack           0.11-1     2019-12-14 [2]
##  readr          * 1.4.0      2020-10-05 [1]
##  readxl         * 1.3.1      2019-03-13 [1]
##  remotes          2.3.0      2021-04-01 [1]
##  reprex           2.0.0      2021-04-02 [1]
##  reshape2         1.4.4      2020-04-09 [1]
##  reticulate       1.14       2019-12-17 [2]
##  rio              0.5.16     2018-11-26 [1]
##  rjson            0.2.20     2018-06-08 [1]
##  rlang            0.4.8      2020-10-08 [1]
##  rmarkdown        2.7        2021-02-19 [1]
##  ROCR             1.0-11     2020-05-02 [1]
##  rprojroot        2.0.2      2020-11-15 [1]
##  rstatix          0.6.0      2020-06-18 [1]
##  rstudioapi       0.13       2020-11-12 [1]
##  rsvd             1.0.3      2020-02-17 [1]
##  Rtsne            0.15       2018-11-10 [2]
##  rvest            0.3.6      2020-07-25 [1]
##  sandwich         2.5-1      2019-04-06 [2]
##  sass             0.3.1      2021-01-24 [1]
##  scales           1.1.1      2020-05-11 [1]
##  sctransform      0.2.1      2019-12-17 [2]
##  SDMTools         1.1-221.2  2019-11-30 [2]
##  sessioninfo      1.1.1      2018-11-05 [2]
##  Seurat         * 3.1.2      2019-12-12 [2]
##  shape            1.4.4      2018-02-07 [1]
##  sn               1.5-4      2019-05-14 [2]
##  stringi          1.5.3      2020-09-09 [1]
##  stringr        * 1.4.0      2019-02-10 [1]
##  survival         3.2-10     2021-03-16 [1]
##  testthat         2.3.2      2020-03-02 [1]
##  TFisher          0.2.0      2018-03-21 [2]
##  TH.data          1.0-10     2019-01-21 [2]
##  tibble         * 3.0.4      2020-10-12 [1]
##  tidyr          * 1.1.2      2020-08-27 [1]
##  tidyselect       1.1.0      2020-05-11 [1]
##  tidyverse      * 1.3.0      2019-11-21 [1]
##  tsne             0.1-3      2016-07-15 [2]
##  usethis          1.5.1      2019-07-04 [2]
##  utf8             1.1.4      2018-05-24 [2]
##  uwot             0.1.5      2019-12-04 [2]
##  vctrs            0.3.5      2020-11-17 [1]
##  vipor            0.4.5      2017-03-22 [2]
##  viridis        * 0.6.0      2021-04-15 [1]
##  viridisLite    * 0.4.0      2021-04-13 [1]
##  withr            2.3.0      2020-09-22 [1]
##  xfun             0.22       2021-03-11 [1]
##  xml2             1.3.2      2020-04-23 [1]
##  yaml             2.2.1      2020-02-01 [1]
##  zip              2.1.1      2020-08-27 [1]
##  zoo              1.8-7      2020-01-10 [2]
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
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  Github (clauswilke/colorblindr@1ac3d4d)
##  CRAN (R 3.6.3)                         
##  Bioconductor                           
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
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
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
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.3)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.2)                         
##  CRAN (R 3.6.3)                         
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
```


