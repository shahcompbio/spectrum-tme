---
title: "MSK SPECTRUM data freeze: major cell type embeddings"
author: "Florian Uhlitz"
date: "2021-05-06"
output: 
  html_document: 
    df_print: kable
    number_sections: yes
    toc: yes
    code_folding: hide
editor_options: 
  chunk_output_type: inline
---



# Figure 3 


```r
library(magrittr)
library(tidyverse)
library(Seurat)
library(readxl)
library(cowplot)
# library(colorblindr)
library(viridis)
# library(magick, lib.loc = "/home/uhlitzf/miniconda3/lib/R/library")
library(ggpubr)
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
## load full seurat objects with expression data
# seu_obj_tc <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v5/T.cell_processed_filtered_sub.rds")
# seu_obj_cc <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v5/Ovarian.cancer.cell_processed_filtered.rds"))
seu_obj_cc <- read_rds("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Ovarian.cancer.super_processed_filtered_annotated.rds")
marker_tbl <- read_tsv("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/Ovarian.cancer.super_marker_table_annotated.tsv")

marker_tbl_top <- marker_tbl %>% 
  filter(avg_logFC > 0.5, 
         p_val_adj < 0.01,
         pct.1 > 0.2,
         pct.2 < 0.8,
         !is.na(cluster_label)) %>% 
  group_by(cluster_label) %>% 
  slice(1:50)

myfeatures <- c("UMAP_1", "UMAP_2", "umapharmony_1", "umapharmony_2", "sample", "doublet", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score", "cell_type")

my_subtypes <- names(clrs$cluster_label$Ovarian.cancer.super)
coi <- "Ovarian.cancer.super"

gois <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M")
pois <- paste0(c("JAK.STAT", "NFkB", "TNFa", "TGFb", "Hypoxia"), ".pathway")
```


```r
plot_data <- cbind(cell_id = colnames(seu_obj_cc), FetchData(seu_obj_cc, c("umapharmony_1", "umapharmony_2", "umappca_1", "umappca_2", "RNA_snn_res.0.2", "sample", "cluster_label", gois, grep("pathway", colnames(seu_obj_cc@meta.data), value = T)))) %>% 
  as_tibble() %>% 
  left_join(scrna_meta_tbl, by = "sample") %>% 
  filter(!is.na(consensus_signature)) %>% 
  mutate(cluster_label = ordered(cluster_label, levels = names(clrs$cluster_label$Ovarian.cancer.super)),
         cluster_number = as.numeric(cluster_label))

base_umap <- ggplot(plot_data) +
  coord_fixed() +
  NoAxes() +
  theme(legend.position = c(0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1, 1, 1, 1),
        legend.text = element_text(size = 14, margin = margin(0, 10, 0, 0)),
        legend.spacing.x = unit(0, "npc"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 22))

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

add_umap_coord <- function(gg_obj, x = -0.015, y = -0.015, width = 0.2, height = 0.2) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno, x = x, y = y, width = width, height = height)
  return(p)
}
```


```r
pt.size <- 0.1
pt.size2 <- 0.2
pt.size.mini <- 0.01
pt.alpha <- 0.05
pt.alpha.mini <- 0.02
dpi <- 150

median_tbl <- plot_data %>%
  group_by(patient_id_short) %>% 
  summarise(umappca_1 = median(umappca_1), 
            umappca_2 = median(umappca_2))

median_tbl_cluster <- plot_data %>%
  group_by(cluster_label, cluster_number) %>% 
  summarise(umapharmony_1 = median(umapharmony_1), 
            umapharmony_2 = median(umapharmony_2))

umap_pca_mutsig <- base_umap + 
  ggrastr::geom_point_rast(aes(umappca_1, umappca_2, color = consensus_signature), 
                           size = pt.size, alpha = pt.alpha, raster.dpi = dpi) +
  # geom_text(aes(umappca_1, umappca_2, label = patient_id_short), data = median_tbl) + 
  geom_label(aes(umappca_1, umappca_2, label = patient_id_short), color = "black",
             data = median_tbl, label.size = unit(0, "mm"), label.r = unit(0, "mm"),
             alpha = 0.5) +
  scale_color_manual(values = clrs$consensus_signature) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1, 
                              label.position = "right")) +
  labs(title = "Uncorrected")

umap_mutsig <- base_umap + 
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = consensus_signature), 
                           size = pt.size2, alpha = pt.alpha, raster.dpi = dpi) +
  scale_color_manual(values = clrs$consensus_signature) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1, 
                              label.position = "right")) +
  guides(color = F) + 
  labs(title = "Corrected (harmony)")

umap_pca_cluster <- base_umap + 
  ggrastr::geom_point_rast(aes(umappca_1, umappca_2, color = cluster_label), 
                           size = pt.size, alpha = pt.alpha, raster.dpi = dpi) +
  geom_text(aes(umappca_1, umappca_2, label = patient_id_short), data = median_tbl) + 
  scale_color_manual(values = clrs$cluster_label$Ovarian.cancer.super) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1, 
                              label.position = "right")) +
  labs(title = "Cluster")

umap_cluster <- base_umap + 
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = cluster_label), 
                           size = pt.size2, alpha = pt.alpha, raster.dpi = dpi) +
  scale_color_manual(values = clrs$cluster_label$Ovarian.cancer.super) +
  geom_point(aes(umapharmony_1, umapharmony_2), 
            color = "white", alpha = 0.5, size = 6, 
            data = median_tbl_cluster) +
  geom_text(aes(umapharmony_1, umapharmony_2, label = cluster_number), 
            color = "black",
            data = median_tbl_cluster) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                              ncol = 1, 
                              label.position = "right")) +
  guides(color = F) + 
  labs(title = "Cluster")

cluster_legend <- cowplot::get_legend(umap_pca_cluster)

umap_pca_jak_stat <- base_umap + 
  # geom_point(aes(umappca_1, umappca_2), color = "grey80",
  #            size = pt.size, alpha = pt.alpha, 
  #            data = filter(plot_data, JAK.STAT.pathway <= 0)) +
  ggrastr::geom_point_rast(aes(umappca_1, umappca_2, color = JAK.STAT.pathway), 
                           size = pt.size, alpha = pt.alpha, raster.dpi = dpi, 
             data = mutate(plot_data, JAK.STAT.pathway = ifelse(JAK.STAT.pathway > 4, 4, JAK.STAT.pathway))) +
  geom_text(aes(umappca_1, umappca_2, label = patient_id_short), data = median_tbl) + 
  scale_color_gradientn(colours = viridis(9), breaks = c(0, 2, 4), limits = c(min(plot_data$JAK.STAT.pathway), 4), labels = c(0, 2, "≥4")) +
  labs(title = "JAK-STAT signaling")

umap_jak_stat <- base_umap + 
  # geom_point(aes(umapharmony_1, umapharmony_2), color = "grey80",
  #            size = pt.size, alpha = pt.alpha, 
  #            data = filter(plot_data, JAK.STAT.pathway <= 0)) +
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = JAK.STAT.pathway), 
             size = pt.size2, alpha = pt.alpha, raster.dpi = dpi,
             data = mutate(plot_data, JAK.STAT.pathway = ifelse(JAK.STAT.pathway > 4, 4, JAK.STAT.pathway))) +
  scale_color_gradientn(colours = viridis(9), breaks = c(0, 2, 4), limits = c(min(plot_data$JAK.STAT.pathway), 4), labels = c(0, 2, "≥4")) +
  labs(title = "JAK-STAT signaling")

umap_genes <- plot_data %>% 
  #sample_n(10000) %>% 
  select(umapharmony_1, umapharmony_2, gois) %>% 
  gather(key, value, -umapharmony_1, -umapharmony_2) %>% 
  mutate(value = ifelse(value > 3, 3, value)) %>% 
  filter(key %in% gois) %>% 
  mutate(key = ordered(key, levels = gois)) %>% 
  ggplot() + 
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = value), 
                           size = pt.size2, alpha = pt.alpha, raster.dpi = 70) +
  # scale_color_gradientn(colours = viridis(9)) +
  geom_point(aes(umapharmony_1, umapharmony_2), 
             color = "white", alpha = 0.25, size = 6, 
             data = median_tbl_cluster) +
  geom_text(aes(umapharmony_1, umapharmony_2, label = cluster_number), 
            color = "black",
            data = median_tbl_cluster) +
  scale_color_gradientn(colours = viridis(9), labels = c(0, 1, 2, "≥3"),
                        breaks = c(0, 1, 2, 3), limits = c(0, 3)) +
  remove_xaxis +
  remove_yaxis + 
  theme(aspect.ratio = 1) +
  facet_wrap(~key) +
  labs(color = "Expression")

umap_pathways <- plot_data %>% 
  #sample_n(10000) %>% 
  select(umapharmony_1, umapharmony_2, contains("pathway")) %>% 
  gather(key, value, -umapharmony_1, -umapharmony_2) %>% 
  mutate(value = ifelse(value > 4, 4, value),
         value = ifelse(value < -1, -1, value)) %>% 
  filter(key %in% pois) %>% 
  mutate(key = ordered(str_remove_all(key, ".pathway"), 
                       levels = str_remove_all(pois, ".pathway"))) %>% 
  ggplot() + 
  ggrastr::geom_point_rast(aes(umapharmony_1, umapharmony_2, color = value), 
                           size = pt.size2, alpha = pt.alpha, raster.dpi = 70) +
  # scale_color_gradientn(colours = viridis(9)) +
  geom_point(aes(umapharmony_1, umapharmony_2), 
             color = "white", alpha = 0.25, size = 6, 
             data = median_tbl_cluster) +
  geom_text(aes(umapharmony_1, umapharmony_2, label = cluster_number), 
            color = "black",
            data = median_tbl_cluster) +
  scale_color_gradientn(colours = viridis(9), labels = c("≤-1", 0, 2, "≥4"),
                        breaks = c(-1, 0, 2, 4), limits = c(-1, 4)) +
  remove_xaxis +
  remove_yaxis + 
  theme(aspect.ratio = 1) +
  facet_wrap(~key) +
  labs(color = "PROGENy\nscore")
```

## cluster compositions

### Cancer.cell.1, 2, 3


```r
source("src/comp_plot.R")

cluster_comp <- plot_data %>%
  filter(consensus_signature != "Undetermined") %>% 
  filter(consensus_signature != "HRD") %>% 
  mutate(sort_short_x = str_replace_all(sort_short, "U", "CD45-")) %>% 
  mutate(sample_id = sample) %>% 
  group_by(cluster_label, consensus_signature, sample_id, sort_short_x, tumor_supersite, therapy) %>%
  tally %>%
  group_by(consensus_signature, sample_id, sort_short_x, tumor_supersite, therapy) %>%
  mutate(nrel = n/sum(n)*100) %>%
  ungroup %>%
  mutate(cell_type = "Ovarian.cancer.cell")

plist1 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.1", cluster_label,
                                vec_plot = F, site_box = T, 
                                super_type = "Ovarian.cancer.super")
plist2 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.2", cluster_label,
                                vec_plot = F, site_box = T, yaxis = F,
                                super_type = "Ovarian.cancer.super")
plist3 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.3", cluster_label,
                                vec_plot = F, site_box = T, yaxis = F,
                                super_type = "Ovarian.cancer.super")
plist4 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.4", cluster_label,
                                vec_plot = F, site_box = T, yaxis = T,
                                super_type = "Ovarian.cancer.super")
plist5 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.5", cluster_label,
                                vec_plot = F, site_box = T, yaxis = F,
                                super_type = "Ovarian.cancer.super")
plist6 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                cluster_label, "Cancer.cell.6", cluster_label,
                                vec_plot = F, site_box = T, yaxis = F,
                                super_type = "Ovarian.cancer.super")

pcomp_grid_p1 <- plot_grid(plotlist = plist1,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.18, 0.42))
pcomp_grid_p2 <- plot_grid(plotlist = plist2,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.18, 0.42))
pcomp_grid_p3 <- plot_grid(plotlist = plist3,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.18, 0.42))
pcomp_grid_p4 <- plot_grid(plotlist = plist4,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.18, 0.42))
pcomp_grid_p5 <- plot_grid(plotlist = plist5,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.18, 0.42))
pcomp_grid_p6 <- plot_grid(plotlist = plist6,
                           ncol = 1, align = "v",
                           rel_heights = c(0.2, 0.2, 0.2, 0.38))

pcomp_grid_full_1 <- plot_grid(pcomp_grid_p1, pcomp_grid_p2, ggdraw(),
                             pcomp_grid_p3, ggdraw(), 
                             nrow = 1, 
                             rel_widths = c(0.48, 0.25, 0.01, 0.25, 0.01))

pcomp_grid_full_2 <- plot_grid(pcomp_grid_p4, pcomp_grid_p5, ggdraw(),
                             pcomp_grid_p6, ggdraw(), 
                             nrow = 1, 
                             rel_widths = c(0.48, 0.25, 0.01, 0.25, 0.01))
```


```r
pcomp_grid_full_1
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_064-1.png" width="576" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_comp_plot_1.pdf", pcomp_grid_full_1, width = 6, height = 5.5)
```

### Cancer.cell.4, 5, 6


```r
pcomp_grid_full_2
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_065-1.png" width="576" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_comp_plot_2.pdf", pcomp_grid_full_2, width = 6, height = 5.5)
```

### Cancer.cell.3


```r
plist3_2 <- default_comp_grid_list(filter(cluster_comp, sort_short_x == "CD45-"),
                                   cluster_label, "Cancer.cell.3", cluster_label,
                                   vec_plot = F, site_box = T, yaxis = T,
                                   super_type = "Ovarian.cancer.super")

pcomp_grid_p3_2 <- plot_grid(plotlist = plist3_2,
                             ncol = 1, align = "v",
                             rel_heights = c(0.2, 0.2, 0.18, 0.42))

pcomp_grid_p3_2
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_066-1.png" width="288" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_comp_plot_cluster3.pdf", pcomp_grid_p3_2, width = 3, height = 5.5)
ggsave("figures/320_cancer_cell_plotting_umap/003_comp_plot_cluster3.png", pcomp_grid_p3_2, width = 3, height = 5.5)
```

## composition test statistics heatmap


```r
test_tbl_mutsig <- wilcoxon_tests(filter(cluster_comp, sort_short_x == "CD45-"), cluster_label, names(clrs$cluster_label$Ovarian.cancer.super), consensus_signature) %>% 
  rename(y_var = consensus_signature) %>% 
  mutate(facet_var = "consensus_signature")

test_tbl_site <- wilcoxon_tests(filter(cluster_comp, sort_short_x == "CD45-"), cluster_label, names(clrs$cluster_label$Ovarian.cancer.super), tumor_supersite) %>% 
  rename(y_var = tumor_supersite) %>% 
  mutate(facet_var = "tumor_supersite")

test_tbl <- bind_rows(test_tbl_mutsig, test_tbl_site) %>% 
  mutate(cluster_label = ordered(cluster_label, levels = names(clrs$cluster_label$Ovarian.cancer.super)))

rank_range <- filter(test_tbl, qscore_sig > 0)$median_rank
rank_limits <- c(min(rank_range), 0, max(rank_range))

wilcoxon_test_heat <- ggplot(test_tbl) +
  geom_point(aes(cluster_label, y_var, color = median_rank, size = qscore_sig)) +
  scale_size_continuous(breaks = -log10(c(0.05, 0.01, 0.001)), labels = c("<0.05", "<0.01", "<0.001"), limits = c(1, max(test_tbl$qscore_sig))) +
  #scale_color_gradient2(high = "red", low = "#67A9CF", limits = c(min(test_tbl$median_rank), max(test_tbl$median_rank))) +
  facet_grid(facet_var~., scales = "free", space = "free") + 
  scale_color_gradientn(colors = c("steelblue", "white", "red"), 
                        values = scales::rescale(rank_limits),
                        limits = rank_limits[c(1, 3)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_blank()) +
  labs(x = "Cluster", y = "Group", color = "\nMedian\nscaled\nrank", size = "q-value")

wilcoxon_test_heat
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_067-1.png" width="480" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_comp_test_stat_heat.pdf", width = 5, height = 5)
```

## UMAP panels


```r
cancer_grid_pca <- ggdraw() +
  draw_plot(add_umap_coord(umap_pca_mutsig), 
            x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(add_umap_coord(umap_pca_cluster + guides(color = F)), 
            x = 0.24, y = 0, width = 0.25, height = 1) +
  draw_grob(cluster_legend, x = 0.5, y = -0.15, height = 1) +
  draw_plot(add_umap_coord(umap_pca_jak_stat), 
            x = 0.7, y = 0, width = 0.25, height = 1)

cancer_grid_harmony <- ggdraw() +
  draw_plot(add_umap_coord(umap_mutsig), 
            x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(add_umap_coord(umap_cluster + guides(color = F)), 
            x = 0.24, y = 0, width = 0.25, height = 1) +
  draw_plot(pcomp_grid_p3_2,
            x = 0.5, y = 0, width = 0.17, height = 1.02) +
  draw_plot(add_umap_coord(umap_jak_stat), 
            x = 0.7, y = 0, width = 0.25, height = 1)

cancer_grid <- ggdraw() +
  draw_plot(add_umap_coord(umap_pca_mutsig), 
            x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(add_umap_coord(umap_mutsig + guides(color = F)), 
            x = 0.24, y = 0, width = 0.25, height = 1) +
  draw_plot(add_umap_coord(umap_cluster + guides(color = F)), 
            x = 0.48, y = 0, width = 0.25, height = 1) +
  draw_plot(add_umap_coord(umap_jak_stat), 
            x = 0.73, y = 0, width = 0.25, height = 1)


cancer_grid_pca
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_070-1.png" width="1920" />

```r
cancer_grid_harmony
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_070-2.png" width="1920" />

```r
cancer_grid
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_070-3.png" width="1920" />

```r
# ggsave("figures/320_cancer_cell_plotting_umap/003_umap_grid.pdf", cancer_grid, width = 20, height = 5)
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_grid_pca.png", cancer_grid_pca, width = 20, height = 5)
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_umap_grid_pca.pdf", cancer_grid_pca, width = 20, height = 5)
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_grid_harmony.png", cancer_grid_harmony, width = 20, height = 5)
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_umap_grid_harmony.pdf", cancer_grid_harmony, width = 20, height = 5)
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_grid.png", cancer_grid, width = 20, height = 5)
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_umap_grid.pdf", cancer_grid, width = 20, height = 5)
```


```r
umap_genes_grid <- ggdraw() + 
  draw_plot(add_umap_coord(umap_genes + remove_guides), x = 0, y = 0, width = 1, height = 1) +
  draw_grob(get_legend(umap_genes), x = 0.85, y = -0.25)
umap_genes_grid
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_071-1.png" width="768" />

```r
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_genes.png", umap_genes_grid, width = 6.75, height = 4.5)
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_genes.pdf", umap_genes_grid, width = 6.75, height = 4.5)
```


```r
umap_pathways_grid <- ggdraw() + 
  draw_plot(add_umap_coord(umap_pathways + remove_guides), x = 0, y = 0, width = 1, height = 1) +
  draw_grob(get_legend(umap_pathways), x = 0.7, y = -0.25)
umap_pathways_grid
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_072-1.png" width="768" />

```r
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_pathways.png", umap_pathways_grid, width = 6.75, height = 4.5)
ggsave("figures/320_cancer_cell_plotting_umap/003_umap_pathways.pdf", umap_pathways_grid, width = 6.75, height = 4.5)
```

## Cluster marker heatmap


```r
plot_data_markers <- as_tibble(FetchData(seu_obj_cc, c("cluster_label", myfeatures, unique(marker_tbl_top$gene)))) %>% 
  gather(gene, value, -c(1:(length(myfeatures)+1))) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>% 
  mutate(cluster_label = ordered(cluster_label, levels = names(clrs$cluster_label$Ovarian.cancer.super))) %>% 
  group_by(cluster_label, gene) %>% 
  summarise(value = mean(value, na.rm = T)) %>% 
  group_by(gene) %>% 
  mutate(value = scales::rescale(value)) %>% 
  left_join(select(marker_tbl_top, cluster_label_x = cluster_label, gene), by = "gene") %>% 
  ## reverse names vector to flip row-order in heatmap
  mutate(cluster_label_x = ordered(cluster_label_x, levels = names(clrs$cluster_label$Ovarian.cancer.super))) %>% 
  na.omit()
```


```r
library(ComplexHeatmap)

highlight_genes <- marker_tbl_top %>% 
  group_by(cluster_label) %>% 
  slice(1:2) %>% 
  ## reverse levels also here for row-order flip
  mutate(cluster_label_x = ordered(cluster_label, levels = names(clrs$cluster_label[[coi]]))) %>% 
  ungroup() %>% 
  select(cluster_label_x, gene) %>% 
  na.omit %>% 
  mutate(highlight = T)

plot_data_markers_mat <- plot_data_markers %>% 
  spread(cluster_label, value) %>% 
  left_join(highlight_genes, by = c("gene", "cluster_label_x")) %>% 
  ungroup %>% 
  arrange(desc(cluster_label_x), gene) %>% 
  select(cluster_label_x, gene, highlight, everything())

ha_row <- rowAnnotation(
  `Cell type` = plot_data_markers_mat$cluster_label_x,
  col = list(`Cell type` = clrs$cluster_label[[coi]]),
  show_legend = F,
  annotation_name_side = "top"
)

ha_col <- columnAnnotation(
  `Cell type` = anno_block(gp = gpar(fill = clrs$cluster_label[[coi]], col = NA),
                           labels = as.character(1:10)), 
  show_legend = F,
  annotation_name_side = "left",
  annotation_name_rot = 0
)

gene_idx <- plot_data_markers_mat$highlight == TRUE

ha_genes <- rowAnnotation(
  link = anno_mark(
    at = which(gene_idx),
    labels = plot_data_markers_mat$gene[which(gene_idx)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"),
    side = "left",
    labels_rot = 0
  )
)

marker_heatmap <- Heatmap(
  as.matrix(plot_data_markers_mat[,-c(1:3)]), 
  heatmap_legend_param = list(
    title = "Scaled expression", 
    title_position = "leftcenter-rot"
  ),
  row_order = 1:length(plot_data_markers_mat$cluster_label_x),
  row_split = plot_data_markers_mat$cluster_label_x, 
  column_split = 1:length(colnames(plot_data_markers_mat)[-c(1:3)]),
  column_order = colnames(plot_data_markers_mat)[-c(1:3)], 
  column_names_side = "top",
  right_annotation = ha_row,
  left_annotation = ha_genes,
  top_annotation = ha_col,
  cluster_rows = F, 
  row_title = NULL,
  column_title = NULL,
  col = viridis(9)
)

# marker_heatmap
heatmap_grob <- grid.grabExpr(draw(marker_heatmap), width = 3.75, height = 5)
heat_grid <- ggdraw() + 
  draw_grob(heatmap_grob)

heat_grid
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_080-1.png" width="360" />

```r
ggsave("figures/320_cancer_cell_plotting_umap/003_cluster_marker_heatmap.png", heat_grid, width = 3.75, height = 5)
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_cluster_marker_heatmap.pdf", heat_grid, width = 3.75, height = 5)
```


## Differential pathway expression in cancer cells


```r
set.seed(42)
sampled_cell_ids <- sample(colnames(seu_obj_cc), 10000)
seu_obj_cc_sub <- subset(seu_obj_cc, cells = sampled_cell_ids)

plot_data_pw <- FetchData(seu_obj_cc, c("umapharmony_1", "umapharmony_2", "sample", "cluster_label", grep("pathway|module", colnames(seu_obj_cc@meta.data), value = T))) %>%
  as_tibble() %>%
  gather(pathway, score, -c(1:4)) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>%
  filter(sort_short == "CD45-", therapy == "pre-Rx") %>% 
  mutate(pathway = str_remove_all(pathway, "\\.pathway")) %>% 
  # filter(str_detect(cluster_label, "Cancer|cancer")) %>% 
  mutate(cluster_label = str_remove_all(cluster_label, "\\.cell")) %>% 
  filter(consensus_signature != "Undetermined") %>% 
  filter(consensus_signature != "HRD")

cut_value <- 2
pathway_summary_wrapper <- . %>% 
  summarise(mean_score = mean(score),
            median_score = median(score)) %>% 
  mutate(median_cut = ifelse(median_score > cut_value, cut_value, 
                             ifelse(median_score < -cut_value, -cut_value, 
                                    median_score))) %>% 
  mutate(mean_cut = ifelse(mean_score > cut_value, cut_value, 
                           ifelse(mean_score < -cut_value, -cut_value, 
                                  mean_score)))


# plot_data_summary_patient <- plot_data_pw %>% 
#   group_by(patient_id_short, consensus_signature, cluster_label, 
#            pathway, tumor_supersite) %>% 
#   pathway_summary_wrapper
# 
# plot_data_summary_mutsig <- plot_data_pw %>% 
#   group_by(consensus_signature, pathway) %>% 
#   pathway_summary_wrapper
# 
# plot_data_summary_mutsig_cluster <- plot_data_pw %>% 
#   group_by(consensus_signature, cluster_label, pathway) %>% 
#   pathway_summary_wrapper

plot_data_summary_cluster <- plot_data_pw %>%
  group_by(cluster_label, pathway) %>%
  pathway_summary_wrapper %>% 
  mutate(cluster_label = ordered(cluster_label, levels = str_remove_all(names(clrs$cluster_label$Ovarian.cancer.super), "\\.cell")))

plot_data_summary_mutsig_patient <- plot_data_pw %>% 
  group_by(consensus_signature, patient_id_short, pathway) %>% 
  pathway_summary_wrapper

plot_data_summary_mutsig_sample <- plot_data_pw %>% 
  group_by(consensus_signature, sample, pathway) %>% 
  pathway_summary_wrapper
```


```r
common_heat_layers <- list(
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), 
                       na.value = "grey10", 
                       breaks = c(-cut_value, 0, cut_value), 
                       labels = c(paste0("≤-", cut_value), "0", paste0("≥", cut_value)),
                       limits = c(-cut_value, cut_value)),
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0))
)

# ggplot(plot_data_summary_patient) + 
#   geom_tile(aes(cluster_label, patient_id_short, fill = median_cut)) +
#   facet_grid(consensus_signature~pathway, scales = "free", space = "free") +
#   common_heat_layers
#   
# ggplot(plot_data_summary_mutsig_cluster) + 
#   geom_tile(aes(pathway, cluster_label, fill = median_cut)) +
#   facet_grid(~consensus_signature, scales = "free", space = "free") +
#   common_heat_layers   
# 
# ggplot(plot_data_summary_mutsig) + 
#   geom_tile(aes(consensus_signature, pathway, fill = median_score)) +
#   common_heat_layers

pw_breaks <- c(seq(min(plot_data_summary_cluster$mean_score), 0, length.out = 5)[-5], 
               0, seq(0, max(plot_data_summary_cluster$mean_score), length.out = 5)[-1])

rdbu_clrs <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
rdbu_clrs[c(1:2)] <- rdbu_clrs[3]

pw_plot1 <- plot_data_summary_cluster %>% 
  mutate(facet_helper = "") %>% 
  filter(!str_detect(pathway, "module")) %>% 
  ggplot() +
  geom_tile(aes(cluster_label, pathway, fill = mean_score)) +
  common_heat_layers +
  facet_grid(~facet_helper, scales = "free", space = "free") + 
  labs(x = "Cluster", y = "Pathway", fill = "Mean\nPROGENy\nscore") +
  scale_fill_gradientn(colours = rdbu_clrs, 
                       values = scales::rescale(pw_breaks),
                       # breaks = c(-signif(min(pw_breaks),1)+0.1, 0, signif(max(pw_breaks),1)-0.1),  
                       na.value = "grey10") +
  theme(axis.text.x = element_blank())


pw_plot1_anno <- plot_data_summary_cluster %>% 
  mutate(facet_helper = "") %>% 
  filter(!str_detect(pathway, "module")) %>% 
  distinct(cluster_label, facet_helper, .keep_all = T) %>% 
  ggplot() +
  geom_tile(aes(cluster_label, facet_helper, fill = cluster_label)) +
  geom_point(aes(cluster_label, facet_helper), color = "white", alpha = 0.5, size = 5) +
  geom_text(aes(cluster_label, facet_helper, label = as.numeric(cluster_label))) +
  common_heat_layers + 
  scale_fill_manual(values = clrs$cluster_label$Ovarian.cancer.super %>% setNames(str_remove_all(names(.), "\\.cell"))) +
  facet_grid(~facet_helper, scales = "free", space = "free") +
  theme(axis.text.y = element_blank(),
        strip.text = element_blank()) +
  guides(fill = F)

# pw_plot2 <- ggplot(plot_data_summary_mutsig_patient) + 
#   geom_tile(aes(patient_id_short, pathway, fill = mean_cut)) +
#   facet_grid(~consensus_signature, scales = "free", space = "free") + 
#   common_heat_layers +
#   labs(x = "Patient", y = "", fill = "Mean\nPROGENy\nscore") +
#   theme(axis.text.x = element_blank())

# pw_plot2_anno <- ggplot(mutate(plot_data_summary_mutsig_patient, facet_helper = "")) +
#   geom_tile(aes(patient_id_short, facet_helper, fill = consensus_signature)) +
#   common_heat_layers + 
#   scale_fill_manual(values = clrs$consensus_signature) +
#   facet_grid(~consensus_signature, scales = "free", space = "free") + 
#   theme(axis.text.y = element_blank(),
#         strip.text = element_blank()) +
#   guides(fill = F)

pw_grid_left <- plot_grid(pw_plot1, pw_plot1_anno, ggdraw(),
                          ncol = 1, align = "v", axis = "lrtb",
                          rel_heights = c(0.65, 0.35, 0))

# pw_grid_right <- plot_grid(pw_plot2, pw_plot2_anno, ggdraw(),
#                            ncol = 1, align = "v", axis = "lrtb",
#                            rel_heights = c(0.7, 0.1, 0.2))

comparison_data_pw <- filter(plot_data_summary_mutsig_sample, pathway %in% str_remove_all(pois, "\\.pathway")) %>% 
  mutate(pathway = ordered(pathway, levels = str_remove_all(pois, "\\.pathway"))) %>% 
  rename(score = mean_score) %>% 
  ungroup

# comparison_data_jakstat <- filter(plot_data_summary_mutsig_sample, pathway == "JAK.STAT") %>% rename(score = mean_score) %>% ungroup %>% mutate(consensus_signature = as.numeric(consensus_signature))
# 
# test_result_jakstat <- compare_means(
#   consensus_signature ~ score, ref.group = ".all.", 
#   p.adjust.method = "fdr", method = "wilcox.test", 
#   data = comparison_data_jakstat
# ) %>% 
#   mutate(y.position = 3.5)

pw_boxplot_mutsig_pw <- filter(plot_data_pw, pathway %in% str_remove_all(pois, "\\.pathway")) %>% 
  mutate(pathway = ordered(pathway, levels = str_remove_all(pois, "\\.pathway"))) %>% 
  ggplot(aes(consensus_signature, score)) +
  geom_violin(aes(consensus_signature, score, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
  geom_boxplot(aes(consensus_signature, score, color = consensus_signature),
               width = 0.5, size = 0.75, outlier.shape = NA) +
  geom_boxplot(aes(consensus_signature, score, fill = consensus_signature),
               color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
  stat_compare_means(aes(consensus_signature, score, color = "red"), 
                     ref.group = ".all.", data = comparison_data_pw,
                     label = "p.signif", label.y = 3.5, hide.ns = T) +
  # stat_pvalue_manual(test_result_jakstat, label = "p.adj") +
  stat_compare_means(aes(consensus_signature, score),
                     data = comparison_data_pw,
                     label.y = 4.3, label.x = 1.5, label.sep = "\n\n") +
  facet_wrap(~pathway) +
  scale_color_manual(values = clrs$consensus_signature) +
  scale_fill_manual(values = clrs$consensus_signature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y = "PROGENy score") +
  coord_cartesian(ylim = c(-1, 5))

# comparison_data_module <- filter(plot_data_summary_mutsig_sample, pathway %in% c("IFNg.signaling.module", "ISG.module")) %>% 
#   mutate(pathway = c(`IFNg.signaling.module` = "IFNg", `ISG.module` = "IFNa")[pathway]) %>% 
#   rename(score = mean_score)
# 
# pw_boxplot_mutsig_module <- filter(plot_data, pathway %in% c("IFNg.signaling.module", "ISG.module")) %>% 
#   mutate(pathway = c(`IFNg.signaling.module` = "IFNg", `ISG.module` = "IFNa")[pathway]) %>%  
#   ggplot(aes(consensus_signature, score)) +
#   geom_violin(aes(consensus_signature, score, fill = consensus_signature), color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
#   geom_boxplot(aes(consensus_signature, score, color = consensus_signature),
#                width = 0.5, size = 0.75, outlier.shape = NA) +
#   geom_boxplot(aes(consensus_signature, score, fill = consensus_signature),
#                color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
#   stat_compare_means(aes(consensus_signature, score), 
#                      ref.group = ".all.", data = comparison_data_module,
#                      label = "p.signif", label.y = 0.9, hide.ns = T) +
#   stat_compare_means(aes(consensus_signature, score),
#                      data = comparison_data_module,
#                      label.y = 1.1, label.x = 1.5, label.sep = "\n\n") +
#   facet_wrap(~pathway) +
#   scale_color_manual(values = clrs$consensus_signature) +
#   scale_fill_manual(values = clrs$consensus_signature) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   labs(x = "", y = "Module score") +
#   coord_cartesian(ylim = c(-0.5, 1.2))
# 
# plot_box_grid <- plot_grid(pw_boxplot_mutsig_jakstat + remove_guides, 
#                            pw_boxplot_mutsig_module + remove_guides, 
#                            nrow = 1, align = "h", 
#                            rel_widths = c(c(0.37, 0.63)))


pw_grid_full <- ggdraw() +
  draw_plot(pw_grid_left, x = 0.01, y = 0, width = 0.41, height = 1) +
  draw_plot(pw_boxplot_mutsig_pw + remove_guides, x = 0.46, y = 0, width = 0.4, height = 1)

pw_grid_full
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_100-1.png" width="960" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_pathway_heatmap_boxplot.pdf", pw_grid_full, 
       width = 10, height = 5)
```



```r
pw_patient_boxplot_wrapper <- function(pw) {

  plot_data_pw_ranked <- filter(plot_data_pw, pathway == pw) %>%
    group_by(patient_id_short) %>% 
    mutate(median = median(score)) %>%
    ungroup %>% 
    arrange(median) %>% 
    mutate(patient_id_short = ordered(patient_id_short, levels = unique(patient_id_short)))
  
  ggplot(plot_data_pw_ranked) +
    geom_violin(aes(patient_id_short, score, fill = consensus_signature), 
                color = "white", adjust = 2, alpha = 0.5, width = 1.5) +
    geom_boxplot(aes(patient_id_short, score, color = consensus_signature),
                 width = 0.5, size = 0.75, outlier.shape = NA) +
    geom_boxplot(aes(patient_id_short, score, fill = consensus_signature),
                 color = "white", width = 0.5, outlier.shape = NA, size = 0.5) +
    facet_wrap(~pathway) +
    scale_color_manual(values = clrs$consensus_signature) +
    scale_fill_manual(values = clrs$consensus_signature) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Patient", y = "PROGENy score",
         fill = "Mutational\nsignature",
         color = "Mutational\nsignature")
  
}

pw_patient_plots_list <- lapply(str_remove_all(pois, ".pathway"), pw_patient_boxplot_wrapper) %>% 
  lapply(add, remove_guides) %>% 
  lapply(add, ylim(c(-1, 5)))

plot_grid(plotlist = pw_patient_plots_list, ncol = 1, align = "v")
```

<img src="320_cancer_cell_plotting_umap_files/figure-html/chunk_110-1.png" width="576" />

```r
ggsave_nodb("figures/320_cancer_cell_plotting_umap/003_pathway_heatmap_boxplot_patient_lvl.pdf", 
            width = 6, height = 12)
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
##  package        * version    date       lib source        
##  abind            1.4-5      2016-07-21 [1] CRAN (R 3.6.3)
##  ape              5.3        2019-03-17 [2] CRAN (R 3.6.2)
##  assertthat       0.2.1      2019-03-21 [2] CRAN (R 3.6.2)
##  backports        1.1.10     2020-09-15 [1] CRAN (R 3.6.3)
##  beeswarm         0.2.3      2016-04-25 [2] CRAN (R 3.6.2)
##  bibtex           0.4.2.2    2020-01-02 [2] CRAN (R 3.6.2)
##  Biobase          2.46.0     2019-10-29 [2] Bioconductor  
##  BiocGenerics     0.32.0     2019-10-29 [2] Bioconductor  
##  bookdown         0.22       2021-04-22 [1] CRAN (R 3.6.2)
##  broom            0.7.6      2021-04-05 [1] CRAN (R 3.6.2)
##  bslib            0.2.4      2021-01-25 [1] CRAN (R 3.6.2)
##  Cairo            1.5-12.2   2020-07-07 [1] CRAN (R 3.6.3)
##  callr            3.6.0      2021-03-28 [1] CRAN (R 3.6.3)
##  car              3.0-8      2020-05-21 [1] CRAN (R 3.6.2)
##  carData          3.0-4      2020-05-22 [1] CRAN (R 3.6.2)
##  cellranger       1.1.0      2016-07-27 [1] CRAN (R 3.6.3)
##  circlize         0.4.10     2020-06-15 [1] CRAN (R 3.6.2)
##  cli              2.5.0      2021-04-26 [1] CRAN (R 3.6.2)
##  clue             0.3-57     2019-02-25 [1] CRAN (R 3.6.2)
##  cluster          2.1.1      2021-02-14 [1] CRAN (R 3.6.3)
##  codetools        0.2-18     2020-11-04 [1] CRAN (R 3.6.3)
##  colorspace       2.0-0      2020-11-11 [1] CRAN (R 3.6.3)
##  ComplexHeatmap * 2.2.0      2019-10-29 [1] Bioconductor  
##  cowplot        * 1.1.1      2020-12-30 [1] CRAN (R 3.6.3)
##  crayon           1.4.1      2021-02-08 [1] CRAN (R 3.6.2)
##  curl             4.3        2019-12-02 [2] CRAN (R 3.6.2)
##  data.table       1.14.0     2021-02-21 [1] CRAN (R 3.6.3)
##  DBI              1.1.0      2019-12-15 [2] CRAN (R 3.6.2)
##  dbplyr           2.0.0      2020-11-03 [1] CRAN (R 3.6.2)
##  desc             1.3.0      2021-03-05 [1] CRAN (R 3.6.3)
##  devtools         2.2.1      2019-09-24 [2] CRAN (R 3.6.2)
##  digest           0.6.25     2020-02-23 [1] CRAN (R 3.6.2)
##  dplyr          * 1.0.2      2020-08-18 [1] CRAN (R 3.6.2)
##  ellipsis         0.3.1      2020-05-15 [1] CRAN (R 3.6.3)
##  evaluate         0.14       2019-05-28 [1] CRAN (R 3.6.3)
##  fansi            0.4.1      2020-01-08 [2] CRAN (R 3.6.2)
##  farver           2.0.3      2020-01-16 [1] CRAN (R 3.6.2)
##  fitdistrplus     1.0-14     2019-01-23 [2] CRAN (R 3.6.2)
##  forcats        * 0.5.1      2021-01-27 [1] CRAN (R 3.6.2)
##  foreign          0.8-74     2019-12-26 [3] CRAN (R 3.6.2)
##  fs               1.5.0      2020-07-31 [1] CRAN (R 3.6.3)
##  future           1.15.1     2019-11-25 [2] CRAN (R 3.6.2)
##  future.apply     1.4.0      2020-01-07 [2] CRAN (R 3.6.2)
##  gbRd             0.4-11     2012-10-01 [2] CRAN (R 3.6.2)
##  generics         0.1.0      2020-10-31 [1] CRAN (R 3.6.3)
##  GetoptLong       1.0.2      2020-07-06 [1] CRAN (R 3.6.2)
##  ggbeeswarm       0.6.0      2017-08-07 [2] CRAN (R 3.6.2)
##  ggplot2        * 3.3.3      2020-12-30 [1] CRAN (R 3.6.2)
##  ggpubr         * 0.4.0      2020-06-27 [1] CRAN (R 3.6.2)
##  ggrastr          0.1.9      2020-06-20 [1] CRAN (R 3.6.2)
##  ggrepel          0.9.1      2021-01-15 [1] CRAN (R 3.6.3)
##  ggridges         0.5.2      2020-01-12 [2] CRAN (R 3.6.2)
##  ggsignif         0.6.0      2019-08-08 [1] CRAN (R 3.6.2)
##  GlobalOptions    0.1.2      2020-06-10 [1] CRAN (R 3.6.2)
##  globals          0.12.5     2019-12-07 [2] CRAN (R 3.6.2)
##  glue             1.3.2      2020-03-12 [1] CRAN (R 3.6.2)
##  gridExtra        2.3        2017-09-09 [1] CRAN (R 3.6.3)
##  gtable           0.3.0      2019-03-25 [1] CRAN (R 3.6.3)
##  haven            2.3.1      2020-06-01 [1] CRAN (R 3.6.2)
##  highr            0.8        2019-03-20 [1] CRAN (R 3.6.3)
##  hms              1.0.0      2021-01-13 [1] CRAN (R 3.6.2)
##  htmltools        0.5.1.1    2021-01-22 [1] CRAN (R 3.6.2)
##  htmlwidgets      1.5.3      2020-12-10 [1] CRAN (R 3.6.3)
##  httr             1.4.2      2020-07-20 [1] CRAN (R 3.6.2)
##  ica              1.0-2      2018-05-24 [2] CRAN (R 3.6.2)
##  igraph           1.2.6      2020-10-06 [1] CRAN (R 3.6.3)
##  irlba            2.3.3      2019-02-05 [2] CRAN (R 3.6.2)
##  jquerylib        0.1.3      2020-12-17 [1] CRAN (R 3.6.2)
##  jsonlite         1.7.2      2020-12-09 [1] CRAN (R 3.6.2)
##  KernSmooth       2.23-18    2020-10-29 [1] CRAN (R 3.6.3)
##  knitr            1.31       2021-01-27 [1] CRAN (R 3.6.3)
##  labeling         0.4.2      2020-10-20 [1] CRAN (R 3.6.3)
##  lattice          0.20-41    2020-04-02 [1] CRAN (R 3.6.3)
##  lazyeval         0.2.2      2019-03-15 [1] CRAN (R 3.6.3)
##  leiden           0.3.1      2019-07-23 [2] CRAN (R 3.6.2)
##  lifecycle        0.2.0      2020-03-06 [1] CRAN (R 3.6.2)
##  listenv          0.8.0      2019-12-05 [2] CRAN (R 3.6.2)
##  lmtest           0.9-37     2019-04-30 [2] CRAN (R 3.6.2)
##  lsei             1.2-0      2017-10-23 [2] CRAN (R 3.6.2)
##  lubridate        1.7.10     2021-02-26 [1] CRAN (R 3.6.2)
##  magrittr       * 2.0.1      2020-11-17 [1] CRAN (R 3.6.2)
##  MASS             7.3-53.1   2021-02-12 [1] CRAN (R 3.6.3)
##  Matrix           1.3-2      2021-01-06 [1] CRAN (R 3.6.3)
##  memoise          1.1.0      2017-04-21 [2] CRAN (R 3.6.2)
##  metap            1.2        2019-12-08 [2] CRAN (R 3.6.2)
##  mnormt           1.5-5      2016-10-15 [2] CRAN (R 3.6.2)
##  modelr           0.1.8      2020-05-19 [1] CRAN (R 3.6.2)
##  multcomp         1.4-12     2020-01-10 [2] CRAN (R 3.6.2)
##  multtest         2.42.0     2019-10-29 [2] Bioconductor  
##  munsell          0.5.0      2018-06-12 [1] CRAN (R 3.6.3)
##  mutoss           0.1-12     2017-12-04 [2] CRAN (R 3.6.2)
##  mvtnorm          1.0-12     2020-01-09 [2] CRAN (R 3.6.2)
##  nlme             3.1-152    2021-02-04 [1] CRAN (R 3.6.3)
##  npsurv           0.4-0      2017-10-14 [2] CRAN (R 3.6.2)
##  numDeriv         2016.8-1.1 2019-06-06 [1] CRAN (R 3.6.3)
##  openxlsx         4.1.5      2020-05-06 [1] CRAN (R 3.6.2)
##  pbapply          1.4-2      2019-08-31 [2] CRAN (R 3.6.2)
##  pillar           1.6.0      2021-04-13 [1] CRAN (R 3.6.2)
##  pkgbuild         1.0.6      2019-10-09 [2] CRAN (R 3.6.2)
##  pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 3.6.3)
##  pkgload          1.0.2      2018-10-29 [2] CRAN (R 3.6.2)
##  plotly           4.9.3      2021-01-10 [1] CRAN (R 3.6.3)
##  plotrix          3.7-7      2019-12-05 [2] CRAN (R 3.6.2)
##  plyr             1.8.6      2020-03-03 [1] CRAN (R 3.6.3)
##  png              0.1-7      2013-12-03 [1] CRAN (R 3.6.3)
##  prettyunits      1.1.1      2020-01-24 [1] CRAN (R 3.6.2)
##  processx         3.5.0      2021-03-23 [1] CRAN (R 3.6.3)
##  ps               1.3.2      2020-02-13 [1] CRAN (R 3.6.2)
##  purrr          * 0.3.4      2020-04-17 [1] CRAN (R 3.6.3)
##  R.methodsS3      1.7.1      2016-02-16 [2] CRAN (R 3.6.2)
##  R.oo             1.23.0     2019-11-03 [2] CRAN (R 3.6.2)
##  R.utils          2.9.2      2019-12-08 [2] CRAN (R 3.6.2)
##  R6               2.4.1      2019-11-12 [1] CRAN (R 3.6.3)
##  RANN             2.6.1      2019-01-08 [2] CRAN (R 3.6.2)
##  rappdirs         0.3.3      2021-01-31 [1] CRAN (R 3.6.3)
##  RColorBrewer     1.1-2      2014-12-07 [1] CRAN (R 3.6.3)
##  Rcpp             1.0.4      2020-03-17 [1] CRAN (R 3.6.2)
##  RcppAnnoy        0.0.16     2020-03-08 [1] CRAN (R 3.6.2)
##  RcppParallel     4.4.4      2019-09-27 [2] CRAN (R 3.6.2)
##  Rdpack           0.11-1     2019-12-14 [2] CRAN (R 3.6.2)
##  readr          * 1.4.0      2020-10-05 [1] CRAN (R 3.6.2)
##  readxl         * 1.3.1      2019-03-13 [1] CRAN (R 3.6.2)
##  remotes          2.3.0      2021-04-01 [1] CRAN (R 3.6.3)
##  reprex           2.0.0      2021-04-02 [1] CRAN (R 3.6.2)
##  reshape2         1.4.4      2020-04-09 [1] CRAN (R 3.6.3)
##  reticulate       1.14       2019-12-17 [2] CRAN (R 3.6.2)
##  rio              0.5.16     2018-11-26 [1] CRAN (R 3.6.2)
##  rjson            0.2.20     2018-06-08 [1] CRAN (R 3.6.2)
##  rlang            0.4.8      2020-10-08 [1] CRAN (R 3.6.2)
##  rmarkdown        2.7        2021-02-19 [1] CRAN (R 3.6.2)
##  ROCR             1.0-11     2020-05-02 [1] CRAN (R 3.6.3)
##  rprojroot        2.0.2      2020-11-15 [1] CRAN (R 3.6.3)
##  rstatix          0.6.0      2020-06-18 [1] CRAN (R 3.6.2)
##  rstudioapi       0.13       2020-11-12 [1] CRAN (R 3.6.2)
##  rsvd             1.0.3      2020-02-17 [1] CRAN (R 3.6.2)
##  Rtsne            0.15       2018-11-10 [2] CRAN (R 3.6.2)
##  rvest            0.3.6      2020-07-25 [1] CRAN (R 3.6.2)
##  sandwich         2.5-1      2019-04-06 [2] CRAN (R 3.6.2)
##  sass             0.3.1      2021-01-24 [1] CRAN (R 3.6.2)
##  scales           1.1.1      2020-05-11 [1] CRAN (R 3.6.3)
##  sctransform      0.2.1      2019-12-17 [2] CRAN (R 3.6.2)
##  SDMTools         1.1-221.2  2019-11-30 [2] CRAN (R 3.6.2)
##  sessioninfo      1.1.1      2018-11-05 [2] CRAN (R 3.6.2)
##  Seurat         * 3.1.2      2019-12-12 [2] CRAN (R 3.6.2)
##  shape            1.4.4      2018-02-07 [1] CRAN (R 3.6.2)
##  sn               1.5-4      2019-05-14 [2] CRAN (R 3.6.2)
##  stringi          1.5.3      2020-09-09 [1] CRAN (R 3.6.3)
##  stringr        * 1.4.0      2019-02-10 [1] CRAN (R 3.6.3)
##  survival         3.2-10     2021-03-16 [1] CRAN (R 3.6.3)
##  testthat         2.3.2      2020-03-02 [1] CRAN (R 3.6.2)
##  TFisher          0.2.0      2018-03-21 [2] CRAN (R 3.6.2)
##  TH.data          1.0-10     2019-01-21 [2] CRAN (R 3.6.2)
##  tibble         * 3.0.4      2020-10-12 [1] CRAN (R 3.6.2)
##  tidyr          * 1.1.2      2020-08-27 [1] CRAN (R 3.6.2)
##  tidyselect       1.1.0      2020-05-11 [1] CRAN (R 3.6.2)
##  tidyverse      * 1.3.0      2019-11-21 [1] CRAN (R 3.6.2)
##  tsne             0.1-3      2016-07-15 [2] CRAN (R 3.6.2)
##  usethis          1.5.1      2019-07-04 [2] CRAN (R 3.6.2)
##  utf8             1.1.4      2018-05-24 [2] CRAN (R 3.6.2)
##  uwot             0.1.5      2019-12-04 [2] CRAN (R 3.6.2)
##  vctrs            0.3.5      2020-11-17 [1] CRAN (R 3.6.2)
##  vipor            0.4.5      2017-03-22 [2] CRAN (R 3.6.2)
##  viridis        * 0.6.0      2021-04-15 [1] CRAN (R 3.6.2)
##  viridisLite    * 0.4.0      2021-04-13 [1] CRAN (R 3.6.2)
##  withr            2.3.0      2020-09-22 [1] CRAN (R 3.6.3)
##  xfun             0.22       2021-03-11 [1] CRAN (R 3.6.3)
##  xml2             1.3.2      2020-04-23 [1] CRAN (R 3.6.2)
##  yaml             2.2.1      2020-02-01 [1] CRAN (R 3.6.2)
##  zip              2.1.1      2020-08-27 [1] CRAN (R 3.6.3)
##  zoo              1.8-7      2020-01-10 [2] CRAN (R 3.6.2)
## 
## [1] /home/uhlitzf/R/lib
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
```

