---
title: "MSK SPECTRUM freeze major subset deep dives"
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



# Figure 6 upper


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

coi <- params$cell_type_super
cell_sort <- params$cell_sort
cell_type_major <- params$cell_type_major
louvain_resolution <- params$louvain_resolution
louvain_cluster <- params$louvain_cluster
pcut <- params$pcut

# p1 <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
# vp_outer <- viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.5, angle = 90)
# vp_inner <- viewport(x = 0.5, y = 0.5, width = 1, height = 1, angle = 0)
# grid.newpage()
# pushViewport(vp_outer)
# grid.rect()
# print(p1, vp = vp_inner)
# popViewport(1)
# g1 <- grid.grab()
# ggdraw() + draw_grob(g1)
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
myfeatures <- c("UMAP_1", "UMAP_2", "umapharmony_1", "umapharmony_2", "sample", "doublet", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score", "cell_type", "cluster_label")

seu_obj_sub <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/", coi, "_processed_filtered_annotated.rds"))

seu_obj_sub_mp <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/Macrophages_processed.rds"))

seu_obj_sub_dc <- read_rds(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/DCs_processed.rds"))

marker_tbl <- read_tsv(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/", coi, "_marker_table_annotated.tsv"))

marker_tbl_top <- marker_tbl %>% 
  filter(avg_logFC > 0.5,
         p_val_adj < 0.01,
         pct.1 > 0.2,
         pct.2 < 0.8,
         !is.na(cluster_label),
         !str_detect(gene, "^RPS|^RPL")) %>% 
  group_by(cluster_label) %>% 
  slice(1:50)

marker_sheet <- read_tsv(paste0("/work/shah/uhlitzf/data/SPECTRUM/freeze/v7/supplementary_tables/", coi, "_marker_sheet.tsv"))

my_subtypes <- names(clrs$cluster_label[[coi]])

plot_data_sub <- list(Macrophage = FetchData(seu_obj_sub_mp, c(myfeatures)),
                      DC = FetchData(seu_obj_sub_dc, c(myfeatures))) %>% 
  bind_rows(.id = "cell_type_myeloid") %>% 
  as_tibble() %>% 
  left_join(scrna_meta_tbl, by = "sample") %>% 
  mutate(cell_type_super = cell_type_super_lookup[cell_type]) %>% 
  mutate(sort_short = str_remove_all(sort_parameters, "singlet, live, ")) %>% 
  mutate(sort_short_x = ifelse(sort_short == "U" & cell_type_super == "Immune", 
                               "CD45+", ifelse(sort_short == "U" & cell_type_super == "Stromal", "CD45-", sort_short)),
         cluster_label = ordered(cluster_label, levels = names(clrs$cluster_label[[coi]])))
  
plot_data_sub <- filter(plot_data_sub, sort_short_x == cell_sort, !is.na(tumor_supersite))
```

## UMAP


```r
alpha_lvl <- ifelse(nrow(plot_data_sub) < 20000, 0.2, 0.1)
pt_size <- ifelse(nrow(plot_data_sub) < 20000, 0.2, 0.05)

common_layers_disc <- list(  
  ggrastr::geom_point_rast(size = pt_size, alpha = alpha_lvl, raster.dpi = 150),
  NoAxes(),
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))),
  labs(color = ""),
  theme(aspect.ratio = 1)
)

common_layers_cont <- list(  
  ggrastr::geom_point_rast(size = pt_size, alpha = alpha_lvl, raster.dpi = 150),
  NoAxes(),
  scale_color_gradientn(colors = viridis(9)),
  guides(color = guide_colorbar())
)

umap_coord_anno <- function(size) {
  ggplot(tibble(group = c("UMAP1", "UMAP2"),
                x = c(0, 0), xend = c(1, 0),
                y = c(0, 0), yend = c(0, 1),
                lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                angle = c(0, 90))) +
    geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                 arrow = arrow(angle = 20, type = "closed", length = unit(0.1, "npc")),
                 size = size/4, lineend = "round") +
    geom_text(aes(lx, ly, label = group, angle = angle), size = size) +
    theme_void() +
    coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))
}

add_umap_coord <- function(gg_obj, x = 0, y = 0, width = 0.3, height = 0.3, size = 4.5) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno(size = size), x = x, y = y, width = width, height = height)
  return(p)
}
```


```r
plot_data_label <- plot_data_sub %>% 
  group_by(cluster_label, cell_type_myeloid) %>% 
  summarise(umapharmony_1 = median(umapharmony_1),
            umapharmony_2 = median(umapharmony_2)) %>% 
  mutate(cluster_number = as.numeric(cluster_label))

umap_cell_type <- ggplot(plot_data_sub, aes(umapharmony_1, umapharmony_2, color = cluster_label)) + 
  common_layers_disc +
  facet_wrap(~cell_type_myeloid, scales = "free") +
  # geom_label(aes(umapharmony_1, umapharmony_2, label = cluster_label), color = "black",
  #            data = plot_data_label, label.size = unit(0, "mm"), label.r = unit(0, "mm"),
  #            alpha = 0.5, nudge_y = 1, nudge_x = -1) +
  geom_point(aes(umapharmony_1, umapharmony_2), color = "white",
             data = plot_data_label, alpha = 0.5, size = 6) +
  geom_text(aes(umapharmony_1, umapharmony_2, label = cluster_number), color = "black",
             data = plot_data_label) +
  scale_color_manual(values = clrs$cluster_label[[coi]]) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

umap_site <- ggplot(plot_data_sub, aes(umapharmony_1, umapharmony_2, color = tumor_supersite)) + 
  common_layers_disc +
  facet_wrap(~cell_type_myeloid, scales = "free") +
  scale_color_manual(values = clrs$tumor_supersite)

umap_cell_type_void <- umap_cell_type + guides(color = F)
umap_site_void <- umap_site + guides(color = F)

umap_cell_type_legend <- cowplot::get_legend(umap_cell_type)
umap_site_legend <- cowplot::get_legend(umap_site)
```

## marker heatmap 


```r
plot_data_markers <- as_tibble(FetchData(seu_obj_sub, c(myfeatures, unique(marker_tbl_top$gene)))) %>% 
  gather(gene, value, -c(1:(length(myfeatures)+1))) %>% 
  left_join(scrna_meta_tbl, by = "sample") %>% 
  mutate(cluster_label = ordered(cluster_label, levels = my_subtypes)) %>% 
  group_by(cluster_label, gene) %>% 
  summarise(value = mean(value, na.rm = T)) %>% 
  group_by(gene) %>% 
  mutate(value = scales::rescale(value)) %>% 
  left_join(select(marker_tbl_top, cluster_label_x = cluster_label, gene), by = "gene") %>% 
  mutate(cluster_label_x = ordered(cluster_label_x, levels = rev(names(clrs$cluster_label[[coi]])))) %>% 
  na.omit()
```


```r
highlight_genes <- marker_tbl_top %>% 
  group_by(cluster_label) %>% 
  slice(1:2) %>% 
  mutate(cluster_label_x = ordered(cluster_label, levels = rev(names(clrs$cluster_label[[coi]])))) %>% 
  ungroup() %>% 
  select(cluster_label_x, gene) %>% 
  na.omit %>% 
  mutate(highlight = T)

plot_data_markers_mat <- plot_data_markers %>% 
  spread(cluster_label, value) %>% 
  left_join(highlight_genes, by = c("gene", "cluster_label_x")) %>% 
  arrange(desc(cluster_label_x), gene) %>% 
  select(cluster_label_x, gene, highlight, everything())

ha_row <- rowAnnotation(
  `Cell type` = plot_data_markers_mat$cluster_label_x,
  col = list(`Cell type` = clrs$cluster_label[[coi]]),
  show_legend = F,
  annotation_name_side = "bottom"
)

ha_col <- columnAnnotation(
  `Cell type` = anno_block(gp = gpar(fill = clrs$cluster_label[[coi]], col = NA),
                           labels = as.character(1:(ncol(plot_data_markers_mat)-3)),
                           labels_rot = 90),
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
    labels_rot = 180
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
  column_names_side = "bottom",
  right_annotation = ha_row,
  left_annotation = ha_genes,
  bottom_annotation = ha_col,
  cluster_rows = F, 
  row_title = NULL,
  column_title = NULL,
  col = viridis(9)
)

# marker_heatmap
```

## composition

### per site


```r
source("src/comp_plot.R")

comp_tbl_sample_sort <- plot_data_sub %>%
  mutate(cell_type_super_ml = ifelse(cluster_label %in% grep("DC", names(clrs$cluster_label[[coi]]), value = T), "Dendritic cell", "Macrophage")) %>% 
  filter(therapy == "pre-Rx", !(consensus_signature %in% c("Undetermined", "HRD"))) %>%
  group_by(sample, tumor_subsite, tumor_supersite, tumor_megasite, patient_id_short,
           therapy, sort_short_x, consensus_signature, cluster_label, cell_type_super_ml) %>%
  tally %>%
  group_by(sample, tumor_subsite, tumor_supersite, tumor_megasite, patient_id_short,
           therapy, sort_short_x, consensus_signature) %>%
  mutate(nrel = n/sum(n)*100,
         log10n = log10(n)) %>%
  mutate(sample_id = sample) %>%
  mutate(tumor_supersite = ordered(tumor_supersite, levels = rev(names(clrs$tumor_supersite))))

comp_tbl_sample_sort_mp <- plot_data_sub %>%
  mutate(cell_type_super_ml = ifelse(cluster_label %in% grep("DC", names(clrs$cluster_label[[coi]]), value = T), "Dendritic cell", "Macrophage")) %>% 
  filter(therapy == "pre-Rx", !(consensus_signature %in% c("Undetermined", "HRD")), cell_type_super_ml == "Macrophage", str_detect(cluster_label, "M[1-2]\\.|Cycling.M|Clearing.M")) %>%
  group_by(sample, tumor_subsite, tumor_supersite, tumor_megasite, patient_id_short,
           therapy, sort_short_x, consensus_signature, cluster_label, cell_type_super_ml) %>%
  tally %>%
  group_by(sample, tumor_subsite, tumor_supersite, tumor_megasite, patient_id_short,
           therapy, sort_short_x, consensus_signature) %>%
  mutate(nrel = n/sum(n)*100,
         log10n = log10(n)) %>%
  mutate(sample_id = sample) %>%
  mutate(tumor_supersite = ordered(tumor_supersite, levels = rev(names(clrs$tumor_supersite))))

plist1 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M1.S100A8",
  cluster_label, super_type = "Myeloid.super", vec_plot = F)

# plist2 <- default_comp_grid_list(
#   comp_tbl_sample_sort_mp, cluster_label, "M1.CDKN1C",
#   cluster_label, super_type = "Myeloid.super", vec_plot = F)

plist2 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.CXCL10",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, nmax = 5000)

plist3 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.SELENOP",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, nmax = 5000)

plist4 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.MARCO",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, nmax = 5000)

plist4_no_yaxis <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.MARCO",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, yaxis = F, nmax = 5000)

plist5 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.ECM.1",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, nmax = 5000)

plist6 <- default_comp_grid_list(
  comp_tbl_sample_sort_mp, cluster_label, "M2.ECM.2",
  cluster_label, super_type = "Myeloid.super", vec_plot = F, nmax = 5000)

plist7 <- default_comp_grid_list(
  comp_tbl_sample_sort, cell_type_super_ml, "Dendritic cell",
  cluster_label, super_type = "Myeloid.super", nmax = 5000, highlight = T)

plist8 <- default_comp_grid_list(
  filter(comp_tbl_sample_sort, tumor_supersite != "Ascites"), 
  cell_type_super_ml, "Dendritic cell",
  cluster_label, super_type = "Myeloid.super", nmax = 5000, highlight = T)

comp_grid1 <- plot_grid(plotlist = plist1, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid2 <- plot_grid(plotlist = plist2, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid3 <- plot_grid(plotlist = plist3, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid4 <- plot_grid(plotlist = plist4, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid4_no_yaxis <- plot_grid(plotlist = plist4_no_yaxis, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid5 <- plot_grid(plotlist = plist5, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))

comp_grid6 <- plot_grid(plotlist = plist6, ncol = 1, align = "v", rel_heights = c(0.2, 0.2, 0.18, 0.42))


comp_grid7 <- plot_grid(plotlist = plist7, ncol = 1, align = "v", rel_heights = c(0.13, 0.13, 0.13, 0.17, 0.44))

comp_grid8 <- plot_grid(plotlist = plist8, ncol = 1, align = "v", rel_heights = c(0.13, 0.13, 0.13, 0.17, 0.44))
```


```r
comp_grid1
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-1.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_1.pdf", comp_grid1, width = 3, height = 6)

comp_grid2
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-2.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_2.pdf", comp_grid2, width = 3, height = 6)

comp_grid3
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-3.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_3.pdf", comp_grid3, width = 3, height = 6)

comp_grid4
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-4.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_4.pdf", comp_grid4, width = 3, height = 6)

comp_grid5
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-5.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_5.pdf", comp_grid5, width = 3, height = 6)

comp_grid6
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_090-6.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_6.pdf", comp_grid5, width = 3, height = 6)
```


```r
comp_grid2_4 <- plot_grid(comp_grid2, comp_grid4_no_yaxis, ncol = 2, rel_widths = c(0.64, 0.36))
comp_grid2_4
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_095-1.png" width="432" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Macrophage.cell_comp_2_4.pdf", comp_grid2_4, width = 4.5, height = 6)
```


```r
comp_grid7
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_100-1.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Myeloid.cell_comp.pdf", comp_grid7, width = 3, height = 10)

comp_grid8
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_100-2.png" width="288" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Myeloid.cell_comp_no_ascites.pdf", comp_grid8, width = 3, height = 10)
```


```r
heatmap_grob <- grid.grabExpr(draw(marker_heatmap), width = 4, height = 8)
vp_outer <- viewport(x = 0.5, y = 0.5, width = 0.5, height = 2, angle = 270)
vp_inner <- viewport(x = 0.5, y = 0.5, width = 1, height = 1, angle = 0)
grid.newpage()
pushViewport(vp_outer)
grid.draw(heatmap_grob)
popViewport(1)
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_110-1.png" width="1536" />

```r
heatmap_grob_rot <- grid.grab()
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_110-2.png" width="1536" />

```r
umap_heat_grid <- ggdraw() +
  draw_plot(add_umap_coord(umap_cell_type_void, size = 4.5, x = 0),
            x = 0, y = 0, width = 0.5, height = 1) +
  # draw_plot(add_umap_coord(umap_site_void, size = 3.5, x = 0),
  #           x = 0, y = 0, width = 0.2, height = 0.5) +
  draw_grob(heatmap_grob_rot,
            x = 0.5, y = 0, width = 0.5, height = 1)
  # draw_plot_label(c("A", "B", "C", "D"), x = c(0, 0.63, 0, 0), y = c(0.995, 0.995, 0.72, 0.45))

umap_heat_grid
```

<img src="610_Myeloid_cell_plotting_umap_comp_heatmap_files/figure-html/chunk_110-3.png" width="1536" />

```r
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Myeloid.cell_umap_heat.png", umap_heat_grid, width = 16, height = 4)
ggsave("figures/610_Myeloid_cell_plotting_umap_comp_heatmap/005_Myeloid.cell_umap_heat.pdf", umap_heat_grid, width = 16, height = 4)
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
##  fs               1.5.0      2020-07-31 [1]
##  future           1.15.1     2019-11-25 [2]
##  future.apply     1.4.0      2020-01-07 [2]
##  gbRd             0.4-11     2012-10-01 [2]
##  generics         0.1.0      2020-10-31 [1]
##  GetoptLong       1.0.2      2020-07-06 [1]
##  ggbeeswarm       0.6.0      2017-08-07 [2]
##  ggplot2        * 3.3.3      2020-12-30 [1]
##  ggrastr          0.1.9      2020-06-20 [1]
##  ggrepel          0.9.1      2021-01-15 [1]
##  ggridges         0.5.2      2020-01-12 [2]
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
##  rjson            0.2.20     2018-06-08 [1]
##  rlang            0.4.8      2020-10-08 [1]
##  rmarkdown        2.7        2021-02-19 [1]
##  ROCR             1.0-11     2020-05-02 [1]
##  rprojroot        2.0.2      2020-11-15 [1]
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
##  zoo              1.8-7      2020-01-10 [2]
##  source                                 
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
##  CRAN (R 3.6.2)                         
## 
## [1] /home/uhlitzf/R/lib
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
```


