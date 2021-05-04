## plotting themes --------------------------------

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 16, font_family = "sans", ...) %+replace%
    theme(strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_blank())
}
theme_set(theme_cowplot2())

remove_xaxis <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank())

remove_yaxis <- theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank())

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F)


## umap helpers --------------------------------------

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
    draw_plot(umap_coord_anno, x = -0.015, y = -0.02, width = 0.4, height = 0.4)
  return(p)
}


## cohort marker genes ----------------------------

markers_v7 <- yaml::read_yaml("resources/annotation/hgsc_v7_major.yaml")

helper_markers <- function(x) dplyr::select(unnest(enframe(x, "subtype", "gene"), cols = gene), gene, subtype)
markers_v7_super <- lapply(yaml::read_yaml("resources/annotation/hgsc_v7_super.yaml"), helper_markers)

## load color code --------------------------------

clrs <- yaml::read_yaml("resources/annotation/colors.yaml") %>%
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

clrs$patient_id_short <- clrs$patient_id
names(clrs$patient_id_short) <- str_remove_all(names(clrs$patient_id), "SPECTRUM-OV-")

shps <- yaml::read_yaml("resources/annotation/shapes.yaml") %>% 
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

## load database ----------------------------------

db <- readr::read_rds("resources/db/tme/SPECTRUM.rds")

## define patients included in the study -----------
protocol_patients <- db$consents %>%
  dplyr::filter(patient_consent_irb == "17-182") %>%
  dplyr::pull(patient_id)

included_patients <- db$patients %>%
  dplyr::filter(patient_inclusion_exclusion=="Included") %>%
  dplyr::filter(patient_cohort_version___2=="Checked") %>%
  dplyr::filter(!patient_id %in% protocol_patients) %>%
  dplyr::pull(patient_id)

## load mutational signatures ----------------------

signature_tbl <- db$mutational_signatures %>%
  dplyr::select(patient_id, consensus_signature) %>% 
  # bind_rows(tibble(patient_id = unique(sort(scrna_meta_tbl$patient_id[!(scrna_meta_tbl$patient_id %in% .$patient_id)])), consensus_signature = "Undetermined")) %>% 
  mutate(consensus_signature = ifelse(is.na(consensus_signature), "Undetermined", consensus_signature)) %>% 
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature))) %>% 
  arrange(patient_id)

## load scRNA meta data -----------------------------

scrna_meta_tbl <- db$sequencing_scrna %>% 
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx") %>% 
  rename(sample = isabl_id) %>% 
  distinct(sample, .keep_all = T) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sort_short = str_remove_all(sort_parameters, "singlet, live, "),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>%
  left_join(db$mutational_signatures, by = "patient_id")

## load mpIF meta data -------------------------------

mpif_slide_meta_tbl <- db$mpif_slide %>%
  mutate(slide_id = str_replace_all(pici_id, " ", "_"),
         sample_id = paste0(patient_id, "_", surgery, str_replace_all(toupper(tumor_subsite), " ", "_"))) %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx") %>% 
  left_join(db$mutational_signatures, by = "patient_id") %>%
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature)))

## load H&E meta data -------------------------------

# hne_meta_tbl <- db$he_slide %>%
#   mutate(sample_id = image_hid) %>%
#   mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
#          tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
#   mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
#                                  "Other", tumor_supersite)) %>% 
#   mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
#   filter(patient_id %in% included_patients,
#          therapy == "pre-Rx") %>% 
#   left_join(db$mutational_signatures, by = "patient_id") %>%
#   mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature)))

## load WGS meta data -------------------------------

bulk_dna_meta_tbl <- db$sequencing_bulk_dna %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients) %>% 
  left_join(db$mutational_signatures, by = "patient_id") %>%
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature)))

## cell type sort fraction -------------------------

cell_type_super_lookup <- c(
  B.cell = "Immune",
  Plasma.cell = "Immune",
  T.cell = "Immune", 
  Myeloid.cell = "Immune", 
  Mast.cell = "Immune", 
  Dendritic.cell = "Immune", 
  Endothelial.cell = "Stromal",
  Fibroblast = "Stromal", 
  Ovarian.cancer.cell = "Stromal", 
  Ov.cancer.cell = "Stromal", 
  Other = "Other"
)
