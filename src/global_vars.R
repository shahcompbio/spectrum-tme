## plotting themes --------------------------------

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 16, ...) %+replace%
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

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F, size = F)

## ggsave wrapper suppressing dingbats symbols 
## for adobe illustrator compatibility
ggsave_pdf <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, ...) {
  ggsave(filename = filename, plot = plot, device = cairo_pdf, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, ...)
}

ggsave_png <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, type = "cairo", ...) {
  ggsave(filename = filename, plot = plot, device = device, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, type = type, ...)
}


## umap helpers --------------------------------------

umap_coord_anno <- function(xlab, ylab, linesize, textsize) {
  arrow <- arrow(angle = 20, type = "closed", length = unit(0.1, "npc"))
  ggplot(tibble(group = c("UMAP1", "UMAP2"),
                x = c(0, 0), xend = c(1, 0),
                y = c(0, 0), yend = c(0, 1),
                lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                angle = c(0, 90))) +
    geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                 arrow = arrow, size = linesize, lineend = "round") +
    geom_text(aes(lx, ly, label = group, angle = angle), size = textsize) +
    theme_void() +
    coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))
}

add_umap_coord <- function(gg_obj, x_offset = -0.015, y_offset = -0.02, width = 0.4, height = 0.4,
                           xlab = "UMAP1", ylab = "UMAP2", linesize = 1, textsize = 4) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno(xlab = xlab, ylab = ylab, linesize = linesize, textsize = textsize), 
              x = x_offset, y = y_offset, width = width, height = height)
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

# Define confirmed HGS patients
hgsoc_patients <- db$gyn_diagnosis %>%
  filter(gyn_diagnosis_histology == "HGS") %>%
  pull(patient_id)

# Define patients on clinical trials
protocol_patients <- db$consents %>%
  filter(patient_consent_irb == "17-182") %>%
  pull(patient_id)

# Create list of patients on the SPECTRUM TME study
# - Exclude non-HGSOC patients
# - Exclude patients on clinical trials (e.g. 17-182)
included_patients <- db$patients %>%
  filter(patient_inclusion_exclusion=="Included") %>%
  filter(patient_cohort_version___2=="Checked") %>%
  filter(patient_id %in% hgsoc_patients) %>%
  filter(!patient_id %in% protocol_patients) %>%
  pull(patient_id)

# Define patients included in the study with scRNA data
scrna_patients <- db$sequencing_scrna %>%
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx",
         platform == "10x 3' GE",
         qc_status == "Pass") %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with mpIF data
mpif_patients <- db$mpif_slide %>%
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx",
         qc_status == "Pass") %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with H&E data
hne_patients <- db$he_slide %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         therapy == "pre-Rx") %>%
  pull(patient_id) %>%
  unique

## load mutational signatures ----------------------

signature_tbl <- db$mutational_signatures %>%
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature)),
         consensus_signature_short = ordered(consensus_signature_short, levels = names(clrs$consensus_signature_short))) %>% 
  arrange(patient_id)

## load scRNA meta data -----------------------------

scrna_meta_tbl <- db$sequencing_scrna %>% 
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx") %>% 
  dplyr::rename(sample = isabl_id) %>% 
  distinct(sample, .keep_all = T) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sort_short = str_remove_all(sort_parameters, "singlet, live, "),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Ascites" ~ "Other",
    tumor_megasite == "Other" ~ "Other",
    tumor_megasite == "Unknown" ~ "Unknown",
    TRUE ~ as.character(tumor_megasite)
  )) %>%
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>%
  mutate(
    tumor_supersite_adnexa = case_when(
      tumor_supersite == "Adnexa" & grepl("Adnexa", tumor_subsite) ~ "Adnexa",
      tumor_supersite == "Adnexa" & grepl("Ovary", tumor_subsite) ~ "Ovary",
      tumor_supersite == "Adnexa" & grepl("Fallopian Tube", tumor_subsite) ~ "Fallopian Tube",
      TRUE ~ as.character(tumor_supersite)
    )
  ) %>%
  left_join(signature_tbl, by = "patient_id")

## load H&E meta data -------------------------------

hne_meta_tbl <- db$he_slide %>%
  mutate(sample_id = image_hid) %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>%
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>%
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Ascites" ~ "Other",
    tumor_megasite == "Other" ~ "Other",
    tumor_megasite == "Unknown" ~ "Unknown",
    TRUE ~ as.character(tumor_megasite)
  )) %>%
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>%
  mutate(
    tumor_supersite_adnexa = case_when(
      tumor_supersite == "Adnexa" & grepl("Adnexa", tumor_subsite) ~ "Adnexa",
      tumor_supersite == "Adnexa" & grepl("Ovary", tumor_subsite) ~ "Ovary",
      tumor_supersite == "Adnexa" & grepl("Fallopian Tube", tumor_subsite) ~ "Fallopian Tube",
      TRUE ~ as.character(tumor_supersite)
    )
  ) %>%
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx") %>%
  left_join(signature_tbl, by = "patient_id")

## load mpIF meta data -------------------------------

mpif_slide_meta_tbl <- db$mpif_slide %>%
  # mutate(slide_id = str_replace_all(pici_id, " ", "_"),
  mutate(slide_id = paste0(patient_id, "_", procedure, str_replace_all(toupper(tumor_subsite), " ", "_")),
         patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sample_id_short = str_remove_all(sample_id, "OV-"),
         aliquot_id_short = str_remove_all(aliquot_id, "OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Ascites" ~ "Other",
    tumor_megasite == "Other" ~ "Other",
    tumor_megasite == "Unknown" ~ "Unknown",
    TRUE ~ as.character(tumor_megasite)
  )) %>%
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  mutate(
    tumor_supersite_adnexa = case_when(
      tumor_supersite == "Adnexa" & grepl("Adnexa", tumor_subsite) ~ "Adnexa",
      tumor_supersite == "Adnexa" & grepl("Ovary", tumor_subsite) ~ "Ovary",
      tumor_supersite == "Adnexa" & grepl("Fallopian Tube", tumor_subsite) ~ "Fallopian Tube",
      TRUE ~ as.character(tumor_supersite)
    )
  ) %>%
  filter(patient_id %in% included_patients,
         therapy == "pre-Rx",
         qc_status == "Pass") %>% 
  left_join(signature_tbl, by = "patient_id")

mpif_fov_meta_tbl <- db$mpif_fov

## load WGS meta data -------------------------------

bulk_dna_meta_tbl <- db$sequencing_bulk_dna %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sample_id_short = str_remove_all(sample_id, "OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients) %>% 
  left_join(signature_tbl, by = "patient_id")

## load MSK-IMPACT meta data ------------------------

impact_meta_tbl <- db$sequencing_msk_impact_custom %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients) %>% 
  left_join(signature_tbl, by = "patient_id")

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
