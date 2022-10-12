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

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F, size = F)


## ggsave wrapper suppressing dingbats symbols 
## for adobe illustrator compatibility
ggsave_pdf <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, ...) {
  ggsave(filename = filename, plot = plot, device = cairo_pdf, path = path,
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, ...)
}

ggsave_png <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, type = "cairo", ...) {
  ggsave(filename = filename, plot = plot, device = device, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, type = type, ...)
}


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
db$mutational_signatures <- db$mutational_signatures %>%
  mutate(consensus_signature = wgs_signature)

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
  filter(patient_inclusion_exclusion=="Included",
         patient_cohort_version___2=="Checked",
         patient_id %in% hgsoc_patients,
         !patient_id %in% protocol_patients) %>%
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
  # HACK
  filter(!patient_id %in% c("SPECTRUM-OV-006","SPECTRUM-OV-044","SPECTRUM-OV-046","SPECTRUM-OV-129")) %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with H&E data
hne_patients <- db$he_slide %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         therapy == "pre-Rx",
         qc_status == "Pass") %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with bulk DNA data
bulk_dna_patients <- db$sequencing_bulk_dna %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         qc_status == "Pass") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

bulk_dna_tumor_patients <- db$sequencing_bulk_dna %>%
  mutate(sample_type = ifelse(!is.na(tumor_site), "Tumor", "Normal")) %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         qc_status == "Pass",
         sample_type == "Tumor") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

bulk_dna_normal_patients <- db$sequencing_bulk_dna %>%
  mutate(sample_type = ifelse(!is.na(tumor_site), "Tumor", "Normal")) %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         qc_status == "Pass",
         sample_type == "Normal") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with IMPACT data
impact_patients <- db$sequencing_msk_impact_custom %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scrna_patients, mpif_patients),
         qc_status == "Pass") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

## load mutational signatures ----------------------

signature_tbl <- db$mutational_signatures %>%
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature)),
         consensus_signature_short = ordered(consensus_signature_short, levels = names(clrs$consensus_signature_short)),
         wgs_signature = ordered(wgs_signature, levels = names(clrs$wgs_signature))) %>% 
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
  mutate(
    tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                            "Non-Adnexa", tumor_supersite),
    tumor_megasite = ordered(tumor_megasite, levels = names(clrs$tumor_megasite))) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Ascites" ~ "Ascites",
    tumor_megasite == "Non-Adnexa" ~ "Non-Adnexa",
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
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa"),
                                 "Non-Adnexa", tumor_supersite)) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Non-Adnexa" ~ "Non-Adnexa",
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
         is_adjacent) %>%
         # is_adjacent | is_site_matched) %>%
  group_by(patient_id, tumor_site, therapy) %>% 
  mutate(section_number = row_number()) %>%
  ungroup() %>%
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
                                 "Non-Adnexa", tumor_supersite)) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Non-Adnexa" ~ "Non-Adnexa",
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
  group_by(patient_id, tumor_site, therapy) %>% 
  mutate(section_number = row_number()) %>%
  ungroup() %>%
  left_join(signature_tbl, by = "patient_id")

mpif_fov_meta_tbl <- db$mpif_fov %>%
  dplyr::mutate(isabl_id_dummy = str_replace_all(isabl_id, "-BO", "_BO"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-IO", "_IO"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-LA", "_LA"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-RA", "_RA"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-PP", "_PP"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-LUQ", "_LUQ"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-RUQ", "_RUQ"),
                isabl_id_dummy = str_replace_all(isabl_id_dummy, "-OTHER", "_OTHER")) %>%
  tidyr::separate(isabl_id_dummy, c("patient_id","site_section","markers",NA,"fov"), sep = "_", remove = FALSE) %>%
  tidyr::separate(site_section, c("site","section"), sep = "-", remove = FALSE) %>%
  dplyr::mutate(pici_id = glue::glue("{patient_id} {site_section}")) %>%
  dplyr::select(-c('project','patient_id','site')) %>%
  dplyr::rename(image_id = isabl_id) %>%
  left_join(db$mpif_slide, by = "pici_id") %>%
  mutate(slide_id = paste0(patient_id, "_", procedure, str_replace_all(toupper(tumor_subsite), " ", "_")),
         patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sample_id_short = str_remove_all(sample_id, "OV-"),
         aliquot_id_short = str_remove_all(aliquot_id, "OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa"),
                                 "Non-Adnexa", tumor_supersite)) %>% 
  mutate(tumor_megasite_adnexa = case_when(
    tumor_megasite == "Adnexa" ~ "Adnexa",
    tumor_megasite == "Non-Adnexa" ~ "Non-Adnexa",
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
  # group_by(patient_id, tumor_site, therapy) %>% 
  # mutate(section_number = row_number()) %>%
  # ungroup() %>%
  left_join(signature_tbl, by = "patient_id")

## load WGS meta data -------------------------------

bulk_dna_meta_tbl <- db$sequencing_bulk_dna %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sample_id_short = str_remove_all(sample_id, "OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ"),
         tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite),
         tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients,
         qc_status == "Pass") %>% 
  left_join(signature_tbl, by = "patient_id")

## load MSK-IMPACT meta data ------------------------

impact_meta_tbl <- db$sequencing_msk_impact_custom %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ"),
         tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite),
         tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
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
