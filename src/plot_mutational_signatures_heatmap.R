format_sig_props <- function(props, hc) {
  props <- props[, hc$labels]
  
  return(props)
}

format_labels <- function(labels, hc, clusters) {
  labels <- inner_join(labels, clusters)
  
  rownames(labels) <- labels$sample_id
  # labels <- labels[hc$labels, ]
  
  return(labels)
}

sample2patient <- function(sample_ids) {
  return(gsub("([[:upper:]]+[[:digit:]]+).*", "\\1", sample_ids))
}

annotate_labels_with_counts <- function(labels) {
  levels <- c()
  
  new_labels <- as.character(labels)
  for (label in unique(sort(labels))) {
    mask <- labels == label
    
    n <- sum(mask)
    new_label <- paste0(label, " (", n, ")")
    levels <- c(levels, new_label)
    new_labels[mask] <- new_label
  }
  
  if (is.factor(labels)) {
    new_labels <- factor(new_labels, levels)
  }
  
  return(new_labels)
}

make_signature_palette <- function(labels) {
  signatures <- c("HRD-Dup", "HRD-Del", "FBI", "TD", "Undetermined")
  
  ordered_labels <- c()
  for(h in signatures) {
    matches <- grepl(paste0("^", h), labels)
    if(any(matches)) {
      ordered_labels <- c(ordered_labels, labels[matches][1])
    }
  }
  
  pal <- clrs$consensus_signature[signatures]
  # pal <- brewer.pal(max(3, length(ordered_labels)), "Set1")
  # names(pal) <- ordered_labels
  # pal <- pal[ordered_labels]
  return(pal)
}

make_hr_status_palette <- function(labels) {
  signatures <- c("HRD", "HRP")
  
  ordered_labels <- c()
  for(h in signatures) {
    matches <- grepl(paste0("^", h), labels)
    if(any(matches)) {
      ordered_labels <- c(ordered_labels, labels[matches][1])
    }
  }
  
  pal <- clrs$hrdetect_signature[signatures]
  # pal <- brewer.pal(max(3, length(ordered_labels)), "Set1")
  # names(pal) <- ordered_labels
  # pal <- pal[ordered_labels]
  return(pal)
}

make_chord_signature_palette <- function(labels) {
  signatures <- c("BRCA1-like", "BRCA2-like")
  
  ordered_labels <- c()
  for(h in signatures) {
    matches <- grepl(paste0("^", h), labels)
    if(any(matches)) {
      ordered_labels <- c(ordered_labels, labels[matches][1])
    }
  }
  
  pal <- clrs$chord_signature[signatures]
  # pal <- brewer.pal(max(3, length(ordered_labels)), "Set1")
  # names(pal) <- ordered_labels
  # pal <- pal[ordered_labels]
  return(pal)
}

make_discrete_palette <- function(labels, pal_name) {
  # pal <- brewer.pal(length(unique(labels)), pal_name)
  pal <- carto_pal(length(unique(labels)), pal_name)
  levels <- mixedsort(unique(labels))
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

make_logical_palette <- function() {
  pal <- c("True" = "#333333", "False" = "#FFFFFF")
  return(pal)
}

make_continuous_palette <- function(pal_name, min, max) {
  max_cols = brewer.pal.info[pal_name, ]$maxcolors
  pal <- brewer.pal(max_cols, pal_name)
  
  vals <- seq(min, max, length.out = max_cols)
  return(colorRamp2(vals, rev(pal)))
}

make_column_top_annotation <- function(labels) {
  labels$Cohort <- labels$cohort
  if (!is.null(labels[["wgs_signature"]])) {
    labels$`MMCTM signature` <- annotate_labels_with_counts(labels$wgs_signature)
    labels$`MMCTM signature` <- labels$wgs_signature
  }
  if (!is.null(labels[["hrdetect_signature"]])) {
    labels$`HRDetect signature` <- annotate_labels_with_counts(labels$hrdetect_signature)
    labels$`HRDetect signature` <- labels$hrdetect_signature
  }
  if (!is.null(labels[["chord_signature"]])) {
    labels$`CHORD signature` <- annotate_labels_with_counts(labels$chord_signature)
    labels$`CHORD signature` <- labels$chord_signature
  }
  
  keep_cols <- c()
  if (length(unique(labels$Cohort)) > 1) {
    labels$SPECTRUM <- "False"
    labels[labels$cohort == "SPECTRUM", ]$SPECTRUM <- "True"
    keep_cols <- c(keep_cols, "SPECTRUM")
  }
  if (length(unique(labels$`MMCTM signature`)) > 1) {
    keep_cols <- c(keep_cols, "MMCTM signature")
  }
  if (length(unique(labels$`HRDetect signature`)) > 1) {
    keep_cols <- c(keep_cols, "HRDetect signature")
  }
  if (length(unique(labels$`CHORD signature`)) > 1) {
    keep_cols <- c(keep_cols, "CHORD signature")
  }

  labels <- labels[, keep_cols, drop = FALSE]
  
  annot_colours <- list()
  show_legend <- c()
  legend_params <- list()
  for (label in colnames(labels)) {
    if (label == "MMCTM signature") {
      show_legend <- c(show_legend, TRUE)
      annot_colours$`MMCTM signature` <- make_signature_palette(
        labels$`MMCTM signature`
      )
    } else if (label == "HRDetect signature") {
      show_legend <- c(show_legend, TRUE)
      annot_colours$`HRDetect signature` <- make_hr_status_palette(
        labels$`HRDetect signature`
      )
    } else if (label == "CHORD signature") {
      show_legend <- c(show_legend, TRUE)
      annot_colours$`CHORD signature` <- make_chord_signature_palette(
        labels$`CHORD signature`
      )
    } else if (
      "True" %in% labels[, label] || "False" %in% labels[, label]
    ) {
      show_legend <- c(show_legend, FALSE)
      annot_colours[[label]] <- make_logical_palette()
    }
    legend_params[[label]] <- list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 10),
      at = names(annot_colours[[label]])
    )
  }
  print(annot_colours)
  annot <- HeatmapAnnotation(
    df = labels,
    col = annot_colours,
    show_annotation_name = TRUE,
    simple_anno_size = unit(0.3, "cm"),
    annotation_name_gp = gpar(fontsize = 10),
    show_legend = show_legend,
    annotation_legend_param = legend_params
  )
  
  return(annot)
}

make_column_bottom_annotation <- function(labels) {
  snv_max <- 4e4
  labels[labels$SNVs > snv_max, "SNVs"] <- snv_max
  print("Warning: clipping # SNVs to 4e4")
  
  annot <- HeatmapAnnotation(
    SNVs = anno_barplot(
      labels$SNVs, axis = TRUE, baseline = 0, border = FALSE,
      gp = gpar(col = "#636363", fill = "#636363")#,
      # axis_param = list(
      # at = c(0, 20000, 40000),
      # gp = gpar(fontsize = font_size - 1)
      # )
    ),
    SVs = anno_barplot(
      labels$SVs, axis = TRUE, baseline = 0, border = FALSE,
      gp = gpar(col = "#636363", fill = "#636363")#,
      # axis_param = list(
      # at = c(0, 350, 700),
      # gp = gpar(fontsize = font_size - 1)
      # )
    ),
    show_annotation_name = TRUE,
    height = unit(1.9, "cm"),
    gap = unit(2.5, "mm"),
    annotation_name_gp = gpar(fontsize = font_size)
  )
  
  return(annot)
}

plot_std_prob_heatmap_mmctm <- function(props, hc, clusters, labels, compact, column_order) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 10
  n_sv <- 9
  
  print(dim(as.matrix(props)))
  print(dim(labels))
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Standardized\nprobabilities\n(MMCTM)",
    col = make_continuous_palette("RdBu", -3, 3),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    top_annotation = make_column_top_annotation(clusters),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 10),
    row_split = rep(c("SNV", "SV"), c(n_snv, n_sv)),
    column_split = clusters["wgs_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = FALSE,
    column_names_max_height = unit(6, "cm"),
    column_order = column_order,
    # top_annotation_height = unit(20, "cm"), 
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(-3, 3),
      labels = as.character(seq(-3, 3)),
      title_gp = gpar(fontsize = 10, lineheight = 0.8),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  return(hm)
}


plot_prob_heatmap_mmctm <- function(props, hc, clusters, labels, compact) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 10
  n_sv <- 9
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Probability",
    col = make_continuous_palette("BuPu", 0.5, 0),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    top_annotation = make_column_top_annotation(labels),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = font_size),
    row_names_gp = gpar(fontsize = font_size),
    row_title_gp = gpar(fontsize = font_size),
    row_split = rep(c("SNV", "SV"), c(n_snv, n_sv)),
    column_split = labels["wgs_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = FALSE,
    column_title = "MMCTM",
    column_names_max_height = unit(6, "cm"),
    # top_annotation_height = unit(20, "cm"), 
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(0, 1),
      labels = as.character(seq(0, 1)),
      title_gp = gpar(fontsize = font_size, lineheight = 0.8),
      labels_gp = gpar(fontsize = font_size)
    )
  )
  
  return(hm)
}


plot_prob_heatmap_hrdetect <- function(props, compact) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 2
  n_id <- 1
  n_sv <- 2
  n_scores <- 1
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Probability",
    col = make_continuous_palette("RdBu", -2, 2),
    # col = make_continuous_palette("BuPu", 1, 0),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    # top_annotation = make_column_top_annotation(labels),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = font_size),
    row_names_gp = gpar(fontsize = font_size),
    row_title_gp = gpar(fontsize = font_size),
    row_split = rep(c("SNV", "ID", "SV", ""), c(n_snv, n_id, n_sv, n_scores)),
    # column_split = labels["consensus_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = FALSE,
    column_title = "HRDetect",
    column_names_max_height = unit(6, "cm"),
    # top_annotation_height = unit(20, "cm"), 
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(0, 1),
      labels = as.character(seq(0, 1)),
      title_gp = gpar(fontsize = font_size, lineheight = 0.8),
      labels_gp = gpar(fontsize = font_size)
    )
  )
  
  return(hm)
}

plot_std_prob_heatmap_hrdetect <- function(props, labels, compact, column_order) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 2
  n_id <- 1
  n_sv <- 2
  n_scores <- 1
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Standardized\nprobabilities\n(HRDetect)",
    col = make_continuous_palette("RdBu", -2, 2),
    # col = make_continuous_palette("BuPu", 1, 0),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    top_annotation = make_column_top_annotation(labels),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 10),
    row_split = rep(c("SNV", "ID", "SV", ""), c(n_snv, n_id, n_sv, n_scores)),
    # column_split = labels["consensus_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = FALSE,
    column_title = "HRDetect",
    column_names_max_height = unit(6, "cm"),
    column_order = column_order,
    # top_annotation_height = unit(20, "cm"), 
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(-2, 2),
      labels = as.character(seq(-2, 2)),
      title_gp = gpar(fontsize = 10, lineheight = 0.8),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  return(hm)
}

plot_prob_heatmap_chord <- function(props, labels, compact) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 2
  n_sv <- 1
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Probability",
    col = make_continuous_palette("BuPu", 1, 0),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    # top_annotation = make_column_top_annotation(labels),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = font_size),
    row_names_gp = gpar(fontsize = font_size),
    row_title_gp = gpar(fontsize = font_size),
    row_split = rep(c("SNV/ID/SV", ""), c(n_snv, n_sv)),
    # column_split = labels["consensus_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = TRUE,
    column_title = "CHORD",
    column_names_max_height = unit(6, "cm"),
    # top_annotation_height = unit(20, "cm"), 
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(0, 1),
      labels = as.character(seq(0, 1)),
      title_gp = gpar(fontsize = font_size, lineheight = 0.8),
      labels_gp = gpar(fontsize = font_size)
    )
  )
  
  return(hm)
}

plot_std_prob_heatmap_chord <- function(props, labels, compact) {
  # props <- format_sig_props(props, hc)
  # labels <- format_labels(labels, hc, clusters)
  
  n_snv <- 2
  n_sv <- 1
  
  hm <- Heatmap(
    as.matrix(props),
    name = "Standardized\nprobabilities\n(CHORD)",
    col = make_continuous_palette("RdBu", -1.5, 1.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_column_names = !compact,
    # column_dend_height = unit(1.5, "cm"),
    top_annotation = make_column_top_annotation(labels),
    # bottom_annotation = make_column_bottom_annotation(labels),
    column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    row_title_gp = gpar(fontsize = 10),
    row_split = rep(c("SNV/ID/SV", ""), c(n_snv, n_sv)),
    # column_split = labels["consensus_signature"],
    cluster_row_slices = FALSE, 
    cluster_column_slices = FALSE,
    column_title = "CHORD",
    column_names_max_height = unit(6, "cm"),
    # top_annotation_height = unit(20, "cm"),
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      at = seq(-1.5, 1.5),
      labels = as.character(seq(-1.5, 1.5)),
      title_gp = gpar(fontsize = 10, lineheight = 0.8),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  return(hm)
}

get_plot_size <- function(props, compact) {
  if (!compact) {
    width <- 2 + 0.15 * ncol(props)
    height <- 10.5
    bottom_pad <- 4.5
  } else {
    width <- 6.5
    height <- 4.5
    bottom_pad <- 0.25
  }
  
  return(list(width = width, height = height, bottom_pad = bottom_pad))
}