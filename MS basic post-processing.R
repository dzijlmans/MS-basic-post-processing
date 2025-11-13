# ------------------------------------------------------------
#  General Proteomics Analysis Workflow
#
#  Description:
#  This script performs a complete DEP-based proteomics analysis
#  compatible with outputs from multiple software pipelines:
#    - DIA-NN
#    - MaxQuant
#    - Proteome Discoverer
#
#  The workflow includes:
#    - Import and pre-processing (contaminant and peptide filtering)
#    - Filtering and normalization using the DEP package
#    - Differential enrichment testing
#    - QC and Volcano plotting
#    - Export of results to browsable Excel file
#
#  The script supports experiments with:
#    - 2 conditions  (simple comparison)
#    - 3 conditions  (e.g., additional grouping variable)
#  
#  Outputs:
#    - Output_plots.pdf (QC and volcano plots)
#    - DEP_results.xlsx (results and optionally full data)
#    - .RData workspace (for reproducibility)
# ------------------------------------------------------------


# ------------------------------------------------------------
# Loading libraries and setting parameters
# ------------------------------------------------------------

if (!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, DEP, openxlsx, ggrepel, tcltk, PTXQC, SummarizedExperiment)

setwd(tk_choose.dir())

parameters <- read.xlsx("Analysis_parameters.xlsx")

message(paste0("Processing data from ", parameters$analysis_method))


# ------------------------------------------------------------
#  DIA-NN Pre-processing
# ------------------------------------------------------------

if (parameters$analysis_method == "DiaNN") {
  
  tkmessageBox(message = paste(
    "Make sure that the '*pg_matrix.tsv' and '*pr_matrix.tsv' files are located in your working directory.\n\n",
    "Please ensure that the abundance column names (and thus the .raw files)",
    "contain the name of the mass spectrometer used (e.g. Astral, Exploris)",
    "and the following information, in order and separated by underscores:\n\n",
    "- Condition A (e.g. Cell line)\n",
    "- Condition B (e.g. WT and KO or treated and untreated)\n",
    "- Replicate\n\n",
    "Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw"
  ))
  
  # Import files
  proteinGroups <- read.delim(list.files(pattern = "pg_matrix.tsv"))
  peptides <- read.delim(list.files(pattern = "pr_matrix.tsv"))
  
  # Filter on number of unique peptides
  filtered_proteinGroups <- data.frame(
    Protein.Group = names(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group)),
    Peptide.Count = as.integer(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group))
  ) %>% 
    filter(Peptide.Count > (parameters$unique_peptides - 1))
  
  proteinGroups <- proteinGroups[proteinGroups$Protein.Group %in% filtered_proteinGroups$Protein.Group, ]
  
  unique_peptides <- filtered_proteinGroups %>% dplyr::rename(Protein.IDs = "Protein.Group", peptides = "Peptide.Count")
  
  # Remove contaminants and keratins
  proteinGroups <-
    make_unique(proteinGroups, "Genes", "Protein.Group", delim = ";")
  
  proteinGroups <- proteinGroups %>% filter(!grepl("Cont_", ID) &
                                              !grepl("KRT", name) &
                                              !grepl("cRAP", name) &
                                              !grepl("cRAP", Protein.Group))
  
  # Tidy column names
  proteinGroups <- proteinGroups %>% dplyr::rename(Protein.IDs = "Protein.Group")
  
  value_columns <- which(str_detect(names(proteinGroups), parameters$mass_spec))
  prefix <- LCSn(colnames(proteinGroups[, value_columns]))
  colnames(proteinGroups) <- gsub(".raw", "", gsub(".*" %>% paste0(prefix), "", colnames(proteinGroups)))
  
}


# ------------------------------------------------------------
#  MaxQuant Pre-processing
# ------------------------------------------------------------

if (parameters$analysis_method == "MaxQuant") {
  
  tkmessageBox(message = paste(
    "Make sure that the 'proteinGroups.txt' file is located in your working directory.\n\n",
    "Please ensure that the abundance column names (and thus the .raw files)",
    "contain the following information, in order and separated by underscores:\n\n",
    "- Condition A (e.g. Cell line)\n",
    "- Condition B (e.g. WT and KO or treated and untreated)\n",
    "- Replicate\n\n",
    "Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw"
  ))
  
  # Import files
  proteinGroups <- read.delim(list.files(pattern = "proteinGroups.txt"))
  
  # Filter on number of unique peptides and remove contaminants
  proteinGroups$Only.identified.by.site[is.na(proteinGroups$Only.identified.by.site)] <- ""
  proteinGroups$Reverse[is.na(proteinGroups$Reverse)] <- ""
  proteinGroups$Potential.contaminant[is.na(proteinGroups$Potential.contaminant)] <- ""
  
  proteinGroups <-
    filter(
      proteinGroups,
      Reverse != "+",
      Potential.contaminant != "+",
      Only.identified.by.site != "+",
      Unique.peptides > (parameters$unique_peptides - 1)
    )
  
  unique_peptides <- proteinGroups %>% select(Protein.IDs, Unique.peptides)

  # Tidy column names
  proteinGroups <- proteinGroups %>% dplyr::select(
    Protein.IDs = "Protein.IDs",
    "Gene.names",
    contains("LFQ.intensity."),
  )
  
  proteinGroups <-
    make_unique(proteinGroups, "Gene.names", "Protein.IDs", delim = ";")
  
  value_columns <- which(str_detect(names(proteinGroups), "LFQ.intensity"))
  prefix <- LCSn(colnames(proteinGroups[, value_columns]))
  colnames(proteinGroups) <- gsub(".raw", "", gsub(".*" %>% paste0(prefix), "", colnames(proteinGroups)))
  
}


# ------------------------------------------------------------
#  Proteome Discoverer Pre-processing
# ------------------------------------------------------------

if (parameters$analysis_method == "Proteome Discoverer") {
  
  tkmessageBox(message = paste(
    "Make sure that the '*Proteins.txt' file is located in your working directory.\n\n",
    "Please ensure that the abundance column names (and thus the .raw files)",
    "contain the following information, in order and separated by underscores:\n\n",
    "- Condition A (e.g. Cell line)\n",
    "- Condition B (e.g. WT and KO or treated and untreated)\n",
    "- Replicate\n\n",
    "Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw"
  ))
  
  # Import files
  proteinGroups <- read.delim(list.files(pattern = "Proteins.txt"))
  
  # Filter on number of unique peptides
  proteinGroups <- proteinGroups %>% filter(Number.of.Unique.Peptides > (parameters$unique_peptides - 1))
  
  unique_peptides <- proteinGroups %>% select(Protein.IDS = Accession, Unique.peptides = Number.of.Unique.Peptides)
  
  # Tidy column names
  proteinGroups <- proteinGroups %>% dplyr::select(
    "Accession",
    "Gene.Symbol",
    contains("Abundances") & contains("Normalized")
  )
  
  proteinGroups <-
    make_unique(proteinGroups, "Gene.Symbol", "Accession", delim = ";")
  
  value_columns <- which(str_detect(names(proteinGroups), "Abundances..Normalized"))
  colnames(proteinGroups) <- gsub("\\.\\.", "_", gsub(".*(Sample|Control)\\.\\.(.*)", "\\2", colnames(proteinGroups)))
  
}


# ------------------------------------------------------------
#  Quality Control (QC) Checks
# ------------------------------------------------------------

try(if (max(str_count(colnames(proteinGroups[, value_columns]), "_")) == 1 &
        parameters$conditions != 2 |
        max(str_count(colnames(proteinGroups[, value_columns]), "_")) == 2 &
        parameters$conditions != 3)
  stop("Incorrect number of conditions"))


# ------------------------------------------------------------
#  DEP Analysis (2 conditions)
# ------------------------------------------------------------

if (parameters$conditions == 2) {
  
  results_A <- list()
  data <- proteinGroups
  
  # Create experimental design
  experimental_design <- data.frame(label = colnames(data[, value_columns]))
  experimental_design <- experimental_design %>% 
    extract(label, into = c("condition", "replicate"), regex = "^(.*)_(.*)$", remove = FALSE)
  experimental_design <- experimental_design[, c("label", "condition", "replicate")]
  
  # DEP filtering and normalization
  data_se <- make_se(data, value_columns, experimental_design)
  data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
  data_norm <- normalize_vsn(data_filt)
  if (parameters$filtering_type == "complete") {
    data_imp <- data_norm
  }
  if (parameters$filtering_type == "condition") {
    data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
  }
  
  # Differential enrichment analysis
  if (parameters$comparison == "control") {
    data_diff_manual <- test_diff(data_imp, type = "control", control = parameters$control)
  } else if (parameters$comparison == "manual") {
    manual_contrasts <- str_split(parameters$contrasts, ";")[[1]]
    data_diff_manual <- test_diff(data_imp, type = "manual", test = manual_contrasts)
  } else if (parameters$comparison == "all") {
    data_diff_manual <- test_diff(data_imp, type = "all")
  }
  
  contrast <- sub("_diff", "", str_subset(
    names(data_diff_manual@elementMetadata@listData),
    "_diff"
  ))
  dep <- add_rejections(data_diff_manual, alpha = parameters$p_value, lfc = parameters$log2_FC)
  
  # Generate QC and volcano plots
  pdf("Output_plots.pdf")
  
  print(plot_frequency(data_se))
  print(plot_numbers(data_filt))
  if (parameters$filtering_type == "condition") {
    print(plot_detect(data_filt))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_normalization(data_filt, data_norm, data_imp))
  } else {
    print(plot_normalization(data_filt, data_norm))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_imputation(data_norm, data_imp))
  }
  
  # Generate Volcano plot for each contrast
  for (p in 1:length(contrast)) {
    
    results <- get_results(dep) %>% select(name, ID, contains(contrast[p]) &
                                             (contains("p.val") | contains("ratio")))
    colnames(results) <- sub("ratio", "log2.FC", colnames(results))
    results <- results %>% mutate(!!paste0(contrast[p], "_p.adj") := p.adjust(get(paste0(contrast[p], "_p.val")), method = "BH"))
    results_A[[contrast[p]]] <- results
    
    left <- str_split(contrast[p], pattern = "_vs_")[[1]][1]
    right <- str_split(contrast[p], pattern = "_vs_")[[1]][2]
    p_val_column <- grep("_p.val$", names(results), value = TRUE)
    log2_fc_column <- grep("_log2.FC$", names(results), value = TRUE)
    
    # Define proteins to highlight
    if (parameters$volcano == "protein list") {
      proteins_to_highlight <- str_split(parameters$proteins_to_highlight, ";")[[1]]
      regex_pattern <- paste0("^", proteins_to_highlight, "(\\.|$)", collapse = "|")
      matching_names <- results$name[str_detect(results$name, regex_pattern)]
      proteins_to_highlight <- unique(c(proteins_to_highlight, matching_names))
    }
    
    if (parameters$volcano == "specify significance") {
      proteins_to_highlight <- results$name[results[[p_val_column]] < parameters$p_value & abs(results[[log2_fc_column]]) > parameters$log2_FC]
    }
    
    if (parameters$volcano == "TopN") {
      sorted_results <- results[order(abs(results[[log2_fc_column]]), decreasing = TRUE), ]
      proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < 0.05, ], parameters$TopN)$name
    }
    
    # Volcano plot
    print(
      ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) +
        geom_point(color = "grey60", alpha = 0.5) +
        geom_point(data = results[results$name %in% proteins_to_highlight, ],
                   color = "#E41A1C",
                   alpha = 0.5) +
        geom_text_repel(data = results[results$name %in% proteins_to_highlight, ],
                        aes(label = name),
                        max.overlaps = Inf,
                        size = 3) +
        geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
        theme_bw(base_size = 12) +
        xlab(paste0("log2 FC (", left, " / ", right, ")")) +
        ylab("-log10(p-value)") +
        ggtitle(paste0(left, " vs ", right))
    )
    
    # Optional: highlighting of imputed proteins
    if (parameters$highlight_imputed == TRUE) {
      proteins_imputed <- get_df_long(data_norm) %>%
        group_by(name, condition) %>%
        summarize(NAs = any(is.na(intensity))) %>%
        filter(NAs) %>%
        pull(name) %>%
        unique()
      results$imputed <- ifelse(results$name %in% proteins_imputed, "Yes", "No")
      
      print(
        ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]), color = imputed, shape = imputed, fill = imputed)) +
          geom_point(alpha = 0.5) +
          scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
          scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
          scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
          geom_text_repel(data = results[results$name %in% proteins_to_highlight, ],
                          aes(label = name),
                          color = "black",
                          max.overlaps = Inf,
                          size = 3) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
          geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
          theme_bw(base_size = 12) +
          xlab(paste0("log2 FC (", contrast, ")")) +
          ylab("-log10(p-value)") +
          ggtitle(paste0(left, " vs ", right)) +
          theme(
            legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
            legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
            legend.key = element_blank()
          ) +
          guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                 color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                 shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050")))) 
      )
    }
  }
  
  dev.off()
  
  # Optional: export complete dataset
  if (parameters$complete_output == TRUE) {
    
    datalist <- list()
    
    # Raw
    df <- as.data.frame(assay(data_filt))
    names(df) <- paste0(names(df), "_raw")
    df$name <- row.names(df)
    datalist[["raw"]] <- df
    
    # Normalized
    df <- as.data.frame(assay(data_norm))
    names(df) <- paste0(names(df), "_normalized")
    df$name <- row.names(df)
    datalist[["normalized"]] <- df  
    
    # Imputed
    if (parameters$filtering_type == "condition") {
      df <- as.data.frame(assay(data_imp))
      names(df) <- paste0(names(df), "_imputed")
      df$name <- row.names(df)
      datalist[["imputed"]] <- df
    }
    
    results_A[["complete_data"]] <- datalist %>%
      purrr::reduce(merge, by = "name") %>%          
      merge(get_results(dep)) %>%                    
      select(-contains("significant"),               
             -contains("centered"),
             -contains("p.adj")) %>%
      rename_with(~ sub("ratio", "log2.FC", .x)) %>%
      left_join(proteinGroups[, c("ID", "Protein.IDs")], by = "ID") %>%
      left_join(unique_peptides, by = "Protein.IDs") %>%
      select(-Protein.IDs) 
  }
  
  write.xlsx(results_A, paste0("DEP_results.xlsx"))
  rm(results_A)
}


# ------------------------------------------------------------
#  DEP Analysis (3 conditions)
# ------------------------------------------------------------

if (parameters$conditions == 3) {
  
  results_A <- list()
  
  # Pivot and split data by variable
  proteinGroups <-
    proteinGroups %>% pivot_longer(
      cols = all_of(value_columns),
      names_to = c("condition_B", "condition_A", "replicate"),
      names_sep = "_",
      values_to = "values"
    ) %>% pivot_wider(names_from = c("condition_A", "replicate"),
                      values_from = "values")
  
  datalist_B <- split(proteinGroups, f = ~ proteinGroups$condition_B)
  
  for (i in 1:length(datalist_B)) {
    name_B <- names(datalist_B[i])
    data <- datalist_B[[i]]
    value_columns <- which(sapply(data, is.numeric))
    
    # Create experimental design
    experimental_design <- data.frame(label = colnames(data[, value_columns]))
    experimental_design <- experimental_design %>% 
      extract(label, into = c("condition", "replicate"), regex = "^(.*)_(.*)$", remove = FALSE)
    experimental_design <- experimental_design[, c("label", "condition", "replicate")]
    
    # DEP filtering and normalization
    data_se <- make_se(data, value_columns, experimental_design)
    data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
    data_norm <- normalize_vsn(data_filt)
    if (parameters$filtering_type == "complete") {
      data_imp <- data_norm
    }
    if (parameters$filtering_type == "condition") {
      data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
    }
    
    # Differential enrichment analysis
    if (parameters$comparison == "control") {
      data_diff_manual <- test_diff(data_imp, type = "control", control = parameters$control)
    } else if (parameters$comparison == "manual") {
      manual_contrasts <- str_split(parameters$contrasts, ";")[[1]]
      data_diff_manual <- test_diff(data_imp, type = "manual", test = manual_contrasts)
    } else if (parameters$comparison == "all") {
      data_diff_manual <- test_diff(data_imp, type = "all")
    }
    contrast <- sub("_diff", "", str_subset(
      names(data_diff_manual@elementMetadata@listData),
      "_diff"
    ))
    dep <- add_rejections(data_diff_manual, alpha = parameters$p_value, lfc = parameters$log2_FC)
    
    # Generate QC and volcano plots
    pdf(paste0("Output_plots_", name_B, ".pdf"))
    
    print(plot_frequency(data_se))
    print(plot_numbers(data_filt))
    if (parameters$filtering_type == "condition") {
      print(plot_detect(data_filt))
    }
    if (parameters$filtering_type == "condition") {
      print(plot_normalization(data_filt, data_norm, data_imp))
    } else {
      print(plot_normalization(data_filt, data_norm))
    }
    if (parameters$filtering_type == "condition") {
      print(plot_imputation(data_norm, data_imp))
    }
    
    # Generate Volcano plot for each contrast
    for (p in 1:length(contrast)) {
      results <- get_results(dep) %>% select(name, ID, contains(contrast[p]) & (contains("p.val") | contains("ratio")))
      colnames(results) <- sub("ratio", "log2.FC", colnames(results))
      results <- results %>% mutate(!!paste0(contrast[p], "_p.adj") := p.adjust(get(paste0(contrast[p], "_p.val")), method = "BH"))
      results_A[[paste0(name_B, "_", contrast[p])]] <- results
      
      left <- str_split(contrast[p], pattern = "_vs_")[[1]][1]
      right <- str_split(contrast[p], pattern = "_vs_")[[1]][2]
      p_val_column <- grep("_p.val$", names(results), value = TRUE)
      log2_fc_column <- grep("_log2.FC$", names(results), value = TRUE)
      
      # Define proteins to highlight
      if (parameters$volcano == "protein list") {
        proteins_to_highlight <- str_split(parameters$proteins_to_highlight, ";")[[1]]
      }
      
      if (parameters$volcano == "specify significance") {
        proteins_to_highlight <- results$name[results[[p_val_column]] < parameters$p_value & abs(results[[log2_fc_column]]) > parameters$log2_FC]
      }
      
      if (parameters$volcano == "TopN") {
        sorted_results <- results[order(abs(results[[log2_fc_column]]), decreasing = TRUE), ]
        proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < 0.05, ], parameters$TopN)$name
      }
      
      # Volcano plot
      print(
        ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) +
          geom_point(color = "grey60", alpha = 0.5) +
          geom_point(data = results[results$name %in% proteins_to_highlight, ],
                     color = "#E41A1C",
                     alpha = 0.5) +
          geom_text_repel(data = results[results$name %in% proteins_to_highlight, ],
                          aes(label = name),
                          max.overlaps = Inf,
                          size = 3) +
          geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
          geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
          theme_bw(base_size = 12) +
          xlab(paste0("log2 FC (", left, " / ", right, ")")) +
          ylab("-log10(p-value)") +
          ggtitle(paste0(left, " vs ", right))
      )
      
      # Optional: highlighting of imputed proteins
      if (parameters$highlight_imputed == TRUE) {
        proteins_imputed <- get_df_long(data_norm) %>%
          group_by(name, condition) %>%
          summarize(NAs = any(is.na(intensity))) %>%
          filter(NAs) %>%
          pull(name) %>%
          unique()
        results$imputed <- ifelse(results$name %in% proteins_imputed, "Yes", "No")
        
        print(
          ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]), color = imputed, shape = imputed, fill = imputed)) +
            geom_point(alpha = 0.5) +
            scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
            scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            geom_text_repel(data = results[results$name %in% proteins_to_highlight, ],
                            aes(label = name),
                            max.overlaps = Inf,
                            size = 3) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
            geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
            theme_bw(base_size = 12) +
            xlab(paste0("log2 FC (", contrast, ")")) +
            ylab("-log10(p-value)") +
            ggtitle(paste0(left, " vs ", right)) +
            theme(
              legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
              legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
              legend.key = element_blank()
            ) +
            guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                   color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                   shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050")))) 
          
        )
      }
    }
    
    dev.off()
    
    # Optional: export complete dataset
    if (parameters$complete_output == TRUE) {
      
      datalist <- list()
      
      # Raw
      df <- as.data.frame(assay(data_filt))
      names(df) <- paste0(names(df), "_raw")
      df$name <- row.names(df)
      datalist[["raw"]] <- df
      
      # Normalized
      df <- as.data.frame(assay(data_norm))
      names(df) <- paste0(names(df), "_normalized")
      df$name <- row.names(df)
      datalist[["normalized"]] <- df  
      
      # Imputed
      if (parameters$filtering_type == "condition") {
        df <- as.data.frame(assay(data_imp))
        names(df) <- paste0(names(df), "_imputed")
        df$name <- row.names(df)
        datalist[["imputed"]] <- df
      }
      
      results_A[[paste0(name_B, "_complete_data")]] <- datalist %>%
        purrr::reduce(merge, by = "name") %>%          
        merge(get_results(dep)) %>%                    
        select(-contains("significant"),               
               -contains("centered"),
               -contains("p.adj")) %>%
        rename_with(~ sub("ratio", "log2.FC", .x)) %>%
        left_join(proteinGroups[, c("ID", "Protein.IDs")], by = "ID") %>%
        left_join(unique_peptides, by = "Protein.IDs") %>%
        select(-Protein.IDs)
    }
    
  }
  
  write.xlsx(results_A, paste0("DEP_results.xlsx"))
  rm(results_A)
}


# ------------------------------------------------------------
# Save workspace and session info
# ------------------------------------------------------------

save.image(paste0(getwd(), paste0("/", parameters$analysis_method, " - DEP analysis.RData")))
sessionInfo()
