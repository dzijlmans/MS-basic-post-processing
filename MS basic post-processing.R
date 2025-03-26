#Starting up
library(tidyverse)
library(DEP)
library(openxlsx)
library(ggrepel)
library(tcltk)
library(PTXQC)

setwd(tk_choose.dir())
parameters <- read.xlsx("Analysis_parameters.xlsx")

#DiaNN pre-processing
if (parameters$Analysis.method == "DiaNN") {
  tkmessageBox(message="Make sure that the 'pg_matrix.tsv' and 'pr_matrix.tsv' files are located in your working directory.
  
  Please ensure that the abundance column names (and thus the .raw files) contain the name of the mass spectrometer used (e.g. Astral, Exploris) and the following information, in order and separated by an underscore:

- Condition A (e.g. Cell line)
- Condition B (e.g. WT and KO or treated and untreated)
- Replicate
             
Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw")
  
  proteinGroups <- read.delim(list.files(pattern = "pg_matrix.tsv"))
  peptides <- read.delim(list.files(pattern = "pr_matrix.tsv"))
  
  filtered_proteinGroups <- data.frame(
    Protein.Group = names(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group)),
    Peptide.Count = as.integer(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group))
  ) %>% 
    filter(Peptide.Count > (parameters$unique_peptides - 1))
  
  proteinGroups <- proteinGroups[proteinGroups$Protein.Group %in% filtered_proteinGroups$Protein.Group, ]
  
  proteinGroups <-
    make_unique(proteinGroups, "Genes", "Protein.Group", delim = ";")
  proteinGroups <- proteinGroups %>% filter(!grepl("Cont_", ID) &
                                              !grepl("KRT", name) &
                                              !grepl("cRAP", name))
  
  value_columns <- which(str_detect(names(proteinGroups), parameters$mass_spec))
  prefix <- LCSn(colnames(proteinGroups[, value_columns]))
  colnames(proteinGroups) <- gsub(".raw", "", gsub(".*" %>% paste0(prefix), "", colnames(proteinGroups)))
  
}


#MaxQuant pre-processing
if (parameters$Analysis.method == "MaxQuant") {
  tkmessageBox(message="Make sure that the 'proteinGroups.txt' file is located in your working directory.
  
  Please ensure that the abundance column names (and thus the .raw files) contain the following information, in order and separated by an underscore:

- Condition A (e.g. Cell line)
- Condition B (e.g. WT and KO or treated and untreated)
- Replicate
             
Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw")
  
  proteinGroups <- read.delim(list.files(pattern = "proteinGroups.txt"))
  
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
  
  proteinGroups <- proteinGroups %>% dplyr::select(
    Protein.IDs = "Protein.IDs",
    "Gene.names",
    "Unique.peptides",
    contains("LFQ.intensity."),
  )
  proteinGroups <-
    make_unique(proteinGroups, "Gene.names", "Protein.IDs", delim = ";")
  
  value_columns <- which(str_detect(names(proteinGroups), "LFQ.intensity"))               #specify Mass spec (this should be in the value column names)
  prefix <- LCSn(colnames(proteinGroups[, value_columns]))
  colnames(proteinGroups) <- gsub(".raw", "", gsub(".*" %>% paste0(prefix), "", colnames(proteinGroups)))
  
}


#Proteome Discoverer pre-processing
if (parameters$Analysis.method == "Proteome Discoverer") {
  tkmessageBox(message="Make sure that the 'Proteins.txt' file is located in your working directory.
  
  Please ensure that the abundance column names (and thus the .raw files) contain the following information, in order and separated by an underscore:

- Condition A (e.g. Cell line)
- Condition B (e.g. WT and KO or treated and untreated)
- Replicate
             
Example: 20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw")
  
  proteinGroups <- read.delim(list.files(pattern = "Proteins.txt"))
  
  proteinGroups <- filter(proteinGroups, proteinGroups[, grep("Unique", colnames(proteinGroups))] > (parameters$unique_peptides - 1))
  
  proteinGroups <- proteinGroups %>% dplyr::select(
    "Accession",
    colnames(proteinGroups)[grep("Unique", colnames(proteinGroups))],
    "Gene.Symbol",
    contains("Abundances") & contains("Normalized")
  )
  
  proteinGroups <-
    make_unique(proteinGroups, "Gene.Symbol", "Accession", delim = ";")
  
  value_columns <- which(str_detect(names(proteinGroups), "Abundances..Normalized"))
  # Use gsub to trim everything before "Sample" or "Control"
  colnames(proteinGroups) <- gsub("\\.\\.", "_", gsub(".*(Sample|Control)\\.\\.(.*)", "\\2", colnames(proteinGroups)))
  
}


##QC
try(if (max(str_count(colnames(proteinGroups[, value_columns]), "_")) == 1 &
        parameters$conditions != 2 |
        max(str_count(colnames(proteinGroups[, value_columns]), "_")) == 2 &
        parameters$conditions != 3)
  stop("incorrect number of conditions"))



if (parameters$conditions == 2) {
  results_A <- list()
  data <- proteinGroups
  experimental_design <- data.frame(label = colnames(data[, value_columns]))
  experimental_design <- experimental_design %>% separate_wider_delim(
    label,
    delim = "_",
    names = c("condition", "replicate"),
    cols_remove = FALSE
  )
  experimental_design <- experimental_design[, c("label", "condition", "replicate")]
  
  #Make object and do analysis
  data_se <- make_se(data, value_columns, experimental_design)
  data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
  data_norm <- normalize_vsn(data_filt)
  if (parameters$filtering_type == "complete") {
    data_imp <- data_norm
  }
  if (parameters$filtering_type == "condition") {
    data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
  }
  data_diff_manual <- test_diff(data_imp, type = parameters$comparison, control = parameters$control)
  contrast <- sub("_diff", "", str_subset(
    names(data_diff_manual@elementMetadata@listData),
    "_diff"
  ))
  dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(2))
  
  ####make pdf containing QCs
  pdf("Output_plots.pdf")
  
  print(plot_frequency(data_se))
  print(plot_numbers(data_filt))
  if (parameters$filtering_type == "condition") {
    print(plot_detect(data_filt))
  }
  print(plot_normalization(data_filt, data_norm, data_imp))
  print(plot_imputation(data_norm, data_imp))
  
  #add dataframe for each contrast to list
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
    
    if (parameters$volcano == "protein list") {
      proteins_to_highlight <- str_split(parameters$proteins_to_highlight, "; ")[[1]]
    }
    
    if (parameters$volcano == "specify significance") {
      proteins_to_highlight <- results$name[results[[p_val_column]] < parameters$p_value & abs(results[[log2_fc_column]]) > parameters$log2_FC]
    }
    
    if (parameters$volcano == "TopN") {
      sorted_results <- results[order(abs(results[[log2_fc_column]]), decreasing = TRUE), ]
      proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < 0.05, ], parameters$TopN)$name
    }
    
    print(ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) + 
      geom_point(color = "grey50", alpha = 0.5) +
      geom_point(data = results[results$name%in%proteins_to_highlight, ], aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]])), color = "red", alpha = 0.5) +
      geom_text_repel(data = results[results$name%in%proteins_to_highlight, ], aes(label = name)) +
      geom_hline(yintercept = -log10(parameters$p_value),linetype="dashed", color = "black") +
      geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC),linetype="dashed", color = "black") +
      theme_bw(base_size = 12) +
      xlab(paste0("log2 FC (", left, " / ", right, ")")) +
      ylab("-log10(p-value)") +
      ggtitle(paste0(left, " vs ", right)))
  }
  dev.off()
  write.xlsx(results_A, paste0("DEP_results.xlsx"))
  rm(results_A)
}

if (parameters$conditions == 3) {
  results_A <- list()
  proteinGroups <-
  proteinGroups %>% pivot_longer(
    cols = all_of(value_columns),
    names_to = c("condition_B", "condition_A", "replicate"),                 #change this depending on the amount of different variables used in your experiment
    names_sep = "_",
    values_to = "values"
  )
  proteinGroups <-
    proteinGroups %>% pivot_wider(names_from = c("condition_A", "replicate"),         #change this if necessary
                                  values_from = "values")
  datalist_B <- split(proteinGroups, f = ~ proteinGroups$condition_B)               #specify variable to be split on
  for (i in 1:length(datalist_B)) {
    name_B <- names(datalist_B[i])
    data <- datalist_B[[i]]
    value_columns <- which(sapply(data, is.numeric))
    
    experimental_design <- data.frame(label = colnames(data[, value_columns]))
    experimental_design <- experimental_design %>% separate_wider_delim(
      label,
      delim = "_",
      names = c("condition", "replicate"),
      cols_remove = FALSE
    )
    experimental_design <- experimental_design[, c("label", "condition", "replicate")]
    data_se <- make_se(data, value_columns, experimental_design)
    data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
    data_norm <- normalize_vsn(data_filt)
    if (parameters$filtering_type == "complete") {
      data_imp <- data_norm
    }
    if (parameters$filtering_type == "condition") {
      data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
    }
    data_diff_manual <- test_diff(data_imp, type = parameters$comparison, control = parameters$control)
    contrast <- sub("_diff", "", str_subset(
      names(data_diff_manual@elementMetadata@listData),
      "_diff"
    ))
    dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(2))
    
    ####make pdf containing QCs
    pdf(paste0("Output_plots_", name_B, ".pdf"))
    
    print(plot_frequency(data_se))
    print(plot_numbers(data_filt))
    if (parameters$filtering_type == "condition") {
      print(plot_detect(data_filt))
    }
    print(plot_normalization(data_filt, data_norm, data_imp))
    print(plot_imputation(data_norm, data_imp))
    
    #add all results to list
    results_A[[paste0(name_B)]] <- get_results(dep)
    #add dataframe for each contrast to list
    for (p in 1:length(contrast)) {
      results <- get_results(dep) %>% select(name, ID, contains(contrast[p]) & (contains("p.val") | contains("ratio")))
      colnames(results) <- sub("ratio", "log2.FC", colnames(results))
      results <- results %>% mutate(!!paste0(contrast[p], "_p.adj") := p.adjust(get(paste0(contrast[p], "_p.val")), method = "BH"))
      
      results_A[[paste0(name_B, "_", contrast[p])]] <- results
      
      left <- str_split(contrast[p], pattern = "_vs_")[[1]][1]
      right <- str_split(contrast[p], pattern = "_vs_")[[1]][2]
      p_val_column <- grep("_p.val$", names(results), value = TRUE)
      log2_fc_column <- grep("_log2.FC$", names(results), value = TRUE)
      
      if (parameters$volcano == "protein list") {
        proteins_to_highlight <- str_split(parameters$proteins_to_highlight, "; ")[[1]]
      }
      
      if (parameters$volcano == "specify significance") {
        proteins_to_highlight <- results$name[results[[p_val_column]] < parameters$p_value & abs(results[[log2_fc_column]]) > parameters$log2_FC]
      }
      
      if (parameters$volcano == "TopN") {
        sorted_results <- results[order(abs(results[[log2_fc_column]]), decreasing = TRUE), ]
        proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < 0.05, ], parameters$TopN)$name
      }
      
      print(ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) + 
              geom_point(color = "grey50", alpha = 0.5) +
              geom_point(data = results[results$name%in%proteins_to_highlight, ], aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]])), color = "red", alpha = 0.5) +
              geom_text_repel(data = results[results$name%in%proteins_to_highlight, ], aes(label = name)) +
              geom_hline(yintercept = -log10(parameters$p_value),linetype="dashed", color = "black") +
              geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC),linetype="dashed", color = "black") +
              theme_bw(base_size = 12) +
              xlab(paste0("log2 FC (", left, " / ", right, ")")) +
              ylab("-log10(p-value)") +
              ggtitle(paste0(left, " vs ", right)))
    }
    dev.off()
  }
  write.xlsx(results_A, paste0("DEP_results.xlsx"))
  rm(results_A)
}


save.image(paste0(getwd(), paste0("/", parameters$Analysis.method, " - DEP analysis.RData")))

