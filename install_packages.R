if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

packages <- c(
  "tidyverse",
  "openxlsx",
  "ggrepel",
  "SummarizedExperiment",
  "arnesmits/DEP",
  "cbielow/PTXQC"
)

pak::pkg_install(packages)