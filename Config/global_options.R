knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center",
  fig.width = 9,
  fig.height = 6,
  output.dir = "/Users/bladenmaxwell/Documents/Cincinnati work/MB Diabetes Analysis/Reports"
)

wd <- "/Users/bladenmaxwell/Documents/Cincinnati work/MB Diabetes Analysis/"
installPackages <- FALSE

if (installPackages) {
  install.packages("BiocManager")
  BiocManager::install("mixOmics")
  
  install.packages("UniprotR")
  install.packages("readxl")
  install.packages("data.table")
}
