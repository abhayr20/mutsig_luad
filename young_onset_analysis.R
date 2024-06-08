# Package Installation 
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

packages <- c("SigProfilerMatrixGeneratorR", "SigProfilerAssignmentR",
              "SigProfilerExtractorR", "SigProfilerPlottingR")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    remotes::install_github(paste0("AlexandrovLab/", package))
  }
}

# install('GRCh37', rsync=FALSE, bash=TRUE)

#Load libraries
library(reticulate)
library(devtools)
library(tidyverse)
library(dplyr)
library(maftools)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerAssignmentR)
library(SigProfilerExtractorR)
library(SigProfilerPlottingR)

setwd("/data/abhay/young_onset/")
use_condaenv('mutational_signatures')

metadata = read.csv("clinical_data/combined_df.csv") %>%
  mutate(patient_id = if_else(str_starts(patient_id, "NSLC"),
                              paste0(patient_id, "-T01"),
                              if_else(str_starts(patient_id, "TCGA"),
                                      paste0(patient_id, "-01"),
                                      patient_id)))

young_40_grp <- metadata %>%
  filter(age <= 40) %>%
  pull(patient_id)

young_45_grp <- metadata %>%
  filter(age <= 45) %>%
  pull(patient_id)

old_60_grp <- metadata %>%
  filter(age >= 60) %>%
  pull(patient_id)

# Generate mutational profiles analysis using SigProfilerMatrixGenerator
matrices <- SigProfilerMatrixGeneratorR(project = "combined_LUAD",
                                        genome = "GRCh37",
                                        matrix_path = "./signatures/SPMG",
                                        plot = F,
                                        exome = T)

# plotSBS(matrix_path = 'signatures/SPMG/output/SBS/combined_LUAD.SBS96.exome',
#         output_path = 'signatures/SPMG/output/SBS/',
#         project = 'young_onset_LUAD',
#         plot_type = '96',
#         percentage = FALSE)

#Get average mutational profiles
mut_matrix = matrices[['96']]
relative_mut_matrix = apply(mut_matrix, 2, prop.table)
average_mut_matrix = rowMeans(relative_mut_matrix)
average_mut_matrix = data.frame(Average_combined = average_mut_matrix)

# Add row names as column and print
average_mut_matrix_to_print = cbind(rownames(average_mut_matrix),
                                    average_mut_matrix)
colnames(average_mut_matrix_to_print)[1] = 'MutationType'
write.table(average_mut_matrix_to_print, 'signatures/avg_combined.SBS96.all',
            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles (note the percentage parameter now)
plotSBS(matrix_path = 'signatures/avg_combined.SBS96.all',
        output_path = 'signatures/',
        project = 'avg_combined',
        plot_type = '96',
        percentage = TRUE)


mm_young_40 <- as.data.frame(rowMeans(relative_mut_matrix[,young_40_grp]))
mm_young_45 <- as.data.frame(rowMeans(relative_mut_matrix[,young_45_grp]))
mm_old_60   <- as.data.frame(rowMeans(relative_mut_matrix[,old_60_grp]))
