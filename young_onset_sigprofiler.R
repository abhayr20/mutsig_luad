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


# Generate mutational profiles analysis using SigProfilerMatrixGenerator
matrices <- SigProfilerMatrixGeneratorR(project = "combined_LUAD",
                                        genome = "GRCh37",
                                        matrix_path = "./signatures/SPMG",
                                        plot = F,
                                        exome = T)

# This will plot the 96 channel mutational profile for all 1060 patients in combined LUAD
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


#Subgroup Analysis by age

metadata = read.csv("clinical_data/combined_df.csv") %>%
  mutate(patient_id = if_else(str_starts(patient_id, "NSLC"),
                              paste0(patient_id, "-T01"),
                              if_else(str_starts(patient_id, "TCGA"),
                                      paste0(patient_id, "-01"),
                                      patient_id)))

perform_subset_analysis <- function(metadata, relative_mut_matrix, age_limit, matrix_name) {
  # Filter patient IDs based on age limit and presence in relative_mut_matrix
  patient_grp <- metadata %>%
    filter(age <= age_limit, patient_id %in% colnames(relative_mut_matrix)) %>%
    pull(patient_id)
  
  # Calculate row means for the filtered group
  mm_grp <- as.data.frame(rowMeans(relative_mut_matrix[, patient_grp, drop = FALSE]))
  
  # Prepare the data frame for output
  mm_grp_to_print <- cbind(rownames(mm_grp), mm_grp)
  colnames(mm_grp_to_print) <- c('MutationType', matrix_name)
  
  # Generate file names based on the subgroup name
  output_file <- paste0('signatures/avg_', matrix_name, '.SBS96.all')
  project_name <- paste0('avg_', matrix_name, '.SBS96.all')
  
  # Write the subgroup relative matrix output to a file
  write.table(mm_grp_to_print, output_file, quote = FALSE, row.names = FALSE, sep = '\t')
  
  message("Subset analysis completed and written to ", output_file)
  
  # Plot average mutational profiles
  plotSBS(matrix_path = output_file,
          output_path = 'signatures/',
          project = project_name,
          plot_type = '96',
          percentage = TRUE)
}

# Perform the analysis for different groups
perform_subset_analysis(metadata, relative_mut_matrix, 40, 'young_40')
perform_subset_analysis(metadata, relative_mut_matrix, 45, 'young_45')
perform_subset_analysis(metadata, relative_mut_matrix, 60, 'old_60')

# De-novo mutational signature extraction
start_time1 <- Sys.time()
sigprofilerextractor(input_type = 'matrix',
                     output = 'signatures/SPE/',
                     input_data = 'signatures/SPMG/output/SBS/combined_LUAD.SBS96.exome',
                     nmf_replicates = 100,
                     minimum_signatures = 1,
                     maximum_signatures = 25,
                     exome = T)
end_time1 <- Sys.time()
execution_time1 <- end_time1 - start_time1 # takes 2-3 hours

# Assign reference mutational signatures
start_time2 <- Sys.time()
cosmic_fit(samples = 'signatures/SPMG/output/SBS/combined_LUAD.SBS96.exome',
           output = 'signatures/SPA',
           input_type='matrix',
           exome = T)
end_time2 <- Sys.time()
execution_time2 <- end_time2 - start_time2 # takes 15 mins


