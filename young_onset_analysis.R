# Package Installation 
packages <- c("SigProfilerMatrixGeneratorR", "SigProfilerAssignmentR",
              "SigProfilerExtractorR", "SigProfilerPlottingR")

for (package in packages) {
  if (!require(package, quietly = TRUE)) {
    devtools::install_github(paste0("AlexandrovLab/", package))
  }
}

#Load libraries
library(reticulate)
library(devtools)
library(tidyverse)
library(dplyr)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerAssignmentR)
library(SigProfilerExtractorR)
library(SigProfilerPlottingR)


setwd("/data/rajan_data/young_onset/luad_tcga_pan_can_atlas_2018/")
use_condaenv('mutational_signatures')
#install('GRCh37', rsync=FALSE, bash=TRUE)

#Read data
maf_cbioportal = read.delim('data_mutations.txt')

# Selection of specific columns needed by SigProfiler
maf_sp = maf_cbioportal %>%
  select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome,
         Start_Position, End_Position, Strand, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode)

# Filter for only considering single base substitutions
maf_sp = maf_sp %>%
  filter(Variant_Type == 'SNP')

dir.create('signatures') #for signature analysis results
dir.create('signatures/SPMG/') #for updated MAF file (needed for SigProfilerMatrixGenerator)

# Write updated MAF file
write.table(maf_sp, 'signatures/SPMG/data_mutations.maf', quote = F,
            row.names = F, sep = '\t')

# Generate mutational profiles analysis using SigProfilerMatrixGenerator
matrices <- SigProfilerMatrixGeneratorR(project = "TCGA_LUAD",
                                        genome = "GRCh37",
                                        matrix_path = "./signatures/SPMG",
                                        plot = F,
                                        exome = T)

plotSBS(matrix_path = 'signatures/SPMG/output/SBS/TCGA_LUAD.SBS96.exome',
        output_path = 'signatures/SPMG/output/SBS/',
        project = 'TCGA_LUAD',
        plot_type = '96',
        percentage = FALSE)

#Get average mutational profiles
mut_matrix = matrices[['96']]
relative_mut_matrix = apply(mut_matrix, 2, prop.table)
average_mut_matrix = rowMeans(relative_mut_matrix)
average_mut_matrix = data.frame(Average_TCGA_LUAD = average_mut_matrix)

# Add row names as column and print
average_mut_matrix_to_print = cbind(rownames(average_mut_matrix),
                                    average_mut_matrix)
colnames(average_mut_matrix_to_print)[1] = 'MutationType'
write.table(average_mut_matrix_to_print, 'signatures/avg_TCGA_LUAD.SBS96.all',
            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles (note the percentage parameter now)
plotSBS(matrix_path = 'signatures/avg_TCGA_LUAD.SBS96.all',
        output_path = 'signatures/',
        project = 'avg_TCGA_LUAD',
        plot_type = '96',
        percentage = TRUE)
