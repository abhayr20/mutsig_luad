library(maftools)
library(dplyr)

setwd("/data/abhay/young_onset/")

#Clinical parameters for 3 studies
tcga_luad <- read.delim("tcga_luad/data_clinical_patient.txt", skip=4, header = TRUE) %>% 
  select(PATIENT_ID, AGE, SEX, AJCC_PATHOLOGIC_TUMOR_STAGE, GENETIC_ANCESTRY_LABEL, RACE) %>%
  rename(patient_id = PATIENT_ID, age = AGE,
         sex = SEX, stage = AJCC_PATHOLOGIC_TUMOR_STAGE,
         ancestry = GENETIC_ANCESTRY_LABEL, race = RACE) %>%
  mutate(smoking_status = NA) %>%
  mutate(study = "TCGA_LUAD") %>%
  mutate(stage = case_when(
    grepl("IV", stage) ~ "IV",
    grepl("III", stage) ~ "III",
    grepl("II", stage) ~ "II",
    grepl("I", stage) ~ "I",
    stage == "" ~ NA_character_,
    TRUE ~ stage
  )) %>%
  mutate(ancestry = if_else(race == "White" & ancestry != "EUR", "EUR", ancestry)) %>%
  mutate(ancestry = if_else(race == "Black or African American", "AFR", ancestry)) %>%
  select(patient_id, age, sex, stage, ancestry, race, smoking_status, study)


oncosg = read.delim("oncosg/data_clinical_patient.txt", skip=4, header = TRUE) %>% 
  select(PATIENT_ID, AGE, SEX, STAGE, ETHNICITY, SMOKING_STATUS) %>%
  rename(patient_id = PATIENT_ID, age = AGE,
         sex = SEX, stage = STAGE,
         ancestry = ETHNICITY, smoking_status = SMOKING_STATUS) %>%
  mutate(study = "ONCOSG",
         race = "EAS",
         smoking_status = ifelse(smoking_status == "Yes", "smoker", "non_smoker")) %>%
  select(patient_id, age, sex, stage, ancestry, race, smoking_status, study)


sherlock <- read.delim("sherlock/data_clinical_patient.txt", skip=4, header = TRUE) %>%
  select(PATIENT_ID, AGE, SEX, TUMOR_STAGE, HISTOLOGY) %>%
  filter(HISTOLOGY == "Adenocarcinomas") %>%
  mutate(ancestry = "EUR",
         race = "White",
         smoking_status = 'non_smoker',
         study = "SHERLOCK") %>%
  rename(patient_id = PATIENT_ID, age = AGE,
         sex = SEX, stage = TUMOR_STAGE) %>%
  mutate(age = as.integer(round(age)))%>%
  mutate(stage = case_when(
    grepl("IV", stage) ~ "IV",
    grepl("III", stage) ~ "III",
    grepl("II", stage) ~ "II",
    grepl("I", stage) ~ "I",
    stage == "" ~ NA_character_,
    TRUE ~ stage
  )) %>%
  select(patient_id, age, sex, stage, ancestry, race, smoking_status, study)


#Patients according to age
combined_df <- bind_rows(tcga_luad, oncosg, sherlock) %>%
  mutate_all(~ifelse(. == "", NA, .))

young_40 <- combined_df %>%
  filter(age <= 40)

young_45 <- combined_df %>%
  filter(age <= 45)

old_60 <- combined_df %>% 
  filter(age >= 60)


# Save the dataframes
dir.create("clinical_data")
write.csv(tcga_luad, file = "clinical_data/tcga_luad.csv", row.names = FALSE)
write.csv(oncosg, file = "clinical_data/oncosg.csv", row.names = FALSE)
write.csv(sherlock, file = "clinical_data/sherlock.csv", row.names = FALSE)
write.csv(combined_df, file = "clinical_data/combined_df.csv", row.names = FALSE)
write.csv(young_40, file = "clinical_data/young_40.csv", row.names = FALSE)
write.csv(young_45, file = "clinical_data/young_45.csv", row.names = FALSE)
write.csv(old_60, file = "clinical_data/old_60.csv", row.names = FALSE)


# Read mutational data
maf_tcga_luad <- read.delim('tcga_luad/data_mutations.txt')
maf_oncosg <- read.delim('oncosg/data_mutations.txt', skip=2, header = TRUE)
maf_sherlock <- read.delim('sherlock/data_mutations.txt')

# Columns to keep in maf
cols_to_keep <- c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
                  "Start_Position", "End_Position", "Strand", "Variant_Classification",
                  "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
                  "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode")

# Selection of specific columns needed by SigProfiler
process_maf <- function(maf_df) {
  maf_df %>%
    select(all_of(cols_to_keep)) %>%
    filter(Variant_Type == 'SNP')
}

# Filter and select columns for each MAF file
maf_tcga_luad <- process_maf(maf_tcga_luad)
maf_oncosg <- process_maf(maf_oncosg)
maf_sherlock <- process_maf(maf_sherlock)
combined_maf <- bind_rows(maf_tcga_luad, maf_oncosg, maf_sherlock)

dir.create('signatures')       #for signature analysis results
dir.create('signatures/SPMG/') #for updated MAF file (needed for SigProfilerMatrixGenerator)
write.table(combined_maf, 'signatures/SPMG/data_mutations.maf', quote = F,
            row.names = F, sep = '\t')