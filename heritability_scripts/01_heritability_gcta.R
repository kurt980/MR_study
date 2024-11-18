---
title: "0912_Heritability"
author: "Chengyan Ji"
date: "2024-09-12"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
library(neuroCombat)
#library(snpStats)
library(tidyr)
library(tibble)
library(dplyr)
#library(AER)
library(readxl)
library(ggplot2)
library(readr)
```

**This is for preparing and running Heritability test in GCTA with full ABCD data and harmonized phenotypes**
**09/12/2024: Started editing**
**09/13/2024: Added harm0 and harm1**

# PART I. Preprocessing for data

## a. Need to get phenotypes harmonized.
### Utilities
```{r}
### Utilities
path <- "utilities/"
load(paste0(path, "coupling_cbcl.RData"))
load(paste0(path, "scfc_coupling_matched.RData"))
new_fam <- read.table(paste0(path, "acspsw03.txt"), header = T) # get new family ID
new_fam <- new_fam %>%
  select(rel_family_id, subjectkey, eventname, race_ethnicity)

CFDRs <- demo_cbcl_base_scfc %>%
  select(subjectkey, site_id_l, eventname, interview_age, sex, race, ethn) %>%
  mutate(
    subjectkey = sub("^(.{4})", "\\1_", subjectkey)  # Add '_' after the 4th character
    #ethn = as.factor(ethn),  # Convert ethnicity to a factor
    #race = as.factor(race)   # Convert race to a factor if needed
  ) %>%
  inner_join(new_fam, by = c("subjectkey", "eventname")) %>%
  select(rel_family_id, subjectkey, site_id_l, interview_age, sex, race, ethn, race_ethnicity) %>%
  rename(fid = rel_family_id, iid = subjectkey) %>%
  mutate(  # Add underscore to subjectkey
    race = gsub(" ", "_", race),                      # Replace spaces in 'race'
    ethn = gsub(" ", "_", ethn)                       # Replace spaces in 'ethnicity'
  )
# # WE DO NOT remove NAs this time!!!
# CFDRs <- na.omit(CFDRs) # should have 6208 left

# Fill NAs in race and ethnicity with mode
calculate_mode <- function(x) {
  ux <- unique(na.omit(x))  # Remove NA values before finding the mode
  ux[which.max(tabulate(match(x, ux)))]
}
CFDRs$race <- ifelse(is.na(CFDRs$race), calculate_mode(CFDRs$race), CFDRs$race)
CFDRs$ethn <- ifelse(is.na(CFDRs$ethn), calculate_mode(CFDRs$ethn), CFDRs$ethn)

# write new id list to local
fid_iid_scfc <- CFDRs %>%
  select(fid, iid)
write.table(fid_iid_scfc, paste0("bfile/fid_iid_scfc.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

##===================================================Covariates===================================================
qcovar_age <- CFDRs %>%
  select(fid, iid, interview_age) %>%
  arrange(iid)
covar_sex_race_ethn <- CFDRs %>%
  select(fid, iid, sex, race, ethn) %>%
  arrange(iid)
# Covariate files
write.table(qcovar_age, paste0("covariates/qcovar_age.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(covar_sex_race_ethn, paste0("covariates/covar_sex_race_ethn.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

##===================================================Phenotype====================================================
SCFC_coup_matrix <- as.data.frame(SCFC_coup %>% unlist() %>% matrix(ncol = 379, byrow = T))
SCFC_coup_matrix[is.na(SCFC_coup_matrix)] = 0 # 1. Fill NAs with 0 (for phenotype)
colnames(SCFC_coup_matrix) <- mmp_subcor
iid <- c() # Get family id and id; format by adding "-" to scfc ids
# Loop through each ID in the list
for (id in sub_SC_FC) {
  # print(id)
  id_format <- paste0(substring(id, 1, 4), "_", substring(id, 5))
  iid <- append(iid, id_format) # Append modified ID to the list
}
SCFC_DF <- data.frame(iid)
SCFC_DF <- cbind(SCFC_DF, SCFC_coup_matrix)
ALL_DF <- CFDRs %>%
  inner_join(SCFC_DF, by = "iid") # This is the big big DF that has everything
# Do Fisher Z Transformation
fisher_transform <- function(r, n) {
  z <- 0.5 * log((1 + r) / (1 - r))
  z_scaled <- z * sqrt(n - 3)
  return(z_scaled)
}
SCFC_coup_matrix_transformed <- apply(SCFC_coup_matrix, 2, function(x) fisher_transform(x, n = 379))

#--------------------------------------HARMONIZATION---------------------------------------
scfc_matrix <- t(as.matrix(SCFC_coup_matrix_transformed))
# get covariates matrix
interview_age = ALL_DF$interview_age
sex = ALL_DF$sex
race = ALL_DF$race
ethn = ALL_DF$ethn
batch <- ALL_DF$site_id_l
mod <- model.matrix(~interview_age+sex+race+ethn)
data_trans_harm0 <- scfc_matrix
data_trans_harm1 <- neuroCombat(dat=scfc_matrix, batch=batch)
data_trans_harm2 <- neuroCombat(dat=scfc_matrix, batch=batch, mod=mod)
# harm2 is harmonized for age, sex race and ethnicity information

#====================================================Write to local==========================================
scfc_trans_harm0 <- cbind(fid_iid_scfc, as.data.frame(t(data_trans_harm0)))
write.csv(scfc_trans_harm0, "utilities/scfc_harm0_sd.csv", row.names = FALSE)
scfc_trans_harm1 <- cbind(fid_iid_scfc, as.data.frame(t(data_trans_harm1$dat.combat)))
write.csv(scfc_trans_harm1, "utilities/scfc_harm1_sd.csv", row.names = FALSE)
scfc_trans_harm2 <- cbind(fid_iid_scfc, as.data.frame(t(data_trans_harm2$dat.combat)))
write.csv(scfc_trans_harm2, "utilities/scfc_harm2_sd.csv", row.names = FALSE)

harms = c(0,1,2)
for (harm in harms) {
  dir.create(paste0("phenotype_roi/harm", harm, ""), recursive = F)
}
for (i in 1:379) {
  phenotype_data_0 <- scfc_trans_harm0[, c(1, 2, i + 2)]
  file_name <- paste0("phenotype_roi/harm0/scfc_sd_harm0_family_roi", i, ".phen")
  write.table(phenotype_data_0, file = file_name, row.names = F, col.names = F, quote = F)
  phenotype_data_1 <- scfc_trans_harm1[, c(1, 2, i + 2)]
  file_name <- paste0("phenotype_roi/harm1/scfc_sd_harm1_family_roi", i, ".phen")
  write.table(phenotype_data_1, file = file_name, row.names = F, col.names = F, quote = F)
  phenotype_data_2 <- scfc_trans_harm2[, c(1, 2, i + 2)]
  file_name <- paste0("phenotype_roi/harm2/scfc_sd_harm2_family_roi", i, ".phen")
  write.table(phenotype_data_2, file = file_name, row.names = F, col.names = F, quote = F)
}
```

### Bfile: Need to parse and qc SNP data
这次不用去掉NA了，但是还是要改fid
```{r}
#====================================Get new family ID==================================================
old_fam <- read.table("bfile/ABCD_all_old.fam")
colnames(old_fam) <- c("fid", "iid", "pid", "mid", "sex", "phenotype")
new_fam <- old_fam %>%
  left_join(fid_iid_scfc, by = "iid") %>%
  mutate(fid = ifelse(is.na(fid.y), fid.x, fid.y)) %>%
  select(fid, iid, pid, mid, sex, phenotype)

write.table(new_fam, paste0("bfile/ABCD_all.fam"), row.names = F, col.names = F, quote = F, sep = "\t")
#=====================================Subset+QC==========================================================
dirpath <- '/Users/kurt/Documents/Data_Science_Project/0912_heritability/'
plink_path <- "/Users/kurt/Documents/Data_Science_Project/0912_heritability/"
# Subset by id
system(paste0(plink_path, "plink --bfile ", dirpath, "bfile/ABCD_all --keep ", dirpath, "bfile/fid_iid_scfc.txt --make-bed --out ", dirpath, "bfile/scfc_all"), intern = FALSE)
#-------------------------------------QC-----------------------------------------------------------------
system(paste0(plink_path, "plink --bfile ", dirpath, "bfile/scfc_all --geno 0.1 --maf 0.05 --mind 0.1 --make-bed --out ", dirpath, "bfile/scfc_all_qc"), intern = FALSE)

#=====================================grm and sparse grm=================================================
lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path
system(
  paste0(dirpath, "gcta64 ",
      "--bfile ", dirpath, "bfile/scfc_all_qc ",
      "--autosome  --make-grm  ", 
      "--thread-num 16 ",
      "--out ", dirpath, "bfile/scfc_all_qc_out"), 
      intern = F)
# system(
#   paste0(dirpath, "gcta64 ",
#       "--grm ", dirpath, "bfile/scfc_all_qc_out ",
#       "--make-bK-sparse 0.05 ",
#       "--thread-num 16 ",
#       "--out ", dirpath, "bfile/scfc_all_qc_sp_out"), 
#       intern = F)
```

# PART II. GCTA

## a. 
```{r}
sample <- 379
H2_harm_df <- data.frame(ROI = mmp_subcor)

## Loop through harm 0 1 and 2
for (harm in harms[1:3]) {
  dir.create(paste0(dirpath, "greml_results/harm", harm, ""), recursive = T)
  dir.create(paste0(dirpath, "greml_results/greml_results_csv/harm", harm, ""), recursive = T)
  library(stringr)
  files <- list.files(path = paste0(dirpath, 'phenotype_roi/harm', harm, "/"))
  phens <- files[stringr::str_detect(files, "\\.phen$")]
  file_pairs <- setNames(paste0('greml_', sub("\\.phen$", "", phens)), phens)
  # print(length(file_pairs))
  # print(file_pairs)
  
  lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
  system(lib_path, intern = TRUE)  # Set library path
  
  readme_path <- paste0(dirpath, "greml_results/harm", harm, "/readme.txt")
  readme_conn <- file(readme_path, "a")
  writeLines("Commands run for GREML analysis:", readme_conn)
  start_time <- Sys.time()
  
  for (idx in seq_along(file_pairs)) {
    if (idx <= sample) {
      phenotype_file <- names(file_pairs)[idx]
      out_file <- file_pairs[[idx]]
      
      command <- paste0(dirpath,
                        "gcta64 ",
                        "--grm ", dirpath, "bfile/scfc_all_qc_out ",
                        "--pheno ", dirpath, "phenotype_roi/harm", harm, "/", phenotype_file, " ",
                        "--qcovar ", dirpath, "covariates/qcovar_age.txt ",
                        "--covar ", dirpath, "covariates/covar_sex_race_ethn.txt ",
                        "--reml ",
                        "--thread-num 16 ",
                        "--out ", dirpath, "greml_results/harm", harm, "/", out_file
                 )
      
      system(command, intern = FALSE)
      writeLines(command, readme_conn, useBytes = TRUE)
      
      cat("\nGREML(all covariates) command executed successfully for the phenotype file:", phenotype_file, "\n\n\n")
    }
  }
  
  end_time <- Sys.time()
  total_duration <- difftime(end_time, start_time, units = "mins")
  writeLines(sprintf("Total time to run all commands: %s minutes", total_duration), readme_conn)
  writeLines("End of command list.", readme_conn)
  close(readme_conn)
}

for (harm in harms) {
  #===================================================Read Results========================================================
  readme_path <- paste0(dirpath, "greml_results/greml_results_csv/harm", harm, "/readme.txt")
  readme_conn <- file(readme_path, "w")
  writeLines("GREML using GCTA\n using all covariates\n\n\n", readme_conn)
  writeLines(paste0("Commands run for GREML analysis: ", command, "\n\n"), readme_conn)
  
  start_time <- Sys.time()
  
  h_2_values <- numeric(379)
  p_values <- numeric(379)
  SEs <- numeric(379)
  Ns <- numeric(379)
  missing <- 0
  for (i in 1:sample) {
    filename <- paste0(dirpath, "greml_results/harm", harm, "/greml_scfc_sd_harm", harm, "_family_roi", i, ".hsq")
    
    if (file.exists(filename)) {
      scfc_hsq <- read.table(filename, header = T, fill = T, sep = "")
      h_2_values[i] <- scfc_hsq[4, 2]
      p_values[i] <- scfc_hsq[9, 2]
      SEs[i] <- scfc_hsq[4,3]
      Ns[i] <- scfc_hsq[10,2]
      output_filename <- paste0(dirpath, "greml_results/greml_results_csv/harm", harm, "/greml_scfc_sd_harm", harm, "_family_roi", i, ".csv")
      write.csv(scfc_hsq, output_filename, row.names = FALSE, quote = FALSE)
    } else {
      missing = missing + 1
      writeLines(paste0("Missing file for ROI ", i, ": ", filename), readme_conn)
      cat("File does not exist: ", filename, "\n")  # Print message to the console as well
    }
  }
  print(h_2_values)
  end_time <- Sys.time()
  total_duration <- difftime(end_time, start_time, units = "mins")
  
  writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
  writeLines(paste("Total missing files:", missing), readme_conn)
  close(readme_conn)
  
  H2_harm_df[[paste0("h2_harm",harm)]]=h_2_values
  H2_harm_df[[paste0("p_value_harm",harm)]]=p_values
  H2_harm_df[[paste0("SE_harm",harm)]]=SEs 
  H2_harm_df[[paste0("N_harm",harm)]]=Ns
}
write.csv(H2_harm_df, paste0(dirpath, "greml_results/H2_harm_sd_df.csv"))
```
