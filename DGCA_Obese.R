
library(dplyr)
library(tidyr)
library(ggplot2)
library(MEGENA)
library(tidyverse)
library(DGCA)

## Loading Data

load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

## Formatting data

NAFLD_O <- dplyr::select(NAFLD_O, -one_of(("Group")))
Control_O <- dplyr::select(Control_O, -one_of(("Group"))) 


##### NAFLD ##########
N_input_mat <-
   bind_rows(NAFLD_O, Control_O) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1]

N_design_mat <- 
  data.frame(Control_O = c(rep(0,nrow(NAFLD_O)),rep(1,nrow(Control_O))),
                           NAFLD_O= c(rep(1,nrow(NAFLD_O)),rep(0,nrow(Control_O)))) %>% 
   apply(., 2,as.numeric) %>% 
   `rownames<-`(rownames(N_input_mat))

Full_metadata <-  
   read.csv(file = "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", row.names = 1) 


excluded = c('BMI.body_weight_variables', 'Age.demographic_variable', 'Gender.demographic_variable', 
             'HOMAIR.glucose_insulin', 'Systolic_pressure.cardiac_variables','Group' )

Full_metadata <- 
   Full_metadata %>% 
   dplyr::select(.,-excluded)


Full_metadata <- 
   Full_metadata[which(rownames(Full_metadata) %in% rownames(N_input_mat)),]


Full_metadata = Full_metadata %>% as.data.frame() %>% .[order(rownames(.)),]
N_input_mat = N_input_mat %>% as.data.frame() %>% .[order(rownames(.)),]
N_design_mat = N_design_mat %>% as.data.frame() %>% .[order(rownames(.)),] %>% as.matrix()





library(DGCA)
set.seed(123)

## ALL Correlations for two groups
N_spearman_ddcor_res <- ddcorAll(inputMat = t(N_input_mat), 
                                 design = N_design_mat,
                                 compare = c("Control_O", "NAFLD_O"),
                                 adjust = "perm",
                                 heatmapPlot = FALSE, 
                                 nPerm = 1000, 
                                 corrType = "spearman") %>%
   filter(Classes != "NonSig") %>% 
   filter(empPVals < 0.05) 

## Remove unchanged correlation ###
N_spearman_ddcor_res <- 
   N_spearman_ddcor_res %>% 
   filter(Classes != "0/0")


## Find clusters ###
## Keep correlations with empval < 0.01 ###
N_megena_re_2 = ddMEGENA(N_spearman_ddcor_res %>% filter(empPVals < 0.01), 
                         adjusted = FALSE,
                         evalCompactness = TRUE,
                         pval_gene_thresh = 1,
                         hubPVal = 0.05,
                         modulePVal = 0.05,
                         nCores = 10)

result_list <- list(NAFLD_Cor_Res = N_spearman_ddcor_res %>% filter(empPVals < 0.01),
                    NAFLD_Megena_Res = N_megena_re_2)

# 
# saveRDS(object = result_list, file =
#            "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")




