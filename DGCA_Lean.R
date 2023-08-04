
#### Alpha Diversity ###
rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MEGENA)
library(tidyverse)
library(DGCA)


load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

NAFLD_L <- dplyr::select(NAFLD_L, -one_of(("Group"))) 
Control_L <- dplyr::select(Control_L, -one_of(("Group"))) 



##### Formating data and metadata ##########
N_input_mat <-
   bind_rows(NAFLD_L, Control_L) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1]


N_design_mat <- data.frame(Control_L = c(rep(0,nrow(NAFLD_L)),rep(1,nrow(Control_L))),
                           NAFLD_L= c(rep(1,nrow(NAFLD_L)),rep(0,nrow(Control_L)))) %>% 
   apply(., 2,as.numeric) %>% 
   `rownames<-`(rownames(N_input_mat))



Full_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", row.names = 1) 


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






## Correlations for Lean for two groups###
N_spearman_ddcor_res <- ddcorAll(inputMat = t(N_input_mat), 
                                   design = N_design_mat,
                                   compare = c("Control_L", "NAFLD_L"),
                                   adjust = "perm",
                                   heatmapPlot = FALSE, 
                                   nPerm = 1000, 
                                   corrType = "spearman") %>% 
   filter(Classes != "NonSig") %>% 
   filter(empPVals < 0.05)


## Removing unchanged correlation
N_spearman_ddcor_res <- 
   N_spearman_ddcor_res %>% 
   filter(Classes != "0/0")


## Finding clusters with correlations emppval < 0.01

N_megena_re_2 = ddMEGENA(N_spearman_ddcor_res %>%   filter(empPVals < 0.01), 
                         adjusted = FALSE,
                         evalCompactness = TRUE,
                         pval_gene_thresh = 1,
                         hubPVal = 0.05,
                         modulePVal = 0.05,
                         nCores = 10)

result_list <- list(NAFLD_Cor_Res = N_spearman_ddcor_res %>% filter(empPVals < 0.01),
                    NAFLD_Megena_Res = N_megena_re_2)

saveRDS(object = result_list, file = 
           "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")

