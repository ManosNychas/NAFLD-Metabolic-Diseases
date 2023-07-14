
#### Alpha Diversity ###
rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyr)
library(RANN)
library(ggplot2)
library(ggsignif)
library(vegan)
library(VennDiagram)
library(rcompanion)
#library(DGCA)
library(MEGENA)
library(tidyverse)
library(data.table)
library(lsr)
library(zCompositions)
library(compositions)

load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")



NAFLD_L <- dplyr::select(NAFLD_L, -one_of(("Group"))) #%>% 
# mutate_if(is.numeric, ~round(.*10^6)) %>% 
#.[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 

# zCompositions::cmultRepl() %>%
# compositions::clr()  
Control_L <- dplyr::select(Control_L, -one_of(("Group"))) #%>% 
# mutate_if(is.numeric, ~round(.*10^6)) %>% 
# .[ ,colSums(.!= 0) >= nrow(.)*0.1]#%>%
# zCompositions::cmultRepl() %>%
# compositions::clr()  



##### NAFLD ##########
N_input_mat <-
   bind_rows(NAFLD_L, Control_L) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1]
# zCompositions::cmultRepl() %>%
# compositions::clr() #%>% 
#t(.)

N_design_mat <- data.frame(Control_L = c(rep(0,nrow(NAFLD_L)),rep(1,nrow(Control_L))),
                           NAFLD_L= c(rep(1,nrow(NAFLD_L)),rep(0,nrow(Control_L)))) %>% 
   apply(., 2,as.numeric) %>% 
   `rownames<-`(rownames(N_input_mat))



Full_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", row.names = 1) 


var_adjust <-
   Full_metadata[,c('BMI.body_weight_variables', 'Age.demographic_variable', 'Gender.demographic_variable')] %>% 
   data.frame() 

excluded = c('BMI.body_weight_variables', 'Age.demographic_variable', 'Gender.demographic_variable', 
             'HOMAIR.glucose_insulin', 'Systolic_pressure.cardiac_variables','Group' )

Full_metadata <- 
   Full_metadata %>% 
   dplyr::select(.,-excluded)

# rownames(BMI) <- rownames(allData)
BMI <- var_adjust$BMI.body_weight_variables
age <- var_adjust$Age.demographic_variable
gender <- var_adjust$Gender.demographic_variable

Full_metadata <- 
   Full_metadata[which(rownames(Full_metadata) %in% rownames(N_input_mat)),]


Full_metadata = Full_metadata %>% as.data.frame() %>% .[order(rownames(.)),]
N_input_mat = N_input_mat %>% as.data.frame() %>% .[order(rownames(.)),]
var_adjust = var_adjust %>% as.data.frame() %>% .[order(rownames(.)),]
N_design_mat = N_design_mat %>% as.data.frame() %>% .[order(rownames(.)),] %>% as.matrix()



var_adjust_1 <- 
   var_adjust %>%
   as.data.frame() %>%
   rownames_to_column() 

Partial_mat <- 
   N_input_mat %>%
   as.data.frame() %>%
   rownames_to_column() %>%
   pivot_longer(2:length(.),names_to = "species") %>%
   left_join(var_adjust_1,by = c("rowname" = "rowname")) %>%
   group_by(species) %>%
   mutate(res = residuals(lm(formula = value ~  
                                BMI.body_weight_variables +
                                Age.demographic_variable + 
                                Gender.demographic_variable))) %>%
   dplyr::select(rowname,species,res) %>%
   pivot_wider(names_from = species, values_from = res) %>%
   column_to_rownames() %>%
   as.matrix() %>% 
   t()





library(DGCA)
source("Documents/DGCA/DGCA/R/ddcorAll.R")
source("Documents/DGCA/DGCA/R/getDCorPerm_multi.R")
## NAFLD Version ###
N_spearman_ddcor_res <- ddcorAll(inputMat = t(N_input_mat), 
                                   design = N_design_mat,
                                   compare = c("Control_L", "NAFLD_L"),
                                   adjust = "perm",
                                   heatmapPlot = FALSE, 
                                   nPerm = 1000, 
                                   corrType = "spearman") %>% 
   filter(Classes != "NonSig") %>% 
   filter(empPVals < 0.05)

N_spearman_ddcor_res_ADJ <- ddcorAll(inputMat = Partial_mat, 
                         design = N_design_mat,
                         compare = c("Control_L", "NAFLD_L"),
                         adjust = "perm",
                         heatmapPlot = FALSE, 
                         nPerm = 1000, 
                         corrType = "spearman") %>% 
   filter(Classes != "NonSig") %>% 
   filter(empPVals < 0.05)

N_spearman_ddcor_res <- 
   N_spearman_ddcor_res %>% 
   filter(Classes != "0/0")


N_spearman_ddcor_res_ADJ <- 
   N_spearman_ddcor_res_ADJ %>% 
   filter(Classes != "0/0")



N_megena_res = ddMEGENA(N_spearman_ddcor_res, 
                        adjusted = FALSE,
                        evalCompactness = TRUE,
                        pval_gene_thresh = 1,
                        hubPVal = 0.05,
                        modulePVal = 0.05,
                        nCores = 10)


result_list <- list(NAFLD_Cor_Res = N_spearman_ddcor_res ,
                    NAFLD_Megena_Res = N_megena_res)


saveRDS(object = result_list, file = 
           "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")




N_megena_re_2 = ddMEGENA(N_spearman_ddcor_res %>%    filter(empPVals < 0.01), 
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



N_megena_res_ADJ = ddMEGENA(N_spearman_ddcor_res_ADJ, 
                            adjusted = FALSE,
                            evalCompactness = TRUE,
                            pval_gene_thresh = 1,
                            hubPVal = 0.05,
                            modulePVal = 0.05,
                            nCores = 10)



result_list <- list(NAFLD_Cor_Res = N_spearman_ddcor_res_ADJ,
                    NAFLD_Megena_Res = N_megena_res_ADJ)

saveRDS(object = result_list, file = 
           "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/ADJ_Prev_All_Strict_DGCA_Control_L_vs_NAFLD_L.rds")




N_megena_re_2_ADJ = ddMEGENA(N_spearman_ddcor_res_ADJ %>%   filter(empPVals < 0.01), 
                             adjusted = FALSE,
                             evalCompactness = TRUE,
                             pval_gene_thresh = 1,
                             hubPVal = 0.05,
                             modulePVal = 0.05,
                             nCores = 10)




result_list <- list(NAFLD_Cor_Res = N_spearman_ddcor_res_ADJ %>% filter(empPVals < 0.01),
                    NAFLD_Megena_Res = N_megena_re_2_ADJ)

saveRDS(object = result_list, file = 
           "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/ADJ_Strict_Prev_All_Strict_DGCA_Control_L_vs_NAFLD_L.rds")





library(tibble)
library(ppcor)
'%!in%' <- function(x,y)!('%in%'(x,y))
calulo_corr <- function(cdata_matrix, spe_matrix){
   
   pval_ <- data.frame(matrix(NA, nrow = ncol(spe_matrix), ncol = ncol(cdata_matrix) ))
   colnames(pval_) <- colnames(cdata_matrix)
   rownames(pval_) <- colnames(spe_matrix)
   cor_ <- pval_
   for (i in colnames(cdata_matrix)) {
      for (j in colnames(spe_matrix )) {
         all <- na.omit(data.frame(spe_matrix[, j], cdata_matrix[,i]))
         tryCatch({
            #ctest <- pcor.test(all[, 1], all[,2], all[,3:ncol(all)], method = 'spearman' )
            ctest <- cor.test(all[, 1], all[,2],  method = 'spearman')
            pval <- ctest$p.value
            cor <- ctest$estimate
            pval_[j,i] <- pval
            cor_[j,i] <- cor
         }, error=function(e){})
      }
   }
   return(list(pvalue=pval_, correlation=cor_))
}
calulo_corr_adjusted <- function(var_adjust, cdata_matrix, spe_matrix ){
   
   pval_ <- data.frame(matrix(NA, nrow = ncol(spe_matrix), ncol = ncol(cdata_matrix) ))
   colnames(pval_) <- colnames(cdata_matrix)
   rownames(pval_) <- colnames(spe_matrix)
   cor_ <- pval_
   for (i in colnames(cdata_matrix)) {
      for (j in colnames(spe_matrix )) {
         all <- na.omit(data.frame(spe_matrix[, j], cdata_matrix[,i], var_adjust))
         tryCatch({
            ctest <- pcor.test(all[, 1], all[,2], all[,3:ncol(all)], method = 'spearman' )
            #ctest <- cor.test(all[, 1], all[,2],  method = 'spearman' )
            pval <- ctest$p.value
            cor <- ctest$estimate
            pval_[j,i] <- pval
            cor_[j,i] <- cor
         }, error=function(e){})
      }
   }
   return(list(pvalue=pval_, correlation=cor_))
}
Fisher_Procedure <-function(Network_Object,adjusted,clusters_allowed){
   # # 
   # Network_Object = N_megena_res
   # 
   # adjusted  = TRUE
   
   
   Clusters = data.frame()
   for (i in 1:length(Network_Object$summary$modules)) {
      current = Network_Object$summary$modules[[i]]
      
      Clusters = rbind(Clusters,
                       data.frame(Cluster = names(Network_Object$summary$modules)[i],
                                  Species =current ))
   }
   
   
   ###3 Correlation ###
   Full_metadata <- 
      Full_metadata[which(rownames(Full_metadata) %in% rownames(full_features)),] 
   
   var_adjust <- 
      var_adjust[which(rownames(var_adjust) %in% rownames(full_features)),] 
   
   
   Full_metadata = Full_metadata[order(rownames(Full_metadata)),] %>% 
      dplyr::select(FLI.liver_function,
                    ALT.liver_function,
                    AST.liver_function,
                    Triglyceride.lipid_profiles,
                    GGT.liver_function,
                    Liver_Fat.liver_function)
   full_features = full_features[order(rownames(full_features)),] %>% as.data.frame()
   var_adjust = var_adjust[order(rownames(var_adjust)),]
   
   
   # 
   ####  Spearman Correlation ### 
   corr_pval_data = data.frame()
   if(adjusted == TRUE) {
      corr_pval_data <- calulo_corr_adjusted(var_adjust = var_adjust, cdata_matrix = Full_metadata , spe_matrix = full_features)
   } else {
      corr_pval_data <- calulo_corr(cdata_matrix = Full_metadata , spe_matrix = full_features)
   }
   
   corr_mat <- as.matrix(corr_pval_data$correlation)
   pmat <- as.matrix(corr_pval_data$pvalue)
   
   corr_mat_long <- melt(corr_mat) %>% magrittr::set_colnames(c("Species","Metadata","Correlation")) 
   #corr_mat_long = corr_mat_long[corr_mat_long$Metadata %like% "liver_function",]
   
   pmat_long <- melt(pmat) %>% magrittr::set_colnames(c("Species","Metadata","Pvalue"))
   #pmat_long = pmat_long[pmat_long$Metadata %like% "liver_function",]
   
   
   ##### Significant
   final_mat_pvalue_1 <-
      left_join(pmat_long,Clusters, by = "Species") %>% 
      dplyr::mutate(Cluster = replace_na(Cluster %>% as.character(), "No_Cluster")) %>% 
      na.omit(.) %>% 
      filter(.,Pvalue < 0.05) %>% 
      dplyr::rename(Clinical = "Metadata") %>% 
      dplyr::select(Clinical,Species,Pvalue,Cluster) %>% 
      dplyr::filter(grepl('FLI.liver_function|ALT.liver_function|AST.liver_function|Triglyceride.lipid_profiles|GGT.liver_function|
                        GGT.liver_function|Liver_Fat.liver_function',Clinical))
   
   
   Input_Species_Clinical <- 
      table(final_mat_pvalue_1[,c(1,4)]) %>% 
      as.data.frame(.) %>% 
      .[which(.$Clinical %in% final_mat_pvalue_1$Clinical),]
   
   Cont_Matrix_Clinical_Species = Input_Species_Clinical %>% spread(., (Clinical), Freq) %>% tibble::column_to_rownames(., var = "Cluster")
   Total_1 = rowSums(Cont_Matrix_Clinical_Species)
   Cont_Matrix_Clinical_Species = cbind(Cont_Matrix_Clinical_Species,Total_1)
   
   Total_2 = colSums(Cont_Matrix_Clinical_Species)
   Cont_Matrix_Clinical_Species = rbind(Cont_Matrix_Clinical_Species,Total_2)
   rownames(Cont_Matrix_Clinical_Species)[nrow(Cont_Matrix_Clinical_Species)] = "Total_2"
   
   ### Not Significant
   final_mat_pvalue_1_notsignif <- 
      left_join(pmat_long,Clusters, by = "Species") %>% 
      dplyr::mutate(Cluster = replace_na(Cluster %>% as.character(), "No_Cluster"))%>% 
      na.omit(.) %>%
      filter(.,Pvalue > 0.05) %>%
      dplyr::rename(Clinical = "Metadata") %>% 
      dplyr::select(Clinical,Species,Pvalue,Cluster) %>% 
      dplyr::filter(grepl('FLI.liver_function|ALT.liver_function|AST.liver_function|Triglyceride.lipid_profiles|GGT.liver_function|
                        GGT.liver_function|Liver_Fat.liver_function',Clinical))
   
   
   
   Input_Species_Clinical_not_signif <-
      table(final_mat_pvalue_1_notsignif[,c(1,4)]) %>%
      as.data.frame(.) %>%
      .[which(.$Clinical %in% final_mat_pvalue_1$Clinical),]
   
   Cont_Matrix_Clinical_Species_not_signif = Input_Species_Clinical_not_signif %>% spread(., (Clinical), Freq) %>% tibble::column_to_rownames(., var = "Cluster")
   Total_1 = rowSums(Cont_Matrix_Clinical_Species_not_signif)
   Cont_Matrix_Clinical_Species_not_signif = cbind(Cont_Matrix_Clinical_Species_not_signif,Total_1)
   
   Total_2 = colSums(Cont_Matrix_Clinical_Species_not_signif)
   Cont_Matrix_Clinical_Species_not_signif = rbind(Cont_Matrix_Clinical_Species_not_signif,Total_2)
   rownames(Cont_Matrix_Clinical_Species_not_signif)[nrow(Cont_Matrix_Clinical_Species_not_signif)] = "Total_2"
   
   
   all_results_FT_Clinical = data.frame()
   for (i in 1:(ncol(Cont_Matrix_Clinical_Species)-1)) {
      Current_Clinial_signif = Cont_Matrix_Clinical_Species[,c(i,1)]
      Current_Clinial_not_signif = Cont_Matrix_Clinical_Species_not_signif[,c(i,1)]
      results_FT_Clinical = data.frame()
      
      for (j in 1:(nrow(Cont_Matrix_Clinical_Species)- 1)) {
         
         Current_Clinical = colnames(Cont_Matrix_Clinical_Species)[i]
         
         ### Significant
         current_signif = Current_Clinial_signif %>% .[c(j,nrow(.)),c(1)]
         Current_Clinical_Total_Significant_Correlations = final_mat_pvalue_1[final_mat_pvalue_1$Clinical == (Current_Clinical),]
         
         
         current_signif[2] =
            Current_Clinical_Total_Significant_Correlations %>%  
            nrow(.) - current_signif[1]
         
         
         #### Not Significant 
         current_not_signif = Current_Clinial_not_signif %>% .[c(j,nrow(.)),c(1)]
         Current_Clinical_Total_Not_Significant_Correlations = final_mat_pvalue_1_notsignif[final_mat_pvalue_1_notsignif$Clinical == (Current_Clinical),]
         
         
         
         current_not_signif[2] =
            Current_Clinical_Total_Not_Significant_Correlations %>% 
            nrow(.)- current_not_signif[1]
         
         
         current = cbind(current_signif,
                         current_not_signif)
         
         FT = fisher.test(current)
         results_FT_Clinical = rbind(results_FT_Clinical,
                                     data.frame(Pvalue = FT$p.value, 
                                                Cluster = rownames(Cont_Matrix_Clinical_Species)[j],
                                                Clinical = colnames(Cont_Matrix_Clinical_Species)[i],
                                                Number_significant_in_cluster_vs_total_number_significant = paste(current_signif[1],"/",(current[1,1]+current[2,1]),sep = ""),
                                                Enrichment_Ratio = (current[1,1]/(current[1,1]+current[1,2]))/((current[1,1]+current[2,1])/(sum(current)))
                                     )
         )
      }
      all_results_FT_Clinical = rbind(all_results_FT_Clinical,results_FT_Clinical)
   }
   
   colnames(all_results_FT_Clinical) = c("Pvalue","Cluster","Metadata","Number_significant_in_cluster_vs_total_number_significant","Enrichment_Ratio")
   all_results_FT_Clinical = all_results_FT_Clinical[,c(3,2,1,4,5)]
   
   Significant_Clinical <-
      all_results_FT_Clinical %>% 
      #all_results_FT_Clinical[which(all_results_FT_Clinical$Pvalue < 0.05),] %>% 
      arrange(Cluster)
   
   return(Significant_Clinical)
   
}


Full_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", row.names = 1) 


var_adjust <-
   Full_metadata[,c(
      'BMI.body_weight_variables', 
      'Age.demographic_variable', 
      'Gender.demographic_variable'
   )] %>% 
   data.frame()

Full_metadata <- 
   Full_metadata %>% 
   dplyr::select(.,-excluded)


full_features = bind_rows(NAFLD_L,
                          Control_L)
full_features[is.na(full_features)] <- 0

full_features<-
   full_features %>%
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] 


DGCA_All_Res_1 <- 
   Fisher_Procedure(N_megena_res,
                    adjusted  = TRUE) %>% 
   as_tibble() %>% 
   dplyr::filter(Pvalue < 0.05) 

DGCA_All_Res_Strict <- 
   Fisher_Procedure(N_megena_re_2,
                    adjusted  = TRUE) %>% 
   as_tibble() %>% 
   dplyr::filter(Pvalue < 0.05) 

DGCA_All_Res_ADJ <- 
   Fisher_Procedure(N_megena_res_ADJ,
                    adjusted  = TRUE) %>% 
   as_tibble() %>% 
   dplyr::filter(Pvalue < 0.05) 

DGCA_All_Res_ADJ_Strict <- 
   Fisher_Procedure(N_megena_re_2_ADJ,
                    adjusted  = TRUE) %>% 
   as_tibble() %>% 
   dplyr::filter(Pvalue < 0.05) 
