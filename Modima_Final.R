library(dplyr)
library(tidyr)
library(reshape)
library(tidyverse)
library(vegan)
library(ade4)
library(energy)
library(permute)
library(matrixStats)
source("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/MODIMA-master/modima.R")


load("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


DGCA_NAFLD_L_vs_Control_L <- readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")
DGCA_NAFLD_O_vs_Control_O <- readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")



### Loading the Data 
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(("Group"))) 
Control_L <- dplyr::select(Control_L, -one_of(("Group"))) 

NAFLD_O <- dplyr::select(NAFLD_O, -one_of(("Group"))) 
Control_O <- dplyr::select(Control_O, -one_of(("Group"))) 
##### NAFLD ##########

Obese_Species <-
   bind_rows(NAFLD_O, Control_O) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   .[order(rownames(.)),] 

Lean_Species <-
   bind_rows(NAFLD_L, Control_L) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   .[order(rownames(.)),] 


lean_coms_met <-
   readRDS( "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Lean_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

obese_coms_met <-
   readRDS( "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Obese_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

List_Signif_Metabolites <-
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/update_NAFLD_up_down.csv", sep = "\t") %>%
   as_tibble() %>%
   mutate(Direction = ifelse(Direction == "up",1,0)) %>%
   dplyr::rename(Species = "Sig_dif") %>%
   mutate(.,Species_Direction = ifelse(Direction == 1,paste(Species,"Up",sep = "_"),paste(Species,"Down",sep = "_"))) %>%
   mutate(Comparison = ifelse(Comparison == "O","NAFLD-O vs Control-O", "NAFLD-L vs Control-L")) %>%
   dplyr::left_join(read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/p_value_significant_metabolites_csv.csv", sep = ",") %>%
                       dplyr::rename(Species = "Metabolites") %>%
                       filter(Category == "NAFLD_L" | Category == "NAFLD_O" ) %>%
                       mutate(Category  = gsub("NAFLD_O","NAFLD-O vs Control-O",Category)) %>%
                       mutate(Category  = gsub("NAFLD_L","NAFLD-L vs Control-L",Category)) %>%
                       dplyr::rename(Comparison = "Category"), by = c("Species","Comparison"))

Lean_mets <- List_Signif_Metabolites %>%
   filter(Comparison =="NAFLD-L vs Control-L" )


Obese_mets <- List_Signif_Metabolites %>%
   filter(Comparison =="NAFLD-O vs Control-O" )


Metabolites_NAFLD_L <-
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
   column_to_rownames(var = "X") %>%
   t() %>%
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>%
   magrittr::add(1) %>%
   log2() %>%
   .[which(rownames(.) %in% rownames(Lean_Species)),] %>%
   as.data.frame() %>%
   .[order(rownames(.)),] %>%
   .[,which(colnames(.) %in% Lean_mets$Species)]



Metabolites_NAFLD_O <-
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
   column_to_rownames(var = "X") %>%
   t() %>%
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>%
   magrittr::add(1) %>%
   log2() %>%
   .[which(rownames(.) %in% rownames(Obese_Species)),] %>%
   as.data.frame() %>%
   .[order(rownames(.)),] %>%
   .[,which(colnames(.) %in% Obese_mets$Species)]

Obese_Species <-
   Obese_Species %>%
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_O)),]

Lean_Species <-
   Lean_Species %>%
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_L)),]



Primary_metadata <-  
   read.csv(file = "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(FLI.liver_function,
                 GGT.liver_function,
                 ALT.liver_function,
                 AST.liver_function,
                 Triglyceride.lipid_profiles,
                 Liver_Fat.liver_function) 



##3 Ordering 
Primary_metadata_obese <- 
   Primary_metadata[which(rownames(Primary_metadata) %in% rownames(Obese_Species)),] %>% 
   na.omit()
Obese_Species <-
   Obese_Species[which(rownames(Obese_Species) %in% rownames(Primary_metadata_obese)),]
Metabolites_NAFLD_O <-
   Metabolites_NAFLD_O[which(rownames(Metabolites_NAFLD_O) %in% rownames(Primary_metadata_obese)),]


Primary_metadata_obese = Primary_metadata_obese[order(rownames(Primary_metadata_obese)),]
Obese_Species = Obese_Species[order(rownames(Obese_Species)),]
Metabolites_NAFLD_O = Metabolites_NAFLD_O[order(rownames(Metabolites_NAFLD_O)),]



Primary_metadata_lean <- 
   Primary_metadata[which(rownames(Primary_metadata) %in% rownames(Lean_Species)),] %>% 
   na.omit()
Lean_Species <-
   Lean_Species[which(rownames(Lean_Species) %in% rownames(Primary_metadata_lean)),]
Metabolites_NAFLD_L <-
   Metabolites_NAFLD_L[which(rownames(Metabolites_NAFLD_L) %in% rownames(Primary_metadata_lean)),]

Primary_metadata_lean = Primary_metadata_lean[order(rownames(Primary_metadata_lean)),]
Lean_Species = Lean_Species[order(rownames(Lean_Species)),]
Metabolites_NAFLD_L = Metabolites_NAFLD_L[order(rownames(Metabolites_NAFLD_L)),]


Metabolites_NAFLD_L_N <- Metabolites_NAFLD_L[which(rownames(Metabolites_NAFLD_L) %in% rownames(NAFLD_L)),]
Metabolites_NAFLD_L_C <- Metabolites_NAFLD_L[which(rownames(Metabolites_NAFLD_L) %in% rownames(Control_L)),]

Metabolites_NAFLD_O_N <- Metabolites_NAFLD_O[which(rownames(Metabolites_NAFLD_O) %in% rownames(NAFLD_O)),]
Metabolites_NAFLD_O_C <- Metabolites_NAFLD_O[which(rownames(Metabolites_NAFLD_O) %in% rownames(Control_O)),]





### Separating Clusters
##Species

Obese_Species_c4 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_4") %>% 
                    dplyr::pull(Genes))
Obese_Species_c7 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_7") %>% 
                    dplyr::pull(Genes))



###Lean
Lean_Species_c2 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_2") %>% 
                    dplyr::pull(Genes))

Lean_Species_c5 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_5") %>% 
                    dplyr::pull(Genes))





modima_function <- function(independent,mediator,dependent) {
   
   # independent = Obese_Species_c7
   # mediator = Metabolites_NAFLD_O
   # dependent = Primary_metadata_obese
   
   independent <- independent %>%
      filter(!rowSums(. == 0) == ncol(.))
   
   mediator <- mediator[which(rownames(mediator) %in% rownames(independent)),] 
   dependent <- dependent[which(rownames(dependent) %in% rownames(independent)),] 
   
   
   independent <- vegan::vegdist(independent %>% as.matrix(), method="bray")
   mediator <- vegan::vegdist(mediator %>% as.matrix(), method="euclidean")
   dependent <- vegan::vegdist(dependent %>% as.matrix(), method="euclidean")
   
   
   set.seed(12345)
   modima_results = modima(exposure=independent, 
                           mediator=mediator, 
                           response=dependent, 
                           nrep=9999)
   
   
   
   set.seed(12345)
   ER_M_list <- pdcor.test (independent, dependent, mediator, R=9999)
   set.seed(12345)
   MR_E_list <- pdcor.test(mediator, dependent, independent, R=9999)
   set.seed(12345)
   beta_MR_E <- pdcor(mediator, dependent, independent)  %>% as.numeric()
   set.seed(12345)
   beta_ER <- bcdcor(independent,dependent)  %>% as.numeric()
   set.seed(12345)
   ER_list <- dcor.test(independent,dependent,R = 9999)
   set.seed(12345)
   beta_EM <- bcdcor(independent,mediator)  %>% as.numeric()
   set.seed(12345)
   EM_list <- dcor.test(independent,mediator,R = 9999)
   
   result_Pollution_Taxa_Benthos_list_all <- data.frame(
      Med_statistic = modima_results$statistic %>% as.numeric(),
      Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
      Med_P = modima_results$p.value,
      beta_EM = beta_EM,
      P_EM = EM_list$p.value,
      beta_MR_E = beta_MR_E,
      P_MR_E = MR_E_list$p.value,
      beta_ER = beta_ER,
      P_ER = ER_list$p.value,
      name = "Pollution_Taxa_Benthos_all"
   )
   
   results <- list(modima = modima_results,
                   values = result_Pollution_Taxa_Benthos_list_all)
   
   
   return(results)
   
}

C4 <- modima_function(Obese_Species_c4, Metabolites_NAFLD_O,Primary_metadata_obese)

C4$values 
# C4_s <- modima_function(Obese_Species_c4, Metabolites_NAFLD_O,
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))
# 
# C7 <- modima_function(Obese_Species_c7,Metabolites_NAFLD_O,Primary_metadata_obese)
# C7_s <- modima_function(Obese_Species_c7,Metabolites_NAFLD_O,
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))

C2 <- modima_function(Lean_Species_c2,Metabolites_NAFLD_L,Primary_metadata_lean)
C2$values
# C2_s <- modima_function(Lean_Species_c2,Metabolites_NAFLD_L,
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,Liver_Fat.liver_function))
# 
# C5 <- modima_function(Lean_Species_c5,Metabolites_NAFLD_L,Primary_metadata_lean)
# C5_s <- modima_function(Lean_Species_c5,Metabolites_NAFLD_L,
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,
#                                                          AST.liver_function,Liver_Fat.liver_function))


# 
# C4_u <- modima_function(Obese_Species_c4  , 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_obese)
# 
# C4_u_s <- modima_function(Obese_Species_c4  , 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))
# 
# 
# C4_d <- modima_function(Obese_Species_c4, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_obese)
# C4_d_s <- modima_function(Obese_Species_c4, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))
# 
# 
# 
# 
# C7_u <- modima_function(Obese_Species_c7, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_obese)
# C7_u_s <- modima_function(Obese_Species_c7, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))
# 
# C7_d <- modima_function(Obese_Species_c7, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_obese)
# C7_d_s <- modima_function(Obese_Species_c7, 
#                         Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_obese %>% select(FLI.liver_function,ALT.liver_function,AST.liver_function,
#                                                           GGT.liver_function,Liver_Fat.liver_function))
# 
# 
# 
# C2_u <- modima_function(Lean_Species_c2, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_lean)
# C2_u_s <- modima_function(Lean_Species_c2, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,Liver_Fat.liver_function))
# 
# 
# 
# C2_d <- modima_function(Lean_Species_c2, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_lean)
# C2_d_s <- modima_function(Lean_Species_c2, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,Liver_Fat.liver_function))
# 
# 
# C5_u <- modima_function(Lean_Species_c5, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_lean)
# C5_u_s <- modima_function(Lean_Species_c5, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,
#                                                          AST.liver_function,Liver_Fat.liver_function))
# 
# C5_d <- modima_function(Lean_Species_c5, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_lean)
# C5_d_s <- modima_function(Lean_Species_c5, 
#                         Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
#                         Primary_metadata_lean %>% select(FLI.liver_function,ALT.liver_function,
#                                                          AST.liver_function,Liver_Fat.liver_function))
# 
# 

Grouped_Modima <- list(C4 = C4,
                       C4_s = C4_s,
                       C4_u = C4_u,
                       C4_u_s = C4_u_s,
                       C4_d = C4_d,
                       C4_d_s = C4_d_s,
                       
                       C7 = C7,
                       C7_s = C7_s,
                       C7_u = C7_u,
                       C7_u_s = C7_u_s,
                       C7_d = C7_d,
                       C7_d_s = C7_d_s,
                       
                       C2 = C2,
                       C2_s = C2_s,
                       C2_u = C2_u,
                       C2_u_s = C2_u_s,
                       C2_d = C2_d,
                       C2_d_s = C2_d_s,
                       
                       C5 = C5,
                       C5_s = C5_s,
                       C5_u = C5_u,
                       C5_u_s = C5_u_s,
                       C5_d = C5_d,
                       C5_d_s = C5_d_s)


#### Modules - All DA metabolites - ALL Clinical Data
#### Significant 

Grouped_Modima$C4$values  #signif
Grouped_Modima$C4_s$modima  #signif

Grouped_Modima$C4_d$modima  #signif
Grouped_Modima$C4_d_s$modima #signif


Grouped_Modima$C7$modima #signif
Grouped_Modima$C7_s$modima #signif


Grouped_Modima$C7_d$modima #signif
Grouped_Modima$C7_d_s$modima #signif


Grouped_Modima$C2$values  #signif
Grouped_Modima$C2_s$modima  #signif


Grouped_Modima$C2_u_s$modima #signif




Grouped_Modima$C5$modima #signif
Grouped_Modima$C5_s$modima #signif


Grouped_Modima$C5_u_s$valuesmodima #signif



res <- rbind(Grouped_Modima$C4$values,
      Grouped_Modima$C4_s$values,
      
      Grouped_Modima$C4_d$values,
      Grouped_Modima$C4_d_s$values,
      
      Grouped_Modima$C7$values,
      Grouped_Modima$C7_s$values,
      
      Grouped_Modima$C7_d$values,
      Grouped_Modima$C7_d_s$values,
      
      Grouped_Modima$C2$values,
      Grouped_Modima$C2_s$values, 
      Grouped_Modima$C2_u_s$values, 
      
      Grouped_Modima$C5$values,
      Grouped_Modima$C5_s$values, 
      Grouped_Modima$C5_u_s$values)


res_2 <- data.frame(Analysis = c('C4',
                        'C4 - 5 Clin',
                        'C4 - Down Met',
                        'C4 - Down Met - 5 Clin',
                        'C7',
                        'C7 - 5 Clin',
                        'C7 - Down Met',
                        'C7 - Down Met - 5 Clin',
                        'C2',
                        'C2 - 3 Clin',
                        'C2 - Up Met - 3 Clin',
                        'C5',
                        'C5 - 3 Clin',
                        'C5 - Up Met - 3 Clin'),
           res)


write_rds(Grouped_Modima,
          "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Grouped_Groups.rds")


write.csv(res_2,
          "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Grouped_Groups_REs.csv")


modima_function_single_clin <- function(independent,mediator,dependent) {
   
   # independent = Obese_Species_c4 %>% select(all_of(SCCA_4 %>% filter(SCCA_Cor == "Positive") %>% pull(Features)))
   # mediator = Metabolites_NAFLD_O
   # dependent = Primary_metadata_obese
   
   
   independent <- independent %>%
      filter(!rowSums(. == 0) == ncol(.))
   
   mediator <- mediator[which(rownames(mediator) %in% rownames(independent)),] 
   dependent <- dependent[which(rownames(dependent) %in% rownames(independent)),] 
   
   
   
   
   mediation_rest <- data.frame()
   
   for (i in 1:ncol(dependent)) {
      print(colnames(dependent[i]))
      independent_c <- vegan::vegdist(independent %>% as.matrix(), method="bray")
      mediator_c <- vegan::vegdist(mediator %>% as.matrix(), method="euclidean")
      dependent_c <- vegan::vegdist(dependent[,i] %>% as.matrix(), method="euclidean")
      
      set.seed(12345)
      modima_results = modima(exposure=independent_c, 
                              mediator=mediator_c, 
                              response=dependent_c, 
                              nrep=9999)
      
      set.seed(12345)
      ER_M_list <- pdcor.test (independent_c, dependent_c, mediator_c, R=9999)
      set.seed(12345)
      MR_E_list <- pdcor.test(mediator_c, dependent_c, independent_c, R=9999)
      set.seed(12345)
      beta_MR_E <- pdcor(mediator_c, dependent_c, independent_c)  %>% as.numeric()
      set.seed(12345)
      beta_ER <- bcdcor(independent_c,dependent_c)  %>% as.numeric()
      set.seed(12345)
      ER_list <- dcor.test(independent_c,dependent_c,R = 9999)
      set.seed(12345)
      beta_EM <- bcdcor(independent_c,mediator_c)  %>% as.numeric()
      set.seed(12345)
      EM_list <- dcor.test(independent_c,mediator_c,R = 9999)
      
      
      result_Pollution_Taxa_Benthos_list_all <- data.frame(
         Clinical = colnames(dependent[i]),
         Perm_P = modima_results$p.value,
         Med_statistic = modima_results$statistic %>% as.numeric(),
         Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
         Med_P = modima_results$p.value,
         beta_EM = beta_EM,
         P_EM = EM_list$p.value,
         beta_MR_E = beta_MR_E,
         P_MR_E = MR_E_list$p.value,
         beta_ER = beta_ER,
         P_ER = ER_list$p.value
      )
      
      
      mediation_rest <- rbind(mediation_rest,
                              result_Pollution_Taxa_Benthos_list_all)
      
   }
   return(mediation_rest)
   
}


C4_C <- modima_function_single_clin(Obese_Species_c4, Metabolites_NAFLD_O,Primary_metadata_obese)
C7_C <- modima_function_single_clin(Obese_Species_c7,Metabolites_NAFLD_O,Primary_metadata_obese)
C2_C <- modima_function_single_clin(Lean_Species_c2,Metabolites_NAFLD_L,Primary_metadata_lean)
C5_C <- modima_function_single_clin(Lean_Species_c5,Metabolites_NAFLD_L,Primary_metadata_lean)



C4_C_u <- modima_function_single_clin(Obese_Species_c4  , 
                        Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
                        Primary_metadata_obese)

C4_C_d <- modima_function_single_clin(Obese_Species_c4, 
                        Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
                        Primary_metadata_obese)


C7_C_u <- modima_function_single_clin(Obese_Species_c7, 
                        Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 1) %>% pull(Species))),
                        Primary_metadata_obese)

C7_C_d <- modima_function_single_clin(Obese_Species_c7, 
                        Metabolites_NAFLD_O %>% select(all_of(Obese_mets %>% filter(Direction == 0) %>% pull(Species))),
                        Primary_metadata_obese)

C2_C_u <- modima_function_single_clin(Lean_Species_c2, 
                        Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
                        Primary_metadata_lean)

C2_C_d <- modima_function_single_clin(Lean_Species_c2, 
                        Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
                        Primary_metadata_lean)



C5_C_u <- modima_function_single_clin(Lean_Species_c5, 
                        Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 1) %>% pull(Species))),
                        Primary_metadata_lean)


C5_C_d <- modima_function_single_clin(Lean_Species_c5, 
                        Metabolites_NAFLD_L %>% select(all_of(Lean_mets %>% filter(Direction == 0) %>% pull(Species))),
                        Primary_metadata_lean)



Grouped_Modima_2 <- list(C4 = C4_C,
                       C4_u = C4_C_u,
                       C4_d = C4_C_d,
                       C7 = C7_C,
                       C7_u = C7_C_u,
                       C7_d = C7_C_d,
                       C2 = C2_C,
                       C2_u = C2_C_u,
                       C2_d = C2_C_d,
                       C5 = C5_C,
                       C5_u = C5_C_u,
                       C5_d = C5_C_d)


write_rds(Grouped_Modima_2,
          "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Grouped_Groups_single_clin.rds")



Grouped_Modima_2$C4 %>% filter(Perm_P < 0.05)

Grouped_Modima_2$C4_d %>% filter(Perm_P < 0.05)

Grouped_Modima_2$C7 %>% filter(Perm_P < 0.05)

Grouped_Modima_2$C7_d %>% filter(Perm_P < 0.05)

Grouped_Modima_2$C2 %>% filter(Perm_P < 0.05)
Grouped_Modima_2$C2_u %>% filter(Perm_P < 0.05)


Grouped_Modima_2$C5 %>% filter(Perm_P < 0.05)
Grouped_Modima_2$C5_u %>% filter(Perm_P < 0.05)


res_c <- rbind(Grouped_Modima_2$C4 %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C4_d %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C7 %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C7_d %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C2 %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C2_u %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C5 %>% filter(Perm_P < 0.05),
      Grouped_Modima_2$C5_u %>% filter(Perm_P < 0.05))

res_3 <- data.frame(Analysis = c(rep("C4",nrow(Grouped_Modima_2$C4 %>% filter(Perm_P < 0.05))),
                        rep("C4 - Down Met",nrow(Grouped_Modima_2$C4_d %>% filter(Perm_P < 0.05))),
                        rep("C7",nrow(Grouped_Modima_2$C7 %>% filter(Perm_P < 0.05))),
                        rep("C7 - Down Met",nrow(Grouped_Modima_2$C7_d %>% filter(Perm_P < 0.05))),
                        
                        rep("C2",nrow(Grouped_Modima_2$C2 %>% filter(Perm_P < 0.05))),
                        rep("C2 - Up Met",nrow(Grouped_Modima_2$C2_u %>% filter(Perm_P < 0.05))),
                        rep("C5",nrow(Grouped_Modima_2$C5 %>% filter(Perm_P < 0.05))),
                        rep("C5 - Up Met",nrow(Grouped_Modima_2$C5_u %>% filter(Perm_P < 0.05)))
                        
                        ),
           res_c)

write.csv(res_3,
          "Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Grouped_Groups_single_clin_REs.csv")

