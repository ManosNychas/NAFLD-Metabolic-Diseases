library(dplyr)
library(tidyr)
library(reshape)
library(tidyverse)
library(vegan)
library(ade4)
library(energy)
library(permute)
library(matrixStats)
source("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/MODIMA-master/modima.R")

### Loading the Data 
load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


DGCA_NAFLD_L_vs_Control_L <- readRDS("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")
DGCA_NAFLD_O_vs_Control_O <- readRDS("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")


### Data processing fo Species - Metabolites and Metadata for Obese and Lean


NAFLD_L <- dplyr::select(NAFLD_L, -one_of(("Group"))) 
Control_L <- dplyr::select(Control_L, -one_of(("Group"))) 

NAFLD_O <- dplyr::select(NAFLD_O, -one_of(("Group"))) 
Control_O <- dplyr::select(Control_O, -one_of(("Group"))) 


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
   readRDS( "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Lean_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

obese_coms_met <-
   readRDS( "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Obese_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

List_Signif_Metabolites <-
   read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/update_NAFLD_up_down.csv", sep = "\t") %>%
   as_tibble() %>%
   mutate(Direction = ifelse(Direction == "up",1,0)) %>%
   dplyr::rename(Species = "Sig_dif") %>%
   mutate(.,Species_Direction = ifelse(Direction == 1,paste(Species,"Up",sep = "_"),paste(Species,"Down",sep = "_"))) %>%
   mutate(Comparison = ifelse(Comparison == "O","NAFLD-O vs Control-O", "NAFLD-L vs Control-L")) %>%
   dplyr::left_join(read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/p_value_significant_metabolites_csv.csv", sep = ",") %>%
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
   read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
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
   read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
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

###Lean
Lean_Species_c2 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_2") %>% 
                    dplyr::pull(Genes))

# modima function + calculation of statistic and pvalues

modima_function <- function(independent,mediator,dependent) {
   

   
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

#Obese CLuster NAFLD- Metabolites - Metadata
C4 <- modima_function(Obese_Species_c4, Metabolites_NAFLD_O,Primary_metadata_obese)

C4$values 


#Lean CLuster NAFLD - Metabolites - Metadata

C2 <- modima_function(Lean_Species_c2,Metabolites_NAFLD_L,Primary_metadata_lean)
C2$values



