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

## Functions --------------------------------------------------------------------
func_main <- function(df,a,b,c)
{
   

   
   # df <-  ALl_mat_Obese
   # a <- "Bifidobacterium_adolescentis"
   # b <- "L-Methionine"
   # c <- "FLI.liver_function"
   data_df <- func_df(df,a,b,c)
   func_step1_list <- func_step1(data_df)
   func_step2_list <- func_step2(data_df)
   if (func_step2_list$step2_re$`Pr(>|t|)` <= 0.05) {
      func_step3_list <- func_step3(data_df)
      Step3_P <- func_step3_list$step3_re %>% filter(terms %in% "med") %>% .$`Pr(>|t|)`
      Step3_beta <- func_step3_list$step3_re %>% filter(terms %in% "iv") %>% .$Estimate %>% abs()
      step1_beta <- func_step1_list$step1_re$Estimate %>% abs()
      Delta_beta <- step1_beta - Step3_beta
      if (Step3_P <= 0.05 && Delta_beta > 0) {
         func_step4_df <- func_step4(func_step2_list$fit.mediator,func_step3_list$fit.dv)
         re_df <- data.frame(
            iv = a,
            med = b,
            dv = c,
            Beta_step1_iv = func_step1_list$step1_re$Estimate,
            Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
            Beta_step2_iv = func_step2_list$step2_re$Estimate,
            Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
            Beta_step3_iv = func_step3_list$step3_re %>%
               filter(terms %in% "iv") %>% .$Estimate,
            Beta_step3_iv_p = func_step3_list$step3_re %>%
               filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
            Beta_step3_med = func_step3_list$step3_re %>%
               filter(terms %in% "med") %>% .$Estimate,
            Beta_step3_med_p = func_step3_list$step3_re %>%
               filter(terms %in% "med") %>% .$`Pr(>|t|)`,
            ACME = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ACME") %>%
               .$Estimate,
            ADE = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ADE") %>%
               .$Estimate,
            ACME_lower_95_CI = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ACME") %>%
               .$`95% CI Lower`,
            ACME_upper_95_CI = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ACME") %>%
               .$`95% CI Upper`,
            ACME_P_value = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ACME") %>%
               .$`p-value`,
            ADE_P_value = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "ADE") %>%
               .$`p-value`,
            Percentage_Medi_P_value = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "Prop. Mediated") %>%
               .$`p-value`,
            Percentage_Medi = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "Prop. Mediated") %>%
               .$Estimate,
            Percentage_Medi_lower_95_CI = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "Prop. Mediated") %>%
               .$`95% CI Lower`,
            Percentage_Medi_upper_95_CI = func_step4_df %>% 
               as.data.frame() %>%
               rownames_to_column() %>% 
               filter(rowname %in% "Prop. Mediated") %>%
               .$`95% CI Upper`
         )
      }else
      {
         re_df <- data.frame(
            iv = a,
            med = b,
            dv = c,
            Beta_step1_iv = func_step1_list$step1_re$Estimate,
            Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
            Beta_step2_iv = func_step2_list$step2_re$Estimate,
            Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
            Beta_step3_iv = func_step3_list$step3_re %>%
               filter(terms %in% "iv") %>% .$Estimate,
            Beta_step3_iv_p = func_step3_list$step3_re %>%
               filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
            Beta_step3_med = func_step3_list$step3_re %>%
               filter(terms %in% "med") %>% .$Estimate,
            Beta_step3_med_p = func_step3_list$step3_re %>%
               filter(terms %in% "med") %>% .$`Pr(>|t|)`,
            ACME = "",
            ADE = "",
            ACME_lower_95_CI = "",
            ACME_upper_95_CI = "",
            ACME_P_value = "",
            ADE_P_value = "",
            Percentage_Medi_P_value = "",
            Percentage_Medi = "",
            Percentage_Medi_lower_95_CI = "",
            Percentage_Medi_upper_95_CI = ""
         )
      }
   }else{
      func_step3_list <- func_step3(data_df)
      re_df <- data.frame(
         iv = a,
         med = b,
         dv = c,
         Beta_step1_iv = func_step1_list$step1_re$Estimate,
         Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
         Beta_step2_iv = func_step2_list$step2_re$Estimate,
         Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
         Beta_step3_iv = func_step3_list$step3_re %>%
            filter(terms %in% "iv") %>% .$Estimate,
         Beta_step3_iv_p = func_step3_list$step3_re %>%
            filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
         Beta_step3_med = func_step3_list$step3_re %>%
            filter(terms %in% "med") %>% .$Estimate,
         Beta_step3_med_p = func_step3_list$step3_re %>%
            filter(terms %in% "med") %>% .$`Pr(>|t|)`,
         ACME = "",
         ADE = "",
         ACME_lower_95_CI = "",
         ACME_upper_95_CI = "",
         ACME_P_value = "",
         ADE_P_value = "",
         Percentage_Medi_P_value = "",
         Percentage_Medi = "",
         Percentage_Medi_lower_95_CI = "",
         Percentage_Medi_upper_95_CI = ""
      )
   }
   return(re_df)
}

func_df <- function(df,a,b,c)
{
   # df <-  ALl_mat
   # a <- "Nitrate Nitrogen"
   # b <- "Nitrospina gracilis"
   # c <- "Polynoidae HK03"
   data_df <- df %>%
      dplyr::select(a,b,c) %>%
      `colnames<-`(c("iv","med","dv"))
   return(data_df)
}

func_step1 <- function(df)
{
   #Step 1: The total effect
   fit.totaleffect=lm(dv~iv,df)
   #summary(fit.totaleffect)
   step1_re <- as.data.frame(coef(summary(fit.totaleffect))) %>%
      rownames_to_column("terms") %>%
      filter(terms %in% "iv") %>%
      mutate(model = "step1")
   return_list <- list(
      fit.totaleffect = fit.totaleffect,
      step1_re = step1_re
   )
   return(return_list)
}

func_step2 <- function(df)
{
   #Step 2: The effect of the IV onto the mediator
   fit.mediator=lm(med~iv,df)
   #summary(fit.mediator)
   step2_re <- as.data.frame(coef(summary(fit.mediator))) %>%
      rownames_to_column("terms") %>%
      filter(terms %in% "iv") %>%
      mutate(model = "step2")
   return_list <- list(
      fit.mediator = fit.mediator,
      step2_re = step2_re
   )
   return(return_list)
}

func_step3 <- function(df)
{
   #Step 3: The effect of the mediator on the dependent variable
   fit.dv=lm(dv~iv+med,df)
   #summary(fit.dv)
   step3_re <- as.data.frame(coef(summary(fit.dv))) %>%
      rownames_to_column("terms") %>%
      filter(terms %in% c("iv","med")) %>%
      mutate(model = "step3")
   return_list <- list(
      fit.dv = fit.dv,
      step3_re = step3_re
   )
   return(return_list)
}

func_step4 <- function(fit.mediator,fit.dv)
{
   #Step 4: Causal Mediation Analysis
   results = mediate(fit.mediator, fit.dv, treat='iv', mediator='med', boot=T)
   #summary(results)
   summary_tb <- extract_mediation_summary(summary(results))
   return(summary_tb)
}

## Sub functions ---------------------------------------------------------------
extract_mediation_summary <- function (x) { 
   
   clp <- 100 * x$conf.level
   isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                     (inherits(x$model.y, "glm") && x$model.y$family$family == 
                         "gaussian" && x$model.y$family$link == "identity") || 
                     (inherits(x$model.y, "survreg") && x$model.y$dist == 
                         "gaussian"))
   
   printone <- !x$INT && isLinear.y
   
   if (printone) {
      
      smat <- c(x$d1, x$d1.ci, x$d1.p)
      smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
      smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
      smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
      
      rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
      
   } else {
      smat <- c(x$d0, x$d0.ci, x$d0.p)
      smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
      smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
      smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
      smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
      smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
      smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
      smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
      smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
      smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
      
      rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                          "ADE (control)", "ADE (treated)", "Total Effect", 
                          "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                          "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
      
   }
   
   colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                       paste(clp, "% CI Upper", sep = ""), "p-value")
   smat
   
}



load("~/Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


DGCA_NAFLD_L_vs_Control_L <- readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")
DGCA_NAFLD_O_vs_Control_O <- readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")



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

Lean_Species %>% dim()

lean_coms_met <-
   readRDS( "Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Lean_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

obese_coms_met <-
   readRDS( "Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Omics_Combination/Obese_Metabolites_metadata_Cors.rds") %>%
   as.data.frame()

List_Signif_Metabolites <-
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Metabolomics/update_NAFLD_up_down.csv", sep = "\t") %>%
   as_tibble() %>%
   mutate(Direction = ifelse(Direction == "up",1,0)) %>%
   dplyr::rename(Species = "Sig_dif") %>%
   mutate(.,Species_Direction = ifelse(Direction == 1,paste(Species,"Up",sep = "_"),paste(Species,"Down",sep = "_"))) %>%
   mutate(Comparison = ifelse(Comparison == "O","NAFLD-O vs Control-O", "NAFLD-L vs Control-L")) %>%
   dplyr::left_join(read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Metabolomics/p_value_significant_metabolites_csv.csv", sep = ",") %>%
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
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
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
   read.csv("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>%
   column_to_rownames(var = "X") %>%
   t() %>%
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>%
   magrittr::add(1) %>%
   log2() %>%
   .[which(rownames(.) %in% rownames(Obese_Species)),] %>%
   as.data.frame() %>%
   .[order(rownames(.)),] %>%
   .[,which(colnames(.) %in% Obese_mets$Species)]


Metabolites_NAFLD_O_N <-Metabolites_NAFLD_O[which(rownames(Metabolites_NAFLD_O) %in% rownames(NAFLD_O)),]
Metabolites_NAFLD_O_C <-Metabolites_NAFLD_O[which(rownames(Metabolites_NAFLD_O) %in% rownames(Control_O)),]


Metabolites_NAFLD_L_N <-Metabolites_NAFLD_L[which(rownames(Metabolites_NAFLD_L) %in% rownames(NAFLD_L)),]
Metabolites_NAFLD_L_C <-Metabolites_NAFLD_L[which(rownames(Metabolites_NAFLD_L) %in% rownames(Control_L)),]


fs <- function(m1,m2) {
  ratio <- m2 / m1
  
  # take the logarithm (base 2) of the ratio
  fold_change <- log2(ratio)
  print(fold_change)
  return(fold_change)
}

fs(mean(Metabolites_NAFLD_L_C$`Guanosine 3  phosphate C10H12N5O8P`),
   mean(Metabolites_NAFLD_L_N$`Guanosine 3  phosphate C10H12N5O8P`))

fs(mean(Metabolites_NAFLD_L_C$`2',3'-Cyclic UMP`),
   mean(Metabolites_NAFLD_L_N$`2',3'-Cyclic UMP`))

fs(mean(Metabolites_NAFLD_L_C$`D-Alanyl-D-alanine`),
   mean(Metabolites_NAFLD_L_N$`D-Alanyl-D-alanine`))



Obese_Species <-
   Obese_Species %>%
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_O)),]

Lean_Species <-
   Lean_Species %>%
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_L)),]



Primary_metadata <-  
   read.csv(file = "~/Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
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





library(foreach)
library(doParallel)
library(mediation)
registerDoParallel(22)
ALl_mat_Obese <- cbind(Metabolites_NAFLD_O,Primary_metadata_obese,Obese_Species)
ALl_mat_Lean <- cbind(Metabolites_NAFLD_L,Primary_metadata_lean,Lean_Species)



Obese_Results <- foreach (i = 1:dim(Obese_Species)[2],.combine=rbind) %dopar% {
   #i <- 3
   set.seed(123)
   #print(i)
   a = colnames(Obese_Species)[i]
   mediation_df <- data.frame()
   for (j in 1:length(Metabolites_NAFLD_O)) {
      #j <- 1
      for (k in 1:length(Primary_metadata_obese)) {
         #k <- 17
         b = colnames(Metabolites_NAFLD_O)[j]
         c = colnames(Primary_metadata_obese)[k]
         #str = paste0(a,"_",b,"_",c)
         #print(str)
         mediation_tmp <- func_main(ALl_mat_Obese,a,b,c)
         mediation_df <- rbind(mediation_df,mediation_tmp)
      }
   }
   mediation_df
}
Lean_Results <- foreach (i = 1:dim(Lean_Species)[2],.combine=rbind) %dopar% {
   #i <- 3
   set.seed(i)
   #print(i)
   a = colnames(Lean_Species)[i]
   mediation_df <- data.frame()
   for (j in 1:length(Metabolites_NAFLD_L)) {
      #j <- 1
      for (k in 1:length(Primary_metadata_lean)) {
         #k <- 17
         b = colnames(Metabolites_NAFLD_L)[j]
         c = colnames(Primary_metadata_lean)[k]
         #str = paste0(a,"_",b,"_",c)
         #print(str)
         mediation_tmp <- func_main(ALl_mat_Lean,a,b,c)
         mediation_df <- rbind(mediation_df,mediation_tmp)
      }
   }
   mediation_df
}


# 
# saveRDS(Obese_Results,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Singled_Obese.rds")
# saveRDS(Lean_Results,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Singled_Lean.rds")
# 
# 



library(compositions)
library(magrittr)
library(vegan)
library(tsensembler)
library(ape)
library(pracma)
library(phyloseq)


distance_matrices <- function(Type,Network,Cluster) {
  set.seed(123)
  
   #TREE <- phyloseq::read_tree("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
  Type = Obese_Species
  Network = DGCA_NAFLD_O_vs_Control_O
  Cluster = "c1_4"
  
  Type <- 
    Type %>% 
    dplyr::select(all_of(Network$NAFLD_Megena_Res$modules %>% 
                           filter(Modules == Cluster) %>% 
                           pull(Genes)))
  
  # EU <- Type %>%  vegdist(method = "euclidean") 
  B <- Type %>%  vegdist(method = "bray") 
  # U <- metaphlanToPhyloseq(t(Type)) %>% 
  #   merge_phyloseq(., TREE = TREE) %>% 
  #   phyloseq::distance(., method = "wunifrac")
  
  ### MDS
    MDS_B <- 
    B %>%  metaMDS(k = 2, trymax = 50) %>% extract("points") %>% 
      as.data.frame() %>% dplyr::select(points.MDS1)
  
  #   MDS_U <- U %>% 
  #   metaMDS(k = 2, trymax = 50) %>%  extract("points") %>%
  #   as.data.frame() %>% dplyr::select(points.MDS1)
  # 
  # MDS_E <- EU %>% 
  #   metaMDS(k = 2, trymax = 50) %>% extract("points") %>%
  #   as.data.frame() %>% dplyr::select(points.MDS1)
  
  
  PCOA_B <- B  %>% 
    pcoa() %>%  extract("vectors") %>%
    as.data.frame() %>%  dplyr::select(vectors.Axis.1)
  
  # PCOA_U <- U  %>% 
  #   pcoa() %>%   extract("vectors") %>%
  #   as.data.frame() %>%  dplyr::select(vectors.Axis.1)
  # 
  # 
  # PCOA_E <- EU  %>% 
  #   pcoa() %>% extract("vectors") %>%
  #   as.data.frame() %>%  dplyr::select(vectors.Axis.1)
  
  
  PCA_B <- B %>% 
    prcomp() %>%  extract("x") %>%
    as.data.frame() %>%  dplyr::select(x.PC1)
  
  # PCA_E <- EU %>% 
  #   prcomp() %>% extract("x") %>%
  #   as.data.frame() %>% dplyr::select(x.PC1)
  # 
  # 
  # PCA_U <- U %>% 
  #   prcomp() %>%  extract("x") %>%
  #   as.data.frame() %>% dplyr::select(x.PC1)
  # 
  # PCA <- Type %>% 
  #   prcomp() %>%  extract("x") %>%
  #   as.data.frame() %>% dplyr::select(x.PC1)
  
  
  Results <- data.frame(MDS_B = MDS_B$points.MDS1,
                        # MDS_U = MDS_U$points.MDS1,
                        # MDS_E = MDS_E$points.MDS1,
                        
                        PCOA_B = PCOA_B$vectors.Axis.1)
                        # PCOA_U = PCOA_U$vectors.Axis.1,
                        # PCOA_E = PCOA_E$vectors.Axis.1,
                        
                        # PCA_B = PCA_B$x.PC1,
                        # PCA_U = PCA_U$x.PC1,
                        # PCA_E = PCA_E$x.PC1,
                        # PCA = PCA$x.PC1)
  
  return(Results)
  
}




Distance_C4 <- distance_matrices(Obese_Species,DGCA_NAFLD_O_vs_Control_O,Cluster = "c1_4")
Distance_C7 <- distance_matrices(Obese_Species,DGCA_NAFLD_O_vs_Control_O,Cluster = "c1_7")
Distance_C2 <- distance_matrices(Lean_Species,DGCA_NAFLD_L_vs_Control_L,Cluster = "c1_2")
Distance_C5 <- distance_matrices(Lean_Species,DGCA_NAFLD_L_vs_Control_L,Cluster = "c1_5")


Sum_Obese_4 <- cbind(Metabolites_NAFLD_O,Primary_metadata_obese,Distance_C4)
Sum_Obese_7 <- cbind(Metabolites_NAFLD_O,Primary_metadata_obese,Distance_C7)
Sum_Lean_2 <- cbind(Metabolites_NAFLD_L,Primary_metadata_lean,Distance_C2)

Sum_Lean_5 <- cbind(Metabolites_NAFLD_L,Primary_metadata_lean,Distance_C5)



the_sum_mediation <-function(Species,Metabolites,Metadata,SumMat) {

  
  set.seed(123)
  # Species = Distance_C4
  # Metabolites = Metabolites_NAFLD_O
  # Metadata =  Primary_metadata_obese
  # SumMat = Sum_Obese_4
  
  
  
  
  the_super_function <- function(Species_2,Metabolites,Metadata,SumMat) {
   
     a = colnames(Species_2)
    mediation_df <- data.frame()
    set.seed(123)
    for (j in 1:length(Metabolites)) {
      #j <- 1
      for (k in 1:length(Metadata)) {
        #k <- 17
        b = colnames(Metabolites)[j]
        c = colnames(Metadata)[k]
        #str = paste0(a,"_",b,"_",c)
        #print(str)
        mediation_tmp <- func_main(SumMat,a,b,c)
        mediation_df <- rbind(mediation_df,mediation_tmp)
      }
    }
    mediation_df <- 
      mediation_df %>% 
      filter(Beta_step2_iv_p <= 0.05) %>% 
      filter(Beta_step3_med_p <= 0.05) 
  }

  MDS_B <- the_super_function(Species[1],Metabolites,Metadata,SumMat)
   # MDS_U <- the_super_function(Species[2],Metabolites,Metadata,SumMat)

  PCOA_B <- the_super_function(Species[2],Metabolites,Metadata,SumMat)
  # PCOA_U <- the_super_function(Species[4],Metabolites,Metadata,SumMat)
  # 
  # PCA_B <- the_super_function(Species[5],Metabolites,Metadata,SumMat)
  # PCA_U <- the_super_function(Species[6],Metabolites,Metadata,SumMat)
  # 
  # PCA <- the_super_function(Species[7],Metabolites,Metadata,SumMat)
  
  The_List <- list(MDS_B = MDS_B,
                   # MDS_U = MDS_U,

                   PCOA_B = PCOA_B)
                   # PCOA_U = PCOA_U,

                   # PCA_B = PCA_B,
                   # PCA_U = PCA_U,
                   # 
                   # PCA = PCA)
  
  return(The_List)
  
  
}

MED_Obese_4 <- the_sum_mediation(Distance_C4,Metabolites_NAFLD_O,Primary_metadata_obese,Sum_Obese_4)
MED_Obese_7 <- the_sum_mediation(Distance_C7,Metabolites_NAFLD_O,Primary_metadata_obese,Sum_Obese_7)
MED_Lean_2 <- the_sum_mediation(Distance_C2,Metabolites_NAFLD_L,Primary_metadata_lean,Sum_Lean_2)
MED_Lean_5 <- the_sum_mediation(Distance_C5,Metabolites_NAFLD_L,Primary_metadata_lean,Sum_Lean_5)


MED_Analysis <- list(MED_Obese_4 =MED_Obese_4,
                     MED_Obese_7 =MED_Obese_7,
                     MED_Lean_2 =MED_Lean_2,
                     MED_Lean_5 =MED_Lean_5)





#write_rds(MED_Analysis,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_Mediations.rds")


#MED_Analysis<-readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_Mediations.rds")

the_sum_mediation <-function(Species,Metabolites,Metadata) {
  
  
  set.seed(123)
  # Species = Distance_C4
  # Metabolites = Metabolites_NAFLD_O
  # Metadata =  Primary_metadata_obese
  # SumMat = Sum_Lean_2
  
  Metabolites <- Metabolites %>% clr() %>% as.data.frame()
  # Metadata <- Metadata %>% clr() %>% as.data.frame()
  
  

  SumMat <- cbind(Metabolites,Metadata,Species)
  
  the_super_function <- function(Species_2,Metabolites,Metadata,SumMat) {
    
    a = colnames(Species_2)
    mediation_df <- data.frame()
    set.seed(123)
    for (j in 1:length(Metabolites)) {
      #j <- 1
      for (k in 1:length(Metadata)) {
        #k <- 17
        b = colnames(Metabolites)[j]
        c = colnames(Metadata)[k]
        #str = paste0(a,"_",b,"_",c)
        #print(str)
        mediation_tmp <- func_main(SumMat,a,b,c)
        mediation_df <- rbind(mediation_df,mediation_tmp)
      }
    }
    mediation_df <- 
      mediation_df %>% 
      filter(Beta_step2_iv_p <= 0.05) %>% 
      filter(Beta_step3_med_p <= 0.05) 
  }
  
  MDS_B <- the_super_function(Species[1],Metabolites,Metadata,SumMat)
  MDS_U <- the_super_function(Species[2],Metabolites,Metadata,SumMat)
  MDS_E <- the_super_function(Species[3],Metabolites,Metadata,SumMat)
  PCOA_B <- the_super_function(Species[4],Metabolites,Metadata,SumMat)
  PCOA_U <- the_super_function(Species[5],Metabolites,Metadata,SumMat)
  PCOA_E <- the_super_function(Species[6],Metabolites,Metadata,SumMat)
  PCA_B <- the_super_function(Species[7],Metabolites,Metadata,SumMat)
  PCA_U <- the_super_function(Species[8],Metabolites,Metadata,SumMat)
  PCA_E <- the_super_function(Species[9],Metabolites,Metadata,SumMat)
  PCA <- the_super_function(Species[10],Metabolites,Metadata,SumMat)
  
  The_List <- list(MDS_B = MDS_B,
                   MDS_U = MDS_U,
                   MDS_E = MDS_E,
                   PCOA_B = PCOA_B,
                   PCOA_U = PCOA_U,
                   PCOA_E = PCOA_E,
                   PCA_B = PCA_B,
                   PCA_U = PCA_U,
                   PCA_E = PCA_E,
                   PCA = PCA)
  
  return(The_List)
  
  
}

MED_Obese_4 <- the_sum_mediation(Distance_C4,Metabolites_NAFLD_O,Primary_metadata_obese)
MED_Obese_7 <- the_sum_mediation(Distance_C7,Metabolites_NAFLD_O,Primary_metadata_obese)
MED_Lean_2 <- the_sum_mediation(Distance_C2,Metabolites_NAFLD_L,Primary_metadata_lean)
MED_Lean_5 <- the_sum_mediation(Distance_C5,Metabolites_NAFLD_L,Primary_metadata_lean)


# MED_Obese_4_CLR <- the_sum_mediation(Distance_C4,Metabolites_NAFLD_O,Primary_metadata_obese)
# MED_Obese_7_CLR <- the_sum_mediation(Distance_C7,Metabolites_NAFLD_O,Primary_metadata_obese)
# MED_Lean_2_CLR <- the_sum_mediation(Distance_C2,Metabolites_NAFLD_L,Primary_metadata_lean)
# MED_Lean_5_CLR <- the_sum_mediation(Distance_C5,Metabolites_NAFLD_L,Primary_metadata_lean)
# 
# 
# MED_Obese_4_CLR2 <- the_sum_mediation(Distance_C4,Metabolites_NAFLD_O,Primary_metadata_obese)
# MED_Obese_7_CLR2 <- the_sum_mediation(Distance_C7,Metabolites_NAFLD_O,Primary_metadata_obese)
# MED_Lean_2_CLR2 <- the_sum_mediation(Distance_C2,Metabolites_NAFLD_L,Primary_metadata_lean)
# MED_Lean_5_CLR2 <- the_sum_mediation(Distance_C5,Metabolites_NAFLD_L,Primary_metadata_lean)


MED_Analysis <- list(MED_Obese_4 =MED_Obese_4,
                     MED_Obese_7 =MED_Obese_7,
                     MED_Lean_2 =MED_Lean_2,
                     MED_Lean_5 =MED_Lean_5)


df <- data.frame(PCA = Distance_C2$PCA,
                 UMP = Metabolites_NAFLD_L$`2',3'-Cyclic UMP`,
                 FLI = Primary_metadata_lean$FLI.liver_function)

fit.mediator=lm(UMP~PCA,df)
fit.dv=lm(FLI~PCA+UMP,df)
results = mediate(fit.mediator, fit.dv, treat='PCA', mediator='UMP', boot=T)

summary(results)



MED_Analysis$MED_Obese_4
MED_Analysis$MED_Obese_7$MDS_B
MED_Analysis$MED_Lean_2$PCA_B
MED_Analysis$MED_Lean_5$PCA_B


# 
# MED_Analysis_CLR <- list(MED_Obese_4 =MED_Obese_4_CLR,
#                      MED_Obese_7 =MED_Obese_7_CLR,
#                      MED_Lean_2 =MED_Lean_2_CLR,
#                      MED_Lean_5 =MED_Lean_5_CLR)
# 
# MED_Analysis_CLR2 <- list(MED_Obese_4 =MED_Obese_4_CLR2,
#                      MED_Obese_7 =MED_Obese_7_CLR2,
#                      MED_Lean_2 =MED_Lean_2_CLR2,
#                      MED_Lean_5 =MED_Lean_5_CLR2)


MED_Analysis <-readRDS("Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_Mediations.rds")

write_rds(MED_Analysis,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_Mediations.rds")
write_rds(MED_Analysis_CLR,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_MediationsCLR.rds")
write_rds(MED_Analysis_CLR2,"Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Distance_MediationsCLR2.rds")

MED_Analysis$MED_Obese_4$MDS_B




