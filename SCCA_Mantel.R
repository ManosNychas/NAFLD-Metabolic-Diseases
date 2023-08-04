rm(list = ls(all = TRUE))
library(dplyr)
library(magrittr)
library(reshape)
library(metacoder)
library(tidyverse)
library(PMA)
library(ade4)
library(vegan)


### Loading the Data 
load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


DGCA_NAFLD_L_vs_Control_L <- readRDS("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")
DGCA_NAFLD_O_vs_Control_O <- readRDS("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")



#Formating Obese and Lean for all datasets : Species- Metadata

NAFLD_L <- dplyr::select(NAFLD_L, -one_of(("Group"))) 
Control_L <- dplyr::select(Control_L, -one_of(("Group"))) 

NAFLD_O <- dplyr::select(NAFLD_O, -one_of(("Group"))) 
Control_O <- dplyr::select(Control_O, -one_of(("Group"))) 

Obese_Species <-
   bind_rows(NAFLD_O, Control_O) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   mutate(Groups = c(rep("NAFLD",nrow(NAFLD_O)),
                     rep("Control",nrow(Control_O))))

Lean_Species <-
   bind_rows(NAFLD_L, Control_L) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   mutate(Groups = c(rep("NAFLD",nrow(NAFLD_L)),
                     rep("Control",nrow(Control_L))))

Primary_metadata <-  
   read.csv(file = "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(FLI.liver_function,
                 GGT.liver_function,
                 ALT.liver_function,
                 AST.liver_function,
                 Triglyceride.lipid_profiles,
                 Liver_Fat.liver_function) 

Adj_metadata <-  
   read.csv(file = "~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(Age.demographic_variable,
                 BMI.body_weight_variables,
                 Gender.demographic_variable) 


Primary_metadata_obese <- 
   Primary_metadata[which(rownames(Primary_metadata) %in% rownames(Obese_Species)),] %>% 
   na.omit()

Adj_metadata_obese <- 
   Adj_metadata[which(rownames(Adj_metadata) %in% rownames(Obese_Species)),] %>% 
   na.omit()

Obese_Species <-
   Obese_Species[which(rownames(Obese_Species) %in% rownames(Primary_metadata_obese)),]

Adj_metadata_obese <- 
   Adj_metadata[which(rownames(Adj_metadata) %in% rownames(Primary_metadata_obese)),] %>% 
   na.omit()

Primary_metadata_obese = Primary_metadata_obese[order(rownames(Primary_metadata_obese)),]
Adj_metadata_obese = Adj_metadata_obese[order(rownames(Adj_metadata_obese)),]
Obese_Species = Obese_Species[order(rownames(Obese_Species)),]




Primary_metadata_lean <- 
   Primary_metadata[which(rownames(Primary_metadata) %in% rownames(Lean_Species)),] %>% 
   na.omit()

Adj_metadata_lean <- 
   Adj_metadata[which(rownames(Adj_metadata) %in% rownames(Lean_Species)),] %>% 
   na.omit()

Lean_Species <-
   Lean_Species[which(rownames(Lean_Species) %in% rownames(Primary_metadata_lean)),]

Adj_metadata_lean <- 
   Adj_metadata[which(rownames(Adj_metadata) %in% rownames(Primary_metadata_lean)),] %>% 
   na.omit()


Primary_metadata_lean = Primary_metadata_lean[order(rownames(Primary_metadata_lean)),]
Adj_metadata_lean = Adj_metadata_lean[order(rownames(Adj_metadata_lean)),]
Lean_Species = Lean_Species[order(rownames(Lean_Species)),]


## Ordering 

Groups_Obese <-
   Obese_Species$Groups%>% 
   as.data.frame() %>% 
   set_rownames(rownames(Obese_Species))
Groups_Lean <- 
   Lean_Species$Groups %>% 
   as.data.frame() %>% 
   set_rownames(rownames(Lean_Species))


### Separating Clusters
##Obese
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






set.seed(123)

### Mantel and CCA for Obese
mantel.partial_CCA_Obese_c4 <- mantel.partial(vegdist(Obese_Species_c4 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c4 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c4 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              
                              method = "spearman", permutations = 9999, na.rm = TRUE)


CCA_Obese_c4 <- CCA.permute(Obese_Species_c4,Primary_metadata_obese,
                            typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                            penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                         0.3,0.3,0.3,0.3,0.3,0.3,
                                         0.4,0.4,0.4,0.4,0.4,0.4,
                                         0.5,0.5,0.5,0.5,0.5,0.5,
                                         0.6,0.6,0.6,0.6,0.6,0.6,
                                         .7,.7,.7,.7,.7,.7),
                            penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))


### Mantel and CCA for Lean
mantel.partial_Lean_c2<- mantel.partial(vegdist(Lean_Species_c2 %>% 
                                                  .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                        dist(Primary_metadata_lean
                                             [which(rownames(Primary_metadata_lean) %in% rownames(Lean_Species_c2 %>% 
                                                                                                    .[rowSums(.!= 0) >= ncol(.)*0.1,])),] 
                                        ), 
                                        
                                        
                                        dist(Adj_metadata_lean
                                             [which(rownames(Adj_metadata_lean) %in% rownames(Lean_Species_c2 %>% 
                                                                                                .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                        
                                        method = "spearman", permutations = 9999, na.rm = TRUE)



CCA_Lean_c2 <- CCA.permute(Lean_Species_c2,Primary_metadata_lean,
                           typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                           penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                        0.3,0.3,0.3,0.3,0.3,0.3,
                                        0.4,0.4,0.4,0.4,0.4,0.4,
                                        0.5,0.5,0.5,0.5,0.5,0.5,
                                        0.6,0.6,0.6,0.6,0.6,0.6,
                                        .7,.7,.7,.7,.7,.7),
                           penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))


## Extract Features Function with optmial Penalties : Species and Metadata that contribute to the CCA

extract_features <- function(D1,M1,CCA_Object) {
  

  DF <- 
    data.frame(Penalty_X = CCA_Object$penaltyxs,
               Penalty_Z = CCA_Object$penaltyzs,
               Z_Stat = CCA_Object$zstats,
               Cors = CCA_Object$cors,
               Pval = CCA_Object$pvals) %>% 
    filter(Pval < 0.05) 
  
  
  print(DF$Penalty_X[nrow(DF)])
  print(DF$Penalty_Z[nrow(DF)])
  
  res_SCCA <- CCA(x=D1,
                  z=M1,
                  typex="standard",
                  typez="standard",
                  penaltyx=DF$Penalty_X[nrow(DF)],
                  penaltyz=DF$Penalty_Z[nrow(DF)],
                  niter=10000,
                  trace=TRUE,K=3)
  Species <-
    res_SCCA$u %>% 
    `rownames<-`(colnames(D1)) %>% 
    as.data.frame(.) %>% 
    dplyr::select(V1) %>% 
    filter(V1 != 0) %>% 
    rownames_to_column() %>% 
    dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
    mutate(Type = rep("Species",nrow(.))) %>% 
    mutate(Cor_Type = ifelse(SCCA_Cor > 0, "Positive","Negative"))
  
  Metabolites <-
    res_SCCA$v %>% 
    `rownames<-`(colnames(M1)) %>% 
    as.data.frame(.) %>% 
    dplyr::select(V1) %>% 
    filter(V1 != 0) %>% 
    rownames_to_column() %>% 
    dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
    mutate(Type = rep("Clinical Parameters",nrow(.))) %>% 
    mutate(Cor_Type = ifelse(SCCA_Cor > 0, "Positive","Negative"))
  
  res <- list(Species = Species,
              Clinical_Data = Metabolites)
  
  return(res)
  
}




extract_Obese_4_0.1 <-extract_features(Obese_Species_c4,Primary_metadata_obese,
                                       CCA_Obese_c4)

extract_Lean_2_0.1 <-extract_features(Lean_Species_c2,Primary_metadata_lean,
                                      CCA_Lean_c2)

t1 <- rbind(extract_Obese_4_0.1$Species,
            extract_Obese_4_0.1$Clinical_Data)


t2 <- rbind(extract_Lean_2_0.1$Species,
            extract_Lean_2_0.1$Clinical_Data)

write.csv(t1,"Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Supplementary_Table_6_obese.csv")
write.csv(t2,"Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Supplementary_Table_6_lean.csv")


