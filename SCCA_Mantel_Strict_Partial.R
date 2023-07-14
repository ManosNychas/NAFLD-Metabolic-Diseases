

rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
#library(compositions)
library(ancom)
library(dplyr)
library(tidyr)
library(RANN)
library(magrittr)
library(reshape)
library(plyr)
library(data.table)
library(metacoder)
library(readxl)
library(mice)
library(ggplot2)
library(textshape)
library(ppcor)
library(tidyverse)
load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


DGCA_NAFLD_L_vs_Control_L <- readRDS("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_L_vs_NAFLD_L.rds")
DGCA_NAFLD_O_vs_Control_O <- readRDS("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Strict_Prev_All_Non_Strict_DGCA_Control_O_vs_NAFLD_O.rds")

DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$summary$module.table
DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$summary$module.table


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
   mutate(Groups = c(rep("NAFLD",nrow(NAFLD_O)),
                     rep("Control",nrow(Control_O))))

Lean_Species <-
   bind_rows(NAFLD_L, Control_L) %>%
   replace(is.na(.), 0) %>% 
   .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   mutate(Groups = c(rep("NAFLD",nrow(NAFLD_L)),
                     rep("Control",nrow(Control_L))))

Primary_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(FLI.liver_function,
                 GGT.liver_function,
                 ALT.liver_function,
                 AST.liver_function,
                 Triglyceride.lipid_profiles,
                 Liver_Fat.liver_function) 



Adj_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
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


##3 Ordering 


Groups_Obese <-
   Obese_Species$Groups%>% 
   as.data.frame() %>% 
   set_rownames(rownames(Obese_Species))
Groups_Lean <- 
   Lean_Species$Groups %>% 
   as.data.frame() %>% 
   set_rownames(rownames(Lean_Species))


### Separating Clusters
##Species
Obese_Species_c2 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_2") %>% 
                    dplyr::pull(Genes))


Obese_Species_c3 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_3") %>% 
                    dplyr::pull(Genes))



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




Obese_Species_c10 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_10") %>% 
                    dplyr::pull(Genes))

Obese_Species_c11 <-
   Obese_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_O_vs_Control_O$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_11") %>% 
                    dplyr::pull(Genes))

###Lean

Lean_Species_c2 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_2") %>% 
                    dplyr::pull(Genes))

Lean_Species_c3 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_3") %>% 
                    dplyr::pull(Genes))
Lean_Species_c4 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_4") %>% 
                    dplyr::pull(Genes))


Lean_Species_c5 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_5") %>% 
                    dplyr::pull(Genes))



Lean_Species_c6 <-
   Lean_Species %>% 
   as.data.frame() %>% 
   dplyr::select(DGCA_NAFLD_L_vs_Control_L$NAFLD_Megena_Res$modules %>% 
                    filter(Modules == "c1_6") %>% 
                    dplyr::pull(Genes))





### Metabolites 

library(PMA)
library(ade4)
set.seed(123)
library(phyloseq)
TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")


library(vegan)
library(PMA)
# mantel.partial_for_every <- function(D1,M1) {
#    # D1 = Obese_Species_c2
#    # M1 = Primary_metadata_obese
#    
#    D1 <- D1 %>%  .[rowSums(.!= 0) >= ncol(.)*0.1,]
#    M1 <- M1[which(rownames(M1) %in% rownames(D1)),] %>% 
#       .[rowSums(.!= 0) >= ncol(.)*0.1,]
#    
#    physeq <-  metaphlanToPhyloseq(t(D1))
#    physeq <- merge_phyloseq(physeq, TREE)
#    dist_Wunifrac <- phyloseq::distance(physeq, method = "wunifrac")
#    
#    mantel.partial1 <- mantel.partial(dist_Wunifrac,
#                      dist(M1), 
#                      method = "spearman", permutations = 9999, na.rm = TRUE)
#    return(mantel.partial1)
# }
# 
# new_mantel.partial_CCA_Obese_c2 <-mantel.partial_for_every(Obese_Species_c2,Primary_metadata_obese)
# new_mantel.partial_CCA_Obese_c3 <-mantel.partial_for_every(Obese_Species_c3,Primary_metadata_obese)
# new_mantel.partial_CCA_Obese_c7 <-mantel.partial_for_every(Obese_Species_c7,Primary_metadata_obese)
# new_mantel.partial_CCA_Obese_c10 <-mantel.partial_for_every(Obese_Species_c10,Primary_metadata_obese)
# new_mantel.partial_CCA_Obese_c11 <-mantel.partial_for_every(Obese_Species_c11,Primary_metadata_obese)
# 
# 
# new_mantel.partial_Lean_c2 <- mantel.partial_for_every(Lean_Species_c2,Primary_metadata_lean)
# new_mantel.partial_Lean_c3 <- mantel.partial_for_every(Lean_Species_c3,Primary_metadata_lean)
# new_mantel.partial_Lean_c4 <- mantel.partial_for_every(Lean_Species_c4,Primary_metadata_lean)
# new_mantel.partial_Lean_c5 <- mantel.partial_for_every(Lean_Species_c5,Primary_metadata_lean)
# new_mantel.partial_Lean_c6 <- mantel.partial_for_every(Lean_Species_c6,Primary_metadata_lean)
# 

data.frame(rownames(Obese_Species_c3),
           rownames(Primary_metadata_obese),
           rownames(Adj_metadata_obese))



mantel.partial_CCA_Obese_c2 <- mantel.partial(vegdist(Obese_Species_c2 %>% 
                                                        .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                              dist(Primary_metadata_obese
                                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c2 %>% 
                                                                                                           .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                              dist(Adj_metadata_obese
                                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c2 %>% 
                                                                                                       .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                              method = "spearman", permutations = 9999, na.rm = TRUE)




mantel.partial_CCA_Obese_c3 <- mantel.partial(vegdist(Obese_Species_c3 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c3 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c3 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              method = "spearman", permutations = 9999, na.rm = TRUE)



mantel.partial_CCA_Obese_c4 <- mantel.partial(vegdist(Obese_Species_c4 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c4 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c4 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              
                              method = "spearman", permutations = 9999, na.rm = TRUE)

mantel.partial_CCA_Obese_c7 <- mantel.partial(vegdist(Obese_Species_c7 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c7 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),] ),
                              
                              
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c7 %>% 
                                                                                        .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              method = "spearman", permutations = 9999, na.rm = TRUE)

mantel.partial_CCA_Obese_c10<- mantel.partial(vegdist(Obese_Species_c10 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c10 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),] ),
                              
                              
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c10 %>% 
                                                                                        .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              method = "spearman", permutations = 9999, na.rm = TRUE)

mantel.partial_CCA_Obese_c11<- mantel.partial(vegdist(Obese_Species_c11 %>% 
                                         .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                              dist(Primary_metadata_obese
                                   [which(rownames(Primary_metadata_obese) %in% rownames(Obese_Species_c11 %>% 
                                                                                            .[rowSums(.!= 0) >= ncol(.)*0.1,])),] ),
                              
                              
                              dist(Adj_metadata_obese
                                   [which(rownames(Adj_metadata_obese) %in% rownames(Obese_Species_c11 %>% 
                                                                                        .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                              method = "spearman", permutations = 9999, na.rm = TRUE)




### 1000 perms
CCA_Obese_c4 <- CCA.permute(Obese_Species_c4,Primary_metadata_obese,
                            typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                            penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                         0.3,0.3,0.3,0.3,0.3,0.3,
                                         0.4,0.4,0.4,0.4,0.4,0.4,
                                         0.5,0.5,0.5,0.5,0.5,0.5,
                                         0.6,0.6,0.6,0.6,0.6,0.6,
                                         .7,.7,.7,.7,.7,.7),
                            penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
CCA_Obese_c7<- CCA.permute(Obese_Species_c7,Primary_metadata_obese,
                           typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                           penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                        0.3,0.3,0.3,0.3,0.3,0.3,
                                        0.4,0.4,0.4,0.4,0.4,0.4,
                                        0.5,0.5,0.5,0.5,0.5,0.5,
                                        0.6,0.6,0.6,0.6,0.6,0.6,
                                        .7,.7,.7,.7,.7,.7),
                           penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
CCA_Obese_c10<- CCA.permute(Obese_Species_c10,Primary_metadata_obese,
                            typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                            penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                         0.3,0.3,0.3,0.3,0.3,0.3,
                                         0.4,0.4,0.4,0.4,0.4,0.4,
                                         0.5,0.5,0.5,0.5,0.5,0.5,
                                         0.6,0.6,0.6,0.6,0.6,0.6,
                                         .7,.7,.7,.7,.7,.7),
                            penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
CCA_Obese_c11<- CCA.permute(Obese_Species_c11,Primary_metadata_obese,
                            typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                            penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                         0.3,0.3,0.3,0.3,0.3,0.3,
                                         0.4,0.4,0.4,0.4,0.4,0.4,
                                         0.5,0.5,0.5,0.5,0.5,0.5,
                                         0.6,0.6,0.6,0.6,0.6,0.6,
                                         .7,.7,.7,.7,.7,.7),
                            penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))


library(vegan)



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

mantel.partial_Lean_c3<- mantel.partial(vegdist(Lean_Species_c3 %>% 
                                                  .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                        dist(Primary_metadata_lean
                                             [which(rownames(Primary_metadata_lean) %in% rownames(Lean_Species_c3 %>% 
                                                                                                    .[rowSums(.!= 0) >= ncol(.)*0.1,])),] 
                                        ), 
                                        dist(Adj_metadata_lean
                                             [which(rownames(Adj_metadata_lean) %in% rownames(Lean_Species_c3 %>% 
                                                                                                .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                        
                                        method = "spearman", permutations = 9999, na.rm = TRUE)

mantel.partial_Lean_c4<- mantel.partial(vegdist(Lean_Species_c4 %>% 
                                                  .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                        dist(Primary_metadata_lean
                                             [which(rownames(Primary_metadata_lean) %in% rownames(Lean_Species_c4 %>% 
                                                                                                    .[rowSums(.!= 0) >= ncol(.)*0.1,])),] 
                                        ), 
                                        
                                        dist(Adj_metadata_lean
                                             [which(rownames(Adj_metadata_lean) %in% rownames(Lean_Species_c4 %>% 
                                                                                                .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                        method = "spearman", permutations = 9999, na.rm = TRUE)


mantel.partial_Lean_c5<- mantel.partial(vegdist(Lean_Species_c5 %>% 
                                                  .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                        dist(Primary_metadata_lean
                                             [which(rownames(Primary_metadata_lean) %in% rownames(Lean_Species_c5 %>% 
                                                                                                    .[rowSums(.!= 0) >= ncol(.)*0.1,])),] 
                                        ), 
                                        
                                        dist(Adj_metadata_lean
                                             [which(rownames(Adj_metadata_lean) %in% rownames(Lean_Species_c5 %>% 
                                                                                                .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                        method = "spearman", permutations = 9999, na.rm = TRUE)


mantel.partial_Lean_c6<- mantel.partial(vegdist(Lean_Species_c6 %>% 
                                                  .[rowSums(.!= 0) >= ncol(.)*0.1,], method = "bray"),
                                        dist(Primary_metadata_lean
                                             [which(rownames(Primary_metadata_lean) %in% rownames(Lean_Species_c6 %>% 
                                                                                                    .[rowSums(.!= 0) >= ncol(.)*0.1,])),] 
                                        ), 
                                        
                                        dist(Adj_metadata_lean
                                             [which(rownames(Adj_metadata_lean) %in% rownames(Lean_Species_c6 %>% 
                                                                                                .[rowSums(.!= 0) >= ncol(.)*0.1,])),]),
                                        method = "spearman", permutations = 10000, na.rm = TRUE)


### Lean
CCA_Lean_c2 <- CCA.permute(Lean_Species_c2,Primary_metadata_lean,
                           typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                           penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                        0.3,0.3,0.3,0.3,0.3,0.3,
                                        0.4,0.4,0.4,0.4,0.4,0.4,
                                        0.5,0.5,0.5,0.5,0.5,0.5,
                                        0.6,0.6,0.6,0.6,0.6,0.6,
                                        .7,.7,.7,.7,.7,.7),
                           penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))


CCA_Lean_c5 <- CCA.permute(Lean_Species_c5,Primary_metadata_lean,
                           typex = "standard",typez = "standard",niter=100, nperms = 1000,trace=TRUE,standardize=TRUE,
                           penaltyxs =c(0.2,0.2,0.2,0.2,0.2,0.2,
                                        0.3,0.3,0.3,0.3,0.3,0.3,
                                        0.4,0.4,0.4,0.4,0.4,0.4,
                                        0.5,0.5,0.5,0.5,0.5,0.5,
                                        0.6,0.6,0.6,0.6,0.6,0.6,
                                        .7,.7,.7,.7,.7,.7),
                           penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))



extract_features <- function(D1,M1,CCA_Object) {
  
  # 
  # D1 = Obese_Species_c7
  # M1 = Metabolites_NAFLD_O
  # CCA_Object = Updated_CCA_obese_c7_0.1
  # filter_res = res_mantel_obese_7
  
  
  DF <- 
    data.frame(Penalty_X = CCA_Object$penaltyxs,
               Penalty_Z = CCA_Object$penaltyzs,
               Z_Stat = CCA_Object$zstats,
               Cors = CCA_Object$cors,
               Pval = CCA_Object$pvals) %>% 
    filter(Pval < 0.05) #%>% 
  # arrange(Cors)
  
  
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


t1 <- rbind(extract_Obese_4_0.1$Species,
      extract_Obese_4_0.1$Clinical_Data)


t2 <- rbind(extract_Lean_2_0.1$Species,
            extract_Lean_2_0.1$Clinical_Data)



write.csv(t1,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Supplementary_Table_6_obese.csv")
write.csv(t2,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Supplementary_Table_6_lean.csv")


extract_Obese_4_0.1 <-extract_features(Obese_Species_c4,Primary_metadata_obese,
                                       CCA_Obese_c4)

extract_Obese_7_0.1 <-extract_features(Obese_Species_c7,Primary_metadata_obese,
                                       CCA_Obese_c7)

extract_Lean_2_0.1 <-extract_features(Lean_Species_c2,Primary_metadata_lean,
                                      CCA_Lean_c2)
extract_Lean_5_0.1 <-extract_features(Lean_Species_c5,Primary_metadata_lean,
                                      CCA_Lean_c5)


extract_Obese_4_0.1$Species

write.csv(rbind(extract_Obese_4_0.1$Species,
                extract_Obese_4_0.1$Clinical_Data),
          "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/SCCA_4_Overweight_Results.csv")


test <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/SCCA_4_Overweight_Results.csv")

write.csv(rbind(extract_Obese_7_0.1$Species,
                extract_Obese_7_0.1$Metabolites),
          "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/SCCA_7_Overweight_Results.csv")


write.csv(rbind(extract_Lean_2_0.1$Species,
                extract_Lean_2_0.1$Metabolites),
          "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/SCCA_2_Lean_Results.csv")


write.csv(rbind(extract_Lean_5_0.1$Species,
                extract_Lean_5_0.1$Metabolites),
          "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/SCCA_5_Lean_Results.csv")





cca_CCA_Obese_Species_c6 <- CCA(x=Lean_Species_c2,
                                z=Primary_metadata_lean,
                                typex="standard",
                                typez="standard",
                                penaltyx=0.5,
                                penaltyz=0.6,
                                niter=10000,
                                trace=TRUE,K=3)
lean_c2 <- 
   cca_CCA_Obese_Species_c6$u %>% 
   `rownames<-`(colnames(Lean_Species_c2)) %>% 
   as.data.frame(.) %>% 
   dplyr::select(V1) %>% 
   filter(V1 != 0) %>% 
   rownames_to_column() %>% 
   dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
   mutate(SCCA_Cor = ifelse(SCCA_Cor > 0, "Positive","Negative"))


write.table(lean_c2 %>% 
              mutate(Features = gsub("_"," ",Features)),
            file = "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Prev_All_Non_Adjust_Strict/Lean/SCCA_2.tsv",
            quote = FALSE,sep = "\t",row.names=FALSE)




cca_CCA_Obese_Species_c6 <- CCA(x=Lean_Species_c5,
                                z=Primary_metadata_lean,
                                typex="standard",
                                typez="standard",
                                penaltyx=0.700,
                                penaltyz=0.700,
                                niter=10000,
                                trace=TRUE,K=3)
lean_c5 <- 
   cca_CCA_Obese_Species_c6$u %>% 
   `rownames<-`(colnames(Lean_Species_c5)) %>% 
   as.data.frame(.) %>% 
   dplyr::select(V1) %>% 
   filter(V1 != 0) %>% 
   rownames_to_column() %>% 
   dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
   mutate(SCCA_Cor = ifelse(SCCA_Cor > 0, "Positive","Negative"))


write.table(lean_c5 %>% 
              mutate(Features = gsub("_"," ",Features)),
            file = "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Prev_All_Non_Adjust_Strict/Lean/SCCA_5.tsv",
            quote = FALSE,sep = "\t",row.names=FALSE)





cca_CCA_Obese_Species_c6 <- CCA(x=Obese_Species_c4,
                                z=Primary_metadata_obese,
                                typex="standard",
                                typez="standard",
                                penaltyx=0.7,
                                penaltyz=0.7,
                                niter=10000,
                                trace=TRUE,K=3)
obese_c4 <-
   cca_CCA_Obese_Species_c6$u %>% 
   `rownames<-`(colnames(Obese_Species_c4)) %>% 
   as.data.frame(.) %>% 
   dplyr::select(V1) %>% 
   filter(V1 != 0) %>% 
   rownames_to_column() %>% 
   dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
   mutate(SCCA_Cor = ifelse(SCCA_Cor > 0, "Positive","Negative"))



write.table(obese_c4 %>% 
          mutate(Features = gsub("_"," ",Features)),
            file = "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Prev_All_Non_Adjust_Strict/Obese/SCCA_4.tsv",
            quote = FALSE,sep = "\t",row.names=FALSE)



cca_CCA_Obese_Species_c6 <- CCA(x=Obese_Species_c7,
                                z=Primary_metadata_obese,
                                typex="standard",
                                typez="standard",
                                penaltyx=0.7,
                                penaltyz=0.7,
                                niter=10000,
                                trace=TRUE,K=3)
obese_c7 <-
   cca_CCA_Obese_Species_c6$u %>% 
   `rownames<-`(colnames(Obese_Species_c7)) %>% 
   as.data.frame(.) %>% 
   dplyr::select(V1) %>% 
   filter(V1 != 0) %>% 
   rownames_to_column() %>% 
   dplyr::rename(Features = "rowname", SCCA_Cor = "V1") %>% 
   mutate(SCCA_Cor = ifelse(SCCA_Cor > 0, "Positive","Negative"))



write.table(obese_c7 %>% 
              mutate(Features = gsub("_"," ",Features)),
            file = "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Prev_All_Non_Adjust_Strict/Obese/SCCA_7.tsv",
            quote = FALSE,sep = "\t",row.names=FALSE)

