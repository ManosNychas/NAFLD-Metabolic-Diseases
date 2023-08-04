rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyr)
library(magrittr)
library(vegan)
library(ggplot2)
library(phyloseq)
library(picante)
library(ggrepel)
library(cowplot)
library(ggsignif)
library(textshape)
library(janitor)

##Load Data

load("Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

df <-
   read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>% 
   as.data.frame() %>% 
   tibble::column_to_rownames(var = "X") %>% 
   t() %>% 
   as.data.frame() 

##### Data Prep ####

NAFLD_O <- df[which(rownames(df) %in% rownames(NAFLD_O)),]%>% 
   magrittr::add(1) %>%
   log2()

NAFLD_L <- df[which(rownames(df) %in% rownames(NAFLD_L)),]%>% 
   magrittr::add(1) %>%
   log2()
Control_O <- df[which(rownames(df) %in% rownames(Control_O)),]%>% 
   magrittr::add(1) %>%
   log2()
Control_L <-df[which(rownames(df) %in% rownames(Control_L)),]%>% 
   magrittr::add(1) %>%
   log2()

Overweight <- df[which(rownames(df) %in% rownames(Overweight)),]%>% 
   magrittr::add(1) %>%
   log2()
T2D_Obese <- df[which(rownames(df) %in% rownames(T2D_Obese)),]%>% 
   magrittr::add(1) %>%
   log2()
Prediabetes <- df[which(rownames(df) %in% rownames(Prediabetes)),]%>% 
   magrittr::add(1) %>%
   log2()

Control_Hypertension <- df[which(rownames(df) %in% rownames(Control_Hypertension)),]%>% 
   magrittr::add(1) %>%
   log2()
Hypertension <- df[which(rownames(df) %in% rownames(Hypertension)),]%>% 
   magrittr::add(1) %>%
   log2()
Pre_Hypertension <- df[which(rownames(df) %in% rownames(Pre_Hypertension)),]%>% 
   magrittr::add(1) %>%
   log2()

Lean_Normal<- df[which(rownames(df) %in% rownames(Lean_Normal)),]%>% 
   magrittr::add(1) %>%
   log2()
T2D_Lean <- df[which(rownames(df) %in% rownames(T2D_Lean)),]%>% 
   magrittr::add(1) %>%
   log2()

Atherosclerosis <- df[which(rownames(df) %in% rownames(Atherosclerosis)),]%>% 
   magrittr::add(1) %>%
   log2()
Control_Atherosclerosis <- df[which(rownames(df) %in% rownames(Control_Atherosclerosis)),]%>% 
   magrittr::add(1) %>%
   log2()



#### Preparing Data Frames 


full_features = rbind(NAFLD_O,
                      NAFLD_L,
                      Control_O,
                      Control_L,
                      Overweight,
                      T2D_Obese,
                      Prediabetes,
                      Control_Hypertension,
                      Hypertension,
                      Pre_Hypertension,
                      Lean_Normal,
                      T2D_Lean,
                      Atherosclerosis,
                      Control_Atherosclerosis)

full_features[is.na(full_features)] <- 0






df <- read.csv("~/Documents/Phd/HKI_PC/Metabolic_Diseases//Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv")

full_metadata  <- df  %>% .[which(df$X %in% rownames(full_features)),]
full_metadata = full_metadata[order(match(full_metadata$X,rownames(full_features))),] 
full_metadata <- full_metadata %>%
   dplyr::select(X,
                 Age.demographic_variable,
                 Gender.demographic_variable,
                 BMI.body_weight_variables,
   ) %>% 
   dplyr::rename(IDs = "X") %>% 
   dplyr::rename(Age = "Age.demographic_variable") %>% 
   dplyr::rename(BMI = "BMI.body_weight_variables") %>% 
   dplyr::rename(Gender = "Gender.demographic_variable") %>% 
   dplyr::mutate(Groups = factor(c(rep("NAFLD-O",nrow(NAFLD_O)),
                                   rep("NAFLD-L",nrow(NAFLD_L)),
                                   rep("Control-0",nrow(Control_O)),
                                   rep("Control-L",nrow(Control_L)),
                                   rep("Overweight",nrow(Overweight)),
                                   rep("T2D_Obese",nrow(T2D_Obese)),
                                   rep("Prediabetes",nrow(Prediabetes)),
                                   rep("Control_Hypertension",nrow(Control_Hypertension)),
                                   rep("Hypertension",nrow(Hypertension)),
                                   rep("Pre_Hypertension",nrow(Pre_Hypertension)),
                                   rep("Lean_Normal",nrow(Lean_Normal)),
                                   rep("T2D_Lean",nrow(T2D_Lean)),
                                   rep("Atherosclerosis",nrow(Atherosclerosis)),
                                   rep("Control_Atherosclerosis",nrow(Control_Atherosclerosis))),
                                 levels =c("NAFLD-O",
                                           "NAFLD-L",
                                           "Control-0",
                                           "Control-L",
                                           "Overweight",
                                           "T2D_Obese",
                                           "Prediabetes",
                                           "Control_Hypertension",
                                           "Hypertension",
                                           "Pre_Hypertension",
                                           "Lean_Normal",
                                           "T2D_Lean",
                                           "Atherosclerosis",
                                           "Control_Atherosclerosis"
                                 ))) %>% 
   textshape::column_to_rownames(., "IDs") %>% 
   .[complete.cases(.), ] 


library(mixOmics)

full_features =  full_features[which(rownames(full_features) %in% rownames(full_metadata)),]
#### PLSDA 

###3 Euclidean Distance PCOA
set.seed(123)

full = cbind(full_features,full_metadata)
PCOA_3_BC <- full_features %>% 
   vegdist(.,method = "euclidean", trymax = 50) %>% 
   ape::pcoa()

distance <- full_features %>% 
   vegdist(.,method = "euclidean", trymax = 50)



Adonis_BC <- adonis2(distance ~  BMI + Age + Gender + Groups, data = full, permutations = 999,by = "terms")

##GG PLot

ggplot_df <-
   data.frame(MDS1 = PCOA_3_BC$vectors[,1],MDS2 = PCOA_3_BC$vectors[,2],Groups = full_metadata$Groups) %>% 
   dplyr::group_by(Groups) %>% 
   dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
  dplyr::rename(.,NMDS1 = MDS1_mean,NMDS2 = MDS2_mean,SE_1 = MDS1_se,SE_2 = MDS2_se) %>% 
  dplyr::mutate(Groups = gsub("Control-0","CTRL-NAFLD-O",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Control-L","CTRL-NAFLD-L",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Overweight","CTRL-T2D-O",Groups)) %>% 
  dplyr::mutate(Groups = gsub("T2D_Obese","T2D-O",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Control_Atherosclerosis","CTRL-ATH",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Pre_Hypertension","PRE-HYP",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Control_Hypertension","CTRL-HYP",Groups)) %>% 
  dplyr::mutate(Groups = gsub("T2D_Lean","T2D-L",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Lean_Normal","CTRL-T2D-L",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Prediabetes","PRE-T2D-O",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Atherosclerosis","ATH",Groups)) %>% 
  dplyr::mutate(Groups = gsub("Hypertension","HYP",Groups)) %>% 
  dplyr::mutate(Groups = factor(Groups, levels = c("NAFLD-O",
                                            "NAFLD-L",
                                            "CTRL-NAFLD-O",
                                            "CTRL-NAFLD-L",
                                            "CTRL-T2D-O",
                                            "T2D-O",
                                            "PRE-T2D-O",
                                            "CTRL-HYP",
                                            "HYP",
                                            "PRE-HYP",
                                            "CTRL-T2D-L",
                                            "T2D-L",
                                            "ATH",
                                            "CTRL-ATH")))


PCOA_3_BC$values$Eigenvalues[1]/sum(PCOA_3_BC$values$Eigenvalues)
PCOA_3_BC$values$Eigenvalues[2]/sum(PCOA_3_BC$values$Eigenvalues)

nmds_bc_beta_div <-ggplot(ggplot_df, aes(x=NMDS1, y=NMDS2,colour = Groups, label = Groups)) + 
   geom_errorbar(aes(ymin=NMDS2 - SE_2, ymax=NMDS2 + SE_2),linetype = 1,size =1.5,width = 0.1) +
   geom_errorbar(aes(xmin=NMDS1 - SE_1, xmax=NMDS1 + SE_1),linetype = 1,size =1.5,width = 0.1) +
   geom_label(size = 4) + 
  scale_color_manual(values=c("#d86b6f",
                                       "#6182a8", #NAFLD
                                       
                                       "#d86b6f",
                                       "#6182a8", #Control
                                       
                                       "#FF8C00",
                                       "#FF8C00",	
                                       "#FF8C00",## T2D Obese 
                                       
                                       "#32CD32",
                                       "#32CD32",	
                                       "#32CD32", #Hypertension
                                       
                                       "#9966cc",	
                                       "#9966cc", #T2D Lean
                                       
                                       "#9d7029", 
                                       "#9d7029" #Atheroslcerosis
                                       
  )) +
   theme_classic() +
   theme(plot.title = element_text(size = 14, hjust = 0.5),  text = element_text(size = 15), 
         legend.position="none",
         axis.text.x=element_blank(), #remove x axis labels
         axis.ticks.x=element_blank(), #remove x axis ticks
         axis.text.y=element_blank(),  #remove y axis labels
         axis.ticks.y=element_blank(), 
         panel.border = element_rect(colour = "black", fill=NA, size=2),
         axis.title=element_text(size=20)
   ) + 
   labs(x = "PCoA 1 15%",
        y = "PCoA 2 1%")



##Adonis Function for each single group vs another single group to test beta diversity between them

adonis_new <- function(Metadata) {
  set.seed(123)
  
  
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]
  
  
  new_full = cbind(Features,Metadata)
  
  
  distance <-   Features %>% 
    vegdist(.,method = "euclidean", trymax = 50)
  
  test2 <- as.matrix(distance)
  new_adonis_Wnifrac_treat <- adonis2(distance ~  Age + Gender + BMI + Groups, 
                                      data = new_full, permutations = 999,
                                      by = "terms")
  
  return(new_adonis_Wnifrac_treat)
}





t1  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "T2D_Obese") %>% 
  adonis_new(.)

t2 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Prediabetes") %>% 
  adonis_new(.)
t3 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Hypertension") %>% 
  adonis_new(.)

t4 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Pre_Hypertension") %>% 
  adonis_new(.)

t5 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "T2D_Lean") %>% 
  adonis_new(.)

t6 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Atherosclerosis") %>% 
  adonis_new(.)


l1  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "T2D_Obese") %>% 
  adonis_new(.)

l2 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Prediabetes") %>% 
  adonis_new(.)
l3 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Hypertension") %>% 
  adonis_new(.)

l4 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Pre_Hypertension") %>% 
  adonis_new(.)

l5 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "T2D_Lean") %>% 
  adonis_new(.)

l6 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Atherosclerosis") %>% 
  adonis_new(.)


tm1  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control-0") %>% 
  adonis_new(.)

tm2  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control-L") %>% 
  adonis_new(.)


tm3  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "NAFLD-L") %>% 
  adonis_new(.)



RES <- data.frame(Comparison = c("NAFLD O vs T2D_Obese",
                                 "NAFLD O vs Prediabetes",
                                 "NAFLD O vs Hypertension",
                                 "NAFLD O vs Pre_Hypertension",
                                 "NAFLD O vs T2D_Lean",
                                 "NAFLD O vs Atherosclerosis",
                                 
                                 "NAFLD L vs T2D_Obese",
                                 "NAFLD L vs Prediabetes",
                                 "NAFLD L vs Hypertension",
                                 "NAFLD L vs Pre_Hypertension",
                                 "NAFLD L vs T2D_Lean",
                                 "NAFLD L vs Atherosclerosis",
                                 
                                 "NAFLD O vs Control O",
                                 "NAFLD L vs Control L",
                                 "NAFLD O vs NAFLD L"),

                  Pvalue = c(as.numeric(t1$`Pr(>F)`[4]),
                             as.numeric(t2$`Pr(>F)`[4]),
                             as.numeric(t3$`Pr(>F)`[4]),
                             as.numeric(t4$`Pr(>F)`[4]),
                             as.numeric(t5$`Pr(>F)`[4]),
                             as.numeric(t6$`Pr(>F)`[4]),
                             as.numeric(l1$`Pr(>F)`[4]),
                             as.numeric(l2$`Pr(>F)`[4]),
                             as.numeric(l3$`Pr(>F)`[4]),
                             as.numeric(l4$`Pr(>F)`[4]),
                             as.numeric(l5$`Pr(>F)`[4]),
                             as.numeric(l6$`Pr(>F)`[4]),
                             as.numeric(tm1$`Pr(>F)`[4]),
                             as.numeric(tm2$`Pr(>F)`[4]),
                             as.numeric(tm3$`Pr(>F)`[4]))
)
