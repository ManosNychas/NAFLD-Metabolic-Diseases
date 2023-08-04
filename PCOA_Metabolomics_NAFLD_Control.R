rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyr)
library(RANN)
library(magrittr)
library(vegan)
library(ggplot2)
library(phyloseq)
library(picante)
library(ggrepel)
library(cowplot)
library(ggsignif)
library(fossil)
library(GUniFrac)
library(textshape)
library(janitor)


##Loading Data

load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

df <-
   read.csv("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>% 
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



#### Preparing Data Frames for Phylosrec


full_features = rbind(NAFLD_O,
                      NAFLD_L,
                      Control_O,
                      Control_L
                      )

full_features[is.na(full_features)] <- 0






df <- read.csv("~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv")

full_metadata  <- df  %>% .[which(df$X %in% rownames(full_features)),]
full_metadata = full_metadata[order(match(full_metadata$X,rownames(full_features))),] 
full_metadata <- full_metadata %>%
   dplyr::select(X,
                 Age.demographic_variable,
                 Gender.demographic_variable,
                 BMI.body_weight_variables,
                 HOMAIR.glucose_insulin,
                 Systolic_pressure.cardiac_variables
   ) %>% 
   dplyr::rename(IDs = "X") %>% 
   dplyr::rename(Age = "Age.demographic_variable") %>% 
   dplyr::rename(BMI = "BMI.body_weight_variables") %>% 
   dplyr::rename(Gender = "Gender.demographic_variable") %>% 
   dplyr::rename(HOMA_IR = "HOMAIR.glucose_insulin") %>% 
   dplyr::rename(SBP = "Systolic_pressure.cardiac_variables") %>% 
   dplyr::mutate(Groups = factor(c(rep("NAFLD-O",nrow(NAFLD_O)),
                                   rep("NAFLD-L",nrow(NAFLD_L)),
                                   rep("Control-O",nrow(Control_O)),
                                   rep("Control-L",nrow(Control_L))),
                                 levels =c("NAFLD-O",
                                           "NAFLD-L",
                                           "Control-O",
                                           "Control-L"
                                 ))) %>% 
   textshape::column_to_rownames(., "IDs") %>% 
   .[complete.cases(.), ] 


library(mixOmics)

full_features =  full_features[which(rownames(full_features) %in% rownames(full_metadata)),]
#### PLSDA 
set.seed(123)
full = cbind(full_features,full_metadata)

PCOA_3_BC <- full_features %>% 
   vegdist(.,method = "euclidean", trymax = 50) %>% 
   ape::pcoa()

distance <-   full_features %>% 
   vegdist(.,method = "euclidean", trymax = 50)

adonis_Wnifrac_treat <- adonis2(distance ~  BMI + Age + Gender + HOMA_IR + SBP + Groups, data = full, 
                                permutations = 999,by = "terms")



## GG plot
ggplot_df <-
   data.frame(MDS1 = PCOA_3_BC$vectors[,1],MDS2 = PCOA_3_BC$vectors[,2],Groups = full_metadata$Groups) %>% 
   group_by(Groups) %>% 
   summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
   rename(.,NMDS1 = MDS1_mean,NMDS2 = MDS2_mean,SE_1 = MDS1_se,SE_2 = MDS2_se) %>% 
  mutate(Groups = gsub("Control-O","CTRL-NAFLD-O",Groups)) %>% 
  mutate(Groups = gsub("Control-L","CTRL-NAFLD-L",Groups)) %>% 
  mutate(Groups = factor(Groups, levels = c("NAFLD-O",
                                            "NAFLD-L",
                                            "CTRL-NAFLD-O",
                                            "CTRL-NAFLD-L")))



PCOA_3_BC$values$Eigenvalues[1]/sum(PCOA_3_BC$values$Eigenvalues)
PCOA_3_BC$values$Eigenvalues[2]/sum(PCOA_3_BC$values$Eigenvalues)

nmds_bc_beta_div <-ggplot(ggplot_df, aes(x=NMDS1, y=NMDS2,colour = Groups, label = Groups)) + 
   geom_errorbar(aes(ymin=NMDS2 - SE_2, ymax=NMDS2 + SE_2),linetype = 1,size =1.5,width = 0.1) +
   geom_errorbar(aes(xmin=NMDS1 - SE_1, xmax=NMDS1 + SE_1),linetype = 1,size =1.5,width = 0.1) +
  geom_label(size = 8) + 
  scale_color_manual(values=c("#d86b6f","#6182a8", "#d86b6f", "#6182a8")) +
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
  annotate("text", x = -4, y = -1, label = paste(("P value = "),adonis_Wnifrac_treat$`Pr(>F)`[6]),size=8) +
  annotate("text", x = -4, y = -1.1, label = paste(("R² = 0.016")),size=8) +
  
   labs(x = "PCoA 1 15%",
        y = "PCoA 2 1%") 



##Adonis Function for each single group vs another single group to test beta diversity between them

adonis_new <- function(Metadata) {
  set.seed(123)
  # 
  # Metadata<-  full_metadata %>% 
  #   dplyr::filter(Groups == "NAFLD-O" | Groups == "Control-O")
  
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]


  new_full = cbind(Features,Metadata)
  
  
  distance <-   Features %>% 
    vegdist(.,method = "euclidean", trymax = 50)

  
  new_adonis_Wnifrac_treat <- adonis2(distance ~  Age + Gender + BMI + HOMA_IR + SBP + Groups, 
                                      data = new_full, permutations = 999,
                                      by = "terms")
  
  return(new_adonis_Wnifrac_treat)
}



Nafld_O_vs_L <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "NAFLD-L") %>% 
  adonis_new(.)

Nafld_O_vs_Control_O <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control-O") %>% 
  adonis_new(.)

Nafld_L_vs_Control_L <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control-L") %>% 
  adonis_new(.)


RES <- data.frame(Comparison = c("NAFLD O vs NAFLD L",
                                 "NAFLD O vs Control O",
                                 "NAFLD L vs Control L"),
                  R = c(as.numeric(Nafld_O_vs_L$R2[6]),
                        as.numeric(Nafld_O_vs_Control_O$R2[6]),
                        as.numeric(Nafld_L_vs_Control_L$R2[6])),
                  Pvalue = c(as.numeric(Nafld_O_vs_L$`Pr(>F)`[6]),
                             as.numeric(Nafld_O_vs_Control_O$`Pr(>F)`[6]),
                             as.numeric(Nafld_L_vs_Control_L$`Pr(>F)`[6]))
)




