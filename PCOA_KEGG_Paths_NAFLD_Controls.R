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

##Load data
load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_KEGGPaths_CPM_Clinical_Data.RData")

##### Data Prep ####

NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))
Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
Control_L <- dplyr::select(Control_L, -one_of(c("Group")))



#### Preparing Data Frames 

full_features = bind_rows(NAFLD_O,
                          NAFLD_L,
                          Control_O,
                          Control_L)

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



full_features =  full_features[which(rownames(full_features) %in% rownames(full_metadata)),]


#### BC PCOA ####
set.seed(123)
BC_Comb = vegdist(full_features,method = "bray")
PCOA = ape::pcoa(BC_Comb)
full = cbind(full_features,full_metadata)


PWA_BC_treat_NMDS_2 <- adonis2(BC_Comb ~  Age + Gender + BMI + HOMA_IR + SBP  + Groups, data = full, permutations = 999,by = "terms")



##GG Plot

ggplot_df <-
   data.frame(PCoA_1 = PCOA$vectors[,1],PCoA_2 = PCOA$vectors[,2],Groups = full$Groups) %>% 
   group_by(Groups) %>% 
   summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
   rename(.,PCoA_1 = PCoA_1_mean,PCoA_2 = PCoA_2_mean,SE_1 = PCoA_1_se,SE_2 = PCoA_2_se) %>% 
   mutate(Groups = gsub("Control-O","CTRL-NAFLD-O",Groups)) %>% 
   mutate(Groups = gsub("Control-L","CTRL-NAFLD-L",Groups)) %>% 
   mutate(Groups = factor(Groups, levels = c("NAFLD-O",
                                            "NAFLD-L",
                                            "CTRL-NAFLD-O",
                                            "CTRL-NAFLD-L")))





PCOA$values$Eigenvalues[1]/sum(PCOA$values$Eigenvalues)
PCOA$values$Eigenvalues[2]/sum(PCOA$values$Eigenvalues)


pcoa_bc_beta_div <-ggplot(ggplot_df, aes(x=PCoA_1, y=PCoA_2,colour = Groups, label = Groups)) + 
   geom_errorbar(aes(ymin=PCoA_2 - SE_2, ymax=PCoA_2 + SE_2),linetype = 1,size =1.5,width = 0.001) +
   geom_errorbar(aes(xmin=PCoA_1 - SE_1, xmax=PCoA_1 + SE_1),linetype = 1,size =1.5,width = 0.001) +
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
         axis.title=element_text(size=20)) +
   annotate("text", x = -0.03, y = 0.02, label = paste(("P value = "),PWA_BC_treat_NMDS_2$`Pr(>F)`[6]),size=8) +
   annotate("text", x = -0.03, y = 0.018, label = paste(("R² = 0.006")),size=8) +
   labs(x = "PCoA 1 68%",
        y = "PCoA 2 14%")


##Adonis Function for each single group vs another single group to test beta diversity between them

adonis_new <- function(Metadata) {
  set.seed(123)
  
  
  
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]
  
  new_BC_Comb = vegdist(Features,method = "bray")
  new_PCOA = ape::pcoa(new_BC_Comb)
  new_full = cbind(Features,Metadata)
  
  
  new_adonis_Wnifrac_treat <- adonis2(new_BC_Comb ~  Age + Gender + BMI + HOMA_IR + SBP + Groups, 
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
                  R = c((Nafld_O_vs_L$R2[6]),
                        (Nafld_O_vs_Control_O$R2[6]),
                        (Nafld_L_vs_Control_L$R2[6])),
                  Pvalue = c(as.numeric(Nafld_O_vs_L$`Pr(>F)`[6]),
                             as.numeric(Nafld_O_vs_Control_O$`Pr(>F)`[6]),
                             as.numeric(Nafld_L_vs_Control_L$`Pr(>F)`[6]))
)




