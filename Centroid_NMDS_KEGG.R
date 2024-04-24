rm(list = ls(all = TRUE))
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


load("Users/manosnychas/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_KeggMods_CPM_Clinical_Data.RData")



##### Data Prep ####


NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))
Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
Control_L <- dplyr::select(Control_L, -one_of(c("Group")))

Overweight <- dplyr::select(Overweight, -one_of(c("Group")))
T2D_Obese <- dplyr::select(T2D_Obese, -one_of(c("Group")))
Prediabetes <- dplyr::select(Prediabetes, -one_of(c("Group")))

Control_Hypertension <- dplyr::select(Control_Hypertension, -one_of(c("Group")))
Hypertension <- dplyr::select(Hypertension, -one_of(c("Group")))
Pre_Hypertension <- dplyr::select(Pre_Hypertension, -one_of(c("Group")))

Lean_Normal<- dplyr::select(Lean_Normal, -one_of(c("Group")))
T2D_Lean <- dplyr::select(T2D_Lean, -one_of(c("Group")))

Atherosclerosis <- dplyr::select(Atherosclerosis, -one_of(c("Group")))
Control_Atherosclerosis <- dplyr::select(Control_Atherosclerosis, -one_of(c("Group")))

#### Preparing Data Frames 

full_features = dplyr::bind_rows(NAFLD_O,
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



df <- read.csv("Users/manosnychas/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv")


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

full_features =  full_features[which(rownames(full_features) %in% rownames(full_metadata)),]

#### Distance Calculations #####
beta_distances_dataset <- function(dataset){
  set.seed(123)
  
  
  current_metadata <- full_metadata[which(rownames(full_metadata) %in% rownames(dataset)),]
  
  
  NMDS_3_BC <- dataset %>% 
    vegdist(.,method = "bray") %>% 
    metaMDS(.,k=3) %>% 
    vegan::scores(.) %>% 
    as.data.frame() %>% 
    dplyr::select(1,2)

  
  
  distance_list <- list(
    BC_NMDS = NMDS_3_BC
  )
  
  restults_list <- list(distances = distance_list)
  #adonis = adonis_list)
  
  return(restults_list)
  
}


Distance_Full <- beta_distances_dataset(full_features)


##### Centroid Calculations ####
full <- 
  full_features %>% 
  tibble::rownames_to_column(., "Names") %>% 
  left_join(full_metadata %>% 
              tibble::rownames_to_column(., "Names") %>% 
              dplyr::select(Names,Groups), by = "Names") %>% 
  tibble::column_to_rownames(var = "Names") %>% 
  dplyr::filter(Groups != "Overweight",
                Groups != "Control_Hypertension",
                Groups != "Lean_Normal",
                Groups != "Control_Atherosclerosis",
                Groups != "Control-0",
                Groups != "Control-L") 


distance_dataset <- Distance_Full$distances$BC_NMDS
#### Centroid Calculations ####
new_distance_dataset <- list()

  
  Centroid_Control_L <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Control_L)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean) 
  
  Centroid_Control_O <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Control_O)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean)
  
  Centroid_Overweight <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Overweight)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean)
  
  Centroid_Control_Hypertension <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Control_Hypertension)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean)
  
  Centroid_Lean_Normal <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Lean_Normal)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean)
  
  Centroid_Control_Atherosclerosis <-
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    filter(Names %in% rownames(Control_Atherosclerosis)) %>% 
    dplyr::select(2,3) %>% 
    apply(.,2,mean)
  
  new_distance_dataset <- 
    distance_dataset %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., "Names") %>% 
    mutate(NMDS1  =ifelse(Names %in% rownames(NAFLD_O), .[[2]] - Centroid_Control_O[[1]],
                          ifelse(Names %in% rownames(NAFLD_L),.[[2]] - Centroid_Control_L[[1]],
                                 ifelse(Names %in% rownames(T2D_Obese),.[[2]] - Centroid_Overweight[[1]],
                                        ifelse(Names %in% rownames(Prediabetes),.[[2]] - Centroid_Overweight[[1]],
                                               ifelse(Names %in% rownames(Hypertension),.[[2]] - Centroid_Control_Hypertension[[1]],
                                                      ifelse(Names %in% rownames(Pre_Hypertension),.[[2]] - Centroid_Control_Hypertension[[1]],
                                                             ifelse(Names %in% rownames(T2D_Lean),.[[2]] - Centroid_Lean_Normal[[1]],
                                                                    ifelse(Names %in% rownames(Atherosclerosis),.[[2]] - Centroid_Control_Atherosclerosis[[1]],
                                                                           .[[2]] -0))))))))) %>% 
    mutate(NMDS2  =ifelse(Names %in% rownames(NAFLD_O), .[[3]] - Centroid_Control_O[[2]],
                          ifelse(Names %in% rownames(NAFLD_L),.[[3]] - Centroid_Control_L[[2]],
                                 ifelse(Names %in% rownames(T2D_Obese),.[[3]] - Centroid_Overweight[[2]],
                                        ifelse(Names %in% rownames(Prediabetes),.[[3]] - Centroid_Overweight[[2]],
                                               ifelse(Names %in% rownames(Hypertension),.[[3]] - Centroid_Control_Hypertension[[2]],
                                                      ifelse(Names %in% rownames(Pre_Hypertension),.[[3]] - Centroid_Control_Hypertension[[2]],
                                                             ifelse(Names %in% rownames(T2D_Lean),.[[3]] - Centroid_Lean_Normal[[2]],
                                                                    ifelse(Names %in% rownames(Atherosclerosis),.[[3]] - Centroid_Control_Atherosclerosis[[2]],
                                                                           .[[3]] -0))))))))) %>% 
    filter(Names %in% rownames(full))



#### BC NMDS ####

ggplot_df <-
  data.frame(MDS1 = new_distance_dataset$NMDS1,
             MDS2 = new_distance_dataset$NMDS2,Groups = full$Groups) %>% 
  group_by(Groups) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
  rename(.,NMDS1 = MDS1_mean,NMDS2 = MDS2_mean,SE_1 = MDS1_se,SE_2 = MDS2_se) %>% 
  mutate(Groups = gsub("Control-0","CTRL-NAFLD-O",Groups)) %>% 
  mutate(Groups = gsub("Control-L","CTRL-NAFLD-L",Groups)) %>% 
  mutate(Groups = gsub("Overweight","CTRL-T2D-O",Groups)) %>% 
  mutate(Groups = gsub("T2D_Obese","T2D-O",Groups)) %>% 
  mutate(Groups = gsub("Control_Atherosclerosis","CTRL-ATH",Groups)) %>% 
  mutate(Groups = gsub("Pre_Hypertension","PRE-HYP",Groups)) %>% 
  mutate(Groups = gsub("Control_Hypertension","CTRL-HYP",Groups)) %>% 
  mutate(Groups = gsub("T2D_Lean","T2D-L",Groups)) %>% 
  mutate(Groups = gsub("Lean_Normal","CTRL-T2D-L",Groups)) %>% 
  mutate(Groups = gsub("Prediabetes","PRE-T2D-O",Groups)) %>% 
  mutate(Groups = gsub("Atherosclerosis","ATH",Groups)) %>% 
  mutate(Groups = gsub("Hypertension","HYP",Groups)) %>% 
  mutate(Groups = factor(Groups, levels = c("NAFLD-O",
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



nmds_bc_beta_div <-ggplot(ggplot_df, aes(x=NMDS1, y=NMDS2,colour = Groups, label = Groups)) + 
  geom_errorbar(aes(ymin=NMDS2 - SE_2, ymax=NMDS2 + SE_2),linetype = 1,size =1) +
  geom_errorbar(aes(xmin=NMDS1 - SE_1, xmax=NMDS1 + SE_1),linetype = 1,size =1) +
  geom_label(size = 7) +
  scale_color_manual(values=c("#d86b6f",
                              "#6182a8", #NAFLD
                              
                              "#FF8C00",
                              "#FF8C00",## T2D Obese
                              
                              "#32CD32",
                              "#32CD32", #Hypertension
                              
                              "#9966cc", #T2D Lean
                              "#9d7029" #Atheroslcerosis
                              
  )) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),  text = element_text(size = 15),
        legend.position="none",
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.title=element_text(size=20)) 

nmds_bc_beta_div

