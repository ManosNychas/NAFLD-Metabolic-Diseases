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
library(GUniFrac)
library(textshape)

## Function for metaphlan to phyloseq
metaphlanToPhyloseq <- function(
      tax,
      metadat=NULL,
      simplenames=TRUE,
      roundtointeger=FALSE,
      split="|"){
   ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
   ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
   ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
   ## if roundtointeger=TRUE, values will be rounded to the nearest integer
   xnames = rownames(tax)
   shortnames = gsub(paste0(".+\\", split), "", xnames)
   if(simplenames){
      rownames(tax) = shortnames
   }
   if(roundtointeger){
      tax = round(tax * 1e4)
   }
   x2 = strsplit(xnames, split=split, fixed=TRUE)
   taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
   colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
   rownames(taxmat) = rownames(tax)
   for (i in 1:nrow(taxmat)){
      taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
   }
   taxmat = gsub("[a-z]__", "", taxmat)
   taxmat = phyloseq::tax_table(taxmat)
   otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
   if(is.null(metadat)){
      res = phyloseq::phyloseq(taxmat, otutab)
   }else{
      res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
   }
   return(res)
}


#Loading Data
load("~/Documents/Phd/HKI_PC/PhD/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

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

#### Preparing Data Frames for Phyloseq

full_features = bind_rows(NAFLD_O,
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



df <- read.csv("~/~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv")



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




# creating the phyloseq data
physeq <- metaphlanToPhyloseq(t(full_features))
SAM <-  sample_data(data.frame(full_metadata))
TREE <- read_tree("~/~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
physeq <- merge_phyloseq(physeq, SAM, TREE)

SAM = sample_data(physeq) %>% data.frame(.)

#### Weighted Unifrac NMD ####

set.seed(123)
full = cbind(full_features,full_metadata)
ord_w_unifrac <- ordinate(physeq, method="PCoA", distance="unifrac", weighted=TRUE)
dist_Wunifrac <- phyloseq::distance(physeq, method = "wunifrac")


adonis_Res <- adonis2(dist_Wunifrac ~  BMI + Age + Gender + Groups, data = full, permutations = 999,by = "terms")


## GG plot
ggplot_df <-
   data.frame(PCoA_1 = ord_w_unifrac$vectors[,1],
              PCoA_2 = ord_w_unifrac$vectors[,2],
              Groups = SAM$Groups) %>% 
   group_by(Groups) %>% 
   summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
   rename(.,PCoA_1 = PCoA_1_mean,PCoA_2 = PCoA_2_mean,SE_1 = PCoA_1_se,SE_2 = PCoA_2_se) %>% 
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


ord_w_unifrac$values$Eigenvalues[1]/sum(ord_w_unifrac$values$Eigenvalues)
ord_w_unifrac$values$Eigenvalues[2]/sum(ord_w_unifrac$values$Eigenvalues)


pcoa_wunifrac_beta_div <-ggplot(ggplot_df, aes(x=PCoA_1, y=PCoA_2,colour = Groups, label = Groups)) + 
   geom_errorbar(aes(ymin=PCoA_2 - SE_2, ymax=PCoA_2 + SE_2),linetype = 1,size =1.5) +
   geom_errorbar(aes(xmin=PCoA_1 - SE_1, xmax=PCoA_1 + SE_1),linetype = 1,size =1.5) +
   geom_label(size = 7) + 
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
  labs(x = "PCoA1 43%",
       y = "PCoA2 9%") +
   theme(plot.title = element_text(size = 14, hjust = 0.5),  text = element_text(size = 15), 
         legend.position="none",
         axis.text.x=element_blank(), 
         axis.ticks.x=element_blank(), 
         axis.text.y=element_blank(),  
         axis.ticks.y=element_blank(),  
         panel.border = element_rect(colour = "black", fill=NA, size=2),
         axis.title=element_text(size=20)
         )  





##Adonis Function for each single group vs another single group to test beta diversity between them

adonis_new <- function(Metadata) {
  set.seed(123)
  
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]
  new_physeq <- metaphlanToPhyloseq(t(Features))
  new_SAM <-  sample_data(data.frame(Metadata))
  new_TREE <- read_tree("~/~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
  new_physeq <- merge_phyloseq(new_physeq, new_SAM, new_TREE)
  
  new_SAM = sample_data(new_physeq) %>% data.frame(.)
  
  new_dist_Wunifrac <- phyloseq::distance(new_physeq, method = "wunifrac")
  set.seed(123)
  new_adonis_Wnifrac_treat <- adonis2(new_dist_Wunifrac ~  Age + Gender + BMI + Groups, 
                                      data = new_SAM, permutations = 999,
                                      method = 'wunifrac',
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
  



