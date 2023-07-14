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


load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_KEGGPaths_CPM_Clinical_Data.RData")

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
   
   #describe parent call function 
   ststri <- ifelse(is.null(strata),'Null',strata)
   fostri <- as.character(x)
   #list to store results
   
   #copy model formula
   x1 <- x
   # extract left hand side of formula
   lhs <- x1[[2]]
   # extract factors on right hand side of formula 
   rhs <- x1[[3]]
   # create model.frame matrix  
   x1[[2]] <- NULL   
   rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
   
   # create unique pairwise combination of factors 
   co <- combn(unique(as.character(rhs.frame[,1])),2)
   
   # create names vector   
   nameres <- c('parent_call')
   for (elem in 1:ncol(co)){
      nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
   }
   #create results list  
   res <- vector(mode="list", length=length(nameres))
   names(res) <- nameres
   
   #add parent call to res 
   res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
   
   
   #start iteration trough pairwise combination of factors  
   for(elem in 1:ncol(co)){
      
      #reduce model elements  
      if(inherits(eval(lhs),'dist')){	
         xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                              rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
      }else{
         xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
      }
      
      mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
      
      # redefine formula
      if(length(rhs) == 1){
         xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
      }else{
         xnew <- as.formula(paste('xred' , 
                                  paste(rhs[-1],collapse= as.character(rhs[1])),
                                  sep='~'))}
      
      #pass new formula to adonis
      if(is.null(strata)){
         ad <- adonis2(xnew,data=mdat1, ... )
      }else{
         perm <- how(nperm = nperm)
         setBlocks(perm) <- with(mdat1, mdat1[,ststri])
         ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
      
      res[nameres[elem+1]] <- list(ad[1:5])
   }
   #names(res) <- names  
   class(res) <- c("pwadstrata", "list")
   return(res)
} 

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


#### Preparing Data Frames for Phylosrec

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



df <- read.csv("~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv")



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


#### BC PCOA ####
set.seed(123)
BC_Comb = vegdist(full_features,method = "bray")
PCOA = ape::pcoa(BC_Comb)
full = cbind(full_features,full_metadata)



PWA_BC_treat_NMDS <- pairwise.adonis2(BC_Comb ~  Groups, data = full, permutations = 999, by = "terms", p.adjust = "FDR")

PWA_BC_treat_NMDS_2 <- adonis2(BC_Comb ~  BMI + Age + Gender + Groups, data = full, permutations = 999,by = "terms")


PWA_Results <- 
   as.data.frame(do.call(rbind, PWA_BC_treat_NMDS)) %>% 
   tibble::rownames_to_column("Comparison") %>% 
   dplyr::filter(grepl("Groups",Comparison)) %>% 
   mutate(Comparison = gsub('.Groups',"",Comparison)) %>% 
  dplyr::rename(Pvalue = `Pr(>F)`) %>% 
  filter(Comparison == "NAFLD-O_vs_Control-0" |
           Comparison == "NAFLD-L_vs_Control-0" |
           Comparison == "NAFLD-O_vs_Control-L" |
           Comparison == "NAFLD-L_vs_Control-L" |
           Comparison == "NAFLD-O_vs_Control_Atherosclerosis" |
           Comparison == "NAFLD-L_vs_Control_Atherosclerosis" |
           Comparison == "NAFLD-O_vs_Control_Hypertension" |
           Comparison == "NAFLD-L_vs_Control_Hypertension" |
           Comparison == "NAFLD-O_vs_Overweight" |
           Comparison == "NAFLD-L_vs_Overweight" |
           Comparison == "NAFLD-O_vs_Lean_Normal" |
           Comparison == "NAFLD-L_vs_Lean_Normal" ) %>%
   mutate(FDR = p.adjust(Pvalue, method = "fdr")) #%>% 
   #filter(FDR < 0.05)


PWA_Results <- 
  as.data.frame(do.call(rbind, PWA_BC_treat_NMDS)) %>% 
tibble::rownames_to_column("Comparison") %>% 
  dplyr::filter(grepl("Groups",Comparison)) %>% 
  mutate(Comparison = gsub('.Groups',"",Comparison)) %>% 
  dplyr::rename(Pvalue = `Pr(>F)`) %>% 
  filter(Comparison == "NAFLD-O_vs_T2D_Obese" |
           Comparison == "NAFLD-L_vs_T2D_Obese" |
           Comparison == "NAFLD-O_vs_Prediabetes" |
           Comparison == "NAFLD-L_vs_Prediabetes" |
           Comparison == "NAFLD-O_vs_Atherosclerosis" |
           Comparison == "NAFLD-L_vs_Atherosclerosis" |
           Comparison == "NAFLD-O_vs_Hypertension" |
           Comparison == "NAFLD-L_vs_Hypertension" |
           Comparison == "NAFLD-O_vs_Pre_Hypertension" |
           Comparison == "NAFLD-L_vs_Pre_Hypertension" |
           Comparison == "NAFLD-O_vs_T2D_Lean" |
           Comparison == "NAFLD-L_vs_T2D_Lean" ) %>% 
  mutate(FDR = p.adjust(Pvalue, method = "fdr")) 



write.csv(PWA_Results, "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/PWA_KEGG_ALL.csv")




ggplot_df <-
   data.frame(PCoA_1 = PCOA$vectors[,1],PCoA_2 = PCOA$vectors[,2],Groups = full$Groups) %>% 
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

PCOA$values$Eigenvalues[1]/sum(PCOA$values$Eigenvalues)
PCOA$values$Eigenvalues[2]/sum(PCOA$values$Eigenvalues)


pcoa_bc_beta_div <-ggplot(ggplot_df, aes(x=PCoA_1, y=PCoA_2,colour = Groups, label = Groups)) + 
   geom_errorbar(aes(ymin=PCoA_2 - SE_2, ymax=PCoA_2 + SE_2),linetype = 1,size =1.5,width = 0.001) +
   geom_errorbar(aes(xmin=PCoA_1 - SE_1, xmax=PCoA_1 + SE_1),linetype = 1,size =1.5,width = 0.001) +
   geom_label(size = 7) + 
  scale_color_manual(values=c("#d86b6f",
                                       "#6182a8", #NAFLD
                                       
                                       "#d86b6f",
                                       "#6182a8", #Control
                                       
                                       "#FF8C00",#FFD700
                                       "#FF8C00",#FF8C00	
                                       "#FF8C00",## T2D Obese 
                                       
                                       "#32CD32",
                                       "#32CD32",#006400	
                                       "#32CD32", #Hypertension
                                       
                                       "#9966cc",	
                                       "#9966cc", #T2D Lean
                                       
                                       "#9d7029", #Atheroslcerosis
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
         axis.title=element_text(size=20)) + 
   labs(x = "PCoA 1 64%",
        y = "PCoA 2 11%")

pdf("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/KeggPaths_BC_PCOA_adjusted_Age_Gender_BMI.pdf", width = 10, height = 10)
pcoa_bc_beta_div
dev.off()




adonis_new <- function(Metadata) {
  set.seed(123)

  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]

  new_BC_Comb = vegdist(Features,method = "bray")
  new_PCOA = ape::pcoa(new_BC_Comb)
  new_full = cbind(Features,Metadata)
  

  new_adonis_Wnifrac_treat <- adonis2(new_BC_Comb ~  Age + Gender + BMI + Groups, 
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
                  # R = c(as.numeric(t1$R2[4]),
                  #       as.numeric(t2$R2[4]),
                  #       as.numeric(t3$R2[4]),
                  #       as.numeric(t4$R2[4]),
                  #       as.numeric(t5$R2[4]),
                  #       as.numeric(t6$R2[4]),
                  #       as.numeric(l1$R2[4]),
                  #       as.numeric(l2$R2[4]),
                  #       as.numeric(l3$R2[4]),
                  #       as.numeric(l4$R2[4]),
                  #       as.numeric(l5$R2[4]),
                  #       as.numeric(l6$R2[4]),
                  #       as.numeric(tm1$R2[4]),
                  #       as.numeric(tm2$R2[4]),
                  #       as.numeric(tm3$R2[4])),
                  
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
) %>% 
  mutate(FDR = p.adjust(Pvalue, method = "fdr")) 



write.csv(RES
          ,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Kegg_PWA_Adjusted.csv")


tm1  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control-0") %>% 
  adonis_new(.)

tm2 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control-L") %>% 
  adonis_new(.)
tm3 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "T2D_Obese" | Groups == "Overweight") %>% 
  adonis_new(.)

tm4 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "Prediabetes" | Groups == "Overweight") %>% 
  adonis_new(.)

tm5 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "Hypertension" | Groups == "Control_Hypertension") %>% 
  adonis_new(.)

tm6 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "Pre_Hypertension" | Groups == "Control_Hypertension") %>% 
  adonis_new(.)


tm7 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "T2D_Lean" | Groups == "Lean_Normal") %>% 
  adonis_new(.)
tm8 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "Atherosclerosis" | Groups == "Control_Atherosclerosis") %>% 
  adonis_new(.)



c(tm1$`Pr(>F)`[4],
  tm2$`Pr(>F)`[4],
  tm3$`Pr(>F)`[4],
  tm4$`Pr(>F)`[4],
  tm5$`Pr(>F)`[4],
  tm6$`Pr(>F)`[4],
  tm7$`Pr(>F)`[4],
  tm8$`Pr(>F)`[4]) %>% p.adjust(method ="fdr")









tm1  <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control-0") %>% 
  adonis_new(.)

tm2 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control-0") %>% 
  adonis_new(.)
tm3 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control_Atherosclerosis") %>% 
  adonis_new(.)

tm4 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control_Atherosclerosis") %>% 
  adonis_new(.)

tm5 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Control_Hypertension") %>% 
  adonis_new(.)

tm6 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Control_Hypertension") %>% 
  adonis_new(.)

tm7 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Lean_Normal") %>% 
  adonis_new(.)
tm8 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Lean_Normal") %>% 
  adonis_new(.)


tm9 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "Overweight") %>% 
  adonis_new(.)
tm10 <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "Overweight") %>% 
  adonis_new(.)



RES <- data.frame(Comparison = c("NAFLD O vs Control-0",
                                 "NAFLD L vs Control-0",
                                 "NAFLD O vs Control_Atherosclerosis-0",
                                 "NAFLD L vs Control_Atherosclerosis",
                                 "NAFLD O vs Control_Hypertension",
                                 "NAFLD L vs Control_Hypertension",
                                 
                                 "NAFLD O vs Lean_Normal",
                                 "NAFLD L vs Lean_Normal",
                                 "NAFLD O vs Overweight",
                                 "NAFLD L vs Overweight"),
                  
                  Pvalue = c(as.numeric(tm1$`Pr(>F)`[4]),
                             as.numeric(tm2$`Pr(>F)`[4]),
                             as.numeric(tm3$`Pr(>F)`[4]),
                             as.numeric(tm4$`Pr(>F)`[4]),
                             as.numeric(tm5$`Pr(>F)`[4]),
                             as.numeric(tm6$`Pr(>F)`[4]),
                             as.numeric(tm7$`Pr(>F)`[4]),
                             as.numeric(tm8$`Pr(>F)`[4]),
                             as.numeric(tm9$`Pr(>F)`[4]),
                             as.numeric(tm10$`Pr(>F)`[4]))
) %>% 
  mutate(FDR = p.adjust(Pvalue, method = "fdr")) 


