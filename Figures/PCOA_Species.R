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


load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")
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


##### Data Prep ####



# load("~/Downloads/Metaphlan3_Clinical_Data.RData")


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




# creating the phyloseq data
physeq <- metaphlanToPhyloseq(t(full_features))
SAM <-  sample_data(data.frame(full_metadata))
TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
physeq <- merge_phyloseq(physeq, SAM, TREE)

SAM = sample_data(physeq) %>% data.frame(.)

#### Weighted Unifrac NMD ####

set.seed(123)
full = cbind(full_features,full_metadata)
ord_w_unifrac <- ordinate(physeq, method="PCoA", distance="unifrac", weighted=TRUE)

dist_Wunifrac <- phyloseq::distance(physeq, method = "wunifrac")


PWA_Wnifrac_treat_NMD <- pairwise.adonis2(dist_Wunifrac ~  Groups, data = SAM, permutations = 999,by = "terms")
adonis_Res <- adonis2(dist_Wunifrac ~  BMI + Age + Gender + Groups, data = full, permutations = 999,by = "terms")

PWA_Results <- 
  as.data.frame(do.call(rbind, PWA_Wnifrac_treat_NMD)) %>% 
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
  mutate(FDR = p.adjust(Pvalue, method = "fdr"))


write.csv(PWA_Results, "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/PWA_Species_ALL.csv")
   
   


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
  labs(x = "PCoA1 43%",
       y = "PCoA2 9%") +
   theme(plot.title = element_text(size = 14, hjust = 0.5),  text = element_text(size = 15), 
         legend.position="none",
         axis.text.x=element_blank(), #remove x axis labels
         axis.ticks.x=element_blank(), #remove x axis ticks
         axis.text.y=element_blank(),  #remove y axis labels
         axis.ticks.y=element_blank(),  
         panel.border = element_rect(colour = "black", fill=NA, size=2),
         axis.title=element_text(size=20)
         )  



pdf("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/WUnifrac_PCOA_Species.pdf", width = 10, height = 10)
pcoa_wunifrac_beta_div
dev.off()


adonis_new <- function(Metadata) {
  set.seed(123)
  
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]
  new_physeq <- metaphlanToPhyloseq(t(Features))
  new_SAM <-  sample_data(data.frame(Metadata))
  new_TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
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
                  # R = c(as.numeric(t1$R2[4]),
                  #            as.numeric(t2$R2[4]),
                  #            as.numeric(t3$R2[4]),
                  #            as.numeric(t4$R2[4]),
                  #            as.numeric(t5$R2[4]),
                  #            as.numeric(t6$R2[4]),
                  #            as.numeric(l1$R2[4]),
                  #            as.numeric(l2$R2[4]),
                  #            as.numeric(l3$R2[4]),
                  #            as.numeric(l4$R2[4]),
                  #            as.numeric(l5$R2[4]),
                  #            as.numeric(l6$R2[4]),
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
          ,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/Species_PWA_Adjusted.csv")



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
