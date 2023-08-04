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
library(GUniFrac)
library(textshape)

#Metaphlan to phylosec function
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

#Loading data

load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


##### Data Prep ####


NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))

Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
Control_L <- dplyr::select(Control_L, -one_of(c("Group")))



#### Preparing Data Frames for Phylosrec

full_features = bind_rows(NAFLD_O,
                          NAFLD_L,
                          Control_O,
                          Control_L)

full_features[is.na(full_features)] <- 0



full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]


full_metadata = full_metadata[order(match(rownames(full_metadata),rownames(full_features))),] 

full_metadata$Groups = factor(c(rep("NAFLD-O",nrow(NAFLD_O)),
                                rep("NAFLD-L",nrow(NAFLD_L)),
                                rep("CTRL-NAFLD-O",nrow(Control_O)),
                                rep("CTRL-NAFLD-L",nrow(Control_L))),
                              levels =c("NAFLD-O","NAFLD-L","CTRL-NAFLD-O","CTRL-NAFLD-L"))





full_metadata =  full_metadata[complete.cases(full_metadata), ]
full_features =  full_features[which(rownames(full_features) %in% rownames(full_metadata)),]


# creating the phyloseq data
physeq <- metaphlanToPhyloseq(t(full_features))
SAM <-  sample_data(data.frame(full_metadata))
TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
physeq <- merge_phyloseq(physeq, SAM, TREE)

SAM = sample_data(physeq) %>% data.frame(.)

#### Weighted Unifrac PCOA ####

set.seed(123)
ord_w_unifrac <- ordinate(physeq, method="PCoA", distance="unifrac", weighted=TRUE)
dist_Wunifrac <- phyloseq::distance(physeq, method = "wunifrac")


adonis_Wnifrac_treat <- adonis2(dist_Wunifrac ~  Age + Gender + BMI + HOMAIR + SBP  + Groups, data = SAM, permutations = 999,
                                method = 'wunifrac',by = "terms")
adonis_Wnifrac_treat

## GG plot
ggplot_df <-
   data.frame(PCoA_1 = ord_w_unifrac$vectors[,1],
              PCoA_2 = ord_w_unifrac$vectors[,2],
              Groups = SAM$Groups) %>% 
   group_by(Groups) %>% 
   summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) %>%
   rename(.,PCoA_1 = PCoA_1_mean,PCoA_2 = PCoA_2_mean,SE_1 = PCoA_1_se,SE_2 = PCoA_2_se)

ord_w_unifrac$values$Eigenvalues[1]/sum(ord_w_unifrac$values$Eigenvalues)
ord_w_unifrac$values$Eigenvalues[2]/sum(ord_w_unifrac$values$Eigenvalues)


pcoa_wunifrac_beta_div <-ggplot(ggplot_df, aes(x=PCoA_1, y=PCoA_2,colour = Groups, label = Groups)) + 
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
         axis.title=element_text(size=20)
   ) +
   annotate("text", x = -0.02, y = 0.03, label = paste(("P value = "),adonis_Wnifrac_treat$`Pr(>F)`[6]),size=8) +
   annotate("text", x = -0.02, y = 0.027, label = paste(("R² = 0.04")),size=8,) +
   labs(x = "PCoA 1 44%",
        y = "PCoA 2 13%") 


##Adonis Function for each single group vs another single group to test beta diversity between them

adonis_new <- function(Metadata) {
  
  Metadata = full_metadata %>% 
    dplyr::filter(Groups == "NAFLD-O" | Groups == "NAFLD-L")
  Features =  full_features[which(rownames(full_features) %in% rownames(Metadata)),]
  new_physeq <- metaphlanToPhyloseq(t(Features))
  new_SAM <-  sample_data(data.frame(Metadata))
  new_TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")
  new_physeq <- merge_phyloseq(new_physeq, new_SAM, new_TREE)
  
  new_SAM = sample_data(new_physeq) %>% data.frame(.)
  
  new_dist_Wunifrac <- phyloseq::distance(new_physeq, method = "wunifrac")
  set.seed(123)
  new_adonis_Wnifrac_treat <- adonis2(new_dist_Wunifrac ~  Age + Gender + BMI + HOMAIR + SBP  + Groups, 
                                      data = new_SAM, permutations = 999,
                                      method = 'wunifrac',
                                      by = "terms")
  
  return(new_adonis_Wnifrac_treat)
}





Nafld_O_vs_L <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "NAFLD-L") %>% 
  adonis_new(.)

Nafld_O_vs_Control_O <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "CTRL-NAFLD-O") %>% 
  adonis_new(.)

Nafld_L_vs_Control_L <- 
  full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-L" | Groups == "CTRL-NAFLD-L") %>% 
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



