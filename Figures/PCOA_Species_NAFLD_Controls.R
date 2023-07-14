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
library(pairwiseAdonis)
library(devtools)




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

load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

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



#### Preparing Data Frames for Phylosrec

full_features = bind_rows(NAFLD_O,
                          NAFLD_L,
                          Control_O,
                          Control_L)

full_features[is.na(full_features)] <- 0



# install.packages("Downloads/pairwiseAdonis-master.zip", repos = NULL, type="source")

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
PWA <- pairwise.adonis2(dist_Wunifrac ~  Groups, data = SAM, permutations = 999,
                        method = 'wunifrac',by = "terms")

adonis_Wnifrac_treat <- adonis2(dist_Wunifrac ~  Age + Gender + BMI + HOMAIR + SBP  + Groups, data = SAM, permutations = 999,
                                method = 'wunifrac',by = "terms")
adonis_Wnifrac_treat

PWA_Results <- 
  as.data.frame(do.call(rbind, PWA)) %>% 
  tibble::rownames_to_column("Comparison") %>% 
  dplyr::filter(grepl("Groups",Comparison)) %>% 
  mutate(Comparison = gsub('.Groups',"",Comparison)) %>% 
  dplyr::rename(Pvalue = `Pr(>F)`) %>% 
  mutate(FDR = p.adjust(Pvalue, method = "fdr")) 

#write.csv(PWA_Results, "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/PWA_Species_ALL_NAFLD.csv")



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
   # annotate("text", x = -0.02, y = 0.024, label = paste(("NAFLD-O_vs_NAFLD-L = "),PWA$`NAFLD-O_vs_NAFLD-L`$`Pr(>F)`[1]),size=5) +
   # annotate("text", x = -0.02, y = 0.021, label = paste(("NAFLD-O_vs_Control-O = "),PWA$`NAFLD-O_vs_Control-O`$`Pr(>F)`[1]),size=5) +
   # annotate("text", x = -0.02, y = 0.018, label = paste(("NAFLD-O_vs_Control-L = "),PWA$`NAFLD-O_vs_Control-L`$`Pr(>F)`[1]),size=5) +
   # annotate("text", x = -0.02, y = 0.015, label = paste(("NAFLD-L_vs_Control-O = "),PWA$`NAFLD-L_vs_Control-O`$`Pr(>F)`[1]),size=5) +
   # annotate("text", x = -0.02, y = 0.012, label = paste(("NAFLD-L_vs_Control-L = "),PWA$`NAFLD-L_vs_Control-L`$`Pr(>F)`[1]),size=5) +
   # annotate("text", x = -0.02, y = 0.009, label = paste(("Control-O_vs_Control-L = "),PWA$`Control-O_vs_Control-L`$`Pr(>F)`[1]),size=5) +
   labs(x = "PCoA 1 44%",
        y = "PCoA 2 13%") 

pdf("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/PCOA_Species_Wunifrac_NAFLD_Controls.pdf", width = 10, height = 10)
pcoa_wunifrac_beta_div
dev.off()


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


tt <- full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "CTRL-NAFLD-L") %>% 
  adonis_new(.)

ttt <- full_metadata %>% 
  dplyr::filter(Groups == "CTRL-NAFLD-O" | Groups == "NAFLD-L") %>% 
  adonis_new(.)

tttt <- full_metadata %>% 
  dplyr::filter(Groups == "NAFLD-O" | Groups == "CTRL-NAFLD-L") %>% 
  adonis_new(.)



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

c(as.numeric(Nafld_O_vs_L$`Pr(>F)`[6]),
  as.numeric(Nafld_O_vs_Control_O$`Pr(>F)`[6]),
  as.numeric(Nafld_L_vs_Control_L$`Pr(>F)`[6]),
  as.numeric(tt$`Pr(>F)`[6]),
  as.numeric(ttt$`Pr(>F)`[6]),
  as.numeric(tttt$`Pr(>F)`[6]))%>% p.adjust(method ="fdr")


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


write.csv(RES
  ,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables/PWA_Species_NAFLD.csv")


