
#### Beta Diversity ###
rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyverse)
library(RANN)
library(magrittr)
library(vegan)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(data.table)
library(ggplot2)
library(compositions)
library(nlme)
library(VennDiagram)


########### metadata mean sd#######

load("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

library(abind)
library(tidyverse)

nafld_vs_controls = bind_rows(NAFLD_L,
                              NAFLD_O,
                              Control_L,
                              Control_O,
                              Overweight,
                              T2D_Obese,
                              Prediabetes,
                              Hypertension,
                              Pre_Hypertension,
                              Control_Hypertension,
                              T2D_Lean,
                              Lean_Normal,
                              Atherosclerosis,
                              Control_Atherosclerosis) %>% 
   dplyr::select(-Group)


library(readxl)
df <- read_excel("Downloads/Full_Metadata_Metastudy_DONE.xlsx",col_names = T,
                 na = "NA")

mean_sd <- function(dataset) {
   man<-df[which((df$X) %in% rownames(dataset)),] %>% 
      apply(.,2,as.numeric) %>% apply(.,2,function(x) mean(x,na.rm=T)) 
   
   B<-df[which((df$X) %in% rownames(dataset)),] %>%
      apply(.,2,as.numeric) %>% apply(.,2,function(x) sd(x,na.rm=T)) 
   
   return(round(abind(man,B,along=2),2) %>% apply(.,1,function(x)paste(x[1],"\u00B1",as.character(x[2]))))
   
}

test_ = rbind(mean_sd(NAFLD_L),
              mean_sd(NAFLD_O),
              mean_sd(Control_L),
              mean_sd(Control_O),
              mean_sd(T2D_Obese),
              mean_sd(Prediabetes),
              mean_sd(Overweight),
              mean_sd(Hypertension),
              mean_sd(Pre_Hypertension),
              mean_sd(Control_Hypertension),
              mean_sd(T2D_Lean),
              mean_sd(Lean_Normal),
              mean_sd(Atherosclerosis),
              mean_sd(Control_Atherosclerosis)) %>% 
   as.data.frame()




rownames(test_) = c("NAFLD_L",
                    "NAFLD_O",
                    "Control_L",
                    "Control_O",
                    "T2D_Obese",
                    "Prediabetes",
                    "Overweight",
                    "Hypertension",
                    "Pre_Hypertension",
                    "Control_Hypertension",
                    "T2D_Lean",
                    "Lean_Normal",
                    "Atherosclerosis",
                    "Control_Atherosclerosis")
llel<-df[which((df$X) %in% rownames(NAFLD_O)),] 


se(llel$CK.other,na.rm=T)
   



   
data.frame(llel$X,llel$TF3.other)

test_ <-
   apply(test_,2,function(x) gsub("NaN ± NA","NA",x)) 

write.csv(test_, file = "Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Tables_all_mean_sd.csv")

