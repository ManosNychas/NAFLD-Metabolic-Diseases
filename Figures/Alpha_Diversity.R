#### Alpha Index ###
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


load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

# load("~/Downloads/Metaphlan3_Clinical_Data.RData")




NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))

Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
Control_L <- dplyr::select(Control_L, -one_of(c("Group")))


TREE <- read_tree("~/Documents/Metabolic_Diseases/Updated_Human3/Metaphlan_3_parsed_species_only.nwk")

alpha_diversities <- function(df){
   list = list()
   list[["Shannon"]] = vegan::diversity(df,index = "shannon")
   list[["Simpson"]] = vegan::diversity(df,index = "simpson")
   list[["Chao"]] = apply(df, 1, fossil::chao1)
   list[["PD"]] = pd(df,tree = TREE) %>% dplyr::select(., PD)
   list[["ACE"]] = apply(df,1,ACE) 
   list[["Inv.Simpson"]] = vegan::diversity(df,index = "invsimpson")
   return(list)
}




Alpha_Div_NAFLD_O = alpha_diversities(NAFLD_O)
Alpha_Div_NAFLD_L= alpha_diversities(NAFLD_L)
Alpha_Div_Control_O = alpha_diversities(Control_O)
Alpha_Div_Control_L = alpha_diversities(Control_L)




Shannon = data.frame(Shannon_Index = c(Alpha_Div_NAFLD_O$Shannon,
                                       Alpha_Div_NAFLD_L$Shannon,
                                       Alpha_Div_Control_O$Shannon,
                                       Alpha_Div_Control_L$Shannon),
                     Group = factor(c(rep("NAFLD-O",length(Alpha_Div_NAFLD_O$Shannon)),
                                      rep("NAFLD-L",length(Alpha_Div_NAFLD_L$Shannon)),
                                      rep("CTRL-NAFLD-O",length(Alpha_Div_Control_O$Shannon)),
                                      rep("CTRL-NAFLD-L",length(Alpha_Div_Control_L$Shannon))),
                                    levels = c("NAFLD-O","NAFLD-L","CTRL-NAFLD-O","CTRL-NAFLD-L")))







Shannon_plot<- ggplot(Shannon, aes(x = Group, y = Shannon_Index)) + 
   geom_boxplot( alpha=0.25, aes(fill = Group, colour=Group), outlier.shape = NA) +
   theme_classic()  +
   theme(plot.title = element_text(size = 14, hjust = 0.5),
         text = element_text(size = 20), 
         legend.title = element_blank() , 
         legend.position = "none",
         #axis.text.x=element_blank(),
         axis.ticks = element_blank(),
         legend.key.size = unit(2, "cm"),) + 
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-O", "NAFLD-L")), map_signif_level=T, test.args=list(paired = F), tip_length = 0,textsize = 6,y_position = 4) +
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-O", "CTRL-NAFLD-O")), map_signif_level=T, test.args=list(paired = FALSE),  tip_length = 0, textsize = 6,y_position = 4.2) +
   #geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-O", "Control-L")), map_signif_level=T, test.args=list(paired = F), tip_length = 0 ,textsize = 6,y_position = 4.4) +
   #geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-L", "Control-O")), map_signif_level=T, test.args=list(paired = FALSE), tip_length = 0, textsize =6,y_position = 4.6) +
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-L", "CTRL-NAFLD-L")), map_signif_level=T, test.args=list(paired = FALSE), tip_length = 0, textsize = 6,y_position = 4.8) +
   #geom_signif(test="wilcox.test", comparisons = list(c("Control-O", "Control-L")), map_signif_level=T, test.args=list(paired = FALSE), tip_length = 0, textsize = 6,y_position = 5) +
   xlab(element_blank()) + ylab('Shannon index') +
   geom_dotplot(binaxis='y', stackdir='center', dotsize=0, binwidth = 0.02) + 
   scale_fill_manual(values=c("#d86b6f","#6182a8", "#d86b6f", "#6182a8")) +
   scale_color_manual(values=c("#d86b6f","#6182a8", "#d86b6f", "#6182a8")) +
   ggtitle("Comparison of Shannon Index ") 
Shannon_plot



ggsave("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Supplemental_Figures/Alpha_Shannon_Index.pdf", 
       Shannon_plot, 
       width = 10, height = 10)


