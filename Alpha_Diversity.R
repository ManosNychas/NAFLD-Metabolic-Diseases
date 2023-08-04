#### Alpha Index ###

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
#library(picante)
library(ggrepel)
#library(cowplot)
library(ggsignif)


load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))

Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
Control_L <- dplyr::select(Control_L, -one_of(c("Group")))



alpha_diversities <- function(df){
   list = list()
   list[["Shannon"]] = vegan::diversity(df,index = "shannon")
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
         axis.ticks = element_blank(),
         legend.key.size = unit(2, "cm"),) + 
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-O", "NAFLD-L")), map_signif_level=T, test.args=list(paired = F), tip_length = 0,textsize = 6,y_position = 4) +
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-O", "CTRL-NAFLD-O")), map_signif_level=T, test.args=list(paired = FALSE),  tip_length = 0, textsize = 6,y_position = 4.2) +
   geom_signif(test="wilcox.test", comparisons = list(c("NAFLD-L", "CTRL-NAFLD-L")), map_signif_level=T, test.args=list(paired = FALSE), tip_length = 0, textsize = 6,y_position = 4.8) +
   xlab(element_blank()) + ylab('Shannon index') +
   geom_dotplot(binaxis='y', stackdir='center', dotsize=0, binwidth = 0.02) + 
   scale_fill_manual(values=c("#d86b6f","#6182a8", "#d86b6f", "#6182a8")) +
   scale_color_manual(values=c("#d86b6f","#6182a8", "#d86b6f", "#6182a8")) +
   ggtitle("Comparison of Shannon Index ") 
Shannon_plot




