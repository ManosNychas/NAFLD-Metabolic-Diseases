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
    
    source("~/Documents/Metabolic_Diseases/ANCOM-master/scripts/ancom_v2.1.R")
    load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")
    load("~/Documents/Metabolic_Diseases/Updated_Human3/R_datasets/Sequencing_Information.RData")
    
    
    # load("~/Downloads/Metaphlan3_Clinical_Data.RData")
    
    NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
    NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))
    
    Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
    Control_L <- dplyr::select(Control_L, -one_of(c("Group")))
    
    
    #### Metagenome Seq ####
    
    ANCOM_II_difference_abundance <- function(Dataset_1,Dataset_2,prevCutoff = 0.1) {
       # # 
       # Dataset_1 = NAFLD_O
       # Dataset_2 = Control_O
       # 
       Counts_1 <- 
         colSums(Dataset_1!= 0) %>% 
         as.data.frame() %>% 
         tibble::rownames_to_column() %>% 
         dplyr::rename(taxa_id = "rowname", Count_Dataset_1 = ".")
       
       Counts_2 <- 
         colSums(Dataset_2!= 0) %>% 
         as.data.frame() %>% 
         tibble::rownames_to_column() %>% 
         dplyr::rename(taxa_id = "rowname", Count_Dataset_2 = ".")
         
       
       
      dataset_name_1 = deparse(substitute(Dataset_1))
      dataset_name_2 = deparse(substitute(Dataset_2))
      
      Dataset_1 = Dataset_1[ ,colSums(Dataset_1!= 0) >= nrow(Dataset_1)*0.1]
      Dataset_2 = Dataset_2[ ,colSums(Dataset_2!= 0) >= nrow(Dataset_2)*0.1]

      ### Setting the sets
      full_features = bind_rows(Dataset_1, Dataset_2)
      full_features[is.na(full_features)] <- 0
      
      ## setting the metadata samples that exist in the current comparison
      full_metadata  <-
        Combined_metadata  %>% 
        .[which(rownames(.) %in% rownames(full_features)),] %>% 
        dplyr::select(Age,Gender,BMI,HOMAIR,SBP) %>% 
        mutate(Groups = ifelse(rownames(.) %in% rownames(Dataset_1), dataset_name_1,dataset_name_2))
      
      
      ### Remove NAs
      full_metadata = full_metadata[order(rownames(full_metadata)),] %>% .[complete.cases(.), ]
      full_features = full_features[order(rownames(full_features)),] %>% .[which(rownames(.) %in% rownames(full_metadata)),]
      
      full_metadata<- full_metadata %>% tibble::rownames_to_column(., "Sample.ID")
    
      full_sequencing <- All_Sequencing %>% .[order(.$Samples),] %>%
        .[which(.$Samples %in% rownames(full_features)),]
      
      #### attempting to transform the species data to counts, bye multiplying with a big number
      #full_features = full_features*full_sequencing$Reads
      full_features = round(full_features*10^6)
      full_features <- t(full_features) %>% as.data.frame(.)
      
      feature_table = full_features;
      meta_data = full_metadata;
      sample_var = "Sample.ID";
      group_var = "Groups";
      out_cut = 0.05;
      zero_cut = 0.95;
      lib_cut = 0
      neg_lb = FALSE
      prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var,
                                         out_cut, zero_cut, lib_cut, neg_lb)
      feature_table = prepro$feature_table # Preprocessed feature table
      meta_data = prepro$meta_data # Preprocessed metadata
      struc_zero = NULL # Structural zero info
      main_var = "Groups"
      p_adj_method = "fdr"
      alpha = 0.2
      adj_formula = c("Age","Gender","BMI","HOMAIR","SBP")
      rand_formula = NULL
      
      t_start = Sys.time()
      res = ANCOM(feature_table, meta_data,struc_zero, main_var, p_adj_method,
                  alpha, adj_formula, rand_formula)
      t_end = Sys.time()
      t_run = t_end - t_start
      
     
     clr_table = apply(as.matrix(full_features) + 1, 2, clr)
     eff_size = apply(clr_table, 1, function(y) 
       lm(y ~ x, data = data.frame(y = y, 
                                   x = meta_data %>% pull(main_var),
                                   check.names = FALSE))$coef[-1]) %>% 
       as.data.frame() %>% 
       rownames_to_column() %>% 
       dplyr::rename(taxa_id = "rowname", CLR_Mean = "." ) 
       
     
     # Calculate clr mean difference

      results <- res$out %>% 
        arrange(.,taxa_id) %>% 
        mutate(.,CLR_Mean_Difference = (res$fig$data %>% arrange(.,taxa_id))$x) %>% 
        dplyr::left_join(Counts_1, by = "taxa_id") %>% 
        dplyr::left_join(Counts_2, by = "taxa_id") %>% 
        mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>% 
        dplyr::left_join(eff_size, by = "taxa_id") %>% 
        mutate(CLR_Mean_Difference = CLR_Mean) %>% 
        dplyr::select(-CLR_Mean) %>% 
        mutate(Species_Direction = ifelse(CLR_Mean_Difference > 0,
                                          paste(taxa_id,"Up",sep = "_"),
                                          paste(taxa_id,"Down",sep = "_"))) 
      
      return(results)
    }  
    
    
    O_NvC <- 
      ANCOM_II_difference_abundance(NAFLD_O,Control_O,All_Sequencing)  
    
    L_NvC <- 
      ANCOM_II_difference_abundance(NAFLD_L,Control_L,All_Sequencing) 
    
    
    Results_ANCOM_II <- 
      rbind(O_NvC,L_NvC) %>% 
      cbind(c(rep("NAFLD-O vs Control-O",nrow(O_NvC)),
              rep("NAFLD-L vs Control-L",nrow(L_NvC))),.) %>% 
      dplyr::rename(Species = "taxa_id")
      
      colnames(Results_ANCOM_II)[1] = "Comparison"
      
      saveRDS(Results_ANCOM_II,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Differential_Abundance/Species/New_Adjusts/Detailed_ANCOM_II_Res.rds")

      
      Results_ANCOM_II <- 
        readRDS("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Differential_Abundance/Species/New_Adjusts/Detailed_ANCOM_II_Res.rds") %>% 
        filter(detected_0.6 == "TRUE") %>% 
        mutate(Species = gsub("_"," ",Species))
      
    library(reshape2)
    Pvalue_DF <- 
      dcast(Results_ANCOM_II, Comparison ~ Species, value.var = "detected_0.7") %>% 
      tibble::column_to_rownames("Comparison") %>% 
      set_rownames(c("Lean","Overweight/Obese"))
    
    Group_Impact_DF <- dcast(Results_ANCOM_II, Comparison ~ Species, value.var = "CLR_Mean_Difference") %>% 
      tibble::column_to_rownames("Comparison") %>% 
      as.matrix(.) %>% 
      apply(.,2,as.numeric) %>% 
      set_rownames(c("Lean","Overweight"))
    


library(circlize)
col_fun = colorRamp2(c(-1.5,-0.75, 0, 0.75, 1.5),  c("#245b87", "#58b7e0","#FFFFFF","#C24933","#AF2912"))

library(ComplexHeatmap)

Healthy_heatmap <-Heatmap(t(Group_Impact_DF),
                          name = "NAFLD vs Control Comparisons",
                          row_title = " Singificantly Different Species",
                          column_title = "ANCOM II - adjusted by 
                          Age + Gender + BMI + HOMAIR + SBP",
                          col = col_fun,
                          border = TRUE,
                          column_order = c("Lean","Overweight"),
                          width = ncol(t(Group_Impact_DF))*unit(15, "mm"),
                          height = nrow(t(Group_Impact_DF))*unit(2.4, "mm"),
                          cluster_rows = F,
                          cluster_columns = F, 
                          show_heatmap_legend = FALSE,
                          column_names_gp = gpar(fontsize = 10), # column color 
                          row_names_gp = gpar(fontsize = 7.5, fontface = "italic"),
)
Healthy_heatmap

# fdr legend
# lgd_pval = Legend(pch = 1,type = "points", labels = "Q value < 0.05",
#                   legend_gp = gpar(col = "white",lwd=2),background = "black",title = " ")

# legend correlation
lgd_correlation = Legend(col_fun = col_fun,title = "Group Impact", 
                         border = "black",direction = "horizontal",
                         legend_width = unit(4.8, "cm") ,title_position = "topcenter")

pd <- packLegend(#lgd_pval, 
                 lgd_correlation, 
                 max_width = unit(26, "cm"), 
                 column_gap = unit(8, "mm"), 
                 row_gap = unit(1, "cm"),
                 direction = "horizontal")


pdf('Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Differential_Abundance/Species/New_Adjusts/ANCOM_II_Heatmap.pdf',height = 10,width = 6)
# par(mar=c(5,6,4,1)+.1)
ComplexHeatmap::draw(Healthy_heatmap,
                     annotation_legend_list = pd,
                     heatmap_legend_side = "bottom",
                     merge_legend = TRUE)
dev.off()


library(ggVennDiagram)
library(ggvenn)



O_NvC <- Results_ANCOM_II %>% filter(Comparison == "NAFLD-O vs Control-O")

L_NvC <- Results_ANCOM_II %>% filter(Comparison == "NAFLD-L vs Control-L")



list <- list(Overweight = O_NvC %>% filter(detected_0.6 == TRUE) %>% pull(Species_Direction),
             Lean = L_NvC %>% filter(detected_0.6 == TRUE) %>% pull(Species_Direction))

venn <- Venn(list)
data <- process_data(venn)


new_ven <- ggvenn(
  list, 
  #fill_color = c("#04b18f","#df6868"),
  fill_color = c("#E41C06","#3366cc"),
  stroke_size = 1,
  set_name_size = 8,
  stroke_linetype = "solid",
  show_percentage = FALSE,
  text_size = 7
)


ggsave("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/NewSpecies_ANCOM_II_Ven.pdf",
       new_ven, width = 10, height = 10)


MS_Ven <- ggVennDiagram(list, color = "black", label_alpha = 0, label = "both",
                        label_percent_digit = 2, edge_lty = 0, set_size = 0,
                        label_size = 15) +
  geom_sf(size = 1, lty = "solid", color = "#b5b5b5", data = venn_setedge(data), show.legend = F) +
  scale_fill_gradient(low = "#6182a8", high = "#d86b6f", name = "Compounds") +
  #ggtitle("ANCOM II - adjusted by  Age + Gender + BMI + HOMAIR + SBP") +
  theme_void() + 
  theme(plot.title = element_text(size = 20, hjust = 0.5),  text = element_text(size = 20), 
      legend.title = element_blank() , axis.text.x=element_blank(), axis.ticks = element_blank(),
      legend.key.size = unit(1.5, "cm"),legend.position="bottom") 

ggsave("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Species_ANCOM_II_Ven.pdf",
       MS_Ven, width = 10, height = 10)

