rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
pacman::p_load(devtools,
               tidyverse,
               BiocManager,
               tibble,
               gdata,
               ggfortify,
               ggpubr,
               cowplot,
               ggsci,
               gplots,
               RColorBrewer,
               data.table,
               openxlsx,
               readxl,
               readr,
               Hmisc,
               corrr,
               magrittr,
               circlize, 
               tidyverse, 
               ggsci, 
               readxl,
               openxlsx,
               janitor,
               ppcor,
               compositions,
               install = TRUE)

calulo_corr <- function(cdata_matrix, spe_matrix){
   
   pval_ <- data.frame(matrix(NA, nrow = ncol(spe_matrix), ncol = ncol(cdata_matrix) ))
   colnames(pval_) <- colnames(cdata_matrix)
   rownames(pval_) <- colnames(spe_matrix)
   cor_ <- pval_
   for (i in colnames(cdata_matrix)) {
      for (j in colnames(spe_matrix )) {
         all <- na.omit(data.frame(spe_matrix[, j], cdata_matrix[,i]))
         tryCatch({
            #ctest <- pcor.test(all[, 1], all[,2], all[,3:ncol(all)], method = 'spearman' )
            ctest <- cor.test(all[, 1], all[,2],  method = 'spearman')
            pval <- ctest$p.value
            cor <- ctest$estimate
            pval_[j,i] <- pval
            cor_[j,i] <- cor
         }, error=function(e){})
      }
   }
   return(list(pvalue=pval_, correlation=cor_))
}
load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")

dataset <- bind_rows(NAFLD_L,
                     NAFLD_O,
                     T2D_Obese,
                     Prediabetes,
                     Hypertension,
                     Pre_Hypertension,
                     T2D_Lean,
                     Atherosclerosis)

metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(-c(FLI.liver_function,
                    GGT.liver_function,
                    ALT.liver_function,
                    AST.liver_function,
                    Triglyceride.lipid_profiles,
                    Liver_Fat.liver_function)) %>% 
   dplyr::select(-c(67:95)) %>% 
   data.frame() %>% 
   .[which(rownames(.) %in% rownames(dataset)),] %>% 
   .[order(rownames(.)),] %>%
   as.data.frame() %>% 
  dplyr::select(-contains(".other")) %>% 
  dplyr::select(-contains(".liver_function")) %>% 
  dplyr::select(-contains("Age")) %>% 
  dplyr::select(-c("Height.body_weight_variables",
                   "Weight.body_weight_variables",
                   "Hip.body_weight_variables")) 
  
  
  ## fixing order
metadata_2 <-
   bind_cols( metadata %>% dplyr::select(contains(".body_weight_variables")),
              metadata %>% dplyr::select(contains(".cardiac_variables")),
              metadata %>% dplyr::select(contains(".glucose_insulin")),
              metadata %>% dplyr::select(contains(".lipid_profiles")),
              metadata %>% dplyr::select(contains(".kidney_function"))) %>% 
  setNames(gsub(".body_weight_variables", "", names(.))) %>% 
  setNames(gsub(".cardiac_variables", "", names(.))) %>% 
  setNames(gsub(".glucose_insulin", "", names(.))) %>% 
  setNames(gsub(".lipid_profiles", "", names(.))) %>% 
  setNames(gsub(".kidney_function", "", names(.))) %>% 
  dplyr::rename('Systolic' = Systolic_pressure,
                'Diastolic' = Diastolic_pressure,
                "LPA" = Lpa,
                "CR" = Cr) %>% 
  setNames(gsub("HOMAIR", "HOMA-IR", names(.)))
   


Primary_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(FLI.liver_function,
                 GGT.liver_function,
                 ALT.liver_function,
                 AST.liver_function,
                 Triglyceride.lipid_profiles,
                 Liver_Fat.liver_function) %>% 
   .[which(rownames(.) %in% rownames(dataset)),] %>% 
   .[order(rownames(.)),] %>%
   as.data.frame() %>% 
  setNames(gsub(".lipid_profiles", "", names(.))) %>% 
  setNames(gsub(".liver_function", "", names(.))) %>% 
  dplyr::rename('Liver Fat' = Liver_Fat)
  
  

spe_matrix <- Primary_metadata
cdata_matrix <- metadata_2

pval_ <- data.frame(matrix(NA, nrow = ncol(spe_matrix), ncol = ncol(cdata_matrix) ))
colnames(pval_) <- colnames(cdata_matrix)
rownames(pval_) <- colnames(spe_matrix)
cor_ <- pval_
for (i in colnames(cdata_matrix)) {
   
   for (j in colnames(spe_matrix )) {
      
      all <- na.omit(data.frame(spe_matrix[, j], cdata_matrix[,i]))
      tryCatch({
         #ctest <- pcor.test(all[, 1], all[,2], all[,3:ncol(all)], method = 'spearman' )
         ctest <- cor.test(all[, 1], all[,2],  method = 'spearman')
         pval <- ctest$p.value
         cor <- ctest$estimate
         pval_[j,i] <- pval
         cor_[j,i] <- cor
      }, error=function(e){})
   }
}

adjpval_<- 
   matrix(p.adjust(as.vector(as.matrix(pval_)), method='fdr'),ncol=ncol(pval_)) %>% 
   `rownames<-`(rownames(cor_)) %>% 
   `colnames<-`(colnames(cor_))
   

library(circlize)
col_fun = colorRamp2(c(-1,-0.5, 0, 0.5, 1),  c("#245b87", "#58b7e0","white","#C24933","#AF2912"))

library(ComplexHeatmap)

Healthy_heatmap <-Heatmap((cor_),
                          name = "NAFLD vs Control Comparisons",
                          col = col_fun,
                          border = TRUE,
                          width = ncol((cor_))*unit(6.5, "mm"),
                          height = nrow((cor_))*unit(4, "mm"),
                          # # label the circle of fdr
                          layer_fun = function(j, i, x, y, width, height, fill) {
                             ind_mat = restore_matrix(j, i, x, y)
                             v = pindex((as.matrix(pval_)), i, j) # pvalue is the matrix : qvalue table
                             l = v < 0.05
                             ind_mat_vec = as.matrix(ind_mat)
                             ind = ind_mat_vec[l]
                             grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "gray27",lwd=1))
                             
                             v = pindex((as.matrix(adjpval_)), i, j) # qvalue is the matrix : qvalue table
                             l = v < 0.05
                             ind_mat_vec = as.matrix(ind_mat)
                             ind = ind_mat_vec[l]
                             grid.points(x[ind], y[ind], pch = 16, size = unit(2, "mm"), gp = gpar(col = "gray27",lwd=1))
                          },
                          cluster_rows = F,
                          cluster_columns = F, 
                          show_heatmap_legend = FALSE,
                          column_names_gp = gpar(fontsize = 12,
                                              #fontface = "italic",
                                              col = c(rep("#5e8656", 9),
                                                      rep("#4f3372", 2),
                                                      rep("#a90000", 9),
                                                      rep("#095f95", 10),
                                                      rep("#da813b", 2))),
                          column_names_rot = 45
)
Healthy_heatmap




# fdr legend
lgd_fdr = Legend(pch = 16, type = "points",labels = "FDR Q < 0.05", 
                 legend_gp = gpar(col = "white",lwd=2),background = "black",title = " ")
# pval legend
lgd_pval = Legend(pch = 1,type = "points", labels = "P value < 0.05",
                  legend_gp = gpar(col = "white",lwd=2),background = "black",title = " ")
# legend correlation
lgd_correlation = Legend(col_fun = col_fun,title = "Spearman Correlation \n Coeficient", 
                         border = "black",direction = "horizontal",
                         legend_width = unit(4.8, "cm") ,title_position = "topcenter")

pd <- packLegend(lgd_fdr, lgd_pval, 
                 lgd_correlation, 
                 max_width = unit(26, "cm"), 
                 column_gap = unit(8, "mm"), 
                 row_gap = unit(1, "cm"),
                 direction = "vertical")

pdf('Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure_1/Cor_Metadata_NAFLD_metadata.pdf',height = 10,width = 12)
# par(mar=c(5,6,4,1)+.1)
ComplexHeatmap::draw(Healthy_heatmap,annotation_legend_list = pd,heatmap_legend_side = "top", merge_legend = TRUE)
dev.off()


