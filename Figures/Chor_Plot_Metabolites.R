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

calulo_corr_adjusted <- function(var_adjust, cdata_matrix, spe_matrix ){
   
   pval_ <- data.frame(matrix(NA, nrow = ncol(spe_matrix), ncol = ncol(cdata_matrix) ))
   colnames(pval_) <- colnames(cdata_matrix)
   rownames(pval_) <- colnames(spe_matrix)
   cor_ <- pval_
   for (i in colnames(cdata_matrix)) {
      for (j in colnames(spe_matrix )) {
         all <- na.omit(data.frame(spe_matrix[, j], cdata_matrix[,i], var_adjust))
         tryCatch({
            ctest <- pcor.test(all[, 1], all[,2], all[,3:ncol(all)], method = 'spearman' )
            #ctest <- cor.test(all[, 1], all[,2],  method = 'spearman' )
            pval <- ctest$p.value
            cor <- ctest$estimate
            pval_[j,i] <- pval
            cor_[j,i] <- cor
         }, error=function(e){})
      }
   }
   return(list(pvalue=pval_, correlation=cor_))
}
correlation_format <- function(data_1,data_2,Group_1,Group_2,adjusted) {
   #  
   #  data_1 = KeggMods_NAFLD_O_L
   #  data_2 = Primary_metadata
   #  
   # 
   # adjusted = TRUE
   # 
   
   
   data_2 <- 
      data_2[which(rownames(data_2) %in% rownames(data_1)),]
   
   data_2 = data_2[order(rownames(data_2)),]
   data_1 = data_1[order(rownames(data_1)),]
   # var_adjust = var_adjust[order(rownames(var_adjust)),]
   
   
   
   corr_pval_data = data.frame()
   if(adjusted == TRUE) {
      corr_pval_data <- calulo_corr_adjusted(var_adjust = var_adjust, cdata_matrix = data_2 , spe_matrix = data_1)
   } else {
      corr_pval_data <- calulo_corr(cdata_matrix = data_2 , spe_matrix = data_1)
   }
   
   corr_mat <-
      corr_pval_data$correlation %>% 
      rownames_to_column(var = "rownames") %>% 
      reshape2::melt(.,id.vars = c("rownames")) %>% 
      dplyr::rename(from = rownames, to = variable, correlation = value) 
   
   pmat <- 
      corr_pval_data$pvalue %>% 
      rownames_to_column(var = "rownames") %>% 
      reshape2::melt(.,id.vars = c("rownames")) %>% 
      dplyr::rename(from = rownames, to = variable, pvalue = value) %>% 
      mutate(adj_pvalue = p.adjust(pvalue, method = "fdr"))
   
   new_corr_mat <- 
      corr_mat %>% 
      left_join(pmat, by = c("from","to")) %>% 
      mutate(absolute_correlation = abs(correlation)) %>% 
      #filter(pvalue <= 0.05 &  absolute_correlation > 0) %>% 
      mutate(from_class = rep(Group_1,nrow(.))) %>% 
      mutate(to_class = rep(Group_2,nrow(.))) %>% 
      dplyr::select(from,from_class,to,to_class,correlation,absolute_correlation,pvalue,adj_pvalue)
   
   return(new_corr_mat)
}

##### Loading Species - KeggPaths - Metabolites - Metadata  Data per Group #####

load("~/Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


List_Signif_Metabolites <- 
   read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/update_NAFLD_up_down.csv", sep = "\t") %>% 
   as_tibble() %>% 
   mutate(Direction = ifelse(Direction == "up",1,0)) %>% 
   dplyr::rename(Species = "Sig_dif") %>% 
   mutate(.,Species_Direction = ifelse(Direction == 1,paste(Species,"Up",sep = "_"),paste(Species,"Down",sep = "_"))) %>%
   mutate(Comparison = ifelse(Comparison == "O","NAFLD-O vs Control-O", "NAFLD-L vs Control-L")) %>% 
   dplyr::left_join(read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/p_value_significant_metabolites_csv.csv", sep = ",") %>% 
                       dplyr::rename(Species = "Metabolites") %>% 
                       filter(Category == "NAFLD_L" | Category == "NAFLD_O" ) %>% 
                       mutate(Category  = gsub("NAFLD_O","NAFLD-O vs Control-O",Category)) %>% 
                       mutate(Category  = gsub("NAFLD_L","NAFLD-L vs Control-L",Category)) %>% 
                       dplyr::rename(Comparison = "Category"), by = c("Species","Comparison")) 

Signif_Metabolites_O <- 
   List_Signif_Metabolites %>% 
   filter(Comparison == "NAFLD-O vs Control-O")

Signif_Metabolites_L <- 
   List_Signif_Metabolites %>% 
   filter(Comparison == "NAFLD-L vs Control-L")


Commorbidities_Metabolites <- 
   read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/direction_All_Comp_significant_metabolites_csv.csv") %>% 
   mutate(Direction = ifelse(Direction == "up",1,0)) %>% 
   filter(Comparison != "NAFLD_L" & Comparison != "NAFLD_O") %>%  
   dplyr::rename(from = 'Metabolites') %>% 
  mutate(.,Species_Direction = ifelse(Direction == 1,paste(from,"Up",sep = "_"),paste(from,"Down",sep = "_"))) %>% 
  filter(Comparison != "Diabetes_overweight" & Comparison != "Hypertension")



Commorbidities_Metabolites_1 <- 
  read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/OLR_to_Lasso_significnat_metabolites.csv",sep = "\t") %>% 
  dplyr::rename(Direction = 'direction') %>% 
  mutate(Direction = ifelse(Direction == "up",1,0)) %>% 
  dplyr::rename(from = 'Metabolite') %>% 
  mutate(.,Species_Direction = ifelse(Direction == 1,paste(from,"Up",sep = "_"),paste(from,"Down",sep = "_"))) %>% 
  mutate(Comparison  = gsub("Hypertension","Hypertension vs Control",Comparison)) %>% 
  mutate(Comparison  = gsub("Prehypertension","Pre_Hypertension vs Control",Comparison)) %>% 
  mutate(Comparison  = gsub("Prediabetes","Prediabetes vs Overweight",Comparison)) %>% 
  mutate(Comparison  = gsub("T2D_Overweight","T2D_Obese vs Overweight",Comparison)) 



Commorbidities_Metabolites <-
  bind_rows(Commorbidities_Metabolites,
            Commorbidities_Metabolites_1) 


KeggPaths_NAFLD_O <-
   bind_rows(
      NAFLD_O %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1] ,
             Control_O %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1]) %>% 
   replace(is.na(.), 0)%>% 
   dplyr::select(., -one_of(c("Group"))) %>% 
   mutate_if(is.numeric, ~round(.*10^6)) %>% 
   zCompositions::cmultRepl() %>%
   clr()  %>% 
   .[order(rownames(.)),] 

KeggPaths_NAFLD_L <-
   bind_rows(
      NAFLD_L %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1] ,
      Control_L %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1]) %>% 
   replace(is.na(.), 0)%>% 
   dplyr::select(., -one_of(c("Group"))) %>% 
   mutate_if(is.numeric, ~round(.*10^6)) %>% 
   zCompositions::cmultRepl() %>%
   clr()  %>% 
   .[order(rownames(.)),] 



KeggPaths_NAFLD_All <-
   bind_rows(
      NAFLD_O %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1] ,
      Control_O %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1],
      NAFLD_L %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1] ,
      Control_L %>%  .[ ,colSums(.!= 0) >= nrow(.)*0.1]
      ) %>% 
   replace(is.na(.), 0)%>% 
   dplyr::select(., -one_of(c("Group"))) %>% 
   mutate_if(is.numeric, ~round(.*10^6)) %>% 
   zCompositions::cmultRepl() %>%
   clr()  %>% 
   .[order(rownames(.)),] 



Metabolites_NAFLD_O <-
   read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>% 
   column_to_rownames(var = "X") %>% 
   t() %>%  
   # .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   # magrittr::add(1) %>%
   # log2() %>%
   .[which(rownames(.) %in% rownames(KeggPaths_NAFLD_O)),] %>% 
   as.data.frame() %>% 
   .[order(rownames(.)),] 


Metabolites_NAFLD_L <-
   read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>% 
   column_to_rownames(var = "X") %>% 
   t() %>%  
   # .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   # magrittr::add(1) %>%
   # log2() %>%
   .[which(rownames(.) %in% rownames(KeggPaths_NAFLD_L)),] %>% 
   as.data.frame() %>% 
   .[order(rownames(.)),] 


Metabolites_NAFLD_All <-
   read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/metabolites_imputation.csv") %>% 
   column_to_rownames(var = "X") %>% 
   t() %>%  
   # .[ ,colSums(.!= 0) >= nrow(.)*0.1] %>% 
   # magrittr::add(1) %>%
   # log2() %>%
   .[which(rownames(.) %in% rownames(KeggPaths_NAFLD_All)),] %>% 
   as.data.frame() %>% 
   .[order(rownames(.)),] 


Primary_metadata <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(FLI.liver_function,
                 GGT.liver_function,
                 ALT.liver_function,
                 AST.liver_function,
                 Triglyceride.lipid_profiles,
                 Liver_Fat.liver_function) 


##### Performing Correlations  #######
var_adjust <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(BMI.body_weight_variables,
                 Age.demographic_variable,
                 Gender.demographic_variable) %>% 
   data.frame() %>% 
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_O)),] %>% 
   .[order(rownames(.)),] %>%
   as.data.frame()


Name_Change_metabolites <-
  read.csv("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Metabolomics/id_to_name.csv") %>% 
  dplyr::select(previous,new) %>% 
  dplyr::rename(from = "previous")

Metab_Metadata_Obese <- correlation_format(Metabolites_NAFLD_O,
                                        Primary_metadata,
                                        "Metabolites O",
                                        "Primary Metadata",
                                        adjusted = TRUE) %>% 
   left_join(Signif_Metabolites_O  %>% 
                dplyr::select(Species,Direction) %>% 
                dplyr::rename(from = "Species"), 
             by = "from") %>% 
   na.omit() %>% 
  left_join(Name_Change_metabolites, by = "from") %>% 
  mutate(from = new) %>% 
  dplyr::select(-new)
  
  


var_adjust <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(BMI.body_weight_variables,
                 Age.demographic_variable,
                 Gender.demographic_variable) %>% 
   data.frame() %>% 
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_L)),] %>% 
   .[order(rownames(.)),] %>%
   as.data.frame()


Metab_Metadata_Lean <- correlation_format(Metabolites_NAFLD_L,
                                        Primary_metadata,
                                        "Metabolites L",
                                        "Primary Metadata",
                                        adjusted = TRUE) %>% 
   left_join(Signif_Metabolites_L  %>% 
                dplyr::select(Species,Direction) %>% 
                dplyr::rename(from = "Species"), 
             by = "from") %>% 
  na.omit() %>% 
  left_join(Name_Change_metabolites, by = "from") %>% 
  mutate(from = new) %>% 
  mutate(from = gsub("4hpro","4-Hydroxy-Proline",from)) %>% 
  dplyr::select(-new)

lists <- list(Obese_Met =Metab_Metadata_Obese %>% dplyr::filter(pvalue < 0.05),
              Lean_Met =Metab_Metadata_Lean%>% dplyr::filter(pvalue < 0.05))

write_rds(lists,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Mediation/Metab_Cors.rds")


var_adjust <-  
   read.csv(file = "~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated_FLI_FIB4.csv", 
            row.names = 1) %>% 
   dplyr::select(BMI.body_weight_variables,
                 Age.demographic_variable,
                 Gender.demographic_variable) %>% 
   data.frame() %>% 
   .[which(rownames(.) %in% rownames(Metabolites_NAFLD_All)),] %>% 
   .[order(rownames(.)),] %>%
   as.data.frame() 
  


Metab_Metadata_All <- correlation_format(Metabolites_NAFLD_All,
                                       Primary_metadata,
                                       "Metabolites All",
                                       "Primary Metadata",
                                       adjusted = TRUE) %>% 
  left_join(Name_Change_metabolites, by = "from") %>% 
  mutate(from = new) %>% 
  mutate(from = gsub("4hpro","4-Hydroxy-Proline",from)) %>% 
  dplyr::select(-new)


Metab_Metadata_Hyper <-
   Metab_Metadata_All %>% 
   left_join(Commorbidities_Metabolites %>% 
                filter(Comparison == "Hypertension vs Control") %>% 
                dplyr::select(from,Direction), 
             by = "from") %>% 
   na.omit() %>% 
   mutate(from  = paste("H",from)) %>% 
   mutate(from_class = gsub("Metabolites All","Hypertension",from_class)) 

Metab_Metadata_Pre_Hyper <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Pre_Hypertension vs Control") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  mutate(from  = paste("P",from)) %>% 
  mutate(from_class = gsub("Metabolites All","Pre_Hypertension",from_class))



Metab_Metadata_Athero <-
   Metab_Metadata_All %>% 
   left_join(Commorbidities_Metabolites %>% 
                filter(Comparison == "Atherosclerosis") %>% 
                dplyr::select(from,Direction), 
             by = "from") %>% 
   na.omit() %>% 
   mutate(from  = paste(from,"A")) %>% 
   mutate(from_class = gsub("Metabolites All","Atherosclerosis",from_class))

Metab_Metadata_Diabetes_overweight <-
   Metab_Metadata_All %>% 
   left_join(Commorbidities_Metabolites %>% 
                filter(Comparison == "T2D_Obese vs Overweight") %>% 
                dplyr::select(from,Direction), 
             by = "from") %>% 
   na.omit() %>% 
   mutate(from  = paste("O",from)) %>% 
   mutate(from_class = gsub("Metabolites All","T2D Overweight",from_class))

Metab_Metadata_Prediabetes <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Prediabetes vs Overweight") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  mutate(from  = paste("T",from)) %>% 
  mutate(from_class = gsub("Metabolites All","Prediabetes",from_class))


Metab_Metadata_Diabetes_Lean <-
   Metab_Metadata_All %>% 
   left_join(Commorbidities_Metabolites %>% 
                filter(Comparison == "Diabetes_Lean") %>% 
                dplyr::select(from,Direction), 
             by = "from") %>% 
   na.omit() %>% 
   filter(from %in% (Commorbidities_Metabolites %>% filter(Comparison == "Diabetes_Lean") %>% pull(from))) %>% 
   mutate(from  = paste("L",from)) %>% 
   mutate(from_class = rep("Diabetes_Lean",nrow(.)))


############# Gathering all the Data ####
all_data_Correlation <-
   dplyr::bind_rows(Metab_Metadata_Obese,
                    Metab_Metadata_Lean,
                    Metab_Metadata_Hyper,
                    Metab_Metadata_Pre_Hyper,
                    Metab_Metadata_Athero,
                    Metab_Metadata_Diabetes_overweight,
                    Metab_Metadata_Prediabetes,
                    Metab_Metadata_Diabetes_Lean) %>%
   mutate(from = ifelse(from_class == 'Metabolites O',paste('Obese_',from,sep = ""),paste('Lean_',from,sep = ""))) 
   



extract_meta <- function(dataset,phrase) {
   result <-
      dataset %>% 
      filter(to == phrase) %>% 
      dplyr::select(from,from_class,correlation,pvalue,Direction) %>%
      mutate(phrase = ifelse(pvalue < 0.05 & correlation > 0,1,
                             ifelse(pvalue < 0.05 & correlation < 0,-1,0))) %>%
      dplyr::select(-c("pvalue","correlation"))
   
   return(result)
}

DF_FLI <- all_data_Correlation %>% extract_meta(phrase = 'FLI.liver_function') %>% dplyr::rename(FLI = phrase)
DF_GGT <- all_data_Correlation %>% extract_meta(phrase = 'GGT.liver_function') %>% dplyr::rename(GGT = phrase)
DF_ALT <- all_data_Correlation %>% extract_meta(phrase = 'ALT.liver_function') %>% dplyr::rename(ALT = phrase)
DF_AST <- all_data_Correlation %>% extract_meta(phrase = 'AST.liver_function') %>% dplyr::rename(AST = phrase)
DF_TGs <- all_data_Correlation %>% extract_meta(phrase = 'Triglyceride.lipid_profiles') %>% dplyr::rename(TGs = phrase)
DF_LF <- all_data_Correlation %>% extract_meta(phrase = 'Liver_Fat.liver_function') %>% dplyr::rename(Liver_Fat = phrase)

DF <-  
   DF_FLI %>% 
left_join(DF_GGT) %>% 
   left_join(DF_ALT) %>% 
   left_join(DF_AST) %>% 
   left_join(DF_TGs) %>% 
   left_join(DF_LF) %>% 
   mutate(from = gsub("Obese_","",from)) %>% 
   mutate(from = gsub("Lean_","",from)) %>% 
   mutate(from = ifelse(from  == 'Dihydroxyacetone' & from_class == "Metabolites L" | 
                           from  == 'L-Histidine' & from_class == "Metabolites L"| 
                           from  == 'Potassium' & from_class == "Metabolites L" | 
                           from  == 'Protoheme C34H30FeN4O4' & from_class == "Metabolites L" | 
                           from  == 'Sulfate'& from_class == "Metabolites L",paste(from,"-"),from)) %>% 
   replace(is.na(.), 0) %>% 
   mutate(from = dplyr::recode(from,"Amylose (n=300 repeat units, alpha-1,4-glc)" = "Amylose")) %>% 
   mutate(from = dplyr::recode(from,"O Glycogen, structure 2 (glycogenin-1,6-{7[1,4-Glc], 4[1,4-Glc]})" = "O Glycogen, structure 2")) %>% 
   mutate(from = dplyr::recode(from,"L Glycogen, structure 4 (glycogenin-1,6-{2[1,4-Glc], [1,4-Glc]})" = "L Glycogen, structure 4")) %>% 
   mutate(from = dplyr::recode(from,"Beta-1,3/1,4-glucan (Barley, n=4, Glc beta1->3,4 Glc)" = "Beta-1,3/1,4-glucan")) %>% 
   mutate(from = dplyr::recode(from,"Galactomannan(n=4 repeat units mannose, alpha-1,4 man)" = "Galactomannan 4")) %>% 
   mutate(from = dplyr::recode(from,"Galactomannan (n=600 repeat units mannose, alpha-1,4 man) A" = "Galactomannan 600 A")) %>% 
   mutate(from = dplyr::recode(from,"Guanosine 3  phosphate C10H12N5O8P" = "Guanosine 3  phosphate")) %>% 
   column_to_rownames(var = "from") %>% 
   mutate(Direction = ifelse(Direction == 1,"#E41C06","#3366cc"))
##### The plot ####


Metab_Metadata_Hyper <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Hypertension vs Control") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05) 

Metab_Metadata_Pre_Hyper <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Pre_Hypertension vs Control") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05) 

Metab_Metadata_Athero <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Atherosclerosis") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05)

Metab_Metadata_Diabetes_overweight <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "T2D_Obese vs Overweight") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05)


Metab_Metadata_Prediabetes <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Prediabetes vs Overweight") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05)


Metab_Metadata_Diabetes_Lean <-
  Metab_Metadata_All %>% 
  left_join(Commorbidities_Metabolites %>% 
              filter(Comparison == "Diabetes_Lean") %>% 
              dplyr::select(from,Direction), 
            by = "from") %>% 
  na.omit() %>% 
  filter(pvalue < 0.05)


# Obese_Coms_meta <- c(Metab_Metadata_Hyper$from,
#                      Metab_Metadata_Athero$from,
#                      Metab_Metadata_Diabetes_overweight$from) %>%  unique(.)
# 
# saveRDS(Obese_Coms_meta,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Obese_Coms_Metab.rds")
# 
# Lean_Coms_meta <- c(Metab_Metadata_Diabetes_Lean$from) %>%  unique(.)
#   
# 
# saveRDS(Lean_Coms_meta,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Networks/New_DGCA/Lean_Coms_Metab.rds")
# 
# 





plot_colors <- c(
   "Metabolites L" =  "#339f2c",
   "Metabolites O" = "#ff5500",
   "Atherosclerosis" = "#1f78b4",
   "Hypertension" = "#693c99",
   "Diabetes_Lean" = "#666666",
   "T2D Overweight" = "#fff233",
   "Prediabetes" = "#212f3d",
   "Pre_Hypertension" = "#48c9b0"
)


Category_1.1<- 
   factor(DF$from_class,levels = c(  "T2D Overweight",
                                     "Prediabetes",
                                     "Diabetes_Lean",
                                    "Atherosclerosis",
                                    "Metabolites O",
                                    "Metabolites L",
                                    "Pre_Hypertension",
                                    "Hypertension") ) 

col_meth = colorRamp2(c(-1,0,1), c("#3366cc","white","#E41C06"))
col_meth_2 = colorRamp2(c(0,1), c("white", "purple"))


pdf("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Chor_plot.pdf", width = 15, height = 15)

circos.clear()
circos.heatmap(DF %>% dplyr::select(from_class), 
               split = Category_1.1,
               col = plot_colors,
               bg.border = "black",
               cell.border	= "black",
               rownames.side = "outside",
               rownames.col = DF$Direction,
               rownames.cex = 0.52,
               cell_width = 5,
               track.height = 0.015)


circos.heatmap(DF %>% dplyr::select(FLI), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(GGT), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(ALT), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(AST), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(TGs), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(c(Liver_Fat)), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)


dev.off()



svg("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Figures_2/Figure 4/Chor_plot.svg", width = 15, height = 15)

circos.clear()
circos.heatmap(DF %>% dplyr::select(from_class), 
               split = Category_1.1,
               col = plot_colors,
               bg.border = "black",
               cell.border	= "black",
               rownames.side = "outside",
               rownames.col = DF$Direction,
               rownames.cex = 0.52,
               cell_width = 5,
               track.height = 0.015)


circos.heatmap(DF %>% dplyr::select(FLI), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(GGT), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(ALT), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(AST), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(TGs), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)

circos.heatmap(DF %>% dplyr::select(c(Liver_Fat)), 
               split = Category_1.1, 
               bg.border = "black",
               cell.border	= "black",
               col = col_meth, 
               rownames.font = 8,
               track.height = 0.015)


dev.off()

