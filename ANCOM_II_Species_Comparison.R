
    library(dplyr)
    library(tidyverse)
    library(vegan)
    library(nlme)
    
    
    #Load Datasets
    source("~/Documents/Phd/HKI_PC/Metabolic_Diseases/ANCOM-master/scripts/ancom_v2.1.R")
    load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")
    load("~/Documents/Phd/HKI_PC/Metabolic_Diseases/Updated_Human3/R_datasets/Sequencing_Information.RData")
    
    
    
    NAFLD_O <- dplyr::select(NAFLD_O, -one_of(c("Group")))
    NAFLD_L <- dplyr::select(NAFLD_L, -one_of(c("Group")))
    
    Control_O <- dplyr::select(Control_O, -one_of(c("Group")))
    Control_L <- dplyr::select(Control_L, -one_of(c("Group")))
    
    
    #### ANCOM Function - Correct Direction ####
    
    ANCOM_II_difference_abundance <- function(Dataset_1,Dataset_2,prevCutoff = 0.1) {

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
    
    #NAFLD-O Comparison
    O_NvC <- 
      ANCOM_II_difference_abundance(NAFLD_O,Control_O,All_Sequencing)  
    #NAFLD-L Comparison
    L_NvC <- 
      ANCOM_II_difference_abundance(NAFLD_L,Control_L,All_Sequencing) 
    
    
    Results_ANCOM_II <- 
      rbind(O_NvC,L_NvC) %>% 
      cbind(c(rep("NAFLD-O vs Control-O",nrow(O_NvC)),
              rep("NAFLD-L vs Control-L",nrow(L_NvC))),.) %>% 
      dplyr::rename(Species = "taxa_id")
      
      colnames(Results_ANCOM_II)[1] = "Comparison"
      
      # 
      
      saveRDS(Results_ANCOM_II,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/Differential_Abundance/Species/New_Adjusts/Detailed_ANCOM_II_Res.rds")



