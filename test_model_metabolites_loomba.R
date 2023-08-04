setwd("/media/AMarfil/data/Documents/NAFLD/new_groups/final_groups/metabolites/lessFeatures/new_mambo/")

pacman::p_load(plyr, caret, readxl, varhandle, pROC, cowplot, tidyverse)

# load data and format
model <- read_rds("../metabolites_model_lessFeatures.rds")

control_metabolites <- read.delim("loomba_control_new_metabolites_c1.tsv", check.names = F) %>%
  mutate(., Group = "Control")

control_species <- read.csv("../../../External_cohort/Loomba/clean/loomba_species_control_clean.csv") %>%
  mutate(., Group = "Control") %>% 
  mutate(Sample = gsub("^N\\d+_", "", Sample)) %>% 
  mutate(Sample = gsub("_", "-", Sample))

disease_metabolites <- read.delim("loomba_nafld_new_metabolites_c1.tsv", check.names = F) %>% 
  mutate(., Group = "NAFLD")

disease_species <- read.csv("../../../External_cohort/Loomba/clean/loomba_species_nafld_clean.csv") %>% 
  mutate(., Group = "NAFLD") %>% 
  mutate(Sample = gsub("^N\\d+_", "", Sample)) %>% 
  mutate(Sample = gsub("_", "-", Sample))

## merge datasets
data_metabolites <- rbind.fill(disease_metabolites, control_metabolites)
data_metabolites[] <- apply(data_metabolites, 2, function(x) ifelse(is.na(x), 0, x))

data_species <- rbind.fill(disease_species, control_species)
data_species[] <- apply(data_species, 2, function(x) ifelse(is.na(x), 0, x))

data <- merge(data_metabolites, data_species, by = "Sample")
data$Group.y <- NULL

features <- names(model$trainingData)
data <- select_if(data, names(data) %in% c("Group.x", features))

data[] <- apply(data, 2, function(x) as.numeric(x))

test_grouping <- data$Group.x %>% 
  as.factor(.)

data$Group.x <- NULL

test_data <- apply(data, 2, function(x) ifelse(is.na(x), 0, x))

## predict
confusion_matrix_ranger <- predict(model, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(model, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("NAFLD", "Control"),
                         predictor = probabilities_ranger[, "NAFLD"])