setwd("/media/AMarfil/data/Documents/NAFLD/new_groups/final_groups/metabolites/diseases/")

pacman::p_load(plyr, tidyverse, caret, readxl, varhandle, pROC, cowplot, tidyverse)

# load control data from Atherosclerosis cohort cohort and NAFLD model
real_group <- read.delim("../test_metabolites.tsv")$Group %>% unfactor(.)
model <- read_rds("../metabolites_model.rds")

control_metabolites <- read.delim("../data/control_atherosclerosis_metabolites.tsv", check.names = F) %>%
  mutate(., Group = "Control")

control_species <- read.csv("../../Control-Atherosclerosis_species.csv") %>%
  mutate(., Group = "Control")

control_data <- merge(control_species, control_metabolites, by = "Sample") %>% 
  column_to_rownames(., "Sample") %>% 
  select_if(., colnames(.) %in% c("Group.x", colnames(model$trainingData))) %>% 
  rename(., "Group" = "Group.x") %>% 
  rownames_to_column(., "Sample")

# portability -------------------------------------------------------------

## load disease data from NAFLD cohort 
disease_data <- read.delim("../test_metabolites.tsv", check.names = F) %>% 
  subset(., Group == "NAFLD")

## process data
data <- rbind.fill(control_data, disease_data)

data <- data[, 3:ncol(data)] 
data[] <- apply(data, 2, function(x) as.numeric(x))
data[] <- apply(data, 2, function(x) ifelse(is.na(x), 0, x))
group <- data$Group %>% 
  as.factor(.)

test_data <- apply(data, 2, function(x) ifelse(is.na(x), 0, x))
test_grouping <- group

## predict

confusion_matrix_ranger <- predict(model, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(model, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("NAFLD", "Control"),
                         predictor = probabilities_ranger[, "NAFLD"])

portability <- abs(pROC_ranger$auc-0.5)*2

predictor <- pROC_ranger$predictor %>% 
  as.data.frame(.) %>% 
  add_column(., "Real" = test_grouping) %>% 
  dplyr::rename(., "Predictor" = colnames(.)[1]) %>% 
  add_column(., "Comparison" = "Portability") %>% 
  mutate(., Real = ifelse(Real == "NAFLD", "NAFLD-O", "Control-ASVD"))

# prediction rate ---------------------------------------------------------

## load disease data from Atherosclerosis

disease_metabolites <- read.delim("../data/atherosclerosis_metabolites.tsv", check.names = F) %>% 
  mutate(., Group = "NAFLD")

disease_species <- read.csv("../../../../Data/Atherosclerosis.csv") %>% 
  mutate(., Group = "NAFLD")

## process data
metabolite_data <- rbind.fill(disease_metabolites, control_metabolites)
metabolite_data[] <- apply(metabolite_data, 2, function(x) ifelse(is.na(x), 0, x))

species_data <- rbind.fill(disease_species, control_species)
species_data[] <- apply(species_data, 2, function(x) ifelse(is.na(x), 0, x))

merged_data <- merge(metabolite_data, species_data, by = "Sample")
merged_data$Group.y <- NULL

features <- names(model$trainingData)
merged_data <- select_if(merged_data, names(merged_data) %in% c("Group.x", features))
group <- merged_data$Group.x %>%
  as.factor(.)
merged_data$Group.x <- NULL
data <- merged_data
data[] <- apply(data, 2, function(x) as.numeric(x))
data[] <- apply(data, 2, function(x) ifelse(is.na(x), 0, x))

test_data <- apply(data, 2, function(x) ifelse(is.na(x), 0, x))
test_grouping <- group

## predict

confusion_matrix_ranger <- predict(model, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(model, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("NAFLD", "Control"),
                         predictor = probabilities_ranger[, "NAFLD"])

## calculate detection rate

cutoff <- 0.552 # determine manually

detection_rate_df <- data.frame(Predictor = pROC_ranger$predictor, Real = group, Cutoff = ifelse(pROC_ranger$predictor > cutoff, "NAFLD", "Control")) %>% 
  subset(., Real == "NAFLD")

detection_rate <- subset(detection_rate_df, Cutoff == "NAFLD")

detection_rate <- length(detection_rate$Predictor)/length(detection_rate_df$Predictor)