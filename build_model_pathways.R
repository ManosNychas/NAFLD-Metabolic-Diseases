setwd("/media/AMarfil/data/Documents/NAFLD/new_groups/final_groups/KEGG_paths/")

pacman::p_load(tidyverse, Boruta, caret, caTools, varhandle, plyr)


# load data ---------------------------------------------------------------

disease_paths <- read.csv("NAFLD-O_paths.csv") %>% mutate(., Group = "NAFLD")
control_paths <- read.csv("Control-O_paths.csv") %>% mutate(., Group= "Control")
disease_species <- read.csv("../../../Data/NAFLD-O_species.csv") %>% mutate(., Group = "NAFLD")
control_species <- read.csv("../../../Data/Overweight_noNAFLD.csv") %>% mutate(., Group = "Control")

# merge datasets ----------------------------------------------------------

disease <- merge(disease_paths, disease_species, by = "Sample")
control <- merge(control_paths, control_species, by = "Sample")

disease.control <- rbind.fill(disease, control)
disease.control <- as.data.frame(apply(disease.control, 2, function(x) ifelse(is.na(x), 0, x)))

disease.control <- disease.control[sample(nrow(disease.control)), ]

disease.control$Sample <- NULL
disease.control$Group.y <- NULL
disease.control <- dplyr::rename(disease.control, "Group" = "Group.x")

disease.control <- mutate(disease.control, Group = as.factor(Group))

disease.control[, 2:ncol(disease.control)] <- unfactor(disease.control[, 2:ncol(disease.control)])


# split -------------------------------------------------------------------

set.seed(1)
split <- sample.split(disease.control$Group, SplitRatio = .8)

training_set <- disease.control[split == T, ]

training_data <- training_set[, 1:(ncol(training_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() 

training_grouping <- as.factor(training_set[, ncol(training_set)])

test_set <- disease.control[split == F, ]

test_data <- test_set[, 1:(ncol(test_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame()

test_grouping <- as.factor(test_set[, ncol(test_set)])


# feature selection -------------------------------------------------------

features <- data.frame(Iteration = 1:100,
                       Features = NA)

# run boruta 100 times

seeds <- sample.int(10000, 100)

for (i in 1:100) {
  
  set.seed(seeds[i])
  boruta_res <- Boruta(x = training_data, 
                       y = training_grouping, 
                       doTrace = 0)  
  
  features_i <- getSelectedAttributes(boruta_res)
  
  stats <- attStats(boruta_res) %>% 
    subset(., decision == "Confirmed") %>% 
    .[order(.$meanImp, decreasing = T), ] %>% 
    .[1:20, ]
  
  features_i <- subset(features_i, features_i %in% rownames(stats))
  
  features$Features[i] <- paste(features_i, collapse = ", ")
  
  setTxtProgressBar(pb, i)
  
}

# select top 20 features
features <- features$Features %>% 
  strsplit(", ") %>% 
  unlist() %>% 
  table() %>% 
  as.data.frame()


# train -------------------------------------------------------------------

training_data <- select_if(training_data, names(training_data) %in% features$.[1:20])

trControl <- trainControl(method = "cv",
                          number = 5,
                          search = "random",
                          sampling = "down",
                          verboseIter = F,
                          returnData = T,
                          savePredictions = T,
                          allowParallel = T,
                          classProbs = T,
                          seeds = NULL)

ranger_mod <- train(x = training_data,
                    y = training_grouping,
                    trControl = trControl,
                    method = "rf")

# test --------------------------------------------------------------------

test_data <- select_if(test_data, names(test_data) %in% features$.[1:20]) 

confusion_matrix_ranger <- predict(ranger_mod, test_data) %>% 
  confusionMatrix(., test_grouping)

probabilities_ranger <- predict(ranger_mod, 
                                test_data, 
                                type = "prob")

pROC_ranger <- pROC::roc(response =  test_grouping,
                         levels = c("NAFLD", "Control"),
                         predictor = probabilities_ranger[, "NAFLD"])