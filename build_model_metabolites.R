setwd("/media/AMarfil/data/Documents/NAFLD/new_groups/final_groups/metabolites/lessFeatures/")

pacman::p_load(tidyverse, varhandle, caret, Boruta, caTools, plyr, HTSSIP, cowplot)


# control data ------------------------------------------------------------

controlO_metabolites <- read.delim("../data/controlo_metabolites.tsv", check.names = F) 
controlO_species <- read.csv("../../../../Data/Overweight_noNAFLD.csv", check.names = F) 
controlO <- merge(controlO_species, controlO_metabolites, by = "Sample") %>% mutate(., Group = "Control")


# NAFLD data --------------------------------------------------------------

nafldO_metabolites <- read.delim("../data/nafldo_metabolites.tsv", check.names = F) 
nafldO_species <- read.csv("../../../../Data/NAFLD-O_species.csv", check.names = F)
nafldO <- merge(nafldO_species, nafldO_metabolites, by = "Sample") %>% mutate(., Group = "NAFLD")


# process data ------------------------------------------------------------

obese_data <- rbind.fill(nafldO, controlO) %>% 
  apply(., 2, function(x) ifelse(is.na(x), 0, x)) %>% 
  as.data.frame(.) %>% 
  column_to_rownames(., "Sample") %>% 
  unfactor(.) %>% 
  mutate(., Group = as.factor(Group)) %>% 
  .[sample(nrow(.)), ]


# split -------------------------------------------------------------------

set.seed(1)
split <- sample.split(obese_data$Group, SplitRatio = .8)

training_set <- obese_data[split == T, ]

training_data <- training_set[, 1:(ncol(training_set)-1)] %>% 
  unfactor() %>% 
  as.data.frame() 

training_grouping <- as.factor(training_set[, ncol(training_set)])

test_set <- obese_data[split == F, ]

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