rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dplyr)
library(tidyr)
library(RANN)
library(plyr)
library(magrittr)
library(reshape)
library(data.table)
library(readxl)
library(mice)
#### Improting Data


################################################## RS NAFLD project --> NAFLD MRI Diagnosed ####################################################################

NAFLD_MRI <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/RS_NAFLD/Metaphlan/RS_NAFLD_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(NAFLD_MRI) <- names(NAFLD_MRI) %>% gsub("_metaphlan_bowtie2", "", .)
colnames(NAFLD_MRI) <- paste("RS_NAFLD_Project", colnames(NAFLD_MRI), sep = "_")
NAFLD_MRI = as.data.frame(NAFLD_MRI)
rownames(NAFLD_MRI) = NAFLD_MRI$RS_NAFLD_Project_ID
NAFLD_MRI <- NAFLD_MRI[,-1] %>% t(.) %>% as.data.frame(.) %>% mutate(.,Group = rep("NAFLD_MRI",nrow(.))) %>% select(Group,everything())


## metadata
NAFLD_MRI_metadata <- read_excel("Documents/Metabolic_Diseases/Original_metadata/RS_BestTreat_metadata.xls") %>% as.data.frame(.) %>% .[.$ID %like% "V1",] %>%
  select(.,ID,age,gender,BMI,SBP,DBP,PG0)

NAFLD_MRI_metadata$ID <- paste("RS_NAFLD_Project", NAFLD_MRI_metadata$ID, sep = "_")
rownames(NAFLD_MRI_metadata) <- NAFLD_MRI_metadata$ID
NAFLD_MRI_metadata = NAFLD_MRI_metadata[,-1]

colnames(NAFLD_MRI_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG")

NAFLD_MRI_metadata$Gender = as.factor(NAFLD_MRI_metadata$Gender)


nrow(NAFLD_MRI)
nrow(NAFLD_MRI[rownames(NAFLD_MRI) %in% rownames(NAFLD_MRI_metadata),])


################################################ NASH project --> NASH Biopsy Diagnosed #######################################################################

### Metadata -- Separate Groups (Controls - Steatosis -Borderline - NASH)

sample_group_nash  <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/NASH/Sample_INfo/sample.data.txt",header = T, sep = "\t") %>% as_tibble()
sample_group_nash =  sample_group_nash[order(sample_group_nash$ID),]

### Species ##
NASH_Biopsy_project <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/NASH/Metaphlan/NASH_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(NASH_Biopsy_project) <- names(NASH_Biopsy_project) %>% gsub("_metaphlan_bowtie2", "", .)
colnames(NASH_Biopsy_project) <- paste("NASH_Biopsy_Project", colnames(NASH_Biopsy_project), sep = "_")
NASH_Biopsy_project = as.data.frame(NASH_Biopsy_project)
rownames(NASH_Biopsy_project) = NASH_Biopsy_project$NASH_Biopsy_Project_ID
NASH_Biopsy_project <- NASH_Biopsy_project[,-1] %>% t(.) %>% as.data.frame(.) 

rownames(NASH_Biopsy_project) = paste(rownames(NASH_Biopsy_project),sample_group_nash$Group,sep = "_")
Group <- gsub("NASH_Biopsy_Project_s.*Borderline", "Borderline_Biopsy", rownames(NASH_Biopsy_project)) %>% gsub("NASH_Biopsy_Project_s.*Normal", "Overweight_NAFLD_Free_Biopsy", .) %>% 
  gsub("NASH_Biopsy_Project_s.*NASH", "NASH_Biopsy", .) %>% gsub("NASH_Biopsy_Project_s.*Steatosis", "Steatosis_Biopsy", .)
NASH_Biopsy_project = cbind(Group,NASH_Biopsy_project)

Steatosis_Biopsy = NASH_Biopsy_project[NASH_Biopsy_project$Group == "Steatosis_Biopsy", ]
Bordeline_Biopsy = NASH_Biopsy_project[NASH_Biopsy_project$Group == "Borderline_Biopsy", ]
Overweight_NAFLD_Free_Biopsy = NASH_Biopsy_project[NASH_Biopsy_project$Group == "Overweight_NAFLD_Free_Biopsy", ]
NASH_Biopsy = NASH_Biopsy_project[NASH_Biopsy_project$Group == "NASH_Biopsy", ]



##metadata

#Hba1c

NASH_Biopsy_project_metadata <- sample_group_nash %>% select(.,ID,AgeAtOperation,Sex,BMI,SystolicBP,DiastolicBP,Glu) %>% as.data.frame(.)
NASH_Biopsy_project_metadata$ID <-paste("NASH_Biopsy_Project", NASH_Biopsy_project_metadata$ID,sample_group_nash$Group, sep = "_")
rownames(NASH_Biopsy_project_metadata) <- NASH_Biopsy_project_metadata$ID
NASH_Biopsy_project_metadata = NASH_Biopsy_project_metadata[,-1]

colnames(NASH_Biopsy_project_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG")

NASH_Biopsy_project_metadata$Age = as.double(NASH_Biopsy_project_metadata$Age)
NASH_Biopsy_project_metadata$SBP = as.double(NASH_Biopsy_project_metadata$SBP)
NASH_Biopsy_project_metadata$DBP = as.double(NASH_Biopsy_project_metadata$DBP) 
NASH_Biopsy_project_metadata$Gender = as.factor(NASH_Biopsy_project_metadata$Gender) %>% revalue(., c( male = 1, female =2  )) %>% factor(.,levels = c(1,2))




###################################################### Atherosclerosis ######################################################################################### 
sample_group_athero  <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Atherosclerosis/Metadata/All_IDs.csv",header = T, sep = ",") %>% as_tibble()
sample_group_athero =  sample_group_athero[order(sample_group_athero$run_accession),]
sample_group_athero$run_accession = as.character(sample_group_athero$run_accession)


#### Changing categories to control samples ###
metadata <- read.csv("Documents/Metabolic_Diseases/Original_metadata/Atherosclerosis_Kristiansen.csv",sep = "\t") %>% as.data.frame(.)
metadata$Systolic.Blood.Pressure..mmHg. = as.numeric(metadata$Systolic.Blood.Pressure..mmHg.)
metadata$Diastolic.Blood.Pressure..mmHg. = as.numeric(metadata$Diastolic.Blood.Pressure..mmHg.)
metadata$HbA1c........ = as.numeric(metadata$HbA1c....)
metadata$FBG = as.numeric(metadata$FBG)
metadata$Body.Mass.Index..Body.Mass.Index..BMI.. = as.numeric(metadata$Body.Mass.Index..BMI.)


newmeta <-metadata %>% 
  replace(is.na(.),-1) %>%
  mutate(.,Disease_Group = 
           ifelse(((Systolic.Blood.Pressure..mmHg. >= 140 | Diastolic.Blood.Pressure..mmHg. >= 90) & Body.Mass.Index..BMI. >= 25 & (FBG >= 7 | HbA1c.... >=6.5)),"Hypertension_Overweight_Diabetes",
                  ifelse(Body.Mass.Index..BMI. >= 25 & (FBG >= 7 | HbA1c.... >= 6.5),"Overweight_Diabetes",
                         ifelse((Systolic.Blood.Pressure..mmHg. >= 140 | Diastolic.Blood.Pressure..mmHg. >= 90) & (FBG >= 7 | HbA1c.... >= 6.5),"Hypertension_Diabetes",
                                ifelse((Systolic.Blood.Pressure..mmHg. >= 140 | Diastolic.Blood.Pressure..mmHg. >= 90) &Body.Mass.Index..BMI. >= 25,"Hypertension_Overweight",
                                       ifelse(((Systolic.Blood.Pressure..mmHg. <= 139 & Systolic.Blood.Pressure..mmHg. > 125) | (Diastolic.Blood.Pressure..mmHg. <= 89 & Diastolic.Blood.Pressure..mmHg. > 80)) & (FBG >= 7 | HbA1c.... >= 6.5),"Pre_Hypertension_Diabetes",
                                              ifelse(((Systolic.Blood.Pressure..mmHg. <= 139 & Systolic.Blood.Pressure..mmHg. > 125) | (Diastolic.Blood.Pressure..mmHg. <= 89 & Diastolic.Blood.Pressure..mmHg. > 80)) & (Body.Mass.Index..BMI. >= 25),"Pre_Hypertension_Overweight",
                                                     ifelse(( Systolic.Blood.Pressure..mmHg. >= 140 | Diastolic.Blood.Pressure..mmHg. >= 90),"Hypertension",
                                                            ifelse(((Systolic.Blood.Pressure..mmHg. <= 139 & Systolic.Blood.Pressure..mmHg. > 125) | (Diastolic.Blood.Pressure..mmHg. <= 89 & Diastolic.Blood.Pressure..mmHg. > 80)) ,"Pre_Hypertension",
                                                                   ifelse((FBG >= 7 | HbA1c.... >= 6.5),"Diabetes",
                                                                          ifelse(Body.Mass.Index..BMI. >= 25, "Overweight","Normal"))))))))))) 


newmeta[newmeta == -1] <- NA
colnames(newmeta)[1] = "Clinical_IDs"
sample_info = left_join(sample_group_athero,newmeta[,c(1,ncol(newmeta))],by = "Clinical_IDs")


sample_group_athero_disease = sample_info[sample_info$Clinical_IDs %like% "ZSL", ]
sample_group_athero_controls= sample_info[!sample_info$Clinical_IDs %like% "ZSL", ]

Atherosclerosis_Project <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Atherosclerosis/Metaphlan/Atherosclerosis_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(Atherosclerosis_Project) <- names(Atherosclerosis_Project) %>% gsub("_metaphlan_bowtie2", "", .)
Atherosclerosis_Project = as.data.frame(Atherosclerosis_Project)
rownames(Atherosclerosis_Project) = Atherosclerosis_Project$ID
Atherosclerosis_Project <- Atherosclerosis_Project[,-1] %>% t(.) %>% as.data.frame(.)


###  Combining Metadata IDs -  Diseased Samples ####
sample_group_athero_disease = sample_group_athero_disease[which(sample_group_athero_disease$run_accession %in% rownames(Atherosclerosis_Project) ),]
Atherosclerosis = Atherosclerosis_Project[which(rownames(Atherosclerosis_Project) %in% sample_group_athero_disease$run_accession),]


Atherosclerosis =  Atherosclerosis[order(rownames(Atherosclerosis)),]
rownames(Atherosclerosis) = paste("Atherosclerosis_Project",sample_group_athero_disease$Clinical_IDs,sep = "_")
Group <- gsub("Atherosclerosis_Project.*", "Atherosclerosis", rownames(Atherosclerosis)) 
Atherosclerosis = cbind(Group,Atherosclerosis)


###  Combining Metadata IDs -  Control Samples ####
sample_group_athero_controls = sample_group_athero_controls[which(sample_group_athero_controls$run_accession %in% rownames(Atherosclerosis_Project) ),]
Control_Atherosclerosis = Atherosclerosis_Project[which(rownames(Atherosclerosis_Project) %in% sample_group_athero_controls$run_accession),]


Control_Atherosclerosis =  Control_Atherosclerosis[order(rownames(Control_Atherosclerosis)),]
rownames(Control_Atherosclerosis) = paste("Atherosclerosis_Project",sample_group_athero_controls$Clinical_IDs,sep = "_")
Group <- gsub("Atherosclerosis_Project.*", "Controls", rownames(Control_Atherosclerosis)) 
Control_Atherosclerosis = cbind(Group,Control_Atherosclerosis)
Control_Atherosclerosis$Group = paste(sample_group_athero_controls$Disease_Group,sep = "_")

#### metadata


Atherosclerosis_project_metadata = metadata %>% select(.,Sample.ID,Age..year.,Gender,Body.Mass.Index..BMI.,
                                                       Systolic.Blood.Pressure..mmHg.,Diastolic.Blood.Pressure..mmHg.,
                                                       FBG,Clopidogrel.Hydrogen.Sulphate.Tablets_117,Metoprolol_19)

Atherosclerosis_project_metadata$Sample.ID <- paste("Atherosclerosis_Project",Atherosclerosis_project_metadata$Sample.ID ,sep = "_")
rownames(Atherosclerosis_project_metadata) = Atherosclerosis_project_metadata$Sample.ID

Atherosclerosis_project_metadata = Atherosclerosis_project_metadata[,-1]

colnames(Atherosclerosis_project_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG","Clopidogrel_Hydrogen_Sulphate","Metoprolol")

Atherosclerosis_project_metadata$Age = as.numeric(Atherosclerosis_project_metadata$Age)
Atherosclerosis_project_metadata$Gender = as.factor(Atherosclerosis_project_metadata$Gender) %>% revalue(., c( male = 1, female =2  )) %>% factor(.,levels = c(1,2))
Atherosclerosis_project_metadata$Clopidogrel_Hydrogen_Sulphate = as.factor(Atherosclerosis_project_metadata$Clopidogrel_Hydrogen_Sulphate)
Atherosclerosis_project_metadata$Metoprolol = as.factor(Atherosclerosis_project_metadata$Metoprolol)





###################################################### Prospective ######################################################################################### 


metadata_prospective_nafld <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Prospective_NAFLD/Prospective_NAFLD_metadata2014.csv") %>% as.data.frame(.)


metadata_prospective_nafld <- metadata_prospective_nafld %>% mutate(., BMI = Weight/(Height1/100)^2)

metadata_prospective_nafld$fasting.blood.glucose = as.numeric(metadata_prospective_nafld$fasting.blood.glucose)
metadata_prospective_nafld$SBP1 = as.numeric(metadata_prospective_nafld$SBP1)
metadata_prospective_nafld$DBP1 = as.numeric(metadata_prospective_nafld$DBP1)
metadata_prospective_nafld$BMI = as.numeric(metadata_prospective_nafld$BMI)
metadata_prospective_nafld$HbA1c= as.numeric(metadata_prospective_nafld$HbA1c)

ID_map = read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Prospective_NAFLD/ID_mapping.csv") %>% as.data.frame(.)
colnames(ID_map)[3] = "id"

##### 3 Samples are removed because of bad quality 

metadata_prospective_nafld = left_join(ID_map,metadata_prospective_nafld,by = "id") %>% .[,-c(1,3,4,5)]

newmeta_prospective_nafld <-metadata_prospective_nafld %>% 
  replace(is.na(.),-1) %>%
  mutate(.,Disease_Group = 
           ifelse(((SBP1 >= 140 | DBP1 >= 90) & BMI >= 25 & (fasting.blood.glucose >= 7 | HbA1c >=6.5)),"Hypertension_Overweight_Diabetes",
                  ifelse(BMI >= 25 & (fasting.blood.glucose >= 7 | HbA1c >= 6.5),"Overweight_Diabetes",
                         ifelse((SBP1 >= 140 | DBP1 >= 90) & (fasting.blood.glucose >= 7 | HbA1c >= 6.5),"Hypertension_Diabetes",
                                ifelse((SBP1 >= 140 | DBP1 >= 90) & BMI >= 25,"Hypertension_Overweight",
                                       ifelse(((SBP1 <= 139 & SBP1 > 125) | (DBP1 <= 89 & DBP1 > 80)) & (fasting.blood.glucose >= 7 | HbA1c >= 6.5),"Pre_Hypertension_Diabetes",
                                              ifelse(((SBP1 <= 139 & SBP1 > 125) | (DBP1 <= 89 & DBP1 > 80)) & (BMI >= 25),"Pre_Hypertension_Overweight",
                                                     ifelse(( SBP1 >= 140 | DBP1 >= 90),"Hypertension",
                                                            ifelse(((SBP1 <= 139 & SBP1 > 125) | (DBP1 <= 89 & DBP1 > 80)) ,"Pre_Hypertension",
                                                                   ifelse((fasting.blood.glucose >= 7 | HbA1c >= 6.5),"Diabetes",
                                                                          ifelse(BMI >= 25, "Overweight","Normal"))))))))))) 

newmeta_prospective_nafld[newmeta_prospective_nafld == -1] <- NA


Prospective_NAFLD <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Prospective_NAFLD/Metaphlan/Prospective_NAFLD_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(Prospective_NAFLD) <- names(Prospective_NAFLD) %>% gsub("_metaphlan_bowtie2", "", .)
Prospective_NAFLD = as.data.frame(Prospective_NAFLD)
rownames(Prospective_NAFLD) = Prospective_NAFLD$ID
Prospective_NAFLD <- Prospective_NAFLD[,-1] %>% t(.) %>% as.data.frame(.)

rownames(Prospective_NAFLD)[1:5] = c("160","173","54","69","73")


newmeta_prospective_nafld = newmeta_prospective_nafld[which(newmeta_prospective_nafld$MetagenomicSampleID %in% rownames(Prospective_NAFLD) ),]
newmeta_prospective_nafld =  newmeta_prospective_nafld[order((newmeta_prospective_nafld$MetagenomicSampleID)),]
Prospective_NAFLD = Prospective_NAFLD[which(rownames(Prospective_NAFLD) %in% newmeta_prospective_nafld$MetagenomicSampleID),]
Prospective_NAFLD =  Prospective_NAFLD[order(rownames(Prospective_NAFLD)),]

rownames(Prospective_NAFLD) = paste("Prospective_NAFLD_Project",rownames(Prospective_NAFLD),sep = "_")
Prospective_NAFLD = cbind(newmeta_prospective_nafld$Disease_Group,Prospective_NAFLD)

colnames(Prospective_NAFLD)[1] = "Group"

Prospective_NAFLD$Group = paste(Prospective_NAFLD$Group,"NAFLD_Free_Ultrasound",sep = "_")


Prospective_NAFLD_metadata <- metadata_prospective_nafld %>% select(.,MetagenomicSampleID,age,gender,BMI,SBP1,DBP1,fasting.blood.glucose)
Prospective_NAFLD_metadata$MetagenomicSampleID = paste("Prospective_NAFLD_Project",Prospective_NAFLD_metadata$MetagenomicSampleID,sep = "_")
rownames(Prospective_NAFLD_metadata) = Prospective_NAFLD_metadata$MetagenomicSampleID

Prospective_NAFLD_metadata = Prospective_NAFLD_metadata[,-1]

colnames(Prospective_NAFLD_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG")

Prospective_NAFLD_metadata$Gender = as.factor(Prospective_NAFLD_metadata$Gender)




###################################################### Hypertension ######################################################################################### 


###  Setting the metadata and the IDs for correct categorization 10 samples have different IDs therefore theya re dismissed
sample_group_hypertension  <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Hypertension/Sample_Info/names_ids.csv", sep = ",") %>% as_tibble()
metadata <- read_excel("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Hypertention_ClinicalData-2020.2.18_ORIGINAL.xls") %>% as.data.frame(.)

metadata <- metadata %>% mutate(., BMI = `WEIGHT kg`/(`HEIGHT cm`/100)^2)
metadata$`FBG mmol/L` = as.numeric(metadata$`FBG mmol/L`)
metadata$BMI= as.numeric(metadata$BMI)

newmeta <-metadata %>% 
  replace(is.na(.),-1) %>%
  mutate(.,Disease_Group = 
           ifelse(( `FBG mmol/L` >= 7 & BMI >= 25),"Overweight_Diabetes",
                  ifelse(`FBG mmol/L` >= 7,"Diabetes",
                         ifelse(BMI >= 25,"Overweight","Normal")))) 

newmeta[newmeta == -1] <- NA

colnames(newmeta)[1] = "Clinical_IDs"

sample_info = left_join(sample_group_hypertension,newmeta[,c(1,ncol(newmeta))],by = "Clinical_IDs") %>% na.omit(.)

#### Importing Profiles

Hypertension_Project <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Hypertension/Metaphlan/Hypertension_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(Hypertension_Project) <- names(Hypertension_Project) %>% gsub("_metaphlan_bowtie2", "", .)
Hypertension_Project = as.data.frame(Hypertension_Project)
rownames(Hypertension_Project) = Hypertension_Project$ID
Hypertension_Project <- Hypertension_Project[,-1] %>% t(.) %>% as.data.frame(.)

sample_info_hypertension = sample_info[which(sample_info$run_accession %in% rownames(Hypertension_Project)),] %>% .[order(.$run_accession),]
Hypertension_Project = Hypertension_Project[which(rownames(Hypertension_Project) %in% sample_info$run_accession ),] %>% .[order(rownames(.)),]


rownames(Hypertension_Project) = paste("Hypertension_Project",sample_info_hypertension$Clinical_IDs,sample_info_hypertension$Group,sep = "_")

Group <- gsub("Hypertension_.*_HTN", "Hypertension", rownames(Hypertension_Project)) %>% gsub("Hypertension_.*pHTN", "Pre_Hypertension", .) %>% 
  gsub("Hypertension_.*Control", "Control_NO_HTN", .)
Hypertension_Project = cbind(Group,Hypertension_Project)

###Spliting into Categoreis

Pre_Hypertension_start = Hypertension_Project[Hypertension_Project$Group == "Pre_Hypertension", ]
Hypertension_start= Hypertension_Project[Hypertension_Project$Group == "Hypertension", ]
Control_Hypertension = Hypertension_Project[Hypertension_Project$Group == "Control_NO_HTN", ]


#### categorizing controls 

sample_info_hypertension_controls  = sample_info_hypertension[sample_info_hypertension$Group == "Control",]
Control_Hypertension <- Control_Hypertension %>% mutate(.,Group = sample_info_hypertension_controls$Disease_Group) %>% select(Group,everything())

### metadata

metadata <- metadata %>% dplyr::rename(Height = `HEIGHT cm`,Weight = `WEIGHT kg`,FBG = `FBG mmol/L` , TC = `TC mmol/L`, TG = `TG mmol/L`)

tempData <- mice(metadata,m=5,maxit=50,meth='pmm',seed=500)

Hypertension_project_metadata <- complete(tempData,1) %>%  select(.,SampleID,Age,Gender,BMI,SBP,DBP,FBG) 
colnames(Hypertension_project_metadata)[1] = "ID"

colnames(Hypertension_project_metadata)[1] = "Clinical_IDs"

Hypertension_project_metadata <- left_join(Hypertension_project_metadata,sample_group_hypertension) %>% select(.,-one_of("X","sample_accession","run_accession"))

Hypertension_project_metadata$Clinical_IDs = paste("Hypertension_Project",Hypertension_project_metadata$Clinical_IDs,Hypertension_project_metadata$Group,sep = "_")

rownames(Hypertension_project_metadata) = Hypertension_project_metadata$Clinical_IDs
Hypertension_project_metadata = Hypertension_project_metadata[,-c(1,8)]


colnames(Hypertension_project_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG")


Hypertension_project_metadata$Gender = revalue(Hypertension_project_metadata$Gender, c( male = 1, female =2  )) %>% as.factor(.)





###################################################### Diabesity ######################################################################################### 
metadata_diabesity <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Diabesity/DiabesityProject_samples182_MetaData_Disease_new20190130.txt") %>% as.data.frame(.)


#Create BMI Column

metadata_diabesity$systolic_pressure = as.numeric(metadata_diabesity$systolic_pressure)
metadata_diabesity$diastolic_pressure = as.numeric(metadata_diabesity$diastolic_pressure)
metadata_diabesity$BMI= as.numeric(metadata_diabesity$BMI)


newmeta_diabesity <-metadata_diabesity %>% 
  replace(is.na(.),-1) %>%
  mutate(.,Disease_Group = 
           ifelse(( (systolic_pressure >= 140 |diastolic_pressure >= 90) & BMI >= 25),"Hypertension_Overweight",
                  ifelse((( (systolic_pressure <= 139 & systolic_pressure > 125 ) | (diastolic_pressure <= 89 & diastolic_pressure > 80 )) & BMI >= 25),"Pre_Hypertension_Overweight",
                         ifelse((( (systolic_pressure <= 139 & systolic_pressure > 125 ) | (diastolic_pressure <= 89 & diastolic_pressure > 80 ))),"Pre_Hypertension",
                                ifelse(systolic_pressure >= 140 |diastolic_pressure >= 90,"Hypertension",
                                       ifelse(BMI >= 25,"Overweight","Normal")))))) 

newmeta_diabesity[newmeta_diabesity == -1] <- NA

newmeta_diabesity =  newmeta_diabesity[order(newmeta_diabesity$ID),]
newmeta_diabesity = newmeta_diabesity[!(newmeta_diabesity$ID) %like% "_160", ]
newmeta_diabesity = newmeta_diabesity[!(newmeta_diabesity$ID)  %like% "_173", ] ## Doesnt exist?
newmeta_diabesity = newmeta_diabesity[!(newmeta_diabesity$ID)  %like% "_54", ]
newmeta_diabesity = newmeta_diabesity[!(newmeta_diabesity$ID)  %like% "_69", ]
newmeta_diabesity = newmeta_diabesity[!(newmeta_diabesity$ID)  %like% "_73", ]


NL = newmeta_diabesity[(grepl ("NL_.*",  newmeta_diabesity$ID)) ,]
TL = newmeta_diabesity[(grepl ("TL_.*",  newmeta_diabesity$ID)) ,]
NO = newmeta_diabesity[(grepl ("NO_.*",  newmeta_diabesity$ID)) ,]
TO = newmeta_diabesity[(grepl ("TO_.*",  newmeta_diabesity$ID)) ,]



Diabesity <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Diabesity/Metaphlan/Diabesity_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(Diabesity) <- names(Diabesity) %>% gsub("_metaphlan_bowtie2", "", .)
Diabesity = as.data.frame(Diabesity)
rownames(Diabesity) = Diabesity$ID
Diabesity <- Diabesity[,-1] %>% t(.) %>% as.data.frame(.)




Diabesity =  Diabesity[order(rownames(Diabesity)),]
Group <- gsub("NL_.*", "Lean_Normal", rownames(Diabesity)) %>% gsub("NO_.*", "Overweight_Normal", .) %>% 
  gsub("TO_.*", "Overweight_T2D", .) %>% gsub("TL_.*", "Lean_T2D", .)
Diabesity = cbind(Group,Diabesity)

rownames(Diabesity) = paste("Diabesity_Project",rownames(Diabesity),sep = "_")

Diabesity = Diabesity[!rownames(Diabesity) %like% "_160", ]
Diabesity = Diabesity[!rownames(Diabesity) %like% "_173", ] ## Doesnt exist?
Diabesity = Diabesity[!rownames(Diabesity) %like% "_54", ]
Diabesity = Diabesity[!rownames(Diabesity) %like% "_69", ]
Diabesity = Diabesity[!rownames(Diabesity) %like% "_73", ]

Diabesity =  Diabesity[order(rownames(Diabesity)),]

Lean_T2D = Diabesity[Diabesity$Group == "Lean_T2D", ] %>% mutate(.,Group = TL$Disease_Group) 
Overweight_T2D = Diabesity[Diabesity$Group == "Overweight_T2D", ]  %>% mutate(.,Group = TO$Disease_Group) 
Overweight_Normal = Diabesity[Diabesity$Group == "Overweight_Normal", ]  %>% mutate(.,Group = NO$Disease_Group) 
Lean_Normal = Diabesity[Diabesity$Group == "Lean_Normal", ]  %>% mutate(.,Group = NL$Disease_Group) 



Diabesity_project_metadata <- metadata_diabesity %>% select(.,ID,Age,gender_1M_2F_,BMI,systolic_pressure,diastolic_pressure)


Additional_metadata_diabesity <-read_excel("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Diabesity/Diabesity_resuced.xlsx") %>% as.data.frame(.) %>%
  select(.,Raw_ID,`fating glucose`)

metadata_map <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/Diabesity/DiabesityProject_samples182_MetaData_General_new_withFGF19.txt",sep = ",") %>% as.data.frame(.) %>%
  select(.,ID,Raw_ID)

Additional_metadata_diabesity <- left_join(Additional_metadata_diabesity,metadata_map)

Diabesity_project_metadata <- left_join(Diabesity_project_metadata,Additional_metadata_diabesity)


Diabesity_project_metadata$ID = paste("Diabesity_Project",Diabesity_project_metadata$ID,sep = "_")

rownames(Diabesity_project_metadata) = Diabesity_project_metadata$ID
Diabesity_project_metadata = Diabesity_project_metadata[,-c(1,7)]


colnames(Diabesity_project_metadata) = c("Age","Gender","BMI","SBP","DBP","FBG")
Diabesity_project_metadata$Gender  = as.factor(Diabesity_project_metadata$Gender)

library(mice)
###################################################### Prediabetes ######################################################################################### 
Prediabetes <- read.csv("Documents/Metabolic_Diseases/Updated_Human3/Data/Prediabetes/Metaphlan/Prediabetes_metaphlan_merged_abundance_table_species.txt", sep = "\t") %>% as_tibble()
names(Prediabetes) <- names(Prediabetes) %>% gsub("_metaphlan_bowtie2", "", .)
colnames(Prediabetes) <- paste("Prediabetes_Project", colnames(Prediabetes), sep = "_")
Prediabetes = as.data.frame(Prediabetes)
rownames(Prediabetes) = Prediabetes$Prediabetes_Project_ID
Prediabetes <- Prediabetes[,-1] %>% t(.) %>% as.data.frame(.) %>% mutate(.,Group = rep("Prediabetes",nrow(.))) %>% select(Group,everything())

rownames(Prediabetes) <- rownames(Prediabetes) %>% gsub(".0w","",.) %>% gsub("X","R10",.)

###metadata

Prediabetes_Project_metadata <- read_excel("Documents/Metabolic_Diseases/Updated_Human3/Data/Metadata/ExerciseProject_Clinical data.xlsx") %>% as.data.frame(.) %>%
  dplyr::rename(Weight = `weight (kg)`,Fat_Mass = `fat mass`,Lean_mass = `lean mass` , HOMA_IR = `HOMA-IR`)

Prediabetes_Project_metadata$No. = as.factor(Prediabetes_Project_metadata$No.)
Prediabetes_Project_metadata$Note = as.factor(Prediabetes_Project_metadata$Note)


### Imputations

tempData <- mice(Prediabetes_Project_metadata,m=5,maxit=50,meth='pmm',seed=500)

Prediabetes_Project_metadata <- complete(tempData,1) %>% select(.,No.,age,BMI,SBP,DBP,glucose) %>% mutate(.,Gender = rep(1,nrow(.)))
colnames(Prediabetes_Project_metadata)[1] = "ID"

Prediabetes_Project_metadata$ID <-paste("Prediabetes_Project", Prediabetes_Project_metadata$ID, sep = "_")
rownames(Prediabetes_Project_metadata) = Prediabetes_Project_metadata$ID

Prediabetes_Project_metadata = Prediabetes_Project_metadata[,-1]

colnames(Prediabetes_Project_metadata) = c("Age","BMI","SBP","DBP","FBG","Gender")

Prediabetes_Project_metadata$Gender = as.factor(Prediabetes_Project_metadata$Gender) %>% factor(.,levels = c(1,2))



Prediabetes_Project_metadata <- Prediabetes_Project_metadata %>% relocate(Gender, .after = Age)




######################################################## Ultimate Grouping######################################################## 


### Metadata 
Combined_metadata = rbind.fill(NAFLD_MRI_metadata,
                               NASH_Biopsy_project_metadata,
                               Prospective_NAFLD_metadata,
                               Atherosclerosis_project_metadata,
                               Prediabetes_Project_metadata,
                               Diabesity_project_metadata,
                               Hypertension_project_metadata)

rownames(Combined_metadata) = c(rownames(NAFLD_MRI_metadata),
                                rownames(NASH_Biopsy_project_metadata),
                                rownames(Prospective_NAFLD_metadata),
                                rownames(Atherosclerosis_project_metadata),
                                rownames(Prediabetes_Project_metadata),
                                rownames(Diabesity_project_metadata),
                                rownames(Hypertension_project_metadata))


Combined_metadata$Clopidogrel_Hydrogen_Sulphate[is.na(Combined_metadata$Clopidogrel_Hydrogen_Sulphate)] <- 0
Combined_metadata$Metoprolol[is.na(Combined_metadata$Metoprolol)] <- 0

Combined_metadata$Clopidogrel_Hydrogen_Sulphate = as.factor(Combined_metadata$Clopidogrel_Hydrogen_Sulphate)
Combined_metadata$Metoprolol = as.factor(Combined_metadata$Metoprolol)


df <- read.csv("~/Documents/Metabolic_Diseases/Updated_Human3/Common_Metadata_Sara/FINAL_metadata_all_samples_Updated.csv")
Combined_metadata <- Combined_metadata[which(rownames(Combined_metadata) %in% df$X),]

df <- df[order(match(df$X,rownames(Combined_metadata))),]

Combined_metadata <- 
  cbind(Combined_metadata,df$HOMAIR.glucose_insulin) %>% 
  dplyr::rename(HOMAIR = `df$HOMAIR.glucose_insulin`) %>% 
  select(Age,Gender,BMI,HOMAIR,SBP,Clopidogrel_Hydrogen_Sulphate,Metoprolol)

######################################################## NAFLD Subjects  ######################################################## 
all_nafld <- bind_rows(NAFLD_MRI,
                       Steatosis_Biopsy,
                       Bordeline_Biopsy,
                       NASH_Biopsy)
all_nafld[is.na(all_nafld)] <- 0

nafld_metadata <- Combined_metadata[which(rownames(Combined_metadata)%in% rownames(all_nafld)),]

all_nafld <-
  data.frame(BMI = nafld_metadata$BMI,
             all_nafld)

NAFLD_O <-
  all_nafld %>% 
  filter(BMI > 25) %>% 
  dplyr::select(., -one_of("BMI"))

NAFLD_L <-
  all_nafld %>% 
  filter(BMI < 25) %>% 
  dplyr::select(., -one_of("BMI"))

######################################################## Control - NAFLD Subjects  ######################################################## 

Control_L  = Prospective_NAFLD[Prospective_NAFLD$Group %like% "Normal",]
rownames(Control_L) = rownames(Prospective_NAFLD[Prospective_NAFLD$Group %like% "Normal",])


Control_O = rbind.fill(Overweight_NAFLD_Free_Biopsy,
                                   Prospective_NAFLD[Prospective_NAFLD$Group %like% "Overweight",])
Control_O[is.na(Control_O)] <- 0
rownames(Control_O) = c(rownames(Overweight_NAFLD_Free_Biopsy),
                                    rownames(Prospective_NAFLD[Prospective_NAFLD$Group %like% "Overweight",]))

######################################################## Overweight T2D Groups######################################################## 

Overweight = rbind.fill(Overweight_Normal[!(Overweight_Normal$Group %like% "Normal"),],
                        Lean_Normal[Lean_Normal$Group == "Overweight",])
Overweight[is.na(Overweight)] <- 0
rownames(Overweight) = c(rownames(Overweight_Normal[!(Overweight_Normal$Group %like% "Normal"),],),
                         rownames(Lean_Normal[Lean_Normal$Group == "Overweight",]))

T2D_Obese = rbind.fill(Overweight_T2D[Overweight_T2D$Group %like% "Overweight",],
                       Lean_T2D[(Lean_T2D$Group %like% "Overweight"),])
T2D_Obese[is.na(T2D_Obese)] <- 0

rownames(T2D_Obese) = c(rownames(Overweight_T2D[Overweight_T2D$Group %like% "Overweight",]),
                        rownames(Lean_T2D[(Lean_T2D$Group %like% "Overweight"),]))

#Prediabetes#

######################################################## Hypertension Groups ######################################################## 

Hypertension = rbind.fill(Hypertension_start,
                          (Prospective_NAFLD[Prospective_NAFLD$Group %like% "Hypertension",] %>% .[!(.$Group %like% "Hypertension_Overweight"),] %>% .[!(.$Group %like% "Hypertension_Diabetes"),] %>% .[!(.$Group %like% "Pre_Hypertension"),]))
Hypertension[is.na(Hypertension)] <- 0
rownames(Hypertension) = c(rownames(Hypertension_start),
                           rownames(Prospective_NAFLD[Prospective_NAFLD$Group %like% "Hypertension",] %>% .[!(.$Group %like% "Hypertension_Overweight"),] %>% .[!(.$Group %like% "Hypertension_Diabetes"),] %>% .[!(.$Group %like% "Pre_Hypertension"),]))



Pre_Hypertension = rbind.fill(Pre_Hypertension_start,
                              Prospective_NAFLD[Prospective_NAFLD$Group %like% "Pre_Hypertension",]) %>% .[!(.$Group %like% "Pre_Hypertension_Diabetes"),] %>% .[!(.$Group %like% "Pre_Hypertension_Overweight_"),]
Pre_Hypertension[is.na(Pre_Hypertension)] <- 0
rownames(Pre_Hypertension) = c(rownames(Pre_Hypertension_start),
                               rownames(Prospective_NAFLD[Prospective_NAFLD$Group %like% "Pre_Hypertension",] %>% .[!(.$Group %like% "Pre_Hypertension_Diabetes"),] %>% .[!(.$Group %like% "Pre_Hypertension_Overweight_"),]) )


list_of_duplicated_Hypertension <- c("Hypertension_Project_nHF611710_Control")
Control_Hypertension = Control_Hypertension[-grep(paste(list_of_duplicated_Hypertension,collapse = "|"), rownames(Control_Hypertension)),]


######################################################## Lean T2D Groups ######################################################## 


T2D_Lean = rbind.fill(Lean_T2D[!(Lean_T2D$Group %like% "Overweight"),],
                      Overweight_T2D[!(Overweight_T2D$Group %like% "Overweight"),],
                      Prospective_NAFLD[Prospective_NAFLD$Group == "Diabetes_NAFLD_Free_Ultrasound",],
                      Prospective_NAFLD[Prospective_NAFLD$Group == "Pre_Hypertension_Diabetes_NAFLD_Free_Ultrasound",],
                      Prospective_NAFLD[Prospective_NAFLD$Group == "Hypertension_Diabetes_NAFLD_Free_Ultrasound",])
T2D_Lean[is.na(T2D_Lean)] <- 0

rownames(T2D_Lean) = c(rownames(Lean_T2D[!(Lean_T2D$Group %like% "Overweight"),]),
                       rownames(Overweight_T2D[!(Overweight_T2D$Group %like% "Overweight"),]),
                       rownames(Prospective_NAFLD[Prospective_NAFLD$Group == "Diabetes_NAFLD_Free_Ultrasound",]),
                       rownames(Prospective_NAFLD[Prospective_NAFLD$Group == "Pre_Hypertension_Diabetes_NAFLD_Free_Ultrasound",]),
                       rownames(Prospective_NAFLD[Prospective_NAFLD$Group == "Hypertension_Diabetes_NAFLD_Free_Ultrasound",]))

#Lean_Normal#

######################################################## Atherosclerosis Groups ######################################################## 

#Atherosclerosis#

list_of_anorexic_ <- c("Atherosclerosis_Project_N1266",
                       "Atherosclerosis_Project_N064",
                       "Atherosclerosis_Project_N201",
                       "Atherosclerosis_Project_N152")
Control_Atherosclerosis  = Control_Atherosclerosis[-grep(paste(list_of_anorexic_,collapse = "|"), rownames(Control_Atherosclerosis)),]


#save.image("Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/R_Datasets/Updated_NAFLD_O_L_Metaphlan3_Clinical_Data.RData")


list_samples <-
  data.frame(Sample_Number =
               c(nrow(NAFLD_L),
                 nrow(NAFLD_O),
                 nrow(Control_L),
                 nrow(Control_O),
                 nrow(Overweight),
                 nrow(T2D_Obese),
                 nrow(Prediabetes),
                 nrow(Hypertension),
                 nrow(Pre_Hypertension),
                 nrow(Control_Hypertension),
                 nrow(T2D_Lean),
                 nrow(Lean_Normal),
                 nrow(Atherosclerosis),
                 nrow(Control_Atherosclerosis)
                 ),
             Dataset = c("NAFLD-L",
                         "NAFLD-O",
                         "Control-L",
                         "Control-O",
                         "Overweight",
                         "T2D Overweight",
                         "Prediabetes",
                         "Hypertension",
                         "Pre_Hypertension",
                         "Control_Hypertension",
                         "T2D_Lean",
                         "Lean_Normal",
                         "Atherosclerosis",
                         "Control_Atherosclerosis"))
#write.csv(list_samples,"Documents/Metabolic_Diseases/Updated_Human3/NAFLD_O_L/List_of_Samples.csv")
