

setwd("D:/ALUS_GLOBI/lep_test")

#Lepidoptera test data


#lead various packages
#it may be useful to check to see if some of these aren't being used and eliminate them

library(dplyr)
library(tidyverse)
library(stringr)
library(janitor)
library(data.table)
library(ape)
library(phylotools)
library(RCurl)
library(Information)
library(mice)
library(stats)
library(doParallel)
library(doRNG)
library(lmtest)
library(matrixStats)
library(Metrics)
library(missForest)
library(mltools)
library(naniar)
library(nnet)
library(pastecs)
library(plotrix)
library(PVR)
library(rsample)
library(simputation)  
library(VIM)
library(viridis)
library(purrr)
library(picante)
#https://github.com/Matgend/TDIP
library(TDIP)
library(writexl)
source("Functions/DataHandling_Functions.R")
source("Functions/Imputation_Functions.R")
source("Functions/Phylo_Functions.R")
source("Functions/Simpute_Functions.R")
source("Functions/Additional_Functions.R")



#read in the complete trait data that I prepared
dfCC<-read.delim("Complete_Leptraits_for_analysis_no_dummy.tsv", sep ="\t") #875 rows

# The original data with the dummy encoded data 
# Just use to compare with the new data dfCC (Revised by JD 20251018)
dfCC1<-read.delim("Complete_Leptraits_for_analysis.tsv", sep ="\t") #875 rows

calculate_na_percentage(dfCC) #confirm there is no missing data already
#read in the tree
# treeName<-choose.files() #Leptree.newick
CCLeptree<-read.tree('Leptree.newick')
#Species in both tree and data frame have genus and species separated by underscore



# Look at column names of dfCC.
colnames(dfCC)
# Create vector of column names that contain taxonomic information.
taxCols <- colnames(dfCC)[1:4]
# Extract trait names.
traits <- setdiff(colnames(dfCC), taxCols)




# Calculate class balance  
# All variables except for the first four variables (revised by JD 20251020)
class_balance <- dfCC %>%
  dplyr::select(all_of(traits)) %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
  dplyr::count(Trait, Value, name = "Freq") %>%
  dplyr::group_by(Trait) %>%
  dplyr::mutate(Proportion = Freq / sum(Freq)) %>%
  dplyr::ungroup()

# Display the results
print(class_balance)


# Plot the class balance
ggplot(class_balance, aes(x = Trait, y = Proportion, fill = as.factor(Value))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Trait", y = "Proportion", fill = "Class (0 or 1)",
       title = "Class Balance for Categorical Traits") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



## Settings

# Specify parameters
missing_type <- "MCAR"  # MCAR or "MAR" or "MNAR", different ways of simulating missingness in data

#MCAR missing completely at random
#Vmissing at random
#MNAR missing not at random

start_missing <- 0.05 #start at 5%   # Revised by JD 20250817
increment <- 0.05 #go up by 5%%     # Revised by JD 20250817
end_missing <- 0.5 #end at 50%

cols <- 5:ncol(dfCC) #specify the columns to simulate missing data

## Run 5 repeats to calculate the average performance (Revised by JD 20250817)
set.seed(1)  # Each repeat uses a different seed (Revised by JD 20250817)

###----Seed (1)------###

## Step 1: Data preprocessing

# Simulate missingness
datasets_with_missing_data <- simulate_missingness(
  missing_type, dfCC, cols, start_missing, increment, end_missing, varying = TRUE
)

# Confirm missingness levels
missing_percentage_list <- lapply(datasets_with_missing_data, calculate_missing_percentage)
missing_percentage_list

# # Convert numeric traits to factor
l_missingData1 <- convert_all_numeric_to_factor(datasets_with_missing_data, traits)
# l_missingData1 <- convert_numeric_to_factor(datasets_with_missing_data, traits)


# # Select predictors
l_predictors_list1 <- purrr::map(l_missingData1, ~ apply_SelectPredictors(.x, traits))

# Prune tree to match dataset
ultraTree <- drop.tip(CCLeptree, CCLeptree$tip.label[!CCLeptree$tip.label %in% dfCC$Species])

# # Append phylogenetic eigenvectors
l_evs_list <- AppendEigenvectorsPVR_List(
  data_list = l_missingData1, vars = traits, tree = ultraTree, predictors_list = l_predictors_list1
)

# Extract dfRaw with eigenvectors
l_dfRawEV <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_dfRawEV[[i]] <- l_evs_list[[i]][[1]]
}

# Extract updated predictor lists with eigenvectors
l_EVPredictors_list <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_EVPredictors_list[[i]] <- l_evs_list[[i]][[2]]
}

# Sort dataframes by Species
l_missingData1 <- lapply(l_missingData1, function(df) df[order(df$Species), ])

m_list <- generate_m_list(l_missingData1)


## Step 2: CASE 1 (taxpaint,  MCAR, 5-50%)

#read in the average data for various genera
AverageLepTraits<-read.delim("AverageLepTraits_new.txt", sep="\t")

taxPaint <- replaceMissingValuesList(l_missingData1, AverageLepTraits, traits)

lapply(taxPaint, function(df) colSums(is.na(df[traits])))

calculate_na_percentage_list(taxPaint)

groups_to_check <- list(
  group1 = c("specialist", "generalist", "mixed")
)

taxPaintviolations<-check_column_groups(taxPaint, groups_to_check)
# #resolve violations
taxPaintcorrected <- resolve_violations(taxPaint, groups_to_check)
# #check again
check_column_groups(taxPaintcorrected, groups_to_check) #no violations



f1_results_1 <- calculate_macro_f1(dfCC, taxPaintcorrected, traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_1))
f1_df_1 <- convert_f1_list_to_df(f1_results_1, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_taxpaint <- f1_df_1



## Step 3:  CASE 2  (MissForest, no phylo, just using the package function,  MCAR, 5-50%)

#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_nophylo <- impute_missing_data_MissF(l_missingData1, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo

f1_results_2 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_nophylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_2))
f1_df_2 <- convert_f1_list_to_df(f1_results_2, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_nophylo <- f1_df_2



## Step 4 : CASE 3 (MissF, phylo, just using the package function,  MCAR, 5-50%)

#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_phylo <- impute_missing_data_MissF_phylo(l_dfRawEV, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo


f1_results_3 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_phylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_3))
f1_df_3 <- convert_f1_list_to_df(f1_results_3, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_phylo <- f1_df_3


write.csv(df_results_taxpaint, "Result/Seed_1_method1.csv", row.names = FALSE)
write.csv(df_results_MissF_nophylo, "Result/Seed_1_method2.csv", row.names = FALSE)
write.csv(df_results_MissF_phylo, "Result/Seed_1_method3.csv", row.names = FALSE)


set.seed(2)

###----Seed (2)------###

## Step 1: Data preprocessing

# Simulate missingness
datasets_with_missing_data <- simulate_missingness(
  missing_type, dfCC, cols, start_missing, increment, end_missing, varying = TRUE
)

# Confirm missingness levels
missing_percentage_list <- lapply(datasets_with_missing_data, calculate_missing_percentage)
missing_percentage_list

# # Convert numeric traits to factor
l_missingData1 <- convert_all_numeric_to_factor(datasets_with_missing_data, traits)
# l_missingData1 <- convert_numeric_to_factor(datasets_with_missing_data, traits)


# # Select predictors
l_predictors_list1 <- purrr::map(l_missingData1, ~ apply_SelectPredictors(.x, traits))

# Prune tree to match dataset
ultraTree <- drop.tip(CCLeptree, CCLeptree$tip.label[!CCLeptree$tip.label %in% dfCC$Species])

# # Append phylogenetic eigenvectors
l_evs_list <- AppendEigenvectorsPVR_List(
  data_list = l_missingData1, vars = traits, tree = ultraTree, predictors_list = l_predictors_list1
)

# Extract dfRaw with eigenvectors
l_dfRawEV <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_dfRawEV[[i]] <- l_evs_list[[i]][[1]]
}

# Extract updated predictor lists with eigenvectors
l_EVPredictors_list <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_EVPredictors_list[[i]] <- l_evs_list[[i]][[2]]
}

# Sort dataframes by Species
l_missingData1 <- lapply(l_missingData1, function(df) df[order(df$Species), ])

m_list <- generate_m_list(l_missingData1)

## Step 2: CASE 1 (TaxPaint, MCAR, 5–50%)

#read in the average data for various genera
AverageLepTraits<-read.delim("AverageLepTraits_new.txt", sep="\t")

taxPaint <- replaceMissingValuesList (l_missingData1, AverageLepTraits, traits)
calculate_na_percentage_list(taxPaint)

groups_to_check <- list(
  group1 = c("specialist", "generalist", "mixed")
)

taxPaintviolations<-check_column_groups(taxPaint, groups_to_check)
#resolve violations
taxPaintcorrected <- resolve_violations(taxPaint, groups_to_check)
#check again
check_column_groups(taxPaintcorrected, groups_to_check) #no violations


f1_results_1 <- calculate_macro_f1(dfCC, taxPaintcorrected, traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_1))
f1_df_1 <- convert_f1_list_to_df(f1_results_1, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_taxpaint <- f1_df_1


## Step 3: CASE 2 (MissForest, no phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_nophylo <- impute_missing_data_MissF(l_missingData1, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo

f1_results_2 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_nophylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_2))
f1_df_2 <- convert_f1_list_to_df(f1_results_2, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_nophylo <- f1_df_2


## Step 4: CASE 3 (MissForest, with phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_phylo <- impute_missing_data_MissF_phylo(l_dfRawEV, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo


f1_results_3 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_phylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_3))
f1_df_3 <- convert_f1_list_to_df(f1_results_3, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_phylo <- f1_df_3


write.csv(df_results_taxpaint, "Result/Seed_2_method1.csv", row.names = FALSE)
write.csv(df_results_MissF_nophylo, "Result/Seed_2_method2.csv", row.names = FALSE)
write.csv(df_results_MissF_phylo, "Result/Seed_2_method3.csv", row.names = FALSE)


set.seed(3)

###----Seed (3)------###

## Step 1: Data preprocessing

# Simulate missingness
datasets_with_missing_data <- simulate_missingness(
  missing_type, dfCC, cols, start_missing, increment, end_missing, varying = TRUE
)

# Confirm missingness levels
missing_percentage_list <- lapply(datasets_with_missing_data, calculate_missing_percentage)
missing_percentage_list

# # Convert numeric traits to factor
l_missingData1 <- convert_all_numeric_to_factor(datasets_with_missing_data, traits)
# l_missingData1 <- convert_numeric_to_factor(datasets_with_missing_data, traits)


# # Select predictors
l_predictors_list1 <- purrr::map(l_missingData1, ~ apply_SelectPredictors(.x, traits))

# Prune tree to match dataset
ultraTree <- drop.tip(CCLeptree, CCLeptree$tip.label[!CCLeptree$tip.label %in% dfCC$Species])

# # Append phylogenetic eigenvectors
l_evs_list <- AppendEigenvectorsPVR_List(
  data_list = l_missingData1, vars = traits, tree = ultraTree, predictors_list = l_predictors_list1
)

# Extract dfRaw with eigenvectors
l_dfRawEV <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_dfRawEV[[i]] <- l_evs_list[[i]][[1]]
}

# Extract updated predictor lists with eigenvectors
l_EVPredictors_list <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_EVPredictors_list[[i]] <- l_evs_list[[i]][[2]]
}

# Sort dataframes by Species
l_missingData1 <- lapply(l_missingData1, function(df) df[order(df$Species), ])


m_list <- generate_m_list(l_missingData1)


## Step 2: CASE 1 (TaxPaint, MCAR, 5–50%)

#read in the average data for various genera
AverageLepTraits<-read.delim("AverageLepTraits_new.txt", sep="\t")

taxPaint <- replaceMissingValuesList (l_missingData1, AverageLepTraits, traits)
calculate_na_percentage_list(taxPaint)

groups_to_check <- list(
  group1 = c("specialist", "generalist", "mixed")
)


taxPaintviolations<-check_column_groups(taxPaint, groups_to_check)
#resolve violations
taxPaintcorrected <- resolve_violations(taxPaint, groups_to_check)
#check again
check_column_groups(taxPaintcorrected, groups_to_check) #no violations

# Factor to numeric (Added by JD 20251020)
taxPaintcorrected <- map(taxPaintcorrected, ~ .x %>%
                           mutate(across(-c(Species, Family, Genus, verbatimSpecies),
                                         ~ as.integer(as.character(.)))
                           )
)


f1_results_1 <- calculate_macro_f1(dfCC, taxPaintcorrected, traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_1))
f1_df_1 <- convert_f1_list_to_df(f1_results_1, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_taxpaint <- f1_df_2


## Step 3: CASE 2 (MissForest, no phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_nophylo <- impute_missing_data_MissF(l_missingData1, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo

f1_results_2 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_nophylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_2))
f1_df_2 <- convert_f1_list_to_df(f1_results_2, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_nophylo <- f1_df_2



## Step 4: CASE 3 (MissForest, with phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_phylo <- impute_missing_data_MissF_phylo(l_dfRawEV, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo


f1_results_3 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_phylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_3))
f1_df_3 <- convert_f1_list_to_df(f1_results_3, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_phylo <- f1_df_3

write.csv(df_results_taxpaint, "Result/Seed_3_method1.csv", row.names = FALSE)
write.csv(df_results_MissF_nophylo, "Result/Seed_3_method2.csv", row.names = FALSE)
write.csv(df_results_MissF_phylo, "Result/Seed_3_method3.csv", row.names = FALSE)



set.seed(4)
###----Seed (4)------###

## Step 1: Data preprocessing

# Simulate missingness
datasets_with_missing_data <- simulate_missingness(
  missing_type, dfCC, cols, start_missing, increment, end_missing, varying = TRUE
)

# Confirm missingness levels
missing_percentage_list <- lapply(datasets_with_missing_data, calculate_missing_percentage)
missing_percentage_list

# # Convert numeric traits to factor
l_missingData1 <- convert_all_numeric_to_factor(datasets_with_missing_data, traits)
# l_missingData1 <- convert_numeric_to_factor(datasets_with_missing_data, traits)



# # Select predictors
l_predictors_list1 <- purrr::map(l_missingData1, ~ apply_SelectPredictors(.x, traits))

# Prune tree to match dataset
ultraTree <- drop.tip(CCLeptree, CCLeptree$tip.label[!CCLeptree$tip.label %in% dfCC$Species])

# # Append phylogenetic eigenvectors
l_evs_list <- AppendEigenvectorsPVR_List(
  data_list = l_missingData1, vars = traits, tree = ultraTree, predictors_list = l_predictors_list1
)

# Extract dfRaw with eigenvectors
l_dfRawEV <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_dfRawEV[[i]] <- l_evs_list[[i]][[1]]
}

# Extract updated predictor lists with eigenvectors
l_EVPredictors_list <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_EVPredictors_list[[i]] <- l_evs_list[[i]][[2]]
}

# Sort dataframes by Species
l_missingData1 <- lapply(l_missingData1, function(df) df[order(df$Species), ])


m_list <- generate_m_list(l_missingData1)


## Step 2: CASE 1 (TaxPaint, MCAR, 5–50%)

#read in the average data for various genera
AverageLepTraits<-read.delim("AverageLepTraits_new.txt", sep="\t")

taxPaint <- replaceMissingValuesList (l_missingData1, AverageLepTraits, traits)
calculate_na_percentage_list(taxPaint)

groups_to_check <- list(
  group1 = c("specialist", "generalist", "mixed")
)


taxPaintviolations<-check_column_groups(taxPaint, groups_to_check)
#resolve violations
taxPaintcorrected <- resolve_violations(taxPaint, groups_to_check)
#check again
check_column_groups(taxPaintcorrected, groups_to_check) #no violations


f1_results_1 <- calculate_macro_f1(dfCC, taxPaintcorrected, traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_1))
f1_df_1 <- convert_f1_list_to_df(f1_results_1, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_taxpaint <- f1_df_1


## Step 3: CASE 2 (MissForest, no phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_nophylo <- impute_missing_data_MissF(l_missingData1, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo

f1_results_2 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_nophylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_2))
f1_df_2 <- convert_f1_list_to_df(f1_results_2, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_nophylo <- f1_df_2


## Step 4: CASE 3 (MissForest, with phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_phylo <- impute_missing_data_MissF_phylo(l_dfRawEV, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo


f1_results_3 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_phylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_3))
f1_df_3 <- convert_f1_list_to_df(f1_results_3, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_phylo <- f1_df_3


write.csv(df_results_taxpaint, "Result/Seed_4_method1.csv", row.names = FALSE)
write.csv(df_results_MissF_nophylo, "Result/Seed_4_method2.csv", row.names = FALSE)
write.csv(df_results_MissF_phylo, "Result/Seed_4_method3.csv", row.names = FALSE)



set.seed(5)

###----Seed (5)------###

## Step 1: Data preprocessing

# Simulate missingness
datasets_with_missing_data <- simulate_missingness(
  missing_type, dfCC, cols, start_missing, increment, end_missing, varying = TRUE
)

# Confirm missingness levels
missing_percentage_list <- lapply(datasets_with_missing_data, calculate_missing_percentage)
missing_percentage_list

# # Convert numeric traits to factor
l_missingData1 <- convert_all_numeric_to_factor(datasets_with_missing_data, traits)
# l_missingData1 <- convert_numeric_to_factor(datasets_with_missing_data, traits)


# # Select predictors
l_predictors_list1 <- purrr::map(l_missingData1, ~ apply_SelectPredictors(.x, traits))

# Prune tree to match dataset
ultraTree <- drop.tip(CCLeptree, CCLeptree$tip.label[!CCLeptree$tip.label %in% dfCC$Species])

# # Append phylogenetic eigenvectors
l_evs_list <- AppendEigenvectorsPVR_List(
  data_list = l_missingData1, vars = traits, tree = ultraTree, predictors_list = l_predictors_list1
)

# Extract dfRaw with eigenvectors
l_dfRawEV <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_dfRawEV[[i]] <- l_evs_list[[i]][[1]]
}

# Extract updated predictor lists with eigenvectors
l_EVPredictors_list <- CreateNamedList(listLength = length(l_evs_list), elementNames = 1:length(l_evs_list))
for (i in seq_along(l_evs_list)) {
  l_EVPredictors_list[[i]] <- l_evs_list[[i]][[2]]
}

# Sort dataframes by Species
l_missingData1 <- lapply(l_missingData1, function(df) df[order(df$Species), ])

m_list <- generate_m_list(l_missingData1)

## Step 2: CASE 1 (TaxPaint, MCAR, 5–50%)

#read in the average data for various genera
AverageLepTraits<-read.delim("AverageLepTraits_new.txt", sep="\t")

taxPaint <- replaceMissingValuesList (l_missingData1, AverageLepTraits, traits)
calculate_na_percentage_list(taxPaint)

groups_to_check <- list(
  group1 = c("specialist", "generalist", "mixed")
)


taxPaintviolations<-check_column_groups(taxPaint, groups_to_check)
#resolve violations
taxPaintcorrected <- resolve_violations(taxPaint, groups_to_check)
#check again
check_column_groups(taxPaintcorrected, groups_to_check) #no violations



f1_results_1 <- calculate_macro_f1(dfCC, taxPaintcorrected, traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_1))
f1_df_1 <- convert_f1_list_to_df(f1_results_1, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_taxpaint <- f1_df_1


## Step 3: CASE 2 (MissForest, no phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_nophylo <- impute_missing_data_MissF(l_missingData1, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo

f1_results_2 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_nophylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_2))
f1_df_2 <- convert_f1_list_to_df(f1_results_2, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_nophylo <- f1_df_2



## Step 4: CASE 3 (MissForest, with phylogeny, package function, MCAR, 5–50%)
#will need more succinct code for this

#the issue here is that impute_missing_data_RF imputes and calculates the summary stats
#we need so there isn't any time to intervene and change FlightDuration

# Apply the function to your list of data frames
imputed_results_MissF_phylo <- impute_missing_data_MissF_phylo(l_dfRawEV, dfCC)
#lets not worry about converting from the individual month data to FlightDuration for now

#the results of this function makes a list for each data frame
#ximp holds the dataframe with imputed data
#OOBerror holds the error data for each column, although these aren't named
#error holds the

#lets parse the data out of imputed_results_RF_nophyo


f1_results_3 <- calculate_macro_f1_missforest(dfCC, imputed_results_MissF_phylo , traits, m_list, impNum=1)

# Convert results to dataframe
dataset_names <- paste0("Dataset_", seq_along(f1_results_3))
f1_df_3 <- convert_f1_list_to_df(f1_results_3, dataset_names)

# Final evaluation results (Macro-F1 only)
df_results_MissF_phylo <- f1_df_3


write.csv(df_results_taxpaint, "Result/Seed_5_method1.csv", row.names = FALSE)
write.csv(df_results_MissF_nophylo, "Result/Seed_5_method2.csv", row.names = FALSE)
write.csv(df_results_MissF_phylo, "Result/Seed_5_method3.csv", row.names = FALSE)



