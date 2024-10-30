# Load Packages
library(dplyr)
library(SurrogateRank)
library(readr)
library(parallel)
library(pbmcapply)
library(caTools)
library(ggplot2)
library(Metrics)
library(xtable)
set.seed(18062024)

# set parameters
alpha = 0.05 # Nominal significance level
beta = 0.2 # (1 - beta)*100% = desired power 

# Load PREVAC data
prevac_df = get(load(file = './data/processed/prevac_merged.RData')) 

# Select gene expression at day 7 following injection with rVSV or placebo + clinical variables
df = prevac_df %>% mutate(treatment = recode_factor(arm_short, "rVSV" = 1, "rVSV2" = 1, "Ad26MVA" = 1, "placebo" = 0)) %>% 
  filter(is.na(Ab_M12) == F, Timepoint == "J7", arm_short %in% (c("placebo", "rVSV"))) %>% 
  select(PID, Timepoint, arm_short, treatment, Ab_M12, A1BG:ZZZ3)

# Total number of genes in the data
n_genes = ncol(df %>% select(A1BG:ZZZ3))

# Gene names
gene_names <- names(df[,-c(1:5)])

# Sample splitting - 66% for screening, balanced treatment groups
train_index = sample.split(df$treatment, SplitRatio = 0.66, group = NULL)
df_train = df[train_index,]
df_test = df[!train_index,]

# save the indices for later use to ensure we validate on the correct data
save(train_index, file = "./output/train_index_final.RData")

source("./R/rise_screen.R")

X = df_train %>% select(A1BG:ZZZ3) %>% as.matrix()
Y = df_train %>% select(Ab_M12) %>% as.matrix()
A = df_train %>% select(treatment) %>% as.matrix()


##--TABLE S2 : impact of epsilon (beta) on validation--##

beta_grid = seq(0.95,0.5,-0.05)
res_screen_sensitivity = list()

for (i in 1:length(beta_grid)){
  res_screen_sensitivity[[i]] = rise_screen(Y = Y, X = X, A = A,
                                     reference = "0", alpha = 0.05,
                                     power_desired = beta_grid[i], p_correction = "BH",
                                     cores = 10)
}

# save the results
save(res_screen_sensitivity, file = './output/application_results_final/RISE_screening_supplementary.Rds')
res = get(load('./output/application_results_final/RISE_screening_supplementary.Rds'))

significant_genes = list()
n_significant_genes = rep(0,length(beta_grid))
for (i in 1:length(beta_grid)){
  significant_genes[[i]] = res_screen_sensitivity[[i]]$significant_markers
  n_significant_genes[i] = length(significant_genes[[i]])
}


# validation
source("./R/rise_validate.R")
validate_results_sensitivity = list()
for (i in 1:length(beta_grid)){
  
  markers = significant_genes[[i]]
  X_validate = df_test %>% select(all_of(markers)) %>% as.matrix()
  Y_validate = df_test %>% select(Ab_M12) %>% as.matrix()
  A_validate = df_test %>% select(treatment) %>% as.matrix()
  weights = res_screen_sensitivity[[i]][["weights"]]
  validate_results_sensitivity[[i]] = rise_validate(
                                 Y_validate = Y_validate,
                                 X_validate = X_validate,
                                 A_validate = A_validate,
                                 individual = T,
                                 power_desired = 0.8,
                                 plot = T,
                                 markers = markers, 
                                 weights = weights)
  print(validate_results_sensitivity[[i]][["validation_results"]])
  print(beta_grid[i])
  print(validate_results_sensitivity[[i]]$plot)
  
}

output = bind_cols("power_screening" = beta_grid, "power_validation" = 0.8, "n_predictors" =  n_significant_genes, bind_rows(lapply(validate_results_sensitivity, function(x) x$validation_results))) %>% 
  as.data.frame() %>% 
  mutate(delta_ci = paste0(round(delta, 3), " (", round(ci_lower, 3), ", ", round(ci_upper, 3), ")"),
         sd_round = round(sd, 3)) %>% 
  select(-marker)

output

results_sensitivity_latex = output %>% select(power_screening, power_validation, n_predictors, delta_ci, sd_round, p_unadjusted)
myxtable <- function(x, ...) xtable(apply(x, 2, as.character), ...)
print(myxtable(results_sensitivity_latex), include.rownames = FALSE)
