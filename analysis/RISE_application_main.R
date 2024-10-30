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

res_screen = rise_screen(Y = Y, X = X, A = A,
                         reference = "0", alpha = 0.05,
                         power_desired = 0.8, p_correction = "BH",
                         cores = 10)


# save the results
save(res_screen, file = './output/application_results_final/RISE_screening_final.Rds')
res_screen = get(load('./output/application_results_final/RISE_screening_final.Rds'))

# View the top genes
res_screen$results %>% filter(p_adjusted < 0.05) %>% arrange(p_adjusted)
# what are the significant genes after multiple testing correction
genes = res_screen$significant_markers

##--TABLE 1 : screening results--##
# Paste directly into latex
res_latex = res_screen$results %>% mutate(delta_ci = paste0(round(delta, 3), " (", round(ci_lower, 3), ", ", round(ci_upper, 3), ")"),
                                          sd_round = round(sd, 3))
paste_results = res_latex %>% filter(p_adjusted < 0.05) %>% arrange(p_adjusted) %>% 
  select(marker, delta_ci, sd_round, p_unadjusted,p_adjusted)
myxtable <- function(x, ...) xtable(apply(x, 2, as.character), ...)
print(myxtable(paste_results), include.rownames = FALSE)

# load the training indices if necessary to ensure we validate on the correct data
train_index = get(load('./output/train_index_final.RData'))
df_train = df[train_index,]
df_test = df[!train_index,]

##--VALIDATION : standardised, weighted sum of predictors--##

markers = genes
X_validate = df_test %>% select(all_of(genes)) %>% as.matrix()
Y_validate = df_test %>% select(Ab_M12) %>% as.matrix()
A_validate = df_test %>% select(treatment) %>% as.matrix()
individual = T
power_desired = 0.8
epsilon = NULL
plot = TRUE
weights = res_screen[["weights"]]

source("./R/rise_validate.R")

markers = genes
X_validate = df_test %>% select(all_of(genes)) %>% as.matrix()
Y_validate = df_test %>% select(Ab_M12) %>% as.matrix()
A_validate = df_test %>% select(treatment) %>% as.matrix()
individual = T
power_desired = 0.8
epsilon = NULL
plot = TRUE
weights = res_screen[["weights"]]

validate_results = rise_validate(Y_validate = Y_validate,
                                 X_validate = X_validate,
                                 A_validate = A_validate,
                                 individual = T,
                                 power_desired = 0.8,
                                 plot = T,
                                 markers = markers, 
                                 weights = res_screen[["weights"]])


validate_results$plot

fig4 = validate_results$plot +
  xlab("Month 12 EBOV-GP Binding Antibodies (rank)") +
  ylab(expression(gamma[S]~"(rank)")) +
  ggtitle("Ranks of primary outcome vs new surrogate in validation data") +
  guides(col = guide_legend(title="Treatment")) +
  scale_color_manual(labels = c("Placebo","Vaccine"), values = c("#1D8A99","#c1121f")) +
  theme_minimal(base_size = 22) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white'))

fig4

ggsave(plot = fig4, 
       filename = "./output/figures/main/figure4.pdf", 
       units = "cm", width = 35, height = 25, dpi = 500)

##--TABLE 2 : validation of individual markers--##
# get the indices of the significant genes
gene_ind = which(colnames(df_test) %in%  genes)
# apply the test to the individual genes in the validation data
results_validate = list()
# gene_ind = which(colnames(df_test) %in%  genes)
# gene_order = colnames(df_test)[gene_ind]

for (gene in genes){
  yone = df_test %>% filter(treatment == 1) %>% select(Ab_M12) %>% sapply(as.numeric)
  yzero = df_test %>% filter(treatment == 0) %>% select(Ab_M12) %>% sapply(as.numeric)
  sone = df_test %>% filter(treatment == 1) %>% select(gene) %>% sapply(as.numeric)
  szero = df_test %>% filter(treatment == 0) %>% select(gene) %>% sapply(as.numeric)
  
  ss.test = test.surrogate(yone = yone, yzero = yzero, sone = sone, szero = szero, power.want.s = 1-beta)
  # ss.test = test.surrogate(yone = yone, yzero = yzero, sone = sone, szero = szero, epsilon = get.u - 0.5)
  p = pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta)
  
  res = c(ss.test$delta.estimate, ss.test$sd.u.s, ss.test$epsilon.used, p)
  results_validate = append(results_validate,list(res))
  print(gene)
}

# manipulate the results
results_vec_validate = cbind(genes, do.call(rbind, results_validate)) %>% as.data.frame()
results_vec_validate[,c(2:5)] = results_vec_validate[,c(2:5)] %>% sapply(as.numeric) 
colnames(results_vec_validate) = c("gene", "delta","sd", "epsilon", "p_unadjusted")

# Compute p values and confidence intervals
results_vec_validate = results_vec_validate %>%
  mutate(
    z_alpha = qnorm(1 - alpha),
    ci_upper = delta + z_alpha * sd,
    ci_lower = -1,
    p_bonf = p.adjust(p_unadjusted, method = "bonferroni"),
    p_bh = p.adjust(p_unadjusted, method = "BH"),
    delta_ci = paste0(round(delta, 3), " (", round(ci_lower, 3), ", ", round(ci_upper, 3), ")"),
    sd_round = round(sd, 3)
  ) %>%
  select(-z_alpha) %>% as.data.frame()

results_vec_validate  %>% arrange(ci_upper) %>% select(gene, delta,ci_upper, ci_lower, sd, epsilon, p_unadjusted, p_bonf, p_bh)
save(results_vec_validate, file = './output/application_results_final/RISE_validation_final.Rds')


# Paste directly into latex for table 2
paste_results = results_vec_validate %>% select(gene, delta_ci, sd_round, p_unadjusted,p_bh)
myxtable <- function(x, ...) xtable(apply(x, 2, as.character), ...)
print(myxtable(paste_results), include.rownames = FALSE)


##--Figure 5 : gene plot--##
## We plot the log of delta against the log of delta divided by its standard deviation to examine visually
## if any groups appear
## Only take the first 100 genes for visual clarity
results_mat = results_vec %>% arrange(log(delta)) %>% head(100)

fig5 = results_mat %>% ggplot(aes(x = log(delta/sd), y = -log(delta), label = gene)) + 
  geom_point() + 
  theme_minimal() +
  geom_text(hjust=-0.1, vjust=0, size = 4) +
  ylab(expression(-log(delta))) +
  xlab(expression(log(delta/sigma))) + 
  ggtitle(expression("-log("*delta*") vs log("*delta*"/"*sigma*")"))+
  theme(title = element_text(size = 27, face = "bold"),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 25)) +
  theme_minimal(base_size = 22) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white'))

fig5 

ggsave(plot = fig5,"./output/figures/main/figure5.pdf", units = "cm", width = 20, height = 15)


