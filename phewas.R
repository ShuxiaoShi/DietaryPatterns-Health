
library(data.table)
library(dplyr)
library(myPheWAS)

args = commandArgs(trailingOnly=TRUE)

icd4 <- fread("icd10.csv")
icd4$code <- as.character(icd4$code)
phenotypes <- createPhenotypes(icd4,min.code.count=2)


cov <- fread("cov.csv")
#cov <- cov %>% rename(id = n_eid)

cov$sex<-as.factor(cov$sex)
cov$edu_level<-as.factor(cov$edu_level)
#cov$race_group<-as.factor(cov$race_group)
cov$smoking<-as.factor(cov$smoking)
cov$drinking<-as.factor(cov$drinking)
#cov$townsend_score_group<-as.factor(cov$townsend_score_group)


#diet_pca <- read.csv("diet_index_quintile.csv")
diet_pca <- read.csv("tradpca_quintile.csv")
#diet_pca <- diet_pca %>% rename(id = eid)

diet_pca2 <-diet_pca %>% select(n_eid,args[1])

diet_cov <-diet_pca %>% select(n_eid, pc3_q_2, pc3_q_3, pc3_q_4) 

diet_cov2 <- left_join(diet_pca2, diet_cov)

data <- left_join(left_join(phenotypes,diet_cov2),cov)


write.csv(data, file = "intermediate.csv", row.names = FALSE)

results <- phewas(phenotypes=names(phenotypes)[-1], genotypes=names(diet_pca2)[-1], data = data, covariates=c("pc3_q_2","pc3_q_3","pc3_q_4","mean_energy","sex","age_attend","edu_level","bmi","townsend_score","drinking","smoking","sum_physical_met"),alpha = 1, cores=20)
results_d <- addPhecodeInfo(results)
write.csv(results_d, file = args[2], row.names = FALSE)

#filter
filter_results <-  results_d %>% filter(p < 7.352941e-05 )
ffq_rank2 <-arrange(filter_results,p)
write.csv(ffq_rank2, file = args[3], row.names = FALSE)

###fdr
results_d$p_fdr <- p.adjust(results_d$p)

diet_filter <- results_d %>% filter(p_fdr < 0.05)
write.csv(ffq_rank2, file = args[4], row.names = FALSE, quote = FALSE)
