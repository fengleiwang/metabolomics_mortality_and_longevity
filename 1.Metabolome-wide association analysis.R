library(tidyverse)
library(dplyr)
library(survival)

load("~/Longevity/readin20221218.RData")
gdata::keep(f, sure=TRUE)

nhs1 <- read.csv("~/Longevity/n1cox.csv", header = T)
nhs2 <- read.csv("~/Longevity/n2cox.csv", header = T)
hpfs <- read.csv("~/Longevity/hpcox.csv", header = T)

nhs1$cohort="nh"
nhs2$cohort="n2"
hpfs$cohort="hp"

pool=bind_rows(nhs1, nhs2)
pool=bind_rows(pool, hpfs)

check=data.frame(var=colnames(pool[,1:396]), per=round(colMeans(is.na(pool[,1:396])), 4))
sum(check$per<0.25) #243

met243=as.character(check[check$per<0.25,"var"])
temp=pool[, (colnames(pool) %in% met243)]

# excluded metabolites
f_153=f[((rownames(f) %in% check$var) & !(rownames(f) %in% met243)), ] 
table(f_153$mean_icc<0.3) # 8
table(f[rownames(f) %in% met243,]$mean_icc<0.3) # 0

# first do the log-z score transformation within run_id, then impute
temp=pool[, colnames(pool) %in% met243]

# calculate skewness
library(moments)
metskewness=data.frame(hmdb=colnames(temp), skew=NA)
for (i in 1:243){
  metskewness[i,2]=skewness(temp[,i], na.rm=T)
}

# only log transform -2<skew<2, PMID: 23495371
sum(abs(metskewness$skew)>2) # 138
metnormal=as.character(metskewness %>% filter(abs(skew)<2) %>% pull(hmdb))
metskew=as.character(metskewness %>% filter(abs(skew)>2) %>% pull(hmdb))
temp_normal=temp %>% select(all_of(metnormal))
temp_skew=temp %>% select(all_of(metskew))

temp_skew_ln=log(temp_skew)
temp_new=cbind(temp_normal, temp_skew_ln)
temp_new$run_id=pool$run_id

# transform to z score
library(gtools)
z_score = function(metab) {
  if (length(metab)<1) {
    return(metab)
  } else {
    
    avg = mean(metab, na.rm = TRUE)
    std = sd(metab, na.rm = TRUE)
    
    for (i in 1:length(metab)) {
      if (!is.na(metab[i])) {
        metab[i] = (metab[i] - avg)/std
      } else {
        metab[i] = NA
      }
    }
    return(metab)
  }
}

z_tranb = function(dataset, metab_col) {
  z.data = dataset
  for (j in 1:metab_col) {
    # the ave function makes the z_score function operate on column j, stratifed by run_id
    z.data[,j] = ave(dataset[,j], dataset$run_id, FUN = z_score)
  }
  return(z.data)
}

temp_new = z_tranb(temp_new, 243)

# random forest imputaion
library(missRanger)
z=missRanger(temp_new[, 1:243], num.trees = 100, seed=1234)

save.image("~/Longevity/202212/mortality_met_1.RData")

pool2=cbind(z, pool[, 397:455])
check2=data.frame(var=colnames(pool2), per=round(colMeans(is.na(pool2)), 4))

factor=c("white", "fast", "asp", "mvt", "diab", "htn", "hcl", "hbpmed", "cholmed", "dbmed",
         "mifh", "dbfh", "smk", "bmicat", "cacoall", "endpoint")
for (i in 1:length(factor)){
  pool2[, factor[i]]=as.factor(pool2[, factor[i]])
}

formu_model1 <- list()
fit_model1 <- list()
formu_model2 <- list()
fit_model2 <- list()
formu_model3 <- list()
fit_model3 <- list()

for (i in 1:length(met243)){
  formu_model1[[i]] <- as.formula(paste0("Surv(tdead, dead) ~ ", met243[i],
                                         " + agemo +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model1[[i]] <- coxph(formu_model1[[i]], data = pool2)
  
  formu_model2[[i]] <- as.formula(paste0("Surv(tdead, dead) ~ ", met243[i],
                                         " + agemo + white + fast + mvt + diab + htn + hcl + hbpmed +
                                         cholmed + smk + bmi + act + calor + alco +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model2[[i]] <- coxph(formu_model2[[i]], data = pool2)
  
  formu_model3[[i]] <- as.formula(paste0("Surv(tdead, dead) ~ ", met243[i],
                                         " + agemo + white + fast + mvt + diab + htn + hcl + hbpmed +
                                         cholmed + smk + bmi + act + calor + alco + ahei_noal +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model3[[i]] <- coxph(formu_model3[[i]], data = pool2)
  
}

#GET RESULTS for death
results_model1 <- data.frame(matrix(NA, length(met243), 6))
rownames(results_model1) <- met243
colnames(results_model1)<-c("beta", "HR", "LCI", "UCI", "HR (95%CI)", "P_value")
results_model2=results_model1
results_model3=results_model1

for (i in 1:length(met243)){
  get_data_sum <- summary(fit_model1[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_model1[i,1] <- round(get_data_sum$coefficients[1,1], 10)
  results_model1[i,2] <- round(hr_data[1], 2)
  results_model1[i,3] <- round(hr_data[2], 2)
  results_model1[i,4] <- round(hr_data[3], 2)
  results_model1[i,6] <- round(get_data_sum$coefficients[1,5], 10)  # Pvalue
  results_model1[i,5] <- paste(format(round(results_model1[i,2], 2), nsmall=2), " (",
                               format(round(results_model1[i,3], 2), nsmall=2), ", ",
                               format(round(results_model1[i,4], 2), nsmall=2), ")", sep="")
  
  get_data_sum <- summary(fit_model2[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_model2[i,1] <- round(get_data_sum$coefficients[1,1], 10)
  results_model2[i,2] <- round(hr_data[1], 2)
  results_model2[i,3] <- round(hr_data[2], 2)
  results_model2[i,4] <- round(hr_data[3], 2)
  results_model2[i,6] <- round(get_data_sum$coefficients[1,5], 10)  # Pvalue
  results_model2[i,5] <- paste(format(round(results_model2[i,2], 2), nsmall=2), " (",
                               format(round(results_model2[i,3], 2), nsmall=2), ", ",
                               format(round(results_model2[i,4], 2), nsmall=2), ")", sep="")
  get_data_sum <- summary(fit_model3[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_model3[i,1] <- round(get_data_sum$coefficients[1,1], 10)
  results_model3[i,2] <- round(hr_data[1], 2)
  results_model3[i,3] <- round(hr_data[2], 2)
  results_model3[i,4] <- round(hr_data[3], 2)
  results_model3[i,6] <- round(get_data_sum$coefficients[1,5], 10)  # Pvalue
  results_model3[i,5] <- paste(format(round(results_model3[i,2], 2), nsmall=2), " (",
                               format(round(results_model3[i,3], 2), nsmall=2), ", ",
                               format(round(results_model3[i,4], 2), nsmall=2), ")", sep="")
}

met_name=read_csv("~/Longevity/met_name.csv")
model1=as.data.frame(results_model1[,c(1,5,6)])
colnames(model1)=c("mod1_beta", "mod1_HR", "mod1_p")
model1$mod1_p_fdr=p.adjust(model1$mod1_p, method="BH", n=length(met243))
model1$mod1_p_bon=p.adjust(model1$mod1_p, method="bonferroni", n=length(met243))
model1$HMDB=met243

model2=as.data.frame(results_model2[,c(1,5,6)])
colnames(model2)=c("mod2_beta", "mod2_HR", "mod2_p")
model2$mod2_p_fdr=p.adjust(model2$mod2_p, method="BH", n=length(met243))
model2$mod2_p_bon=p.adjust(model2$mod2_p, method="bonferroni", n=length(met243))
model2$HMDB=met243

model3=as.data.frame(results_model3[,c(1,5,6)])
colnames(model3)=c("mod3_beta", "mod3_HR", "mod3_p")
model3$mod3_p_fdr=p.adjust(model3$mod3_p, method="BH", n=length(met243))
model3$mod3_p_bon=p.adjust(model3$mod3_p, method="bonferroni", n=length(met243))
model3$HMDB=met243

results=merge(met_name, model1, by="HMDB")
results=merge(results, model2, by="HMDB")
results=merge(results, model3, by="HMDB")
for (i in c("mod1_p", "mod1_p_fdr", "mod1_p_bon",
            "mod2_p", "mod2_p_fdr", "mod2_p_bon",
            "mod3_p", "mod3_p_fdr", "mod3_p_bon")){
  results[,i]=round(results[,i], 3)
}
