library(readr)
load("~/Longevity/202212/mortality_met_2.RData")
met_name=read_csv("~/Longevity/met_name.csv")
met=pool2[, met_name$HMDB]
gdata::keep(pool2, met_name, met, sure=T)

library(WGCNA)
library(doParallel)
set.seed(1111)
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(met, dataIsExpr = TRUE, powerVector = powers, corFnc = cor, corOptions = list(use = 'p'),
                      networkType = "signed",RsquaredCut = 0.8)
# Plot the results
par(mfrow = c(1,2), oma = c(1,1,1,1), mar = c(4,4,4,1))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n", 
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = sft$powerEstimate; softPower
adj=adjacency(met,type = "signed", power = softPower)
TOM=TOMsimilarityFromExpr(met, networkType = "signed", TOMType = "signed", 
                          power = softPower, corType  = "pearson")

colnames(TOM) = met_name$HMDB
rownames(TOM) = met_name$HMDB
dissTOM=1-TOM

library(flashClust)
geneTree = hclust(as.dist(dissTOM), method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3)

# Set the minimum module size
minModuleSize = 10

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,  
                            distM = dissTOM,
                            deepSplit = 2,
                            pamRespectsDendro = F,
                            method="hybrid", 
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# blue     brown     green turquoise    yellow 
# 49        24        10       142        18 

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# plot the Topological overlap map of metaboltie modules
TOMna <- TOM
dissTOMna <- dissTOM
diag(TOMna) = NA
diag(dissTOMna) = NA

hhh=data.frame(HMDB=met_name$HMDB, module=dynamicMods, colour=dynamicColors)
hhh=merge(hhh, met_name, by="HMDB")
datME=moduleEigengenes(met, dynamicColors)$eigengenes

# standardize the modules
modules=colnames(datME)
for (i in 1:length(modules)) {
  temp=modules[i]
  avg=mean(datME[,temp])
  sd=sd(datME[,temp])
  datME[, temp]=(datME[, temp]-avg)/sd
}

pool3=cbind(pool2, datME)

# Cox
formu_model1 <- list()
fit_model1 <- list()
formu_model2 <- list()
fit_model2 <- list()
formu_model3 <- list()
fit_model3 <- list()

library(survival)
for (i in 1:length(modules)){
  formu_model1[[i]] <- as.formula(paste0("Surv(tdead, dead) ~ ", modules[i],
                                         " + agemo + white + fast + mvt + diab + htn + hcl + hbpmed + 
                                         cholmed + smk + bmi + act + calor + alco + ahei_noal +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model1[[i]] <- coxph(formu_model1[[i]], data = pool3)
  
  formu_model2[[i]] <- as.formula(paste0("Surv(tdead, dead_cvd) ~ ", modules[i],
                                         " + agemo + white + fast + mvt + diab + htn + hcl + hbpmed + 
                                         cholmed + smk + bmi + act + calor + alco + ahei_noal +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model2[[i]] <- coxph(formu_model2[[i]], data = pool3)
  
  formu_model3[[i]] <- as.formula(paste0("Surv(tdead, dead_can) ~ ", modules[i],
                                         " + agemo + white + fast + mvt + diab + htn + hcl + hbpmed + 
                                         cholmed + smk + bmi + act + calor + alco + ahei_noal +
                                         strata(cacoall, cohort, endpoint)"))
  fit_model3[[i]] <- coxph(formu_model3[[i]], data = pool3)
  
}

# GET RESULTS for death
results_all <- data.frame(matrix(NA, length(modules), 5))
rownames(results_all) <- modules
colnames(results_all)<-c("HR", "LCI", "UCI", "HR (95%CI)", "P_value")
results_cvd=results_all
results_can=results_all

for (i in 1:length(modules)){
  get_data_sum <- summary(fit_model1[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_all[i,1] <- format(round(hr_data[1], 2), nsmall=2)
  results_all[i,2] <- format(round(hr_data[2], 2), nsmall=2)
  results_all[i,3] <- format(round(hr_data[3], 2), nsmall=2)
  results_all[i,5] <- format(round(get_data_sum$coefficients[1,5], 3), nsmall=3)  # Pvalue
  results_all[i,4] <- paste(results_all[i,1], " (", results_all[i,2], ", ", results_all[i,3], ")", sep="")
  
  get_data_sum <- summary(fit_model2[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_cvd[i,1] <- format(round(hr_data[1], 2), nsmall=2)
  results_cvd[i,2] <- format(round(hr_data[2], 2), nsmall=2)
  results_cvd[i,3] <- format(round(hr_data[3], 2), nsmall=2)
  results_cvd[i,5] <- format(round(get_data_sum$coefficients[1,5], 3), nsmall=3)  # Pvalue
  results_cvd[i,4] <- paste(results_cvd[i,1], " (", results_cvd[i,2], ", ", results_cvd[i,3], ")", sep="")
  
  get_data_sum <- summary(fit_model3[[i]])
  hr_data <- get_data_sum$conf.int[1,c("exp(coef)", "lower .95", "upper .95")]
  results_can[i,1] <- format(round(hr_data[1], 2), nsmall=2)
  results_can[i,2] <- format(round(hr_data[2], 2), nsmall=2)
  results_can[i,3] <- format(round(hr_data[3], 2), nsmall=2)
  results_can[i,5] <- format(round(get_data_sum$coefficients[1,5], 3), nsmall=3)  # Pvalue
  results_can[i,4] <- paste(results_can[i,1], " (", results_can[i,2], ", ", results_can[i,3], ")", sep="")
}
