load("~/Longevity/202212/Fig2_bcd.RData")
gdata::keep(mortality_all, mortality_cvd, mortality_can, 
            longevity, sure=T)

met_name=read_csv("~/Longevity/met_name.csv")
results_all=merge(mortality_all[,c("HMDB", "mod3_beta")], met_name[,c("HMDB", "CATE2")],all.x=T)
cate=unique(results_all$CATE2)
cate=cate[-1]

library(fgsea)
library(tidyverse)
##########################################################
#               All-cause mortality
##########################################################
cate_list_all=list()
for (i in 1:length(cate)){
  temp=cate[i]
  cate_list_all[[i]]=results_all %>% filter(CATE2 == temp) %>% pull(HMDB)
}
names(cate_list_all)=cate
estimate_all=results_all$mod3_beta
names(estimate_all)=results_all$HMDB
fgsea_all <- fgsea(pathways = cate_list_all, 
                   stats = estimate_all, 
                   minSize=5, maxSize=500) %>% arrange(NES) 

##########################################################
#                     CVD mortality
##########################################################
results_cvd=merge(mortality_cvd[,c("HMDB", "mod3_beta")], met_name[,c("HMDB", "CATE2")],all.x=T)
cate_list_cvd=list()
for (i in 1:length(cate)){
  temp=cate[i]
  cate_list_cvd[[i]]=results_cvd %>% filter(CATE2 == temp) %>% pull(HMDB)
}
names(cate_list_cvd)=cate
estimate_cvd=results_cvd$mod3_beta
names(estimate_cvd)=results_cvd$HMDB
fgsea_cvd <- fgsea(pathways = cate_list_cvd, 
                   stats = estimate_cvd, 
                   minSize=5, maxSize=500) %>% arrange(NES) 

##########################################################
#                   Cancer mortality
##########################################################
results_can=merge(mortality_can[,c("HMDB", "mod3_beta")], met_name[,c("HMDB", "CATE2")],all.x=T)
cate_list_can=list()
for (i in 1:length(cate)){
  temp=cate[i]
  cate_list_can[[i]]=results_can %>% filter(CATE2 == temp) %>% pull(HMDB)
}
names(cate_list_can)=cate
estimate_can=results_can$mod3_beta
names(estimate_can)=results_can$HMDB
fgsea_can <- fgsea(pathways = cate_list_can, 
                   stats = estimate_can, 
                   minSize=5, maxSize=500) %>% arrange(NES) 

##########################################################
#                     Longevity
##########################################################
results_lon=merge(longevity[,c("HMDB", "mod3_beta")], met_name[,c("HMDB", "CATE2")],all.x=T)
cate_list_lon=list()
for (i in 1:length(cate)){
  temp=cate[i]
  cate_list_lon[[i]]=results_lon %>% filter(CATE2 == temp) %>% pull(HMDB)
}
names(cate_list_lon)=cate
estimate_lon=results_lon$mod3_beta
names(estimate_lon)=results_lon$HMDB
fgsea_lon <- fgsea(pathways = cate_list_lon, 
                   stats = estimate_lon, 
                   minSize=5, maxSize=500) %>% arrange(NES) 
