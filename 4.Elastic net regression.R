load("~/Longevity/202212/mortality_met_2.RData")
all=pool2
gdata::keep(all, sure=T)

### training set and test set (70/30)   
set.seed(1234)
n=floor(nrow(all)*0.7)
train_ind=sample(seq_len(nrow(all)), size = n)
train=all[train_ind,]
test=all[-train_ind,]

library(glmnet)
library(survival)

record_dead = data.frame(repNo=NA,NoMetabs=NA,C_index=NA)
met_coef_dead = data.frame(HMDB=colnames(train[,c(1:243)]))

for (inrep in 1:100) { 
  x=as.matrix(train[,c(1:243)])
  y=as.matrix(data.frame(time=train$tdead, status=train$dead))
  
  # Harell C index, a higher C index means better prediction performance
  Training_CV=cv.glmnet(x, y, family = "cox", alpha=0.5, type.measure = "C", nfolds = 10)
  lambda_1se_10F = Training_CV$lambda.1se
  
  met_coef_dead[,dim(met_coef_dead)[2]+1]=data.frame(coef(Training_CV, s=lambda_1se_10F)[,])
  names(met_coef_dead)[dim(met_coef_dead)[2]] = paste("inrep",inrep,sep='')
  
  lambda_c=data.frame(lambda=Training_CV$lambda, c=Training_CV$cvm)
  record_dead[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record_dead[inrep,"NoMetabs"] = sum(met_coef_dead[,dim(met_coef_dead)[2]]!=0)
  record_dead[inrep,"C_index"] = lambda_c[lambda_c$lambda==lambda_1se_10F,]$c
  
}

### find the score for dead
table(record_dead$NoMetabs)
# 58 62 67 75 79 80 84 87 
# 2  5 17 42 22  9  2  1

# calculate metabolite score, inrep1
identical(colnames(test)[1:243], met_coef_dead$HMDB) # TRUE
test0=test
test0$metscore_dead=colSums(t(test[, 1:243])*met_coef_dead$inrep1)

### leave-one-out in the training set
for (i in 1:n) { 
  temp=train[-i, ]
  x=as.matrix(temp[,c(1:243)])
  y=as.matrix(data.frame(time=temp$tdead, status=temp$dead))
  
  # Harell C index, a higher C index means better prediction performance
  Training_CV=cv.glmnet(x, y, family = "cox", alpha=0.5, type.measure = "C", nfolds = 10)
  lambda_1se_10F = Training_CV$lambda.1se
  
  met_coef_dead[,dim(met_coef_dead)[2]+1]=data.frame(coef(Training_CV, s=lambda_1se_10F)[,])
  names(met_coef_dead)[dim(met_coef_dead)[2]] = paste("coef_",i,sep='')
  
  record_dead[i,"NoMetabs"] = sum(met_coef_dead[,dim(met_coef_dead)[2]]!=0)
}
