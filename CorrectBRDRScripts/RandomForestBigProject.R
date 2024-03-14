library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
#install.packages("randomForest")
library(randomForest)
library(knitr)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(ggpubr)

##Now try the complete data
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/MetadataCorrectBeef/")
metadata <- read.csv("BRDBigqPCRBeefTotalNamesRow.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
str(metadata)
#metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies) * 100
#metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) * 100
#metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) * 100
#metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) * 100
metadata2= metadata

rownames(metadata2) = metadata2$ID
metadata2 = metadata2[-c(1,2,4:7,11,13,16,19,21)]
metadata2$HSPM = (metadata2$HScopies + metadata2$PMcopies)
metadata2$HSMB = (metadata2$HScopies + metadata2$MBcopies)
metadata2$HSMH = (metadata2$HScopies + metadata2$MHcopies)
metadata2$PMMB = (metadata2$PMcopies + metadata2$MBcopies)
metadata2$PMMH = (metadata2$PMcopies + metadata2$MHcopies)
metadata2$MBMH = (metadata2$MBcopies + metadata2$MHcopies)
metadata2$HSPMMB = (metadata2$HScopies + metadata2$PMcopies + metadata2$MBcopies)
metadata2$HSPMMH = (metadata2$HScopies + metadata2$PMcopies + metadata2$MHcopies)
metadata2$HSMBMH = (metadata2$HScopies + metadata2$MBcopies + metadata2$MHcopies)
metadata2$PMMBMH = (metadata2$PMcopies + metadata2$MBcopies + metadata2$MHcopies)
metadata2$Bact4 = (metadata2$HScopies + metadata2$PMcopies + metadata2$MBcopies + metadata2$MHcopies)

metadata2<- mutate(metadata2, HSPM_log = log10(HSPM + 1))
metadata2<- mutate(metadata2, HSMB_log = log10(HSMB + 1))
metadata2<- mutate(metadata2, HSMH_log = log10(HSMH + 1))
metadata2<- mutate(metadata2, PMMB_log = log10(PMMB + 1))
metadata2<- mutate(metadata2, PMMH_log = log10(PMMH + 1))
metadata2<- mutate(metadata2, MBMH_log = log10(MBMH + 1))
metadata2<- mutate(metadata2, HSPMMB_log = log10(HSPMMB + 1))
metadata2<- mutate(metadata2, HSPMMH_log = log10(HSPMMH + 1))
metadata2<- mutate(metadata2, HSMBMH_log = log10(HSMBMH + 1))
metadata2<- mutate(metadata2, PMMBMH_log = log10(PMMBMH + 1))
metadata2<- mutate(metadata2, Bact4_log = log10(Bact4 + 1))

metadata2 = metadata2[metadata2$X16S.copies != "230", ]
#removing NA's
metadata2 <- na.omit(metadata2)

##Prepare for Random Forest
rm(train)
rm(test)
rm(sample)
set.seed(147562) #test different set seed
sample <- sample.int(n = nrow(metadata2), size = floor(.60*nrow(metadata2)), replace = F)
train <- metadata2[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata2[-sample, ]
test %>% tally()
test %>% count(Status)

### Run random forest
rm(qpcr1)
set.seed(147562)
qpcr1 <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr1)
#rfcv：Random Forest Cross Validation
features <- train[-c(1)]
result = rfcv(features, train$Status, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#Factor importance
varImpPlot(qpcr1, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr1$importance)))
impo = as.data.frame(qpcr1[["importance"]])
write.csv(impo, "impoRFset.seed147562Beef.csv")

#testing set
rm(predQPCR1)
set.seed(147562)
predQPCR1 = predict(qpcr1, newdata=test) #predictions
predQPCR1
predQPCR1 <- as.data.frame(predQPCR1)
predQPCR1$ID <- rownames(predQPCR1)
predQPCR1 %>% count(predQPCR1)
aPCR1 <- table(test$Status, predQPCR1$predQPCR)
Apqpcr1 <- as.data.frame(aPCR1)

#Confusion matrix
matrixqpcr1 = Apqpcr1
matrixqpcr1 #0.37 error rate

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR1$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

### removing HS and PM

metadata5= metadata2[-c(3,4,9:15,17:19,21:26,28:30,32)]

##Prepare for Random Forest
rm(sample)
rm(train)
rm(test)
set.seed(938535) #test different set seed
sample <- sample.int(n = nrow(metadata5), size = floor(.60*nrow(metadata5)), replace = F)
train <- metadata5[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata5[-sample, ]
test %>% tally()
test %>% count(Status)

### Run random forest
rm(qpcr3)
set.seed(938535)
qpcr3 <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr3) # if i remove ID samples, accuracy 28.26%

#rfcv：Random Forest Cross Validation
features <- train[-c(1)]
result = rfcv(features, train$Status, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#Factor importance
varImpPlot(qpcr3, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr3$importance)))
impo2 = as.data.frame(qpcr3[["importance"]])
write.csv(impo2, "impoRFset.seed938535BeefNoHSPM.csv")

#testing set
rm(predQPCR3)
set.seed(938535)
predQPCR3 = predict(qpcr3, newdata=test) #predictions
predQPCR3
predQPCR3 <- as.data.frame(predQPCR3)
predQPCR3$ID <- rownames(predQPCR3)
predQPCR3 %>% count(predQPCR3)
aPCR3 <- table(test$Status, predQPCR3$predQPCR)
Apqpcr3 <- as.data.frame(aPCR3)

#Confusion matrix
matrixqpcr3 = Apqpcr3
matrixqpcr3

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR3$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred)) #0.3571

### removing MB
metadata6= metadata2[-c(5,6,12,14,16,17,19,20,21,23,25,27,28,30:32)]

##Prepare for Random Forest
rm(sample)
rm(train)
rm(test)
set.seed(676463) #test different set seed
sample <- sample.int(n = nrow(metadata6), size = floor(.60*nrow(metadata6)), replace = F)
train <- metadata6[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata6[-sample, ]
test %>% tally()
test %>% count(Status)

### Run random forest
rm(qpcr3)
set.seed(676463)
qpcr3 <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr3) # if i remove ID samples, accuracy 28.26%

#rfcv：Random Forest Cross Validation
features <- train[-c(1)]
result = rfcv(features, train$Status, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#Factor importance
varImpPlot(qpcr3, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr3$importance)))
impo2 = as.data.frame(qpcr3[["importance"]])
write.csv(impo2, "impoRFset.seed676463BeefNoMB.csv")

#testing set
rm(predQPCR3)
set.seed(676463)
predQPCR3 = predict(qpcr3, newdata=test) #predictions
predQPCR3
predQPCR3 <- as.data.frame(predQPCR3)
predQPCR3$ID <- rownames(predQPCR3)
predQPCR3 %>% count(predQPCR3)
aPCR3 <- table(test$Status, predQPCR3$predQPCR)
Apqpcr3 <- as.data.frame(aPCR3)

#Confusion matrix
matrixqpcr3 = Apqpcr3
matrixqpcr3

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR3$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred)) #0.3571

### only 16S rNA
metadata1= metadata2[c(1,2)]

##Prepare for Random Forest
rm(sample)
rm(train)
rm(test)
set.seed(023984) #test different set seed
sample <- sample.int(n = nrow(metadata1), size = floor(.60*nrow(metadata1)), replace = F)
train <- metadata1[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata1[-sample, ]
test %>% tally()
test %>% count(Status)

### Run random forest
rm(qpcr3)
set.seed(023984)
qpcr3 <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr3) # if i remove ID samples, accuracy 28.26%

#rfcv：Random Forest Cross Validation
features <- train[-c(1)]
result = rfcv(features, train$Status, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#Factor importance
varImpPlot(qpcr3, sort=TRUE,scale =FALSE, n.var=min(1,nrow(qpcr3$importance)))
impo2 = as.data.frame(qpcr3[["importance"]])
write.csv(impo2, "impoRFset.seed023984BeefOnly16S.csv")

#testing set
rm(predQPCR3)
set.seed(023984)
predQPCR3 = predict(qpcr3, newdata=test) #predictions
predQPCR3
predQPCR3 <- as.data.frame(predQPCR3)
predQPCR3$ID <- rownames(predQPCR3)
predQPCR3 %>% count(predQPCR3)
aPCR3 <- table(test$Status, predQPCR3$predQPCR)
Apqpcr3 <- as.data.frame(aPCR3)

#Confusion matrix
matrixqpcr3 = Apqpcr3
matrixqpcr3

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR3$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred)) #0.3571

### removing MB and PM and 16S
metadata3= metadata2[c(1,7:10,15,26)]

##Prepare for Random Forest
rm(sample)
rm(train)
rm(test)
set.seed(257375) #test different set seed
sample <- sample.int(n = nrow(metadata3), size = floor(.60*nrow(metadata3)), replace = F)
train <- metadata3[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata3[-sample, ]
test %>% tally()
test %>% count(Status)

### Run random forest
rm(qpcr3)
set.seed(257375)
qpcr3 <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr3) # if i remove ID samples, accuracy 28.26%

#rfcv：Random Forest Cross Validation
features <- train[-c(1)]
result = rfcv(features, train$Status, cv.fold=10)
result$error.cv #corresponding vector of error rates or MSEs at each step
#the mean squared error (MSE) or mean squared deviation (MSD) of an estimator 
#(of a procedure for estimating an unobserved quantity) measures the average of the squares of the errors—that is, the average squared difference between the estimated values and the actual value.
#plot
#10% of the data is used to test the function
# Draw validation, results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

#Factor importance
varImpPlot(qpcr3, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr3$importance)))
impo2 = as.data.frame(qpcr3[["importance"]])
write.csv(impo2, "impoRFset.seed257375BeefNoMBHS16S.csv")

#testing set
rm(predQPCR3)
set.seed(257375)
predQPCR3 = predict(qpcr3, newdata=test) #predictions
predQPCR3
predQPCR3 <- as.data.frame(predQPCR3)
predQPCR3$ID <- rownames(predQPCR3)
predQPCR3 %>% count(predQPCR3)
aPCR3 <- table(test$Status, predQPCR3$predQPCR)
Apqpcr3 <- as.data.frame(aPCR3)

#Confusion matrix
matrixqpcr3 = Apqpcr3
matrixqpcr3

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR3$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred)) #0.3571

### get the mean and confidence intervals for the error rates
errors <- read.csv("RFBeeferros.csv", na.strings = c("","NA"), header=TRUE)

summary2 <-errors %>%
  group_by(Test) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(error),
    SD = sd(error),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(summary2, "RFsummary2Beef.csv")

#generate confidence interval for each error 
training = subset(errors, Test =="Training")
testing = subset(errors, Test =="Testing")
sample.mean <- mean(testing$error)
print(sample.mean)

sample.n <- length(testing$error)
sample.sd <- sd(testing$error)
sample.se <- sample.sd/sqrt(sample.n)
print(sample.se)

alpha = 0.05
degrees.freedom = sample.n - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
print(t.score)

margin.error <- t.score * sample.se
lower.bound <- sample.mean - margin.error
upper.bound <- sample.mean + margin.error

TrainCI = c(lower.bound,upper.bound)
TrainCI
TestCI = c(lower.bound,upper.bound)
TestCI

### Importance RF values
imporFeatures <- read.csv("ImportanceRFBeef.csv", na.strings = c("","NA"), header=TRUE)
str(imporFeatures)

summary2 <-imporFeatures  %>%
  group_by(Feature) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(MeanDecreaseAccuracy),
    SD = sd(MeanDecreaseAccuracy),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(summary2, "RFsummary2BeefImportantFeature.csv")


