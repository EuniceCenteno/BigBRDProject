## Logistic regresssion
#Logistic regression in R Programming is a classification algorithm used to find the probability of event success and event failure. Logistic regression is used when the dependent variable is binary(0/1, True/False, Yes/No) in nature. Logit function is used as a link function in a binomial distribution.

#Test following 
#https://www.geeksforgeeks.org/logistic-regression-in-r-programming/
# Loading package
library(dplyr)
# Loading package
library(caTools)
library(ROCR) 

#AUC ranges in value from 0 to 1. A model whose predictions are 100% wrong has an AUC of 0.0; one whose predictions are 100% correct has an AUC of 1.0.
#AUC represents the probability that a random positive (green) example is positioned to the right of a random negative (red) example.
# test how well the model is

#first you need to build the ROC curve
#An ROC curve (receiver operating characteristic curve) is a graph showing the performance of a classification model at all classification thresholds.

### Notw we need to test the dairy and beef samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/") 

metadata <- read.csv("BRDBigqPCRDairyTotal_RF_CqThre.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)
metadata$Status <- as.factor(metadata$Status)
#now let's select for the variables we will test
metadata <- metadata[-5]
rownames(metadata) <- metadata$Name


#removinf NA's
metadata <- na.omit(metadata) #to remove misisng values 406 samples to 313 (347 without alpha diversity and date of collection) 
#when I don't include the date of collection | beef 528 to 439
str(metadata)
metadata$State <- as.factor(metadata$State)
levels(metadata$State)
levels(metadata$Status) <- list("0"="Healthy", "1"="BRD") #we make it binomial

select <- c("CA", "NY", "IN")
select2 <- c("TX")
train <- metadata[metadata$State %in% select,]

meanHS <-train %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    HSmean = mean(HScopies),
  )
train <- merge(train, meanHS[-2], by.x = "Status", by.y = "Status")
train$HSDistMean <- train$HSmean - train$HScopies

meanPM <-train %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    PMmean = mean(PMcopies),
  )
train <- merge(train, meanPM[-2], by.x = "Status", by.y = "Status")
train$PMDistMean <- train$PMmean - train$PMcopies

meanMB <-train %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    MBmean = mean(MBcopies),
  )
train <- merge(train, meanMB[-2], by.x = "Status", by.y = "Status")
train$MBDistMean <- train$MBmean - train$MBcopies

meanMH <-train %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    MHmean = mean(MHcopies),
  )
train <- merge(train, meanMH[-2], by.x = "Status", by.y = "Status")
train$MHDistMean <- train$MHmean - train$MHcopies
rownames(train) <- train$Name
#Dairy
train = train[,-c(2,31,33,35,37)] # to remove distance mean columns
train2 = train
train = train2[c(1,30:33)] #only dist
train = train2[c(1,6,10,14,18)] #bacteria 
train = train2[c(1,7,11,15,19)] #bacteria copy
#Beef 
train = train[,-c(2,32,34,36,38)] # to remove mean columns
train2 = train
train = train2[c(1,31:34)] #only dist
train = train2[c(1,7,11,15,19)] #bacteria 
train = train2[c(1,8,13,17,21)] #bacteria copy


test <- metadata[metadata$State %in% select2,]
test <- merge(test, meanHS[-2], by.x = "Status", by.y = "Status")
test$HSDistMean <- test$HSmean - test$HScopies

test <- merge(test, meanPM[-2], by.x = "Status", by.y = "Status")
test$PMDistMean <- test$PMmean - test$PMcopies

test <- merge(test, meanMB[-2], by.x = "Status", by.y = "Status")
test$MBDistMean <- test$MBmean - test$MBcopies

test <- merge(test, meanMH[-2], by.x = "Status", by.y = "Status")
test$MHDistMean <- test$MHmean - test$MHcopies

rownames(test) <- test[,2]
#dairy
test= test[,-c(2,31,33,35,37)]
test2 = test
test = test2[c(1,30:33)] #only dist
test = test2[c(1,6,10,14,18)] #bact abund
test = test2[c(1,7,11,15,19)] #bact copt
#beef
test = test[,-c(2,32,34,36,38)] # to remove mean columns
test2 = test
test = test2[c(1,31:34)] #only dist
test = test2[c(1,7,11,15,19)] #bacteria 
test = test2[c(1,8,13,17,21)] #bacteria copy

# Training model
#dependent variables is binary 
#Logistic regression is used when the dependent variable is binary(0/1, True/False, Yes/No) in nature. 
str(train)
model <- glm(Status ~ PMDistMean:HSDistMean + PMDistMean:MBDistMean +PMDistMean:MHDistMean +
               HSDistMean:MBDistMean +HSDistMean:MHDistMean +MBDistMean:MHDistMean, 
             data = train, 
             family = "binomial")
model
summary(model)

model <- glm(Status ~PMcopy_log:Hscopy_log + PMcopy_log:MBcopy_log + PMcopy_log:MHcopy_log +
             Hscopy_log:MBcopy_log + Hscopy_log:MHcopy_log + MBcopy_log:MHcopy_log,
                      data = train, 
                      family = "binomial")
model
summary(model)

str(train)
model <- glm(Status ~ HScopies + MBcopies + MHcopies + PMcopies,
             data = train, 
             family = "binomial")
model
summary(model)

model <- glm(Status ~ PMcopies:HScopies + PMcopies:MBcopies + PMcopies:MHcopies +
               HScopies:MBcopies + HScopies:MHcopies + MBcopies:MHcopies,
             data = train, 
             family = "binomial")
model
summary(model)

# Predict test data based on model
predict = predict(model, newdata=test, type = "response") #predictions
predict  

# Changing probabilities
predQPCR <- ifelse(predict >0.5, 1, 0)
predQPCR
predQPCR2 <- as.data.frame(predQPCR)
predQPCR2$ID <- rownames(predQPCR2)
predQPCR2 %>% count(predQPCR)

#Misclasification rate
test %>% count(Status)
test%>%
  mutate(lda.pred = (predQPCR2$predQPCR)) %>%
  summarise(lda.error = mean(Status != lda.pred))

# Evaluating model accuracy
# using confusion matrix
missing_classerr <- mean(predQPCR!= test$Status)
print(paste('Accuracy =', 1 - missing_classerr)) # this is the accuracy of the model

# ROC-AUC Curve to check the sensitivity and specificity
ROCPred <- prediction(predQPCR, test$Status) # using the predectec variable
ROCPer <- performance(ROCPred, measure = "tpr", 
                      x.measure = "fpr")

auc <- performance(ROCPred, measure = "auc")
auc <- auc@y.values[[1]]
auc

# Plotting curve
plot(ROCPer)
plot(ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC CURVE")
abline(a = 0, b = 1)

auc <- round(auc, 4)
legend(.6, .4, auc, title = "AUC", cex = 1)

##test regression
##Now try the complete data
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/")
metadata <- read.csv("BRDBigqPCRDairyTotal_RF_CqThre.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
str(metadata)
#metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies) * 100
#metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) * 100
#metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) * 100
#metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) * 100
#attach(metadata)
#metadata <- metadata[order(State),]

# Converting date of collection to numeric values
#str(metadata)
#metadata$date <- metadata$Date.of.colleciton
#metadata$Date.of.clleciton <- as.Date(metadata$Date.of.colleciton, "%m/%d/%y")
#d<- as.Date('02/10/22', "%m/%d/%y") #use to calculate the days
#metadata$Date.of.colleciton <- as.Date(d) -as.Date(metadata$Date.of.colleciton) 
#metadata$Date.of.colleciton <- as.numeric(metadata$Date.of.colleciton)
#str(metadata$Date.of.colleciton)
metadata = metadata[,-c(5)]

#removinf NA's
#metadata3 <- na.omit(metadata) #to remove misisng values 406 samples to 313 (347 without alpha diversity and date of collection) 
#when I don't include the date of collection | beef 528 to 439
str(metadata)
metadata$State <- as.factor(metadata$State)
levels(metadata$State)
metadata2 = metadata
#metadata2 = metadata[c(1:3,9,14,18,22)] #beef data fram
metadata2 = metadata[c(1:3,8,12,16,20)] #dairy data fram
summary(metadata2)
metadata2$HScopiesNA <- metadata2$HScopies
metadata2$HScopiesNA <- as.numeric(metadata2$HScopiesNA)
summary(metadata2)
metadata2$HScopiesNA[metadata2$HScopiesNA == '190'] <- ''
metadata2$HScopiesNA <- as.numeric(metadata2$HScopiesNA)

summary(metadata2)
metadata2$MBcopiesNA <- metadata2$MBcopies
metadata2$MBcopiesNA <- as.numeric(metadata2$MBcopiesNA)
summary(metadata2)
metadata2$MBcopiesNA[metadata2$MBcopiesNA == '670'] <- ''
metadata2$MBcopiesNA <- as.numeric(metadata2$MBcopiesNA)

summary(metadata2)
metadata2$MHcopiesNA <- metadata2$MHcopies
metadata2$MHcopiesNA <- as.numeric(metadata2$MHcopiesNA)
summary(metadata2)
metadata2$MHcopiesNA[metadata2$MHcopiesNA == '170'] <- ''
metadata2$MHcopiesNA <- as.numeric(metadata2$MHcopiesNA)

summary(metadata2)
metadata2$PMcopiesNA <- metadata2$PMcopies
metadata2$PMcopiesNA <- as.numeric(metadata2$PMcopiesNA)
summary(metadata2)
metadata2$PMcopiesNA[metadata2$PMcopiesNA == '190'] <- ''
metadata2$PMcopiesNA <- as.numeric(metadata2$PMcopiesNA)
levels(metadata2$State)

select <- c("TX", "CA", "IN")
select2 <- c("NY")
train <- metadata2[metadata2$State %in% select,]
rownames(train) <- train$Name
train <- train[-c(1,3)]
summary(train)
train$Status <- as.factor(train$Status)
train$ID <- rownames(train)

meanHS <-train[,-c(2:5,7:10)] %>%
  na.omit(HScopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(HScopiesNA), list(HSmean = mean))
meanHS
meanHS = meanHS[-c(1),]
#train <- merge(train, meanHS, by.x = "Status", by.y = "Status")
train = cbind(train, meanHS[,-c(1)])
train$HSDistMean <- train$HSmean - train$HScopies

meanPM <-train[,-c(2:8, 10:12)] %>%
  na.omit(PMcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(PMcopiesNA), list(PMmean = mean))
meanPM
meanPM = meanPM[-c(1),]
train = cbind(train, meanPM[,-c(1)])
#train <- merge(train, meanPM, by.x = "Status", by.y = "Status")
train$PMDistMean <- train$PMmean - train$PMcopies

meanMB <-train[,-c(2:6, 8:14)] %>%
  na.omit(MBcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(MBcopiesNA), list(MBmean = mean))
meanMB
meanMB = meanMB[-c(1),]
train = cbind(train, meanMB[,-c(1)])
#train <- merge(train, meanMB, by.x = "Status", by.y = "Status")
train$MBDistMean <- train$MBmean - train$MBcopies

meanMH <-train[,-c(2:7, 9:16)] %>%
  na.omit(MHcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(MHcopiesNA), list(MHmean = mean))
meanMH
meanMH = meanMH[-c(1),]
train = cbind(train, meanMH[,-c(1)])
#train <- merge(train, meanMH, by.x = "Status", by.y = "Status")
train$MHDistMean <- train$MHmean - train$MHcopies
trainB = train
rownames(train) <- train$ID
train = train[,-c(2:11,13,15,17)]
train2 <- merge(metadata,train,  by.x = "Name", by.y = 0)
rownames(train2) = train2$Name
train2 <- train2[,c(2,32:35)]
train2 <- na.omit(train2) # from 433

#train = train[,-c(2,31,33,35,37)]
#train = train2[c(1,34:37)] #only dist
#train = train2[c(1,6,10,14,18)] #bacteria 

test <- metadata2[metadata2$State %in% select2,]
test  <- cbind(test,meanHS[,-c(1)])
#test <- merge(test, meanHS, by.x = "Status", by.y = "Status")
test$HSDistMean <- test$HSmean - test$HScopies

#test <- merge(test, meanPM, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanPM[,-c(1)])
test$PMDistMean <- test$PMmean - test$PMcopies

#test <- merge(test, meanMB, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanMB[,-c(1)])
test$MBDistMean <- test$MBmean - test$MBcopies

#test <- merge(test, meanMH, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanMH[,-c(1)])
test$MHDistMean <- test$MHmean - test$MHcopies
rownames(test) <- test$Name
test2 = test
test = test[,-c(2:12,14,16,18)]
test = test[-1]
test2 = merge(metadata, test, by.x = "Name", by.y = 0)
rownames(test2) <- test2$Name
test2 = test2[,c(2, 31:34)]
test2 <- na.omit(test2)
#test= test[,-c(1,35,37,39,41)]
#test= test[,-c(2,32,34,36,38)]
#test = test2[c(1,34:37)] #only dist
#test = test2[c(1,7,11,15,19)] #bact abund

str(train)
model <- glm(Status.x ~ PMDistMean + HSDistMean + MBDistMean +MHDistMean,
             data = train2, 
             family = "binomial")
model
summary(model)

model <- glm(Status ~ PMDistMean:HSDistMean + PMDistMean:MBDistMean +PMDistMean:MHDistMean +
               HSDistMean:MBDistMean +HSDistMean:MHDistMean +MBDistMean:MHDistMean, 
             data = train, 
             family = "binomial")
model
summary(model)

model <- glm(Status ~PMcopy_log:Hscopy_log + PMcopy_log:MBcopy_log + PMcopy_log:MHcopy_log +
               Hscopy_log:MBcopy_log + Hscopy_log:MHcopy_log + MBcopy_log:MHcopy_log,
             data = train, 
             family = "binomial")
model
summary(model)

str(train)
model <- glm(Status ~ HScopies + MBcopies + MHcopies + PMcopies,
             data = train, 
             family = "binomial")
model
summary(model)

model <- glm(Status ~ PMcopies:HScopies + PMcopies:MBcopies + PMcopies:MHcopies +
               HScopies:MBcopies + HScopies:MHcopies + MBcopies:MHcopies,
             data = train, 
             family = "binomial")
model
summary(model)

# Predict test data based on model
predict = predict(model, newdata=test2, type = "response") #predictions
predict  

# Changing probabilities
predQPCR <- ifelse(predict >0.5, 1, 0)
predQPCR
predQPCR2 <- as.data.frame(predQPCR)
predQPCR2$ID <- rownames(predQPCR2)
predQPCR2 %>% count(predQPCR)

#Misclasification rate
test2 %>% count(Status)
test2%>%
  mutate(lda.pred = (predQPCR2$predQPCR)) %>%
  summarise(lda.error = mean(Status != lda.pred))

# Evaluating model accuracy
# using confusion matrix
missing_classerr <- mean(predQPCR!= test$Status)
print(paste('Accuracy =', 1 - missing_classerr)) # this is the accuracy of the model

# ROC-AUC Curve to check the sensitivity and specificity
ROCPred <- prediction(predQPCR, test$Status) # using the predectec variable
ROCPer <- performance(ROCPred, measure = "tpr", 
                      x.measure = "fpr")

auc <- performance(ROCPred, measure = "auc")
auc <- auc@y.values[[1]]
auc

# Plotting curve
plot(ROCPer)
plot(ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC CURVE")
abline(a = 0, b = 1)

auc <- round(auc, 4)
legend(.6, .4, auc, title = "AUC", cex = 1)

