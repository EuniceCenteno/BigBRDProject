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
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/")
metadata <- read.csv("BRDBigqPCRDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
str(metadata)
metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies) * 100
metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) * 100
metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) * 100
metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) * 100

# Converting date of collection to numeric values
#str(metadata)
#metadata$date <- metadata$Date.of.colleciton
#metadata$Date.of.clleciton <- as.Date(metadata$Date.of.colleciton, "%m/%d/%y")
#d<- as.Date('02/10/22', "%m/%d/%y") #use to calculate the days
#metadata$Date.of.colleciton <- as.Date(d) -as.Date(metadata$Date.of.colleciton) 
#metadata$Date.of.colleciton <- as.numeric(metadata$Date.of.colleciton)
#str(metadata$Date.of.colleciton)
metadata = metadata[,-c(4:6)] # remove farm, date of collection and rectal temperature

#removing NA's
metadata2 <- na.omit(metadata)

##Prepare for Random Forest
set.seed(789) #test different set seed
sample <- sample.int(n = nrow(metadata2), size = floor(.60*nrow(metadata2)), replace = F)
train <- metadata2[sample, ]
str(train)
train %>% tally()
train %>% count(Status)
rownames(train) = train$Name
train = train[-1]

test  <- metadata2[-sample, ]
test %>% tally()
test %>% count(Status)
rownames(test) = test$Name
test = test[-1]

### Run random forest
set.seed(123)
qpcr <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr)

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
varImpPlot(qpcr, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr$importance)))

#testing set
rm(predQPCR)
set.seed(123)
predQPCR = predict(qpcr, newdata=test) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(test$Status, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)

#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

### Now separate the state and run RF for each state
IN = subset(metadata2, State=="IN")
CA = subset(metadata2, State=="CA") #change to CO for the beef samples, change to CA for dairy
NY = subset(metadata2, State=="NY")#change to ID for the beef samples, change to NY for dairy
TX = subset(metadata2, State=="TX")

##Prepare for Random Forest
meta2 = TX
set.seed(789) #test different set seed
sample <- sample.int(n = nrow(meta2), size = floor(.60*nrow(meta2)), replace = F)
train <- meta2[sample, ]
str(train)
train %>% tally()
train %>% count(Status)
rownames(train) = train$Name
train = train[-1]

test  <- meta2[-sample, ]
test %>% tally()
test %>% count(Status)
rownames(test) = test$Name
test = test[-1]

### Run random forest
set.seed(123)
qpcr <- randomForest(Status~ .,data=train, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr)

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
varImpPlot(qpcr, sort=TRUE,scale =FALSE, n.var=min(40,nrow(qpcr$importance)))

#testing set
rm(predQPCR)
set.seed(123)
predQPCR = predict(qpcr, newdata=test) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(test$Status, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)

#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
test%>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

