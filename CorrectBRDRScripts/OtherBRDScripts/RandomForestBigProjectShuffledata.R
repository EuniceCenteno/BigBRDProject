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
metadata <- read.csv("BRDBigqPCRDairyTotal_RF_CqThre.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
str(metadata)
metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies) * 100
metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) * 100
metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) * 100
metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) * 100
attach(metadata)
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

## create other fearures: Bacteria bacteria intereaction
str(metadata)
#addition
metadata$InteMbMh <- (metadata$MBcopies + metadata$MHcopies)
metadata$InteMbHS <- (metadata$MBcopies + metadata$HScopies)
metadata$InteMbPM <- (metadata$MBcopies + metadata$PMcopies)
metadata$InteMhPM <- (metadata$MHcopies + metadata$PMcopies)
metadata$InteMhHs <- (metadata$MHcopies + metadata$HScopies)
metadata$IntePmHs <- (metadata$PMcopies + metadata$HScopies)
#multiplication
metadata$InteMbMhMult <- (metadata$MBcopies * metadata$MHcopies)
metadata$InteMbHSMult <- (metadata$MBcopies * metadata$HScopies)
metadata$InteMbPMMult <- (metadata$MBcopies * metadata$PMcopies)
metadata$InteMhPMMult <- (metadata$MHcopies * metadata$PMcopies)
metadata$InteMhHsMult <- (metadata$MHcopies * metadata$HScopies)
metadata$IntePmHsMult <- (metadata$PMcopies * metadata$HScopies)

#addition 3 bacteria
metadata$InteMbMhHs <- (metadata$MBcopies + metadata$MHcopies + metadata$HScopies)
metadata$InteMbHSPm <- (metadata$MBcopies + metadata$HScopies + metadata$PMcopies)
metadata$InteMbPMMh <- (metadata$MBcopies + metadata$PMcopies + metadata$MHcopies)
metadata$InteMhPMHs <- (metadata$MHcopies + metadata$PMcopies + metadata$HScopies)

#multiplication 3 bacteria
metadata$InteMbMhHsMult <- (metadata$MBcopies * metadata$MHcopies * metadata$HScopies)
metadata$InteMbHSPmMult <- (metadata$MBcopies * metadata$HScopies * metadata$PMcopies)
metadata$InteMbPMMhMult <- (metadata$MBcopies * metadata$PMcopies * metadata$MHcopies)
metadata$InteMhPMHsMult <- (metadata$MHcopies * metadata$PMcopies * metadata$HScopies)

#addition 4 bacteria
metadata$bacteria4 <- (metadata$MBcopies + metadata$MHcopies + metadata$HScopies + metadata$PMcopies)
metadata$bacteria4Mult <- (metadata$MHcopies * metadata$PMcopies * metadata$HScopies * metadata$MBcopies)
metadata$TotalBact <- (metadata$MHcopies + metadata$PMcopies + metadata$HScopies + metadata$MBcopies)

metadata$HSRelBactTotal <- (metadata$HScopies / metadata$TotalBact) * 100
metadata$PMRelBactTotal <- (metadata$PMcopies / metadata$TotalBact) * 100
metadata$MBRelBactTotal <- (metadata$MBcopies / metadata$TotalBact) * 100
metadata$MHRelBactTotal <- (metadata$MHcopies / metadata$TotalBact) * 100

#write.csv(metadata, "RFDairyMetadata_morefeatures.csv")

#import new data
metadata2 <- read.csv("RFBeefMetadata_morefeatures.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
str(metadata)
metadata = metadata[,-c(5,6)] #for dairy farms

#attach(metadata)
#metadata <- metadata[order(State),]

str(metadata)
metadata$State <- as.factor(metadata$State)
levels(metadata$State)
metadata <- na.omit(metadata)
metadata2 = metadata
metadata2 = metadata[c(1:3,7,12,16,20)] #beef data fram
#metadata2 = metadata[c(1:3,7,11,15,19)] #dairy data fram
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

#all the information
metadata <- merge(metadata, metadata2[-c(2:7)], by.x = "Name", by.y = "Name")

str(metadata)
rm(train)
rm(test)
set.seed(789)
sample <- sample.int(n = nrow(metadata), size = floor(.60*nrow(metadata)), replace = F)
train <- metadata[sample, ]
str(train)
train %>% tally()
train %>% count(Status)

test  <- metadata[-sample, ]
test %>% tally()
test %>% count(Status)
str(test)

rownames(train) <- train$Name
train <- train[-c(1)]
summary(train)
train$Status <- as.factor(train$Status)
train$ID <- rownames(train)

meanHS <-train[,c(1:2, 60:63)] %>%
  na.omit(HScopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(HScopiesNA), list(HSmean = mean))
meanHS
meanHS = meanHS[-c(1),]
#train <- merge(train, meanHS, by.x = "Status", by.y = "Status")
rownames(train) = train$ID
train = cbind(train,meanHS[,-c(1)])
train$HSDistMean <- train$HSmean - train$HScopies
str(train)
train$HSDistMeanNew <- train$HSDistMean
str(train)
train$HSDistMeanNew <- as.factor(train$HSDistMeanNew)
t1 = train[c(1,2,6,67)]
levels(train$HSDistMeanNew)[levels(train$HSDistMeanNew)=='7095244.66603147'] <- '190'
train$HSDistMeanNew <- as.numeric(as.character(train$HSDistMeanNew))
str(train)

meanPM <-train[,c(1:2, 60:63)] %>%
  na.omit(PMcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(PMcopiesNA), list(PMmean = mean))
meanPM
meanPM = meanPM[-c(1),]
train = cbind(train, meanPM[,-c(1)])
#train <- merge(train, meanPM, by.x = "Status", by.y = "Status")
train$PMDistMean <- train$PMmean - train$PMcopies
str(train)
train$PMDistMeanNew <- train$PMDistMean
str(train)
train$PMDistMeanNew <- as.factor(train$PMDistMeanNew)
t2 = train[c(1,2,19,70)]
levels(train$PMDistMeanNew)[levels(train$PMDistMeanNew)=='796799.831205083'] <- '190'
train$PMDistMeanNew <- as.numeric(as.character(train$PMDistMeanNew))

meanMB <-train[,c(1:2, 60:63)] %>%
  na.omit(MBcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(MBcopiesNA), list(MBmean = mean))
meanMB
meanMB = meanMB[-c(1),]
train = cbind(train, meanMB[,-c(1)])
#train <- merge(train, meanMB, by.x = "Status", by.y = "Status")
train$MBDistMean <- train$MBmean - train$MBcopies
str(train)
train$MBDistMeanNew <- train$MBDistMean
str(train)
train$MBDistMeanNew <- as.factor(train$MBDistMeanNew)
t3 = train[c(1,2,11,73)]
levels(train$MBDistMeanNew)[levels(train$MBDistMeanNew)=='2558435.39200886'] <- '670'
train$MBDistMeanNew <- as.numeric(as.character(train$MBDistMeanNew))

meanMH <-train[,c(1:2, 60:63)] %>%
  na.omit(MHcopiesNA) %>%
  group_by(Status) %>%
  summarise_at(vars(MHcopiesNA), list(MHmean = mean))
meanMH
meanMH = meanMH[-c(1),]
train = cbind(train, meanMH[,-c(1)])
#train <- merge(train, meanMH, by.x = "Status", by.y = "Status")
train$MHDistMean <- train$MHmean - train$MHcopies
str(train)
train$MHDistMeanNew <- train$MHDistMean
str(train)
train$MHDistMeanNew <- as.factor(train$MHDistMeanNew)
t4 = train[c(1,2,15,76)]
levels(train$MHDistMeanNew)[levels(train$MHDistMeanNew)=='123900.915885056'] <- '170'
train$MHDistMeanNew <- as.numeric(as.character(train$MHDistMeanNew))

str(train)
#trainB = train
rownames(train) <- train$ID
train2 = train[,-c(3,60:66,68:69,71:72,74:75)] #beef samples
#train2 = train[,-c(3,60:66,68,69,71,72,74,75)]
train2 %>% count(Status)
train3 <- train2[,-c(2)] #remove state and farm features
train4 <- train2[,-c(2,20:27)] #only features related with bacteria abundance

test
rownames(test) <- test$Name
test <- test[-c(1)]
summary(test)
test$Status <- as.factor(test$Status)
test$ID <- rownames(test)

test  <- cbind(test,meanHS[,-c(1)])
#test <- merge(test, meanHS, by.x = "Status", by.y = "Status")
test$HSDistMean <- test$HSmean - test$HScopies
test$HSDistMeanNew <- test$HSDistMean
str(test)
test$HSDistMeanNew <- as.factor(test$HSDistMeanNew)
levels(test$HSDistMeanNew)[levels(test$HSDistMeanNew)=='7095244.66603147'] <- '190'
test$HSDistMeanNew <- as.numeric(as.character(test$HSDistMeanNew))

#test <- merge(test, meanPM, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanPM[,-c(1)])
test$PMDistMean <- test$PMmean - test$PMcopies
test$PMDistMeanNew <- test$PMDistMean
str(test)
test$PMDistMeanNew <- as.factor(test$PMDistMeanNew)
levels(test$PMDistMeanNew)[levels(test$PMDistMeanNew)=='796799.831205083'] <- '190'
test$PMDistMeanNew <- as.numeric(as.character(test$PMDistMeanNew))

#test <- merge(test, meanMB, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanMB[,-c(1)])
test$MBDistMean <- test$MBmean - test$MBcopies
test$MBDistMeanNew <- test$MBDistMean
str(test)
test$MBDistMeanNew <- as.factor(test$MBDistMeanNew)
levels(test$MBDistMeanNew)[levels(test$MBDistMeanNew)=='2558435.39200886'] <- '670'
test$MBDistMeanNew <- as.numeric(as.character(test$MBDistMeanNew))

#test <- merge(test, meanMH, by.x = "Status", by.y = "Status")
test  <- cbind(test,meanMH[,-c(1)])
test$MHDistMean <- test$MHmean - test$MHcopies
test$MHDistMeanNew <- test$MHDistMean
str(test)
test$MHDistMeanNew <- as.factor(test$MHDistMeanNew)
levels(test$MHDistMeanNew)[levels(test$MHDistMeanNew)=='123900.915885056'] <- '170'
test$MHDistMeanNew <- as.numeric(as.character(test$MHDistMeanNew))
str(test)

rownames(test) <- test$ID
test2 = test[,-c(3,60:66,68:69,71:72,74:75)] #for beef
#test2 = test[,-c(3,60:66,68,69,71,72,74,75)] #for dairy
str(test2)
test2 %>% count(Status)
test3 <- test2[,-c(2)] #without state
test4 <- test2[,-c(2,20:27)] #only features related with bacteria abundance

rm(qpcr)
set.seed(123)
qpcr <- randomForest(Status~ .,data=train2, importance=TRUE)
qpcr <- randomForest(Status~ .,data=train3, importance=TRUE)
qpcr <- randomForest(Status~ .,data=train4, importance=TRUE)
# For numeric variables, NAs are replaced with column medians. 
print(qpcr)
plotqpcr

#rfcv：Random Forest Cross Validation
features <- train2[-c(1)]
features <- train3[-c(1)]
features <- train4[-c(1)]

result = rfcv(features, train2$Status, cv.fold=10)
result = rfcv(features, train3$Status, cv.fold=10)
result = rfcv(features, train4$Status, cv.fold=10)
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
predQPCR = predict(qpcr, newdata=test2) #predictions
predQPCR = predict(qpcr, newdata=test3) #predictions wihtout state
predQPCR = predict(qpcr, newdata=test4) #predictions only including bacteria features
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(test2$Status, predQPCR$predQPCR)
aPCR <- table(test3$Status, predQPCR$predQPCR)#no state feature
aPCR <- table(test4$Status, predQPCR$predQPCR)#bacteria features
Apqpcr <- as.data.frame(aPCR)

#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
test2%>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

test3%>% #no state
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

test4%>% #bacteria features
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/") 
##Try RF with 16S data
metadata <- read.csv("BRDBigMetadataDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 
metadata2= metadata[c(1,6)]
metadata2= merge(metadata2, metadata[,c(3:4),], by.x=0, by.y=0)
row.names(metadata2) <- metadata2$ID
metadata2= metadata2[,-c(1:2)]

ASVs <- read_qza("table-no-pseudo.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data)
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

#change names on the ASV table 
ASV_table <- as.data.frame(ASV_table)
ASV_table <- merge(metadata, ASV_table, by.x = "cattle", by.y = 0)
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:7)]
ASV_table <- as.matrix(ASV_table)

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomyNoPseudo.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)

#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub("d__*", "", tax.final[,1])
tax.final$Phylum <- sub("p__*", "", tax.final[,1])
tax.final$Phylum <- sub("c__*", "", tax.final[,1])
tax.final$Class <- sub("d__*", "", tax.final[,2])
tax.final$Class <- sub("p__*", "", tax.final[,2])
tax.final$Class <- sub("c__*", "", tax.final[,2])
tax.final$Class <- sub("o__*", "", tax.final[,2])
tax.final$Order <- sub("d__*", "", tax.final[,3])
tax.final$Order <- sub("p__*", "", tax.final[,3])
tax.final$Order <- sub("c__*", "", tax.final[,3])
tax.final$Order <- sub("o__*", "", tax.final[,3])
tax.final$Order <- sub("f__*", "", tax.final[,3])
tax.final$Family <- sub("d__*", "", tax.final[,4])
tax.final$Family <- sub("p__*", "", tax.final[,4])
tax.final$Family <- sub("c__*", "", tax.final[,4])
tax.final$Family <- sub("o__*", "", tax.final[,4])
tax.final$Family <- sub("f__*", "", tax.final[,4])
tax.final$Family <- sub("g__*", "", tax.final[,4])
tax.final$Family <- sub("D_9__*", "", tax.final[,4])
tax.final$Genus <- sub("d__*", "", tax.final[,5])
tax.final$Genus <- sub("p__*", "", tax.final[,5])
tax.final$Genus <- sub("c__*", "", tax.final[,5])
tax.final$Genus <- sub("o__*", "", tax.final[,5])
tax.final$Genus <- sub("f__*", "", tax.final[,5])
tax.final$Genus <- sub("g__*", "", tax.final[,5])
tax.final$Genus <- sub("s__*", "", tax.final[,5])
tax.final$Species <- sub("d__*", "", tax.final[,6])
tax.final$Species <- sub("p__*", "", tax.final[,6])
tax.final$Species <- sub("c__*", "", tax.final[,6])
tax.final$Species <- sub("o__*", "", tax.final[,6])
tax.final$Species <- sub("f__*", "", tax.final[,6])
tax.final$Species <- sub("g__*", "", tax.final[,6])
tax.final$Species <- sub("s__*", "", tax.final[,6])

TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.csv(TaxASV,"taxonomy.csv", row.names = TRUE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 12051 taxa and 402 samples ] dairy
#[ 26496 taxa and 498 samples ] beef

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != " Mitochondria" & 
                Genus != " Mitochondria" &
                Species != " Mitochondria" &
                Order != " Chloroplast" &
                Family != " Chloroplast" &
                Genus != " Chloroplast" &
                Species != " Chloroplast") -> physeq_deseq
physeq_deseq #[ 11992 taxa and 402 samples ] Dairy
#[ 26365 taxa and 498 samples ] beef

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

#only include the ASV with >50 sequences
otu_abund <- colSums(NewASVTable) ##the abundance of each ASV across all samples
otu_abund = subset(otu_abund, otu_abund >=50)
otu_abund = as.data.frame(otu_abund) #from 11992 to 3466

NewASV = t(NewASVTable)
NewASV <- t(NewASV[rownames(NewASV) %in% rownames(otu_abund),]) 


#all the information
metadata2 <- merge(metadata, NewASVTable, by.x = 0, by.y = 0)
row.names(metadata2) <- metadata2[,1]
metadata2 = metadata2[,-c(1)]

#only the ASV >=50
metadata2 <- merge(metadata, NewASV, by.x = 0, by.y = 0)
row.names(metadata2) <- metadata2[,1]
metadata2 = metadata2[,-c(1)]

str(metadata)
rm(train)
rm(test)
set.seed(789)
sample <- sample.int(n = nrow(metadata2), size = floor(.60*nrow(metadata2)), replace = F)
train <- metadata2[sample, ]
str(train)
train %>% tally()
train %>% count(Status)
train2=train

test  <- metadata2[-sample, ]
test %>% tally()
test %>% count(Status)
str(test)
test2=test

summary(train)
train$Status <- as.factor(train$Status)

rm(qpcr)
set.seed(123)
qpcr <- randomForest(Status~ .,data=train2, importance=TRUE)
print(qpcr)
plotqpcr

#rfcv：Random Forest Cross Validation
features <- train2[-c(1)]

result = rfcv(features, train2$Status, cv.fold=10)
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
predQPCR = predict(qpcr, newdata=test2) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(test2$Status, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)

#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
test2%>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

NewTax2 = NewTax
NewTax2$ASV = rownames(NewTax2)

##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/") 
###beef samples
metadata <- read.csv("BRDBigMetadataBeef.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 
metadata = metadata[,-c(1,2,6,7)] #for dairy farms

ASVs <- read_qza("table-no-pseudo.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data)
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

#change names on the ASV table 
ASV_table <- as.data.frame(ASV_table)
ASV_table <- merge(metadata, ASV_table, by.x = "ID", by.y = 0)
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:7)]
ASV_table <- as.matrix(ASV_table)

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomy2.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)

#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub("d__*", "", tax.final[,1])
tax.final$Phylum <- sub("p__*", "", tax.final[,1])
tax.final$Phylum <- sub("c__*", "", tax.final[,1])
tax.final$Class <- sub("d__*", "", tax.final[,2])
tax.final$Class <- sub("p__*", "", tax.final[,2])
tax.final$Class <- sub("c__*", "", tax.final[,2])
tax.final$Class <- sub("o__*", "", tax.final[,2])
tax.final$Order <- sub("d__*", "", tax.final[,3])
tax.final$Order <- sub("p__*", "", tax.final[,3])
tax.final$Order <- sub("c__*", "", tax.final[,3])
tax.final$Order <- sub("o__*", "", tax.final[,3])
tax.final$Order <- sub("f__*", "", tax.final[,3])
tax.final$Family <- sub("d__*", "", tax.final[,4])
tax.final$Family <- sub("p__*", "", tax.final[,4])
tax.final$Family <- sub("c__*", "", tax.final[,4])
tax.final$Family <- sub("o__*", "", tax.final[,4])
tax.final$Family <- sub("f__*", "", tax.final[,4])
tax.final$Family <- sub("g__*", "", tax.final[,4])
tax.final$Family <- sub("D_9__*", "", tax.final[,4])
tax.final$Genus <- sub("d__*", "", tax.final[,5])
tax.final$Genus <- sub("p__*", "", tax.final[,5])
tax.final$Genus <- sub("c__*", "", tax.final[,5])
tax.final$Genus <- sub("o__*", "", tax.final[,5])
tax.final$Genus <- sub("f__*", "", tax.final[,5])
tax.final$Genus <- sub("g__*", "", tax.final[,5])
tax.final$Genus <- sub("s__*", "", tax.final[,5])
tax.final$Species <- sub("d__*", "", tax.final[,6])
tax.final$Species <- sub("p__*", "", tax.final[,6])
tax.final$Species <- sub("c__*", "", tax.final[,6])
tax.final$Species <- sub("o__*", "", tax.final[,6])
tax.final$Species <- sub("f__*", "", tax.final[,6])
tax.final$Species <- sub("g__*", "", tax.final[,6])
tax.final$Species <- sub("s__*", "", tax.final[,6])

TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.csv(TaxASV,"taxonomy.csv", row.names = TRUE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 26496 taxa and 498 samples ] beef

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != " Mitochondria" & 
                Genus != " Mitochondria" &
                Species != " Mitochondria" &
                Order != " Chloroplast" &
                Family != " Chloroplast" &
                Genus != " Chloroplast" &
                Species != " Chloroplast") -> physeq_deseq
physeq_deseq #[ 26365 taxa and 498 samples ]

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

#only include the ASV with >50 sequences
otu_abund <- colSums(NewASVTable) ##the abundance of each ASV across all samples
otu_abund = subset(otu_abund, otu_abund >=50)
otu_abund = as.data.frame(otu_abund) #from 11992 to 3466

NewASV = t(NewASVTable)
NewASV <- t(NewASV[rownames(NewASV) %in% rownames(otu_abund),]) 


#all the information
metadata2 <- merge(metadata, NewASVTable, by.x = 0, by.y = 0)
row.names(metadata2) <- metadata2[,1]
metadata2 = metadata2[,-c(1)]

#only the ASV >=50
metadata2 <- merge(metadata, NewASV, by.x = 0, by.y = 0)
row.names(metadata2) <- metadata2[,1]
metadata2 = metadata2[,-c(1)]

str(metadata)
rm(train)
rm(test)
set.seed(789)
sample <- sample.int(n = nrow(metadata2), size = floor(.60*nrow(metadata2)), replace = F)
train <- metadata2[sample, ]
str(train)
train %>% tally()
train %>% count(Status)
train2=train

test  <- metadata2[-sample, ]
test %>% tally()
test %>% count(Status)
str(test)
test2=test

summary(train)
train$Status <- as.factor(train$Status)

rm(qpcr)
set.seed(123)
qpcr <- randomForest(Status~ .,data=train2, importance=TRUE)
print(qpcr)
plotqpcr

#rfcv：Random Forest Cross Validation
features <- train2[-c(1)]

result = rfcv(features, train2$Status, cv.fold=10)
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
predQPCR = predict(qpcr, newdata=test2) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(test2$Status, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)

#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
test2%>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(Status != lda.pred))

NewTax2 = NewTax
NewTax2$ASV = rownames(NewTax2)

#Function
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}

 
