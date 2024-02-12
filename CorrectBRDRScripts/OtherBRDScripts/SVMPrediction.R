#Support Vector Machines Algorith
#https://www.datacamp.com/tutorial/support-vector-machines-r
library(e1071)
library(mda)
library(ggplot2)

##now test with the actual data
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/")
metadata <- read.csv("RFDairyMetadata_morefeatures.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("0"="Healthy", "1"="BRD")
str(metadata)
str(metadata)
models = metadata[c(1,2,8,12,16,20)] #bacterial data
#models= metadata[c(1,2,32,33,35,41)] #features that were important in RF
#models= metadata[c(1,2,9,13,17,21)] #bacteria copy number
is.na(models)
models <-models[order(models$Status, decreasing=TRUE),]
M1 = models[c(2,3,4)]
#M1 <- na.omit(M1) 
M2 = models[c(2,3,5)]
#M2 <- na.omit(M2) 
M3 = models[c(2,3,6)]
#M3 <- na.omit(M3)
M4 = models[c(2,4,5)]
#M4 <- na.omit(M4)
M5 = models[c(2,4,6)]
#M5 <- na.omit(M5)
M6 = models[c(2,5,6)]
#M6 <- na.omit(M6)

##create svm
fit=svm(factor(Status)~.,data=M6,scale=FALSE,kernel="linear",cost=5)
fit
summary(fit)
fit$index #the total number of supported vector (points close to the decision boundary or on the wrong side)

#create grid from scratch
str(M6)
x = M6[,-c(1)]
y <- as.numeric(as.character(M$Status))
plot(x, col = y+1) #red is BRD and black is healthy

grange=apply(x,2,range) #this step gives the range of values for each of the predictors
x1=seq(from=grange[1,1],to=grange[2,1],length=203) #creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
x2=seq(from=grange[1,2],to=grange[2,2],length=203)# Same but for x2, creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
xgrid=expand.grid(MHcopy_log=x1,PMcopy_log=x2) #creates a data frame with all the combinations of the supplied vectors
ygrid=predict(fit,xgrid) #we predict the classification based on the model and the grids we generated
plot(fit, M6)

test = cbind(xgrid, ygrid)
grange2=apply(test,2,range) 
#this plot creates the decision boundaries for the data set

#add the contour of the data
func = predict(fit,xgrid,decision.values = TRUE)
func = attributes(func)$decision
plot(xgrid)
plot(xgrid,col=c("black","red")[as.numeric(ygrid)],pch=15,cex=.2)
points(x,col=y+1,pch=19) #x is where the actual data is and y is the color
points(x[fit$index,],pch=5,cex=2)#black is healthy and red is BRD
contour(x1,x2, matrix(func,203,203),level=0, add=TRUE)
contour(x1,x2, matrix(func,203,203),level=0.5, add=TRUE, col='blue') #add the bayes decision boundary
legend("topright", legend = c("BRD", "Healthy"),
       lwd = 1, col = c("red", "black"))

##test
xgridtest = xgrid[c(1:25),]
plot(xgridtest)
ygridtest = ygrid[c(1:25)]
plot(xgridtest,col=c("black","red")[as.numeric(ygridtest)],pch=15,cex=.2)
ygridtest
#you make predictions for every grid


#extract the samples 
samplesIndex = x[fit$index,]#black is healthy and red is BRD, samples that are the support vector 
a = x
b = y

#subset to see what was the values used to make the boundary
class1 = subset(test, ygrid =="1")
class1 = class1[,-c(3)]
class1range=apply(class1,2,range) #to be classified as BRD

class0 = subset(test, ygrid =="0")
class0 = class0[,-c(3)]
class0range=apply(class0,2,range) #to be classified as BRD

#include the Cq to make sure if they are better predictors
cqvalues <- read.csv("DairyBacteriaCqValues.csv", na.strings = c("","NA"), header=TRUE)
cqvalues = merge(cqvalues, metadata, by.x = "ID", by.y = "Name")
cqvalues = cqvalues[,-c(7:66)]

ggplot(cqvalues, aes(x=Status, y=MhCq)) + 
  geom_boxplot()

ggplot(cqvalues, aes(x=Status, y=MbCq)) + 
  geom_boxplot()

ggplot(cqvalues, aes(x=Status, y=HsCq)) + 
  geom_boxplot()

ggplot(cqvalues, aes(x=Status, y=PmCq)) + 
  geom_boxplot()

#test the classifier
M1 = cqvalues[c(2,3,6)]
#M1 <- na.omit(M1) 
M2 = cqvalues[c(2,4,6)]
#M2 <- na.omit(M2) 
M3 = cqvalues[c(2,5,6)]
#M3 <- na.omit(M3)
M4 =cqvalues[c(3,4,6)]
#M4 <- na.omit(M4)
M5 = cqvalues[c(3,5,6)]
#M5 <- na.omit(M5)
M6 = cqvalues[c(4,5,6)]
#M6 <- na.omit(M6)

##create svm
fit=svm(factor(Status)~.,data=M1,scale=FALSE,kernel="linear",cost=5)
fit
summary(fit)
fit$index #the total number of supported vector (points close to the decision boundary or on the wrong side)

#create grid from scratch
str(M1)
x = M1[,-c(1)]
y <- as.numeric(as.character(M1$Status))
plot(x, col = y+1) #red is BRD and black is healthy

grange=apply(x,2,range) #this step gives the range of values for each of the predictors
x1=seq(from=grange[1,1],to=grange[2,1],length=100) #creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
x2=seq(from=grange[1,2],to=grange[2,2],length=100)# Same but for x2, creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
xgrid=expand.grid(Hscopy_log=x1,MBcopy_log=x2) #creates a data frame with all the combinations of the supplied vectors
ygrid=predict(fit,xgrid) #we predict the classification based on the model and the grids we generated
plot(fit, M1)

test = cbind(xgrid, ygrid)
grange2=apply(test,2,range) 
#this plot creates the decision boundaries for the data set

#add the contour of the data
func = predict(fit,xgrid,decision.values = TRUE)
func = attributes(func)$decision
plot(xgrid)
plot(xgrid,col=c("black","red")[as.numeric(ygrid)],pch=15,cex=.2)
points(x,col=y+1,pch=19) #x is where the actual data is and y is the color
points(x[fit$index,],pch=5,cex=2)#black is healthy and red is BRD
contour(x1,x2, matrix(func,100,100),level=0, add=TRUE)
contour(x1,x2, matrix(func,100,100),level=0.5, add=TRUE, col='blue') #add the bayes decision boundary
legend("topright", legend = c("BRD", "Healthy"),
       lwd = 1, col = c("red", "black"))


##3-------- now the beef samples
metadata <- read.csv("RFBeefMetadata_morefeatures.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("0"="Healthy", "1"="BRD")
str(metadata)
#models = metadata[c(1,2,8,12,16,20)]#bacteria
models = metadata[c(1,2,9,14,18,22)]#bacteria
models <-models[order(models$Status, decreasing=TRUE),]
M1 = models[c(2,3,4)]
dplyr::count(M1,Status)
M2 = models[c(2,3,5)]
M3 = models[c(2,3,6)]
M4 = models[c(2,4,5)]
M5 = models[c(2,4,6)]
M6 = models[c(2,5,6)]

##create svm
fit=svm(factor(Status)~.,data=M6,scale=FALSE,kernel="linear",cost=5)
fit
summary(fit)
fit$index #the total number of supported vector (points close to the decision boundary or on the wrong side)

#create grid from scratch
str(M6)
x = M6[,-c(1)]
y <- as.numeric(as.character(M6$Status))
plot(x, col = y+1) #red is BRD and black is healthy

grange=apply(x,2,range) #this step gives the range of values for each of the predictors
x1=seq(from=grange[1,1],to=grange[2,1],length=264) #creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
x2=seq(from=grange[1,2],to=grange[2,2],length=264)# Same but for x2, creates a sequence of the values from the lowest range value on X1 to the highest range and creates 100 data points
xgrid=expand.grid(MHcopies=x1,PMcopies=x2) #creates a data frame with all the combinations of the supplied vectors
ygrid=predict(fit,xgrid) #we predict the classification based on the model and the grids we generated
plot(fit, M6)

test = cbind(xgrid, ygrid)
#this plot creates the decision boundaries for the data set

#add the contour of the data
func = predict(fit,xgrid,decision.values = TRUE)
func = attributes(func)$decision
plot(xgrid,col=c("black","red")[as.numeric(ygrid)],pch=15,cex=.2)
points(x,col=y+1,pch=15) #x is where the actual data is and y is the color
points(x[fit$index,],pch=5,cex=2)#black is healthy and red is BRD
contour(x1,x2, matrix(func,264,264),level=0, add=TRUE)
contour(x1,x2, matrix(func,264,264),level=0.5, add=TRUE, col='blue') #add the bayes decision boundary
legend("topright", legend = c("BRD", "Healthy"),
       lwd = 1, col = c("red", "black"))

#extract the samples 
samplesIndex = x[fit$index,]#black is healthy and red is BRD, samples that are the support vector 
a = x
b = y

#subset to see what was the values used to make the boundary
class1 = subset(test, ygrid =="1")
class1 = class1[,-c(3)]
class1range=apply(class1,2,range) #to be classified as BRD

class0 = subset(test, ygrid =="0")
class0 = class0[,-c(3)]
class0range=apply(class0,2,range) #to be classified as BRD
