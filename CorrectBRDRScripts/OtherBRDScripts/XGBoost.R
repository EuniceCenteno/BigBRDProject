##xboost another classification method
#Boosting is a sequential ensemble technique in which the model is improved using the information from previously grown weaker models. This process is continued for multiple iterations until a final model is built which will predict a more accurate outcome. There are 3 types of boosting techniques: 1. Adaboost 2. Gradient Descent. 3. Xgboost Xgboost (extreme gradient boosting) is an advanced version of the gradient descent boosting technique, which is used for increasing the speed and efficiency of computation of the algorithm.
#this model works with iterations

library(xgboost) # for fitting the xgboost model
library(caret)   # for general data preparation and model fitting
library(e1071) 

data <- iris               # reads the dataset
head(data)           # head() returns the top 6 rows of the dataframe
summary(data)       # returns the statistical summary of the data columns
dim(data)

#train the data
# createDataPartition() function from the caret package to split the original dataset into a training and testing set and split data into training (80%) and testing set (20%)
parts = createDataPartition(data$Species, p = 0.8, list = F)
train = data[parts, ]
test = data[-parts, ]

#training and testing are matrices
X_train = data.matrix(train[,-5]) # independent variables for train, contain the predictos
y_train = train[,5]  # dependent variables for train, contains the classes

X_test = data.matrix(test[,-5]) # independent variables for test, contain the predictors
y_test = test[,5]  # dependent variables for test,  contains the classes

# convert the train and test data into xgboost matrix type.
xgboost_train = xgb.DMatrix(data=X_train, label=y_train) #merging the features and classes
xgboost_test = xgb.DMatrix(data=X_test, label=y_test)

#we create the model
# train a model using our training data
model <- xgboost(data = xgboost_train,          # the data   
                 max.depth=3, ,                 # max depth 
                 nrounds=50)                   # max number of boosting iterations
summary(model)

#we make predictions
#use model to make predictions on test data
pred_test = predict(model, xgboost_test) #using the testing set
pred_test

#convert predictor to factors
pred_test[(pred_test>3)] = 3 #everything that's bigger than 3, it's converted to 3
pred_y = as.factor((levels(y_test))[round(pred_test)]) #we now transform the values into the classification levels
print(pred_y)

#create matrix
conf_mat = confusionMatrix(y_test, pred_y)
print(conf_mat)

### now test it with the data
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/")
metadata <- read.csv("RFDairyMetadata_morefeatures.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
rownames(metadata) = metadata$Name
metadata = metadata[-1]
str(metadata)
models = metadata[c(1,7,11,15,19)] #bacterial data
str(models)
#models= metadata[c(1,2,32,33,35,41)] #features that were important in RF

#prepare the training and testing
M1 = createDataPartition(models$Status, p = 0.6, list = F)
train = models[M1, ]
test = models[-M1, ]

#training and testing are matrices
X_train1 = data.matrix(train[,-1]) # independent variables for train, contain the predictos
y_train1 = train[,1]  # dependent variables for train, contains the classes

X_test1 = data.matrix(test[,-1]) # independent variables for test, contain the predictors
y_test1 = test[,1]  # dependent variables for test,  contains the classes

# convert the train and test data into xgboost matrix type.
xgboost_train1 = xgb.DMatrix(data=X_train1, label=y_train1) #merging the features and classes
xgboost_test1 = xgb.DMatrix(data=X_test1, label=y_test1)

#create the model
class <- xgboost(data = xgboost_train1, max.depth = 2, eta = 1, nthread = 2, nrounds = 5, verbose = 1)
summary(class)

#we make predictions
#use model to make predictions on test data
pred_test1 = predict(class, xgboost_test1) #using the testing set
pred_test1

#convert predictor to factors
pred_test1[(pred_test1>2)] = 2 #everything that's bigger than 2, it's converted to 2
pred_y1 = as.factor((levels(y_test1))[round(pred_test1)]) #we now transform the values into the classification levels
print(pred_y1)

#create matrix
conf_mat1 = confusionMatrix(y_test1, pred_y1)
print(conf_mat1)


#### multisensi
#helps to identify the sensitivity of the classification
install.packages('multisensi')
library(multisensi)

# NOT RUN {
## Test case : the Winter Wheat Dynamic Models (WWDM)
#  input factors design
data(biomasseX) #includes all the factors
str(biomasseX)
summary(biomasseX)
# input climate variable
data(Climat)
# output variables (precalculated to speed up the example)
data(biomasseY)

# to do dynsi process
# argument reduction=NULL
resD <- multisensi(design=biomasseX, model=biomasseY, reduction=NULL,
                   dimension=NULL, analysis=analysis.anoasg,
                   analysis.args=list(formula=2,keep.outputs = FALSE))
summary(resD)


# to do gsi process
#------------
# with dimension reduction by PCA
# argument reduction=basis.ACP
resG1 <- multisensi(design=biomasseX, model=biomasseY, reduction=basis.ACP,
                    dimension=0.95, analysis=analysis.anoasg,
                    analysis.args=list(formula=2,keep.outputs = FALSE))
summary(resG1)

plot(x=resG1, beside=FALSE)

#------------
# with dimension reduction by o-splines basis
# arguments reduction=basis.osplines
# and basis.args=list(knots= ... , mdegree= ... )

resG2 <- multisensi(design=biomasseX, model=biomasseY, reduction=basis.osplines,
                    dimension=NULL, center=FALSE, scale=FALSE,
                    basis.args=list(knots=11, mdegree=3), analysis=analysis.anoasg, 
                    analysis.args=list(formula=2,keep.outputs = FALSE))
summary(resG2)

#------------
library(sensitivity) # to use fast99

# with dimension reduction by o-splines basis
# and sensitivity analysis with sensitivity:fast99
resG3 <- multisensi(design=fast99, model=biomasse,
                    analysis=analysis.sensitivity, 
                    design.args=list(factors = names(biomasseX), n = 100,
                                     q = "qunif", q.arg = list(list(min = 0.9, max = 2.8),
                                                               list(min = 0.9, max = 0.99), list(min = 0.6, max = 0.8), 
                                                               list(min = 3, max = 12), list(min = 0.0035, max = 0.01),
                                                               list(min = 0.0011, max = 0.0025),
                                                               list(min = 700, max = 1100))), climdata=Climat,
                    reduction=basis.osplines,
                    basis.args=list(knots=7, mdegree=3),
                    center=FALSE,scale=FALSE,dimension=NULL)
summary(resG3)

# }

