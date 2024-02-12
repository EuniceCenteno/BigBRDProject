library(tidyverse)
library(ggpubr)
library(rstatix)
library(afex)
library(emmeans)
library(lme4)
library(lattice)
library("latticeExtra")
library(dplyr)
library(car)
library(scales)
library(reshape2)
library(tidyr)

#to identify for each animal if they tested positive or negative.

setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/BigRF/")
metadata <- read.csv("BRDBigqPCRDairyTotal_RF_CqThre.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
metadata$State <- as.factor(metadata$State)
metadata$Name <- as.factor(metadata$Name)
str(metadata)

#only prevalence data
prev = metadata[,c(1,10,14,18,22)]

## CA samples 
CAdata<- read.csv("GoldStarID.csv", na.strings = c("","NA"), header=TRUE)
CAdata$LabID <- as.factor(CAdata$LabID)
CAdata = merge(prev, CAdata, by.x = "Name", by.y = "LabID")
CAdata$HSpre_abs = as.factor(CAdata$HSpre_abs)
levels(CAdata$HSpre_abs) <- list("Positive"="1", "Negative"="0")
CAdata$MBpre_abs = as.factor(CAdata$MBpre_abs)
levels(CAdata$MBpre_abs) <- list("Positive"="1", "Negative"="0")
CAdata$MHpre_abs = as.factor(CAdata$MHpre_abs)
levels(CAdata$MHpre_abs) <- list("Positive"="1", "Negative"="0")
CAdata$PMpre_abs = as.factor(CAdata$PMpre_abs)
levels(CAdata$PMpre_abs) <- list("Positive"="1", "Negative"="0")
#write.csv(CAdata, "CAdata.csv")

str(CAdata)
M1<- CAdata %>%
  group_by(Status, PMpre_abs) %>%
  summarise(,
            count = n())



