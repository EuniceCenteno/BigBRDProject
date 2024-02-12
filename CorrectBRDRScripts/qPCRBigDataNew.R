#qPCR data
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
library(dplyr)
library(sjstats)
library(ggplot2)

### Dairy samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/MetadataCorrectDairy/")
metadata <- read.csv("BRDBigqPCRDairyTotalOneRow.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$State <- as.factor(metadata$State)
str(metadata)
metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies)
metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) 
metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) 
metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) 

metadata$BactABund = (metadata$HScopies + metadata$PMcopies + metadata$MBcopies + metadata$MHcopies)
metadata$HSRelBact <- (metadata$HScopies / metadata$BactABund)
metadata$PMRelBact <- (metadata$PMcopies / metadata$BactABund) 
metadata$MBRelBact <- (metadata$MBcopies / metadata$BactABund) 
metadata$MHRelBact <- (metadata$MHcopies / metadata$BactABund) 
metadata <- mutate(metadata, X16S_log = log10(X16S.copies + 1))
write.csv(metadata, "CompleteDairyqPCRDatawithRelativeAbundance.csv")

meta1 = metadata[metadata$X16S.copies != "230", ]
write.csv(meta1, "qPCRDataNo16Scopies230.csv")

#generate model
str(metadata)
M1 <- mixed(Hscopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$Hscopy_log)
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

M2 <- mixed(MBcopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)
shapiro.test(metadata$MBcopy_log)
plot(M2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

M3 <- mixed(MHcopy_log ~Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
shapiro.test(metadata$MHcopy_log)
plot(M3$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

M4 <- mixed(PMcopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)
shapiro.test(metadata$PMcopy_log)
plot(M4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

## statistical analysis for 16S rRNA abundance 
#generate model
str(metadata)
M1 <- mixed(X16S.copies ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$X16S.copies)
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model)) # not normaly distributed

#transform the 16S rRNA copies to log values
M2 <- mixed(X16S_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M2
anova(M2)
shapiro.test(metadata$X16S_log)
plot(M2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model)) # not normaly distributed

## Identifying overall effect between BRD and healthy and the abundance of BRD-pathobionts
str(metadata)

meta2=meta1
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

model_16T <- map(.x = params_list, 
                .f = ~ mwu(meta2, X16S.copies,Status))

M16T = model_16T[["1"]][["tab.df"]]
M16T$Bacteria = "16S rRNA"

model_16L <- map(.x = params_list, 
                 .f = ~ mwu(meta2, X16S_log,Status))

M16L = model_16L[["1"]][["tab.df"]]
M16L$Bacteria = "16S rRNA log10"

model_MBT <- map(.x = params_list, 
                .f = ~ mwu(meta2, MBcopy_log,Status))

MMBT = model_MBT[["1"]][["tab.df"]]
MMBT$Bacteria = "M.bovis"

model_MHT <- map(.x = params_list, 
                 .f = ~ mwu(meta2, MHcopy_log,Status))

MMHT = model_MHT[["1"]][["tab.df"]]
MMHT$Bacteria = "M.haemolytica"

model_PMT <- map(.x = params_list, 
                 .f = ~ mwu(meta2, PMcopy_log,Status))

MPMT = model_PMT[["1"]][["tab.df"]]
MPMT$Bacteria = "P.multocida"

model_HST <- map(.x = params_list, 
                 .f = ~ mwu(meta2, Hscopy_log,Status))

MHST = model_HST[["1"]][["tab.df"]]
MHST$Bacteria = "H.somni"

MannWuT_pvals <- do.call(rbind, list(M16T,M16L,MMBT, MMHT, MPMT, MHST))
MannWuT_pvals$p.adj<-p.adjust(MannWuT_pvals$p, method = "BH")
write.csv(MannWuT_pvals, "OverallBacteriaabundanceonlyStatusDairy.csv")
# no significant effect on the overall BRD pathobiont abundance between healthy and BRD

### identify the difference of 16S rRNA abundance among the different states
my_comparisonsD <- list( c("CA", "IN"), c("CA", "NY"), c("CA", "TX"),
                        c("IN", "NY"), c("IN", "TX"), c("NY", "TX")) ## comparison for dairy samples
summary2 <-meta1 %>%
  group_by(State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(X16S_log),
    SD = sd(X16S_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary2$Bacteria = c("16S log")
SumTotalS = summary2
SumTotalS = rbind(SumTotalS, summary2)

X16S = subset(SumTotalS, Bacteria == "16S copies")
X16L = subset(SumTotalS, Bacteria == "16S log")

my_color2 <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)# for dairy samples

ggplot(metadata, aes(x=State, y=X16S.copies)) +
  geom_jitter(aes(x=State, y=X16S.copies,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16S, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= X16S, aes(State, y=Mean, color=State)) +
  ylab("16S rRNA copy number") +xlab ("State") +
  scale_color_manual(values = my_color2) +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisonsD, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))

v = ggplot(metadata, aes(x=State, y=X16S_log)) +
  geom_jitter(aes(x=State, y=X16S_log,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16L, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= X16L, aes(State, y=Mean, color=State)) +
  ylab("16S rRNA log10") +xlab ("State") +
  scale_color_manual(values = my_color2) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisonsD, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))
v


##prepare the data for the loop, to identify the effect within state
IN= subset(meta1, State=="IN")
NY= subset(meta1, State=="NY") # change this to NY for dairy, ID for beef
CA= subset(meta1, State=="CA")# change this to CA for dairy, CO for beef
TX= subset(meta1, State=="TX")

meta2=NY
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_HS <- map(.x = params_list, 
                .f = ~ mwu(meta2, Hscopy_log,Status))
MHS = model_HS[["1"]][["tab.df"]]
MHS$Bacteria = "H.somni"
model_MB <- map(.x = params_list, 
                .f = ~ mwu(meta2, MBcopy_log,Status))
MMB = model_MB[["1"]][["tab.df"]]
MMB$Bacteria = "M.bovis"
model_MH <- map(.x = params_list, 
                .f = ~ mwu(meta2, MHcopy_log,Status))
MMH = model_MH[["1"]][["tab.df"]]
MMH$Bacteria = "M.heameolytica"
model_PM <- map(.x = params_list, 
                .f = ~ mwu(meta2, PMcopy_log,Status))
MPM = model_PM[["1"]][["tab.df"]]
MPM$Bacteria = "P.multocida"

model_16 <- map(.x = params_list, 
                .f = ~ mwu(meta2, X16S.copies,Status))
M16 = model_16[["1"]][["tab.df"]]
M16$Bacteria = "16S rRNA"

model_16Log <- map(.x = params_list, 
                .f = ~ mwu(meta2, X16S_log,Status))
M16log = model_16Log [["1"]][["tab.df"]]
M16log$Bacteria = "16S rRNA log"

## bacteria relative abundance 16
model_HSRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, HSRel,Status))
MHS_Rel = model_HSRel[["1"]][["tab.df"]]
MHS_Rel$Bacteria = "H.somni Rel"

model_MBRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, MBRel,Status))
MMB_Rel = model_MBRel[["1"]][["tab.df"]]
MMB_Rel$Bacteria = "M.bovis Rel"

model_MHRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, MHRel,Status))
MMH_Rel = model_MHRel[["1"]][["tab.df"]]
MMH_Rel$Bacteria = "M.haemolytica Rel"

model_PMRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, PMRel,Status))
MPM_Rel = model_PMRel[["1"]][["tab.df"]]
MPM_Rel$Bacteria = "P.multocida Rel"

## abundance based on BRD-pathobiont abundance
model_HSBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, HSRelBact,Status))
MHS_RelBact = model_HSBact[["1"]][["tab.df"]]
MHS_RelBact$Bacteria = "H.somni RelTotal"

model_MBBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, MBRelBact,Status))
MMB_RelBact = model_MBBact[["1"]][["tab.df"]]
MMB_RelBact$Bacteria = "M.bovis RelTotal"

model_MHBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, MHRelBact,Status))
MMH_RelBact = model_MHBact[["1"]][["tab.df"]]
MMH_RelBact$Bacteria = "M.haemolytica RelTotal"

model_PMBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, PMRelBact,Status))
MPM_RelBact = model_PMBact[["1"]][["tab.df"]]
MPM_RelBact$Bacteria = "P.multocida RelTotal"

MannWu_pvals <- do.call(rbind, list(M16, M16log, MHS, MMB, MMH, MPM, MHS_Rel, MMB_Rel, MMH_Rel,
                                    MPM_Rel, MHS_RelBact, MMB_RelBact, MMH_RelBact, MPM_RelBact))
MannWu_pvals$State <- c("NY")
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResultsTotal = MannWu_pvals
ResultsTotal =rbind(ResultsTotal,MannWu_pvals)
write.csv(ResultsTotal, "ResultsTotalDairy.csv") #only CA- M. haemolytica abundance was significant

#plots
metadata2 <- read.csv("BRDBigqPCRDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

### calculate what is teh average copy log per pathobiont
str()
sum <-metadata2 %>%
  group_by(Bacteria, State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(sum, "BRDpathobiontAbundance.csv")

#get the mean of the relative abundance and create plots for the significant results
str(metadata)
sumRel <-meta1 %>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(PMRel),
    SD = sd(PMRel),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumRel$Bacteria = c("P. multocida 16S Rel") #NY, PM, HS, MH
#sumRelT = sumRel
sumRelT = rbind(sumRelT, sumRel)
write.csv(sumRelT, "summaryBactRelDairy.csv")

## plot relative abundance significant results
sumRelT = read.csv("summaryBactRelDairy.csv", na.strings = c("","NA"), header=TRUE)
sumRelNY = sumRelT[c(1:2),]
sumRelNY$Status = as.factor(sumRelNY$Status)
levels(sumRelNY$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

sumRelNY2 = sumRelT[c(3:4),]
sumRelNY2$Status = as.factor(sumRelNY2$Status)
levels(sumRelNY2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

sumRelCA = sumRelT[c(5:6),]
sumRelCA$Status = as.factor(sumRelCA$Status)
levels(sumRelCA$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

my_colors <- c("dodgerblue3","goldenrod3")

ggplot(sumRelNY, aes(x=Bacteria, y=Mean,color=Status)) +
  #geom_point(aes(x=Bacteria, y=Mean, color=Status),position = position_jitterdodge(0.5), alpha=0.2, size=1) +
  geom_errorbar(data=sumRelNY, aes(x=Bacteria, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_grid(Bacteria~.) + 
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelNY, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean Relative Abundance") + ggtitle("NY") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

ggplot(sumRelNY2, aes(x=Bacteria, y=Mean,color=Status)) +
  #geom_point(aes(x=Bacteria, y=Mean, color=Status),position = position_jitterdodge(0.5), alpha=0.2, size=1) +
  geom_errorbar(data=sumRelNY2, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_grid(Bacteria~.) + 
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelNY2, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean Relative Abundance") + ggtitle("NY") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))


ggplot(sumRelCA, aes(x=Bacteria, y=Mean,color=Status)) +
  geom_errorbar(data=sumRelCA, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelCA, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean Relative Abundance") + ggtitle("CA") + facet_grid(Bacteria~.) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

ggplot(sumRelIN, aes(x=Bacteria, y=Mean,color=Status)) +
  geom_errorbar(data=sumRelIN, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelIN, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean Relative Abundance") + ggtitle("IN") + facet_grid(Bacteria~.) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

ggarrange(a,ncol =1,labels = c("A"),
          font.label = list(size = 18))

ggarrange(b,c,ncol =2,labels = c("B", "C"),
          font.label = list(size = 18))

## to plot the significance of M. haemolytica in CA samples
meta2 = metadata2[metadata2$X16S.copies != "230", ]

summary <-meta2%>%
  group_by(State, Bacteria, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$comparison = paste(summary$State, summary$Bacteria, summary$Status, sep = "-")
my_colors <- c("dodgerblue3","goldenrod3")

CAMh= subset(summary, State == 'CA' & Bacteria == 'M.haemolytica')
CAMh2= subset(metadata2, State == 'CA' & Bacteria == 'M.haemolytica')

ggplot(CAMh2, aes(x=Bacteria, y=copy_log,color=Status)) +
  geom_point(aes(x=Bacteria, y=copy_log, color=Status),position = position_jitterdodge(0.4), alpha=0.2, size=1) +
  geom_errorbar(data=CAMh, aes(x=Bacteria, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_grid(State~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= CAMh, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Copy number (log10)") +xlab ("M. haemolytica") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12, face = "italic"), axis.text.y =
          element_text(color = "black", size = 12))

##### ---------- Bacteria prevalence------------------
metadata2
#assign numerical values to factors
str(metadata2)
metadata2$Bacteria <- as.factor(metadata2$Bacteria)
levels(metadata2$Bacteria) <- list("H.somni"="H. somni", "P.multocida"="P.multocida", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica")
metadata2$pre_abs = as.factor(metadata2$pre_abs)
levels(metadata2$pre_abs) <- list("Negative"="0", "Positive"="1")

Bact2<- metadata2 %>% 
  group_by(Bacteria, State, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Bacteria, State, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup
str(Bact2)
write.csv(Bact2, "BacteriaPrevalenceDairy.csv")
Bact2$Status = as.factor(Bact2$Status)
levels(Bact2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
Bact2$Bacteria = as.factor(Bact2$Bacteria)
levels(Bact2$Bacteria) <- list("H.somni"="H.somni", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica", "P.multocida"="P.multocida")
Bact2$Samples = paste(Bact2$State,Bact2$Bacteria, Bact2$Status, Bact2$pre_abs , sep = "-")
Bact2P = subset(Bact2, pre_abs== "Positive")
write.csv(Bact2P, "BacteriaPrevalencePositiveDairy.csv")

ggplot(data=Bact2P, aes(x=Bacteria, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive %") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size= 11, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

## Check if there is a difference in the bacteria abundance in the Warren and Southview samples
#because some samples were collected from cows and calves from the same farm

warren = subset(metadata2, Farm=="Warren")
Southview = subset(metadata2, Farm=="Southview")
Animal = rbind(warren, Southview)

sumAni <-Animal %>%
  group_by(Farm,animal, Status,Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(sumAni, "BacteriaAbunSummaryCowCalfWarranSouthview.csv")

#3 Run stats
##prepare the data for the loop, to identify the effect within state
meta3=IN
meta3$animal = as.factor(meta3$animal)
combins <- combn(levels(meta3$animal), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_AnimalHS <- map(.x = params_list, 
                .f = ~ mwu(meta3, Hscopy_log,animal))
ManiHS = model_AnimalHS[["1"]][["tab.df"]]
ManiHS$Bacteria = "H.somni"

model_AnimalPM <- map(.x = params_list, 
                      .f = ~ mwu(meta3, PMcopy_log,animal))
ManiPM = model_AnimalPM[["1"]][["tab.df"]]
ManiPM$Bacteria = "P.multocida"

model_AnimalMB <- map(.x = params_list, 
                      .f = ~ mwu(meta3, MBcopy_log,animal))
ManiMB = model_AnimalMB[["1"]][["tab.df"]]
ManiMB$Bacteria = "M.bovis"

model_AnimalMH <- map(.x = params_list, 
                      .f = ~ mwu(meta3, MHcopy_log,animal))
ManiMH = model_AnimalMH[["1"]][["tab.df"]]
ManiMH$Bacteria = "M.heamolytica"
MAnimal <- do.call(rbind, list(ManiHS, ManiPM, ManiMB, ManiMH))
MAnimal$State <- c("IN")
MAnimal$p.adj<-p.adjust(MAnimal$p, method = "BH")
#ResultsAnimal = MAnimal
ResultsAnimal =rbind(ResultsAnimal,MAnimal)
write.csv(ResultsAnimal, "ResultsAnimalNYIN.csv")


###########################################################################################
### -------------------------Beef samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata")
metadata <- read.csv("./BigRF/BRDBigqPCRBeefTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$State <- as.factor(metadata$State)
str(metadata)
metadata$HSRel <- (metadata$HScopies / metadata$X16S.copies)
metadata$PMRel <- (metadata$PMcopies / metadata$X16S.copies) 
metadata$MBRel <- (metadata$MBcopies / metadata$X16S.copies) 
metadata$MHRel <- (metadata$MHcopies / metadata$X16S.copies) 

metadata$BactABund = (metadata$HScopies + metadata$PMcopies + metadata$MBcopies + metadata$MHcopies)
metadata$HSRelBact <- (metadata$HScopies / metadata$BactABund)
metadata$PMRelBact <- (metadata$PMcopies / metadata$BactABund) 
metadata$MBRelBact <- (metadata$MBcopies / metadata$BactABund) 
metadata$MHRelBact <- (metadata$MHcopies / metadata$BactABund) 
metadata <- mutate(metadata, X16S_log = log10(X16S.copies + 1))
metadata= metadata[,-c(5:6)]
write.csv(metadata, "BeefqPCRRelativeAbundance.csv")

#generate model
str(metadata)
M1 <- mixed(Hscopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$Hscopy_log)
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

M2 <- mixed(MBcopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)
shapiro.test(metadata$MBcopy_log)
plot(M2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

M3 <- mixed(MHcopy_log ~Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
shapiro.test(metadata$MHcopy_log)
plot(M3$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

M4 <- mixed(PMcopy_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)
shapiro.test(metadata$PMcopy_log)
plot(M4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

## statistical analysis for 16S rRNA abundance 
#generate model
str(metadata)
M1 <- mixed(X16S.copies ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$X16S.copies)
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model)) # not normaly distributed

#transform the 16S rRNA copies to log values
M2 <- mixed(X16S_log ~ Status + Farm +  (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M2
anova(M2)
shapiro.test(metadata$X16S_log)
plot(M2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model)) # not normaly distributed

## Identifying overall effect between BRD and healthy and the abundance of BRD-pathobionts
str(metadata)

meta2=metadata
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

model_16T <- map(.x = params_list, 
                 .f = ~ mwu(meta2, X16S.copies,Status))

M16T = model_16T[["1"]][["tab.df"]]
M16T$Bacteria = "16S rRNA"

model_16L <- map(.x = params_list, 
                 .f = ~ mwu(meta2, X16S_log,Status))

M16L = model_16L[["1"]][["tab.df"]]
M16L$Bacteria = "16S rRNA log10"

model_MBT <- map(.x = params_list, 
                 .f = ~ mwu(meta2, MBcopy_log,Status))

MMBT = model_MBT[["1"]][["tab.df"]]
MMBT$Bacteria = "M.bovis"

model_MHT <- map(.x = params_list, 
                 .f = ~ mwu(meta2, MHcopy_log,Status))

MMHT = model_MHT[["1"]][["tab.df"]]
MMHT$Bacteria = "M.haemolytica"

model_PMT <- map(.x = params_list, 
                 .f = ~ mwu(meta2, PMcopy_log,Status))

MPMT = model_PMT[["1"]][["tab.df"]]
MPMT$Bacteria = "P.multocida"

model_HST <- map(.x = params_list, 
                 .f = ~ mwu(meta2, Hscopy_log,Status))

MHST = model_HST[["1"]][["tab.df"]]
MHST$Bacteria = "H.somni"

MannWuT_pvals <- do.call(rbind, list(M16T,M16L,MMBT, MMHT, MPMT, MHST))
MannWuT_pvals$p.adj<-p.adjust(MannWuT_pvals$p, method = "BH")
write.csv(MannWuT_pvals, "OverallBacteriaabundanceonlyStatusBeef.csv")
#overall effect: 4 bacteria significant between BRD and Healthy

metadata2 <- read.csv("BRDBigqPCRBeefTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata2 = metadata2[,-c(5,6)]

## PLot the overall BRD-pathobionts abundance between BRD and healthy animals
# PLOT ONLY FOR Beef samples
sumRel <-metadata %>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(MHRelBact),
    SD = sd(MHRelBact),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumRel$Bacteria = c("MHRelBact")
#sumRelT = sumRel
sumRelT = rbind(sumRelT, sumRel)
write.csv(sumRelT, "summaryBactRelBeef.csv")

##PLot significant relative abudance results
sumRelT = read.csv("summaryBactRelBeef.csv", na.strings = c("","NA"), header=TRUE)
sumRelT <- sumRelT[sumRelT$State %in% c("CO", "IN", "TX"), ]
sumRelT <- sumRelT[sumRelT$Bacteria %in% c("HSRel", "MBRel", "MHRelBact", "HSRelBact", "PMRelBact"), ]

sumRelCO = sumRelT[c(1,2,7,8),]
sumRelCO$Status = as.factor(sumRelCO$Status)
levels(sumRelCO$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

sumRelIN = sumRelT[c(15,16,21,22),]
sumRelIN$Status = as.factor(sumRelIN$Status)
levels(sumRelIN$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

sumRelTX = sumRelT[c(17,18,23,24,29,30),]
sumRelTX$Status = as.factor(sumRelTX$Status)
levels(sumRelTX$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(sumRelIN, aes(x=Bacteria, y=Mean,color=Status)) +
  #geom_point(aes(x=Bacteria, y=Mean, color=Status),position = position_jitterdodge(0.5), alpha=0.2, size=1) +
  geom_errorbar(data=sumRelIN, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Bacteria~.) + 
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelIN, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean RelAbundance~ BRD pathobionts") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 12))

ggplot(sumRelCO, aes(x=Bacteria, y=Mean,color=Status)) +
  #geom_point(aes(x=Bacteria, y=Mean, color=Status),position = position_jitterdodge(0.5), alpha=0.2, size=1) +
  geom_errorbar(data=sumRelCO, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Bacteria~.) + 
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelCO, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean Rel Abundance") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 12))

ggplot(sumRelTX, aes(x=Bacteria, y=Mean,color=Status)) +
  #geom_point(aes(x=Bacteria, y=Mean, color=Status),position = position_jitterdodge(0.5), alpha=0.2, size=1) +
  geom_errorbar(data=sumRelTX, aes(x=Bacteria, ymin=Mean - SE,
                                   ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Bacteria~.) + 
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Mean RelAbundance~ BRD pathobionts") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 12))

## summary for the big qPCR plots
summary <-metadata2%>%
  group_by(Bacteria, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$comparison = paste(summary$Status, summary$Bacteria, sep = "-")
my_colors <- c("dodgerblue3","goldenrod3")

##Plot included the overall abundance of BRD-pathobionts between BRD and heathy
ggplot(metadata2, aes(x=Bacteria, y=copy_log,color=Status)) +
  geom_point(aes(x=Bacteria, y=copy_log, color=Status),position = position_jitterdodge(0.4), alpha=0.2, size=1) +
  geom_errorbar(data=summary, aes(x=Bacteria, ymin=Mean - SE,
                              ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data=summary, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Copy number (log10)")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12, face = "italic"), axis.text.y =
          element_text(color = "black", size = 12))

### identify the difference in BRD abundance between states
my_comparisons <- list( c("CO", "ID"), c("CO", "IN"), c("CO", "TX"),
                        c("ID", "IN"), c("ID", "TX"), c("IN", "TX")) # comparison for beef samples
summary2 <-metadata %>%
  group_by(State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(X16S.copies),
    SD = sd(X16S.copies),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary2$Bacteria = c("16S copies")
#SumTotalS = summary2
SumTotalS = rbind(SumTotalS, summary2)

X16S = subset(SumTotalS, Bacteria == "16S copies")
X16L = subset(SumTotalS, Bacteria == "16S log")

my_color1 <- c(
  "#653936", "#CD9BCD","orange","#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
) #for beef samples 

a = ggplot(metadata, aes(x=State, y=X16S.copies)) +
  geom_jitter(aes(x=State, y=X16S.copies,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16S, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= X16S, aes(State, y=Mean, color=State)) +
  ylab("16S rRNA copy number") +xlab ("State") +
  scale_color_manual(values = my_color2) +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))

b = ggplot(metadata, aes(x=State, y=X16S_log)) +
  geom_jitter(aes(x=State, y=X16S_log,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16L, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= X16L, aes(State, y=Mean, color=State)) +
  ylab("16S rRNA log10") +xlab ("State") +
  scale_color_manual(values = my_color2) +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))

ggarrange(a,b,ncol =2, labels = c("A", "B"),
          font.label = list(size = 11))

##prepare the data for the loop, to identify the effect within state
IN= subset(metadata, State=="IN")
ID= subset(metadata, State=="ID") 
CO= subset(metadata, State=="CO")
TX= subset(metadata, State=="TX")

meta2=CO
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_HS <- map(.x = params_list, 
                .f = ~ mwu(meta2, Hscopy_log,Status))
MHS = model_HS[["1"]][["tab.df"]]
MHS$Bacteria = "H.somni"
model_MB <- map(.x = params_list, 
                .f = ~ mwu(meta2, MBcopy_log,Status))
MMB = model_MB[["1"]][["tab.df"]]
MMB$Bacteria = "M.bovis"
model_MH <- map(.x = params_list, 
                .f = ~ mwu(meta2, MHcopy_log,Status))
MMH = model_MH[["1"]][["tab.df"]]
MMH$Bacteria = "M.heamolytica"
model_PM <- map(.x = params_list, 
                .f = ~ mwu(meta2, PMcopy_log,Status))
MPM = model_PM[["1"]][["tab.df"]]
MPM$Bacteria = "P.multocida"

model_16 <- map(.x = params_list, 
                .f = ~ mwu(meta2, X16S.copies,Status))
M16 = model_16[["1"]][["tab.df"]]
M16$Bacteria = "16S rRNA"

model_16Log <- map(.x = params_list, 
                   .f = ~ mwu(meta2, X16S_log,Status))
M16log = model_16Log [["1"]][["tab.df"]]
M16log$Bacteria = "16S rRNA log"

## bacteria relative abundance 16
model_HSRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, HSRel,Status))
MHS_Rel = model_HSRel[["1"]][["tab.df"]]
MHS_Rel$Bacteria = "H.somni Rel"

model_MBRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, MBRel,Status))
MMB_Rel = model_MBRel[["1"]][["tab.df"]]
MMB_Rel$Bacteria = "M.bovis Rel"

model_MHRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, MHRel,Status))
MMH_Rel = model_MHRel[["1"]][["tab.df"]]
MMH_Rel$Bacteria = "M.haemolytica Rel"

model_PMRel <- map(.x = params_list, 
                   .f = ~ mwu(meta2, PMRel,Status))
MPM_Rel = model_PMRel[["1"]][["tab.df"]]
MPM_Rel$Bacteria = "P.multocida Rel"

## abundance based on BRD-pathobiont abundance
model_HSBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, HSRelBact,Status))
MHS_RelBact = model_HSBact[["1"]][["tab.df"]]
MHS_RelBact$Bacteria = "H.somni RelTotal"

model_MBBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, MBRelBact,Status))
MMB_RelBact = model_MBBact[["1"]][["tab.df"]]
MMB_RelBact$Bacteria = "M.bovis RelTotal"

model_MHBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, MHRelBact,Status))
MMH_RelBact = model_MHBact[["1"]][["tab.df"]]
MMH_RelBact$Bacteria = "M.haemolytica RelTotal"

model_PMBact <- map(.x = params_list, 
                    .f = ~ mwu(meta2, PMRelBact,Status))
MPM_RelBact = model_PMBact[["1"]][["tab.df"]]
MPM_RelBact$Bacteria = "P.multocida RelTotal"

MannWu_pvals <- do.call(rbind, list(M16, M16log, MHS, MMB, MMH, MPM, MHS_Rel, MMB_Rel, MMH_Rel,
                                    MPM_Rel, MHS_RelBact, MMB_RelBact, MMH_RelBact, MPM_RelBact))
MannWu_pvals$State <- c("CO")
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResultsTotal = MannWu_pvals
ResultsTotal =rbind(ResultsTotal,MannWu_pvals)
write.csv(ResultsTotal, "ResultsTotalBeef.csv")

### for beef samples
metadata2
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

summary <-metadata2%>%
  group_by(State, Bacteria, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$comparison = paste(summary$State, summary$Bacteria, summary$Status, sep = "-")
my_colors <- c("dodgerblue3","goldenrod3")

ggplot(metadata2, aes(x=Bacteria, y=copy_log,color=Status)) +
  geom_point(aes(x=Bacteria, y=copy_log, color=Status),position = position_jitterdodge(0.4), alpha=0.2, size=1) +
  geom_errorbar(data=summary, aes(x=Bacteria, ymin=Mean - SE,
                                  ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_grid(State~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= summary, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("Copy number (log10)")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12, face = "italic"), axis.text.y =
          element_text(color = "black", size = 12))

## Additional Beef plot fot 16S rRNA abundance
summary3 <-metadata %>%
  group_by(State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(X16S_log),
    SD = sd(X16S_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary3$Bacteria = c("16S log")
#SumTotal16 = summary3
SumTotal16 = rbind(SumTotal16, summary3)
SumTotal16 = SumTotal16[c(5:8, 13:16),]

X16S = subset(SumTotal16, Bacteria == "16S copies")
X16L = subset(SumTotal16, Bacteria == "16S log")
filtered_data <- metadata[metadata$State %in% c("IN", "TX"), ]

x = ggplot(filtered_data, aes(x=Status, y=X16S.copies)) +
  geom_jitter(aes(x=Status, y=X16S.copies,color=Status), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16S, aes(x=Status, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ facet_grid(State~.) +
  geom_point(data= X16S, aes(Status, y=Mean, color=Status)) +
  ylab("16S rRNA copy number") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))

y = ggplot(filtered_data, aes(x=Status, y=X16S_log)) +
  geom_jitter(aes(x=Status, y=X16S_log,color=Status), width = 0.25, alpha=0.4) +
  geom_errorbar(data=X16L, aes(x=Status, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ facet_grid(State~.) +
  geom_point(data= X16L, aes(Status, y=Mean, color=Status)) +
  ylab("16S rRNA log copy") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 9))

ggarrange(x,y,ncol =2, labels = c("A", "B"),
          font.label = list(size = 11))

##### ---------- Bacteria prevalence------------------
#assign numerical values to factors
str(metadata2)
metadata2$Bacteria <- as.factor(metadata2$Bacteria)
levels(metadata2$Bacteria) <- list("H.somni"="H. somni", "P.multocida"="P.multocida", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica")
metadata2$pre_abs = as.factor(metadata2$pre_abs)
levels(metadata2$pre_abs) <- list("Negative"="0", "Positive"="1")

Bact2<- metadata2 %>% 
  group_by(Bacteria, State, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Bacteria, State, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup
str(Bact2)
Bact2$Status = as.factor(Bact2$Status)
levels(Bact2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
Bact2$Bacteria = as.factor(Bact2$Bacteria)
levels(Bact2$Bacteria) <- list("H.somni"="H.somni", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica", "P.multocida"="P.multocida")
Bact2$Samples = paste(Bact2$State,Bact2$Bacteria, Bact2$Status, Bact2$pre_abs , sep = "-")
Bact2P = subset(Bact2, pre_abs== "Positive")

ggplot(data=Bact2P, aes(x=Bacteria, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive %") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size= 11, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

  