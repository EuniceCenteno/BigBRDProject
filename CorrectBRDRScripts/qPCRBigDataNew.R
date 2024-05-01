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
set_sum_contrasts() # important for afex
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

my_color2 <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)# for dairy samples

x= ggplot(meta1, aes(x=State, y=X16S_log)) +
  geom_boxplot(aes(x=State, y=X16S_log,fill=State),  alpha= 0.9) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color2) + theme_bw() +
  geom_jitter(aes(x=State, y=X16S_log,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("16S rRNA log10") +xlab ("State") + 
  scale_color_manual(values = my_color2) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisonsD, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold"), axis.text.y =
          element_text(color = "black", size = 9))

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

##### plot with bacteria abundance
metadata2 <- read.csv("BRDBigqPCRDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
meta2 = metadata2[metadata2$X16S.copies != "230", ]

summary <-meta2%>%
  group_by(State, Bacteria, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$comparison = paste(summary$State, summary$Bacteria, summary$Status, sep = "-")
my_colors <- c("dodgerblue3","goldenrod3")
meta2$Status = as.factor(meta2$Status)

ggplot(data=meta2, aes(x=Bacteria, y=copy_log)) +
  geom_point(aes(x=Bacteria, y=copy_log, color=Status),position = position_jitterdodge(0.2),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Bacteria, y=copy_log,fill= Status), position=position_dodge(width=0.8), alpha=0.4, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_point(data=summary, aes(x=Bacteria, y=Mean, color=Status), position = "identity", size=3.5, shape=17) +
  facet_grid(State~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Copy Number (log10)") + 
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

## total abundance regardless of state and bacteria
summary <-meta2%>%
  group_by(Status, State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
my_colors <- c("dodgerblue3","goldenrod3")
meta2$Status = as.factor(meta2$Status)

ggplot(data=meta2, aes(x=Status, y=copy_log)) +
  geom_point(aes(x=Status, y=copy_log, color=Status),position = position_jitterdodge(0.2),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Status, y=copy_log,fill= Status), position=position_dodge(width=0.8), alpha=0.4) +
  scale_fill_manual(values = my_colors) +
  geom_point(data=summary, aes(x=Status, y=Mean, color=Status), position = "identity", size=3.5, shape=17) +
  facet_grid(State~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Copy Number (log10)") + 
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

### statistical test for prevalence data
#only work with prevalence of 1 to detect differences
meta2= CA
meta2$HSpre_abs = as.factor(meta2$HSpre_abs)
levels(meta2$HSpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$MBpre_abs = as.factor(meta2$MBpre_abs)
levels(meta2$MBpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$MHpre_abs = as.factor(meta2$MHpre_abs)
levels(meta2$MHpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$PMpre_abs = as.factor(meta2$PMpre_abs)
levels(meta2$PMpre_abs) <- list("Detected"="1", "NotDetected"="0")
str(meta2)
## now we run chi-squere test for the prevalence data
result <- chisq.test(table(meta2$Status, meta2$PMpre_abs))
PMPre = as.data.frame(result[["p.value"]])
PMPre$Bacteria = "P.multocidaPre"

result <- chisq.test(table(meta2$Status, meta2$MBpre_abs))
MBPre = as.data.frame(result[["p.value"]])
MBPre$Bacteria = "M.bovisPre"

result <- chisq.test(table(meta2$Status, meta2$MHpre_abs))
MHPre = as.data.frame(result[["p.value"]])
MHPre$Bacteria = "M.haemolyticaPre"

result <- chisq.test(table(meta2$Status, meta2$HSpre_abs))
HSPre = as.data.frame(result[["p.value"]])
HSPre$Bacteria = "H.somniPre"
Prevalence = rbind(PMPre, MBPre,MHPre,HSPre)
str(Prevalence)
Prevalence$p.adjust<-p.adjust(Prevalence$`result[["p.value"]]`, method = "BH")
Prevalence$State = "TX"
#TotalPreva= Prevalence
TotalPreva = rbind(Prevalence, TotalPreva)

write.csv(TotalPreva, "DairyBactPrevalence.csv")

#now run the wilcoxon test
#bacteria log abundance

## total abundance of status for each bacteria
summary6 <-meta2%>%
  group_by(Status, Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = quantile(copy_log, probs = 0.25),
    `75%` = quantile(copy_log, probs = 0.75)
  )
my_colors <- c("dodgerblue3","goldenrod3")
meta2$Status = as.factor(meta2$Status)

ggplot(data=meta2, aes(x=Status, y=copy_log)) +
  geom_point(aes(x=Status, y=copy_log, color=Status),position = position_jitterdodge(0.2),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Status, y=copy_log,fill= Status), position=position_dodge(width=0.8), alpha=0.4) +
  scale_fill_manual(values = my_colors) +
  geom_point(data=summary6, aes(x=Status, y=Mean, color=Status), position = "identity", size=3.5, shape=17) +
  facet_grid(Bacteria~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Copy Number (log10)") + 
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#threshold 
thre = subset(summary6, Status=="Healthy")
thre$`75%` #75% quantile
hsT = 4.355794
mbT = 2.826075
mhT= 4.007284
pmT= 4.795613

Hsv = subset(meta2, Bacteria =="H.somni")
dplyr::count(Hsv, State, Status)
Hsv = Hsv[Hsv$copy_log > hsT , ]
dplyr::count(Hsv, State, Status)

Mbv = subset(meta2, Bacteria =="M.bovis")
dplyr::count(Mbv, State, Status)
Mbv = Mbv[Mbv$copy_log > mbT , ]
dplyr::count(Mbv, State, Status)

Mhv = subset(meta2, Bacteria =="M.haemolytica")
dplyr::count(Mhv, State, Status)
Mhv = Mhv[Mhv$copy_log > mhT , ]
dplyr::count(Mhv, State, Status)

Pmv = subset(meta2, Bacteria =="P.multocida")
dplyr::count(Pmv, State, Status)
Pmv = Pmv[Pmv$copy_log > pmT , ]
dplyr::count(Pmv, State, Status)

#calculate 60% quantile
summary6 <-meta2%>%
  group_by(Status, Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `60%` = quantile(copy_log, probs = 0.60)
  )
thre = subset(summary6, Status=="Healthy")
thre$`60%` #60% quantile
hsT = 3.721381
mbT = 2.826075
mhT= 2.230449
pmT= 3.575700

Hsv = Hsv[Hsv$copy_log > hsT , ]
dplyr::count(Hsv, State, Status)

Mbv = Mbv[Mbv$copy_log > mbT , ]
dplyr::count(Mbv, State, Status)

Mhv = Mhv[Mhv$copy_log > mhT , ]
dplyr::count(Mhv, State, Status)

Pmv = Pmv[Pmv$copy_log > pmT , ]
dplyr::count(Pmv, State, Status)

#plot for 16S copy number ~ NY samples
my_colors <- c("dodgerblue3","goldenrod3")
str(NY)
y = ggplot(NY, aes(x=Status, y=X16S_log,fill=Status)) +
  geom_boxplot(aes(x=Status, y=X16S_log,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  geom_jitter(aes(x=Status, y=X16S_log,color=Status), width = 0.25, alpha=0.8, shape=19) +
  scale_color_manual(values = my_colors)+
  ylab("16S rRNA copy number (log10)") +xlab("NY") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

ggarrange(x,y,ncol =2,labels = c("a", "b"),
          font.label = list(size = 18))


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
my_colors <- c("dodgerblue3","goldenrod3")

a = ggplot(NY, aes(x=Status, y=MBRel,fill=Status)) +
  geom_boxplot(aes(x=Status, y=MBRel,fill=Status)) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= CO, aes(x=Status, y=HSRel, fill=Status), position=position_dodge(width = .8)) +
  ylab("H.somni 16S rRNA Rel") + ggtitle("NY") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold"), axis.text.y =
          element_text(color = "black", size = 10))
a

b = ggplot(NY, aes(x=Status, y=PMRel,fill=Status)) +
  geom_boxplot(aes(x=Status, y=PMRel,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("P.multocida 16S rRNA Rel") + ggtitle("NY") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold"), axis.text.y =
          element_text(color = "black", size = 10))
b

c = ggplot(CA, aes(x=Status, y=MHRelBact,fill=Status)) +
  geom_boxplot(aes(x=Status, y=MHRelBact,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("M.haemolytica Bact Rel") + ggtitle("CA") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
c
ggarrange(a,b,c,ncol =3,labels = c("a", "b", "c"),
          font.label = list(size = 18))

## to plot the significance of M. haemolytica in CA samples

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
write.csv(summary, "SummaryPathobiontsDairy.csv")
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
meta2
#assign numerical values to factors
str(meta2)
meta2$Bacteria <- as.factor(meta2$Bacteria)
levels(meta2$Bacteria) <- list("H.somni"="H. somni", "P.multocida"="P.multocida", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica")
meta2$pre_abs = as.factor(meta2$pre_abs)
levels(meta2$pre_abs) <- list("Negative"="0", "Positive"="1")

Bact2<- meta2 %>% 
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
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

## Check if there is a difference in the bacteria abundance in the farm3 and farm4 samples
#because some samples were collected from cows and calves from the same farm

farm3 = subset(metadata2, Farm=="Farm3")
farm4 = subset(metadata2, Farm=="Farm4")
Animal = rbind(farm3, farm4)

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
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/MetadataCorrectBeef/")
metadata <- read.csv("BRDBigqPCRBeefTotalNamesRow.csv", na.strings = c("","NA"), header=TRUE)
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
write.csv(metadata, "CompleteBeefqPCRDatawithRelativeAbundance.csv")

meta1 = metadata[metadata$X16S.copies != "230", ]
write.csv(meta1, "qPCRDataNo16Scopies230Beef.csv")

##how many samples we have
sum<- meta1 %>% 
  group_by(State, Status) %>% 
  summarise (n = n())

#generate model
str(metadata)
set_sum_contrasts() # important for afex
M1 <- mixed(Hscopy_log ~ Status + State +  (1|State), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$Hscopy_log) #not normal
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

M2 <- mixed(MBcopy_log ~ Status + State +  (1|State), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)
shapiro.test(metadata$MBcopy_log) #not normal
plot(M2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

M3 <- mixed(MHcopy_log ~Status + State +  (1|State), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
shapiro.test(metadata$MHcopy_log) #not normal
plot(M3$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

M4 <- mixed(PMcopy_log ~ Status + State +  (1|State), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)
shapiro.test(metadata$PMcopy_log) #not normal 
plot(M4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

## statistical analysis for 16S rRNA abundance 
#generate model
str(metadata)
M1 <- mixed(X16S.copies ~ Status + State +  (1|State), data = meta1,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M1
anova(M1)
shapiro.test(metadata$X16S.copies) #not normal
plot(M1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model)) # not normaly distributed

#transform the 16S rRNA copies to log values
M2 <- mixed(X16S_log ~ Status + State +  (1|State), data = meta1,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
M2
anova(M2)
shapiro.test(metadata$X16S_log) #not normal
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
write.csv(MannWuT_pvals, "OverallBacteriaabundanceonlyStatusBeef.csv")
# Mb, Mh and HS are significant 


#to create plots
metadata2 <- read.csv("BRDBigqPCRBeefTotal.csv", na.strings = c("","NA"), header=TRUE)
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
meta1total = metadata2[metadata2$X16S.copies != "230", ]
qPCR = meta1total

my_colors <- c("dodgerblue3","goldenrod3")

ggplot(data=qPCR, aes(x=Status, y=copy_log)) +
  geom_point(aes(x=Status, y=copy_log, color=Status),position = position_jitterdodge(0.8),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Status, y=copy_log,fill= Status), alpha=0.5, outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  scale_fill_manual(values = my_colors) +
  facet_grid(.~Bacteria) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Copy Number (log10)") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.x = element_text(face = "italic", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=12, face="bold"))

str(qPCR)
summary2 <-meta1total %>%
  group_by(Status, Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary2
write.csv(summary2, "OverallPathobiontAbundStatus.csv")

### identify the difference of 16S rRNA abundance among the different states
my_comparisonsD <- list( c("CO", "ID"), c("CO", "IN"), c("CO", "TX"),
                         c("ID", "IN"), c("ID", "TX"), c("IN", "TX")) ## comparison for dairy samples
summary2 <-meta1 %>%
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
summary2$Bacteria = c("16S log")
write.csv(summary2, "16SrRNAAbundanceSummaryBeef.csv")

my_color2 <- c(
  "#653936", "#CD9BCD","orange","#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)# for beef samples

ggplot(meta1, aes(x=State, y=X16S_log)) +
  geom_boxplot(aes(x=State, y=X16S_log,fill=State),  alpha= 0.9) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color2) + theme_bw() +
  geom_jitter(aes(x=State, y=X16S_log,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("16S rRNA log10") +xlab ("State") + 
  scale_color_manual(values = my_color2) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisonsD, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold"), axis.text.y =
          element_text(color = "black", size = 9))

##prepare the data for the loop, to identify the effect within state
IN= subset(meta1, State=="IN")
ID= subset(meta1, State=="ID")  
CO= subset(meta1, State=="CO")
TX= subset(meta1, State=="TX")

meta2=TX
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
MannWu_pvals$State <- c("TX")
MannWu_pvals$p <- gsub("<", "", MannWu_pvals$p)
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResultsTotal = MannWu_pvals
ResultsTotal =rbind(ResultsTotal,MannWu_pvals)
write.csv(ResultsTotal, "ResultsTotalBeef.csv") 
#plots

### calculate what is teh average copy log per pathobiont
str()
sum <-meta1total %>%
  group_by(Bacteria, State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median = median(copy_log), 
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(sum, "BRDpathobiontAbundanceBeef.csv")

#get the mean of the relative abundance and create plots for the significant results
str(metadata)
sumRel <-meta1 %>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(HSRelBact),
    SD = sd(HSRelBact),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumRel$Bacteria = c("H.somni Bact Rel") 
#sumRelT = sumRel
sumRelT = rbind(sumRelT, sumRel)
write.csv(sumRelT, "summaryBactRelBeef.csv")

## plot relative abundance significant results
CO = subset(meta1, State == "CO")
IN = subset(meta1, State == "IN")
ID = subset(meta1, State == "ID")
TX= subset(meta1, State == "TX")

## PLot for 16S rRNA abundance IN samples
my_colors <- c("dodgerblue3","goldenrod3")
str(IN)
ggplot(IN, aes(x=Status, y=X16S_log,fill=Status)) +
  geom_boxplot(aes(x=Status, y=X16S_log,fill=Status)) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= CO, aes(x=Status, y=HSRel, fill=Status), position=position_dodge(width = .8)) +
  ylab("16S rRNA log10") + ggtitle("IN") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))


a = ggplot(CO, aes(x=Status, y=HSRel,fill=Status)) +
  geom_boxplot(aes(x=Status, y=HSRel,fill=Status)) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= CO, aes(x=Status, y=HSRel, fill=Status), position=position_dodge(width = .8)) +
  ylab("H.somni 16S rRNA Rel") + ggtitle("CO") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
a

b = ggplot(TX, aes(x=Status, y=HSRel,fill=Status)) +
  geom_boxplot(aes(x=Status, y=HSRel,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("H.somni 16S rRNA Rel") + ggtitle("TX") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
b

c = ggplot(TX, aes(x=Status, y=HSRelBact,fill=Status)) +
  geom_boxplot(aes(x=Status, y=HSRelBact,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("H.somni Bact Total") + ggtitle("TX") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
c

d = ggplot(IN, aes(x=Status, y=HSRelBact,fill=Status)) +
  geom_boxplot(aes(x=Status, y=HSRelBact,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("H.somni Bact Total") + ggtitle("IN") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
d

e = ggplot(IN, aes(x=Status, y=PMRelBact,fill=Status)) +
  geom_boxplot(aes(x=Status, y=PMRelBact,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("P.multocida Bact Total") + ggtitle("IN") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
e

f = ggplot(TX, aes(x=Status, y=PMRelBact,fill=Status)) +
  geom_boxplot(aes(x=Status, y=PMRelBact,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("P.multocida Bact Total") + ggtitle("TX") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
f

g = ggplot(TX, aes(x=Status, y=MBRel,fill=Status)) +
  geom_boxplot(aes(x=Status, y=MBRel,fill=Status)) +
  scale_fill_manual(values = my_colors) + theme_bw() +
  #geom_point(data= sumRelTX, aes(Bacteria, y=Mean, color=Status), position=position_dodge(width = .8)) +
  ylab("M.bovis 16S rRNA Rel") + ggtitle("TX") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
g

ggarrange(a,b,g, ncol =3,labels = c("a", "b", "c"),
          font.label = list(size = 18))
ggarrange(d,c,e,f, ncol =2,nrow=2, labels = c("a", "b", "c", "d"),
          font.label = list(size = 18))

## to plot the significance of M. haemolytica in CA samples
meta2 = metadata2[metadata2$X16S.copies != "230", ]

summary <-meta2%>%
  group_by(State, Bacteria, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = quantile(copy_log, probs = 0.25),
    `75%` = quantile(copy_log, probs = 0.75)
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$comparison = paste(summary$State, summary$Bacteria, summary$Status, sep = "-")
my_colors <- c("dodgerblue3","goldenrod3")
meta2$Status = as.factor(meta2$Status)

ggplot(data=meta2, aes(x=Bacteria, y=copy_log)) +
  geom_point(aes(x=Bacteria, y=copy_log, color=Status),position = position_jitterdodge(0.2),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Bacteria, y=copy_log,fill= Status), position=position_dodge(width=0.8), alpha=0.4, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_point(data=summary, aes(x=Bacteria, y=Mean, color=Status), position = "identity", size=3.5, shape=17) +
  facet_grid(State~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Copy Number (log10)") + 
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#threshold 
summary6 <-meta2%>%
  group_by(Status,Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = quantile(copy_log, probs = 0.25),
    `75%` = quantile(copy_log, probs = 0.75)
  )
thre = subset(summary6, Status=="Healthy")
thre$`75%` #75% quantile
hsT = 5.883568
mbT = 3.818378
mhT= 4.111656
pmT= 5.599793

Hsv = subset(meta2, Bacteria =="H.somni")
dplyr::count(Hsv, State, Status)
Hsv = Hsv[Hsv$copy_log > hsT , ]
dplyr::count(Hsv, State, Status)

Mbv = subset(meta2, Bacteria =="M.bovis")
dplyr::count(Mbv, State, Status)
Mbv = Mbv[Mbv$copy_log > mbT , ]
dplyr::count(Mbv, State, Status)

Mhv = subset(meta2, Bacteria =="M.haemolytica")
dplyr::count(Mhv, State, Status)
Mhv = Mhv[Mhv$copy_log > mhT , ]
dplyr::count(Mhv, State, Status)

Pmv = subset(meta2, Bacteria =="P.multocida")
dplyr::count(Pmv, State, Status)
Pmv = Pmv[Pmv$copy_log > pmT , ]
dplyr::count(Pmv, State, Status)

#calculate 60% quantile
summary6 <-meta2%>%
  group_by(Status,Bacteria) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    Median= median(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `60%` = quantile(copy_log, probs = 0.60)
  )
thre = subset(summary6, Status=="Healthy")
thre$`60%` #60% quantile
hsT = 5.112249
mbT = 2.826075
mhT= 3.355963
pmT= 4.944655

Hsv = Hsv[Hsv$copy_log > hsT , ]
dplyr::count(Hsv, State, Status)

Mbv = Mbv[Mbv$copy_log > mbT , ]
dplyr::count(Mbv, State, Status)

Mhv = Mhv[Mhv$copy_log > mhT , ]
dplyr::count(Mhv, State, Status)

Pmv = Pmv[Pmv$copy_log > pmT , ]
dplyr::count(Pmv, State, Status)

##### ---------- Bacteria prevalence------------------
meta2
#assign numerical values to factors
str(meta2)
meta2$Bacteria <- as.factor(meta2$Bacteria)
levels(meta2$Bacteria) <- list("H.somni"="H. somni", "P.multocida"="P.multocida", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica")
meta2$pre_abs = as.factor(meta2$pre_abs)
levels(meta2$pre_abs) <- list("Negative"="0", "Positive"="1")

Bact2<- meta2 %>% 
  group_by(Bacteria, State, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Bacteria, State, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup
str(Bact2)
write.csv(Bact2, "BacteriaPrevalenceBeef.csv")
Bact2$Status = as.factor(Bact2$Status)
levels(Bact2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
Bact2$Bacteria = as.factor(Bact2$Bacteria)
levels(Bact2$Bacteria) <- list("H.somni"="H.somni", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica", "P.multocida"="P.multocida")
Bact2$Samples = paste(Bact2$State,Bact2$Bacteria, Bact2$Status, Bact2$pre_abs , sep = "-")
Bact2P = subset(Bact2, pre_abs== "Positive")
write.csv(Bact2P, "BacteriaPrevalencePositiveBeef.csv")

ggplot(data=Bact2P, aes(x=Bacteria, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive %") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

### statistical test for prevalence data
meta2= CO
meta2$HSpre_abs = as.factor(meta2$HSpre_abs)
levels(meta2$HSpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$MBpre_abs = as.factor(meta2$MBpre_abs)
levels(meta2$MBpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$MHpre_abs = as.factor(meta2$MHpre_abs)
levels(meta2$MHpre_abs) <- list("Detected"="1", "NotDetected"="0")
meta2$PMpre_abs = as.factor(meta2$PMpre_abs)
levels(meta2$PMpre_abs) <- list("Detected"="1", "NotDetected"="0")
str(meta2)
## now we run chi-squere test for the prevalence data
result <- chisq.test(table(meta2$Status, meta2$PMpre_abs))
PMPre = as.data.frame(result[["p.value"]])
PMPre$Bacteria = "P.multocidaPre"

result <- chisq.test(table(meta2$Status, meta2$MBpre_abs))
MBPre = as.data.frame(result[["p.value"]])
MBPre$Bacteria = "M.bovisPre"

result <- chisq.test(table(meta2$Status, meta2$MHpre_abs))
MHPre = as.data.frame(result[["p.value"]])
MHPre$Bacteria = "M.haemolyticaPre"

result <- chisq.test(table(meta2$Status, meta2$HSpre_abs))
HSPre = as.data.frame(result[["p.value"]])
HSPre$Bacteria = "H.somniPre"
Prevalence = rbind(PMPre, MBPre,MHPre,HSPre)
str(Prevalence)
Prevalence$p.adjust<-p.adjust(Prevalence$`result[["p.value"]]`, method = "BH")
Prevalence$State = "TX"
#TotalPreva= Prevalence
TotalPreva = rbind(Prevalence, TotalPreva)

write.csv(TotalPreva, "BeefBactPrevalence.csv")
