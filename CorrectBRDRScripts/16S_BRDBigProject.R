library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(ggplot2)
library("cowplot")
theme_set(theme_grey())
#install.packages("sjstats")
library(jtools)
library(ggpubr)
library(sjstats)
library(ggdist)
library(tidyquant)
library(dplyr)
library(tidyverse)
#install.packages("wesanderson")
library(wesanderson)
library(qiime2R)

rm(list = ls ())

setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/CorrectQiimeDairy/") 

metadata <- read.csv("BRDBigMetadataDairy.csv", na.strings = c("","NA"), header=TRUE) #368 samples
str(metadata)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
dplyr::count(metadata, Status) 
metadata$State <- as.factor(metadata$State)
dplyr::count(metadata, State, Status) 

summary <-metadata %>%
  group_by(State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Faith PD"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumDairyAlphaStates.csv")

summary <-metadata %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    Median= median(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Faith PD"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumDairyAlphaStatus.csv")

summary <-metadata %>%
  group_by(State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Faith PD"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumDairyAlphaStatesStatus.csv")

## check if the data is normal distributed or not
set_sum_contrasts() # important for afex

# full model
str(metadata)


M1 <- mixed(pielou_e ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(faith_pd ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otus ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(shannon ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

M5 <- mixed(chao1  ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5)

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata$shannon)
plot(M2$full_model)
shapiro.test(metadata$observed_otus)
plot(M3$full_model)
shapiro.test(metadata$chao1)
plot(M5$full_model)
shapiro.test(metadata$faith_pd)

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

qqnorm(residuals(M5$full_model))
qqline(residuals(M5$full_model))

##test log transformation
metadata <- mutate(metadata, pielou_e_log = log10(pielou_e + 1))
metadata <- mutate(metadata, faith_pdlog = log10(faith_pd + 1))
metadata <- mutate(metadata, observed_otuslog = log10(observed_otus + 1))
metadata <- mutate(metadata, shannonlog = log10(shannon + 1))
metadata <- mutate(metadata, chao1log = log10(chao1 + 1))


M1 <- mixed(pielou_e_log ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) ##NA for farm

M2 <- mixed(faith_pdlog ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otuslog ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(shannonlog ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) #NA for farm

M5 <- mixed(chao1log  ~ Status + Farm + (1|Farm), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5) #NA for farm

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata$shannonlog)
plot(M2$full_model)
shapiro.test(metadata$observed_otuslog)
plot(M3$full_model)
shapiro.test(metadata$chao1log)
plot(M5$full_model)
shapiro.test(metadata$faith_pdlog)

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

qqnorm(residuals(M5$full_model))
qqline(residuals(M5$full_model))

## non-parametric data, we will use 
#Mann-Whitney test to analyze disease status
#kruskal wallis to detect different betweeen farms

##---- Identify the overall effect of disease and location in the alpha diversity
summary <-metadata %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(pielou_e),
    SD = sd(pielou_e),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$Alpha = c("Pielou_e")
#SumTotal = summary
SumTotal = rbind(SumTotal, summary)

chaT = subset(SumTotal, Alpha == "Chao1")
obsT = subset(SumTotal, Alpha == "Observed ASVs")
pieT = subset(SumTotal, Alpha == "Pielou_e")
shaT = subset(SumTotal, Alpha == "Shannon")
faT = subset(SumTotal, Alpha == "Faith_PD")

##now let's run the Mann-Whitney test
meta2=metadata
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

model_Cha <- map(.x = params_list, 
                .f = ~ mwu(meta2, chao1,Status))
Cha = model_Cha[["1"]][["tab.df"]]
Cha$Alpha = "Chao1"

model_Obs <- map(.x = params_list, 
                 .f = ~ mwu(meta2, observed_otus,Status))
Obs = model_Obs[["1"]][["tab.df"]]
Obs$Alpha = "Observed_ASVs"

model_Pie <- map(.x = params_list, 
                 .f = ~ mwu(meta2, pielou_e,Status))
Pie = model_Pie[["1"]][["tab.df"]]
Pie$Alpha = "Pielou_e"

model_Sha <- map(.x = params_list, 
                 .f = ~ mwu(meta2, shannon,Status))
Sha = model_Sha[["1"]][["tab.df"]]
Sha$Alpha = "Shannon"

model_Fa <- map(.x = params_list, 
                 .f = ~ mwu(meta2, faith_pd,Status))
Fa = model_Fa[["1"]][["tab.df"]]
Fa$Alpha = "Faith_PD"

MannWu_pvals <- do.call(rbind, list(Cha, Obs, Pie, Sha, Fa))
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
write.csv(MannWu_pvals, "ResultsTotalDairyAlphaStatus.csv")

#Plots
my_colors <- c("dodgerblue3","goldenrod3")
metadata2 <- read.csv("BRDBigMetadataDairyTotal.csv", na.strings = c("","NA"), header=TRUE) #368 samples
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata2$alpha = as.factor(metadata2$alpha)
levels(metadata2$alpha) <- list("Observed ASVs"="Observed ASVs", "Chao1"="Chao1",
                                "Pielou_e" = "Pielou_e", "Shannon"= "Shannon", "Faith PD"= "Faith PD")
metadata2$State = as.factor(metadata2$State)
# we only need Observed ASVs, Chao1, Pielou, Shannon
metadata3 = metadata2
level_to_remove <- "Faith PD"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]
level_to_remove <- "Shannon"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]
level_to_remove <- "Chao1"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]

ggplot(data=metadata3, aes(x=Status, y=value)) +
  geom_point(aes(x=Status, y=value, color=Status),position = position_jitterdodge(0.8),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Status, y=value,fill= Status), alpha=0.5, outlier.shape = NA) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  scale_fill_manual(values = my_colors) +
  facet_free(alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Alpha Diversity Metrics") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "bold")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#--- now overall State effect
my_comparisons <- list( c("CA", "IN"), c("CA", "NY"), c("CA", "TX"),
c("IN", "NY"), c("IN", "TX"), c("NY", "TX"))

summary2 <-metadata %>%
  group_by(State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary2$Alpha = c("Faith_PD")
#SumTotalS = summary2
SumTotalS = rbind(SumTotalS, summary2)

my_color <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

a = ggplot(metadata, aes(x=State, y=observed_otus)) +
  geom_boxplot(aes(x=State, y=observed_otus,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color) + theme_bw() +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  geom_jitter(aes(x=State, y=observed_otus,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Observed ASVs") +xlab ("Farm") + labs(fill = "Farm") +
  scale_color_manual(values = my_color) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
a

c = ggplot(metadata, aes(x=State, y=pielou_e)) +
  geom_boxplot(aes(x=State, y=pielou_e,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color) + theme_bw() +
  geom_jitter(aes(x=State, y=pielou_e,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Pielou_e") +xlab ("Farm") + labs(fill = "Farm") +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  scale_color_manual(values = my_color) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
c

e = ggplot(metadata, aes(x=State, y=faith_pd)) +
  geom_boxplot(aes(x=State, y=faith_pd,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color) + theme_bw() +
  geom_jitter(aes(x=State, y=faith_pd,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Faith PD") +xlab ("Farm") + labs(fill = "Farm") +
  scale_color_manual(values = my_color) + guides(color="none") +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
e
e

ggarrange(a,c,e,ncol =3)


##prepare the data for the loop
IN= subset(metadata, State=="IN")
NY= subset(metadata, State=="NY") 
CA= subset(metadata, State=="CA")
TX= subset(metadata, State=="TX")

meta2 = CA
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the Mann-Whitney test
model_ChaS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, chao1,Status))
ChaS = model_ChaS[["1"]][["tab.df"]]
ChaS$Alpha = "Chao1"

model_ObsS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, observed_otus,Status))
ObsS = model_ObsS[["1"]][["tab.df"]]
ObsS$Alpha = "Observed_ASVs"
ObsS$p <- gsub("<", "", ObsS$p)

model_PieS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, pielou_e,Status))
PieS = model_PieS[["1"]][["tab.df"]]
PieS$Alpha = "Pielou_e"

model_ShaS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, shannon,Status))
ShaS = model_ShaS[["1"]][["tab.df"]]
ShaS$Alpha = "Shannon"
ShaS$p <- gsub("<", "", ShaS$p)

model_FaS <- map(.x = params_list, 
                .f = ~ mwu(meta2, faith_pd,Status))
FaS = model_FaS[["1"]][["tab.df"]]
FaS$Alpha = "Faith_PD"
FaS$p <- gsub("<", "", FaS$p)

MannWu_pvals <- do.call(rbind, list(ChaS, ObsS, PieS, ShaS, FaS))
MannWu_pvals$State = "CA"
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResulTotalStatus =MannWu_pvals
ResulTotalStatus = rbind(ResulTotalStatus, MannWu_pvals)
write.csv(ResulTotalStatus, "ResultsTotalDairyAlphaDiv.csv")

#plots
metadata2 <- read.csv("BRDBigMetadataDairyTotal.csv", na.strings = c("","NA"), header=TRUE) #368 samples
metadata2$Status = as.factor(metadata2$Status)
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata2$alpha = as.factor(metadata2$alpha)
levels(metadata2$alpha) <- list("Observed ASVs"="Observed ASVs", "Chao1"="Chao1",
                                "Pielou_e" = "Pielou_e", "Shannon"= "Shannon", "Faith PD"= "Faith PD")
metadata2$State = as.factor(metadata2$State)

meanAlDairy <- read.csv("SumDairyAlphaStatesStatus.csv", na.strings = c("","NA"), header=TRUE) #368 samples
meanAlDairy$Status = as.factor(meanAlDairy$Status)
levels(meanAlDairy$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
meanAlDairy$alpha = as.factor(meanAlDairy$alpha)
levels(meanAlDairy$alpha) <- list("Observed ASVs"="Observed ASVs", "Chao1"="Chao1",
                                "Pielou_e" = "Pielou_e", "Shannon"= "Shannon", "Faith PD"= "Faith PD")
level_to_remove <- "Shannon"
metadata2 <- metadata2[metadata2$alpha != level_to_remove, ]
level_to_remove <- "Chao1"
metadata2 <- metadata2[metadata2$alpha != level_to_remove, ]

ggplot(data=metadata2, aes(x=State, y=value)) +
  geom_point(aes(x=State, y=value, color=Status),position = position_jitterdodge(0.4),alpha=0.6, size=1) +
  geom_boxplot(aes(x=State, y=value,fill= Status), alpha=0.5, outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  facet_free(alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Alpha Diversity Metrics") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "bold")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=12, face="bold"))

## Check if there is a difference in the bacteria abundance in the Farm3 and Farm3 samples
#because some samples were collected from cows and calves from the same farm

farm3 = subset(metadata, Farm=="Farm3")
farm4 = subset(metadata, Farm=="Farm4")
Animal = rbind(farm3, farm4)

sumAni <-Animal %>%
  group_by(Farm,animal, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(observed_otus),
    SD = sd(observed_otus),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumAni$alpha = "Observed ASVs"
pie = sumAni
cha = sumAni
sha = sumAni
fai = sumAni
obs = sumAni
AlphaAnimal =rbind(pie, cha, sha, fai, obs)
write.csv(sumAni, "AlphaSummaryCowCalfFarm3-4.csv")

#3 Run stats
##prepare the data for the loop, to identify the effect within state
meta3=IN
meta3$animal = as.factor(meta3$animal)
combins <- combn(levels(meta3$animal), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_pie <- map(.x = params_list, 
                      .f = ~ mwu(meta3, pielou_e,animal))
Mpie = model_pie[["1"]][["tab.df"]]
Mpie$Alpha = "Pielou_e"

model_cha <- map(.x = params_list, 
                 .f = ~ mwu(meta3, chao1,animal))
Mcha = model_cha[["1"]][["tab.df"]]
Mcha$Alpha = "Chao1"
Mcha$p <- gsub("<", "", Mcha$p)

model_sha <- map(.x = params_list, 
                 .f = ~ mwu(meta3, shannon,animal))
Msha = model_sha[["1"]][["tab.df"]]
Msha$Alpha = "Shannon"

model_obs <- map(.x = params_list, 
                 .f = ~ mwu(meta3, observed_otus,animal))
Mobs = model_obs[["1"]][["tab.df"]]
Mobs$Alpha = "Observed ASVs"

model_fai <- map(.x = params_list, 
                 .f = ~ mwu(meta3, faith_pd,animal))
Mfai = model_fai[["1"]][["tab.df"]]
Mfai$Alpha = "Faith PD"

MAnimal <- do.call(rbind, list(Mpie, Mcha, Mobs, Mfai, Msha))
MAnimal$State <- c("IN")
MAnimal$p.adj<-p.adjust(MAnimal$p, method = "BH")
#ResultsAnimal = MAnimal
ResultsAnimal =rbind(ResultsAnimal,MAnimal)
write.csv(ResultsAnimal, "ResultsAlphaAnimalNYIN.csv")

## check if the presence of cows in the NY samples introduce bias/error
columnName <- "animal"  # Replace with the actual column name
valueToRemove <- "Cow"     # Replace with the value to be removed

# Create a new data frame without the specified rows
NY_filtered <- NY[NY[[columnName]] != valueToRemove, ] #from 98 to 79

meta4=NY_filtered
meta4$Status = as.factor(meta4$Status)
combins <- combn(levels(meta4$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

model_cha <- map(.x = params_list, 
                 .f = ~ mwu(meta4, chao1,Status))
Mcha = model_cha[["1"]][["tab.df"]]
Mcha$Alpha = "Chao1"
Mcha$p <- gsub("<", "", Mcha$p)
model_obs <- map(.x = params_list, 
                 .f = ~ mwu(meta4, observed_otus,Status))
Mobs = model_obs[["1"]][["tab.df"]]
Mobs$Alpha = "Observed ASVs"
MAnimal <- do.call(rbind, list(Mcha, Mobs))
MAnimal$p.adj<-p.adjust(MAnimal$p, method = "BH")

#---------------------beef data
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/") 
metadata2 <- read.csv("BRDBigqPCRBeefTotalNamesRow.csv", na.strings = c("","NA"), header=TRUE)
str(metadata2)
metadata2$Status <- as.factor(metadata2$Status)
metadata2$State <- as.factor(metadata2$State)
otu_table <- read.table("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/exported/observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
faith <- read.table("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/exported/faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
evenness <- read.table("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/exported/evenness.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(otu_table, faith, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, evenness, by.x = "Row.names", by.y = 0)
colnames(alpha_diversity) <- c("Row.names", "observed_otus","faith_pd", "pielou_e")
metadata2 <- merge(metadata2, alpha_diversity, by.x = "ID", by.y = "Row.names")
row.names(metadata2) <- metadata2$ID
metadata2 = metadata2[-2]
write.csv(metadata2, "BRDBigBeefMetadataTotal.csv")

levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
dplyr::count(metadata2, Status) 
dplyr::count(metadata2, State, Status) 
# Converting date of collection to numeric values
#the highest day value is the date of the samples collected first 

## check if the data is normal distributed or not
summary <-metadata2 %>%
  group_by(State) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    Median= median(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Faith PD"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumBeefAlphaStates.csv")

## check if the data is normal distributed or not
summary <-metadata2 %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    Median= median(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Faith PD"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumBeefAlphaStatus.csv")

summary <-metadata2 %>%
  group_by(State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(shannon),
    Median= median(shannon),
    SD = sd(shannon),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$alpha = "Shannon"
#SumT = summary
SumT = rbind(SumT, summary)
write.csv(SumT, "SumBeefAlphaStatesStatus.csv")

## check if the data is normal distributed or not
set_sum_contrasts() # important for afex
# full model
str(metadata)

M1 <- mixed(pielou_e ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(faith_pd ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otus ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(shannon ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

M5 <- mixed(chao1  ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5)

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata2$shannon) # NS
plot(M2$full_model)
shapiro.test(metadata2$observed_otus) #NS
plot(M3$full_model)
shapiro.test(metadata2$chao1) #NS
plot(M4$full_model)
shapiro.test(metadata2$pielou_e) #NS
plot(M5$full_model)
shapiro.test(metadata2$faith_pd) #NS

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

qqnorm(residuals(M5$full_model))
qqline(residuals(M5$full_model))

##test log transformation
metadata2 <- mutate(metadata2, pielou_e_log = log10(pielou_e + 1))
metadata2 <- mutate(metadata2, faith_pdlog = log10(faith_pd + 1))
metadata2 <- mutate(metadata2, observed_otuslog = log10(observed_otus + 1))
metadata2 <- mutate(metadata2, shannonlog = log10(shannon + 1))
metadata2 <- mutate(metadata2, chao1log = log10(chao1 + 1))

M1 <- mixed(pielou_e_log ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) ##NA for farm

M2 <- mixed(faith_pdlog ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otuslog ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(shannonlog ~ Status +State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) #NA for farm

M5 <- mixed(chao1log  ~ Status + State + (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5) #NA for farm

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata2$shannonlog) #NS
plot(M2$full_model)
shapiro.test(metadata2$observed_otuslog) #NS
plot(M3$full_model)
shapiro.test(metadata2$chao1log) #NS
plot(M4$full_model)
shapiro.test(metadata2$pielou_e_log) #NS
plot(M5$full_model)
shapiro.test(metadata2$faith_pdlog) #NS

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

qqnorm(residuals(M5$full_model))
qqline(residuals(M5$full_model))

## non-parametric data, we will use 
#Mann-Whitney test to analyze disease status
#kruskal wallis to detect different betweeen farms

##---- Identify the overall effect of disease and location in the alpha diversity
summary <-metadata2 %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    Median= median(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$Alpha = c("faith PD")
#SumTotal = summary
SumTotal = rbind(SumTotal, summary)

chaT = subset(SumTotal, Alpha == "Chao1")
obsT = subset(SumTotal, Alpha == "Observed ASVs")
pieT = subset(SumTotal, Alpha == "Pielou_e")
shaT = subset(SumTotal, Alpha == "Shannon")
faT = subset(SumTotal, Alpha == "faith PD")

##now let's run the Mann-Whitney test
meta2=metadata2
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

model_Cha <- map(.x = params_list, 
                 .f = ~ mwu(meta2, chao1,Status))
Cha = model_Cha[["1"]][["tab.df"]]
Cha$Alpha = "Chao1"

model_Obs <- map(.x = params_list, 
                 .f = ~ mwu(meta2, observed_otus,Status))
Obs = model_Obs[["1"]][["tab.df"]]
Obs$Alpha = "Observed_ASVs"

model_Pie <- map(.x = params_list, 
                 .f = ~ mwu(meta2, pielou_e,Status))
Pie = model_Pie[["1"]][["tab.df"]]
Pie$Alpha = "Pielou_e"

model_Sha <- map(.x = params_list, 
                 .f = ~ mwu(meta2, shannon,Status))
Sha = model_Sha[["1"]][["tab.df"]]
Sha$Alpha = "Shannon"

model_Fa <- map(.x = params_list, 
                .f = ~ mwu(meta2, faith_pd,Status))
Fa = model_Fa[["1"]][["tab.df"]]
Fa$Alpha = "Faith_PD"

MannWu_pvals <- do.call(rbind, list(Cha, Obs, Pie, Sha, Fa))
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
write.csv(MannWu_pvals, "ResultsTotalBeefAlphaStatus.csv")

my_colors <- c("dodgerblue3","goldenrod3")
meta3 <- read.csv("BRDBigBeefMetadataAlphaTotal.csv", na.strings = c("","NA"), header=TRUE) #368 samples
meta3$Status = as.factor(meta3$Status)
levels(meta3$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
meta3$alpha = as.factor(meta3$alpha)
levels(meta3$alpha) <- list("Observed ASVs"="Observed ASVs", "Chao1"="Chao1",
                                "Pielou_e" = "Pielou_e", "Shannon"= "Shannon", "Faith PD"= "Faith PD")
meta3$State = as.factor(meta3$State)
# we only need Observed ASVs, Chao1, Pielou, Shannon
metadata3 = meta3
level_to_remove <- "Shannon"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]
levels(metadata3$alpha)
level_to_remove <- "Pielou_e"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]
levels(metadata3$alpha)
level_to_remove <- "Chao1"
metadata3 <- metadata3[metadata3$alpha != level_to_remove, ]
levels(metadata3$alpha)

ggplot(data=metadata3, aes(x=Status, y=value)) +
  geom_point(aes(x=Status, y=value, color=Status),position = position_jitterdodge(0.8),alpha=0.6, size=1) +
  geom_boxplot(aes(x=Status, y=value,fill= Status), alpha=0.5,outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2.5, color="red", fill="red")+
  facet_free(alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  theme_bw() + ylab("Alpha Diversity Metrics") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(size= 12, face = "bold")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#--- now overall State effect
my_comparisons <- list( c("CO", "ID"), c("CO", "IN"), c("CO", "TX"),
                        c("ID", "IN"), c("ID", "TX"), c("IN", "TX"))

my_color2 <- c(
  "#653936", "#CD9BCD","orange","#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)# for beef samples

a = ggplot(metadata2, aes(x=State, y=observed_otus)) +
  geom_boxplot(aes(x=State, y=observed_otus,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color2) + theme_bw() +
  geom_jitter(aes(x=State, y=observed_otus,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Observed ASVs") +xlab ("Farm") + labs(fill = "Farm") +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  scale_color_manual(values = my_color2) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
a

c = ggplot(metadata2, aes(x=State, y=pielou_e)) +
  geom_boxplot(aes(x=State, y=pielou_e,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color2) + theme_bw() +
  geom_jitter(aes(x=State, y=pielou_e,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Pielou_e") +xlab ("Farm") + labs(fill = "Farm") +
  scale_color_manual(values = my_color2) + guides(color="none") +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
c

e = ggplot(metadata2, aes(x=State, y=faith_pd)) +
  geom_boxplot(aes(x=State, y=faith_pd,fill=State),  alpha= 0.9, outlier.shape = NA) +
  #facet_grid(.~Bacteria) + 
  scale_fill_manual(values = my_color2) + theme_bw() +
  geom_jitter(aes(x=State, y=faith_pd,color=State), width = 0.25, alpha=0.6, shape=19) +
  ylab("Faith PD") +xlab ("Farm") + labs(fill = "Farm") +
  scale_color_manual(values = my_color2) + guides(color="none") +
  stat_summary(fun.y=mean, geom="point", shape=17, size=2, color="red", fill="red")+
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 9))
e
e

ggarrange(a,c,e,ncol =3)

##prepare the data for the loop
IN= subset(metadata2, State=="IN")
ID= subset(metadata2, State=="ID") 
CO= subset(metadata2, State=="CO")
TX= subset(metadata2, State=="TX")

meta2 = TX
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the Mann-Whitney test
model_ChaS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, chao1,Status))
ChaS = model_ChaS[["1"]][["tab.df"]]
ChaS$Alpha = "Chao1"

model_ObsS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, observed_otus,Status))
ObsS = model_ObsS[["1"]][["tab.df"]]
ObsS$Alpha = "Observed_ASVs"

model_PieS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, pielou_e,Status))
PieS = model_PieS[["1"]][["tab.df"]]
PieS$Alpha = "Pielou_e"

model_ShaS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, shannon,Status))
ShaS = model_ShaS[["1"]][["tab.df"]]
ShaS$Alpha = "Shannon"

model_FaS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, faith_pd,Status))
FaS = model_FaS[["1"]][["tab.df"]]
FaS$Alpha = "Faith_PD"

MannWu_pvals <- do.call(rbind, list(ChaS, ObsS, PieS, ShaS, FaS))
MannWu_pvals$State = "TX"
MannWu_pvals$p <- gsub("<", "", MannWu_pvals$p)
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResulTotalStatus =MannWu_pvals
ResulTotalStatus = rbind(ResulTotalStatus, MannWu_pvals)
write.csv(ResulTotalStatus, "ResultsTotalBeefAlphaDivStatus.csv")
