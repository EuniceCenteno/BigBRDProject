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

my_colors <- c("dodgerblue3","goldenrod3")

a = ggplot(metadata, aes(x=Status, y=chao1)) +
  geom_jitter(aes(x=Status, y=chao1,color=Status), width = 0.25, alpha=0.2) +
  geom_errorbar(data=chaT, aes(x=Status, ymin=Mean - SE,
                                     ymax=Mean + SE, y=NULL, color=Status),  width=0.2) +
  theme_bw()+
  geom_point(data= chaT, aes(x=Status, y=Mean, color=Status)) + theme(legend.position = "none") +
  ylab("Chao1") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

b = ggplot(metadata, aes(x=Status, y=observed_otus)) +
  geom_jitter(aes(x=Status, y=observed_otus,color=Status), width = 0.25, alpha=0.2) +
  geom_errorbar(data=obsT, aes(x=Status, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL, color=Status),  width=0.2) +
  theme_bw()+
  geom_point(data= obsT, aes(x=Status, y=Mean, color=Status)) + theme(legend.position = "none") +
  ylab("Observed ASVs") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

c = ggplot(metadata, aes(x=Status, y=pielou_e)) +
  geom_jitter(aes(x=Status, y=pielou_e,color=Status), width = 0.25, alpha=0.2) +
  geom_errorbar(data=pieT, aes(x=Status, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL, color=Status),  width=0.2) +
  theme_bw()+
  geom_point(data= pieT, aes(x=Status, y=Mean, color=Status)) + theme(legend.position = "none") +
  ylab("Pielou_e") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

d = ggplot(metadata, aes(x=Status, y=shannon)) +
  geom_jitter(aes(x=Status, y=shannon,color=Status), width = 0.25, alpha=0.2) +
  geom_errorbar(data=shaT, aes(x=Status, ymin=Mean - SE,
                              ymax=Mean + SE, y=NULL, color=Status),  width=0.2) +
  theme_bw()+
  geom_point(data= shaT, aes(x=Status, y=Mean, color=Status)) + theme(legend.position = "none") +
  ylab("Shannon") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

e = ggplot(metadata, aes(x=Status, y=faith_pd)) +
  geom_jitter(aes(x=Status, y=faith_pd,color=Status), width = 0.25, alpha=0.2) +
  geom_errorbar(data=faT, aes(x=Status, ymin=Mean - SE,
                              ymax=Mean + SE, y=NULL, color=Status),  width=0.2) +
  theme_bw()+
  geom_point(data= faT, aes(x=Status, y=Mean, color=Status)) + theme(legend.position = "none") +
  ylab("Faith_PD") +xlab ("Status") +
  scale_color_manual(values = my_colors) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

ggarrange(a,b,c,d,nrow =4,
          font.label = list(size = 11))

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

chaTS = subset(SumTotalS, Alpha == "Chao1")
obsTS = subset(SumTotalS, Alpha == "Observed ASVs")
pieTS = subset(SumTotalS, Alpha == "Pielou_e")
shaTS = subset(SumTotalS, Alpha == "Shannon")
faTS = subset(SumTotalS, Alpha == "Faith_PD")

my_color <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

w = ggplot(metadata, aes(x=State, y=chao1)) +
  geom_jitter(aes(x=State, y=chao1,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=chaTS, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= chaTS, aes(State, y=Mean, color=State)) +
  ylab("Chao1") +xlab ("State") +
  scale_color_manual(values = my_color) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))
w

v = ggplot(metadata, aes(x=State, y=observed_otus)) +
  geom_jitter(aes(x=State, y=observed_otus,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=obsTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= obsTS, aes(State, y=Mean, color=State)) +
  ylab("Observed ASVs") +xlab ("State") +
  scale_color_manual(values = my_color) + guides(color="none") +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))
v

x = ggplot(metadata, aes(x=State, y=pielou_e)) +
  geom_jitter(aes(x=State, y=pielou_e,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=pieTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= pieTS, aes(State, y=Mean, color=State)) +
  ylab("Pielou_e") +xlab ("State") + guides(color="none") +
  scale_color_manual(values = my_color) +
  #stat_compare_means(label.y =1.20,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))
x

y = ggplot(metadata, aes(x=State, y=shannon)) +
  geom_jitter(aes(x=State, y=shannon,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=shaTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= shaTS, aes(State, y=Mean, color=State)) +
  ylab("Shannon") +xlab ("State") +
  scale_color_manual(values = my_color) + guides(color="none") +
  #stat_compare_means(label.y =11,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))
y

z = ggplot(metadata, aes(x=State, y=faith_pd)) +
  geom_jitter(aes(x=State, y=faith_pd,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=faTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= faTS, aes(State, y=Mean, color=State)) +
  ylab("Faith_PD") +xlab ("State") +
  scale_color_manual(values = my_color) + guides(color="none") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  #stat_compare_means(label.y =62,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))
z
ggarrange(w,v,y,x,z,ncol =3,nrow = 2,
          font.label = list(size = 11))

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
chao1 = metadata[,c(1,2,5,12)]
sha = metadata[,c(1,2,5,11)]
fa = metadata[,c(1,2,5,9)]
obs = metadata[,c(1,2,5,10)]
eve = metadata[,c(1,2,5,8)]

my_colors <- c("dodgerblue3","goldenrod3")

sumCha <-chao1%>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(chao1),
    SD = sd(chao1),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumCha$Status = as.factor(sumCha$Status)
levels(sumCha$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
sumCha$comparison = paste(sumCha$State,  sumCha$Status, sep = "-")
sumCha$Alpha = c("Chao1")

A = ggplot(chao1, aes(x=State, y=chao1,color=Status)) +
  geom_point(aes(x=State, y=chao1, color=Status),position = position_jitterdodge(0.4), alpha=0.3, size=1) +
  geom_errorbar(data=sumCha, aes(x=State, ymin=Mean - SE,
                                  ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Alpha~.) + theme_bw() +  scale_color_manual(values = my_colors) +
  geom_point(data= sumCha, aes(State, y=Mean, color=Status), position=position_dodge(width = .8)) +
   theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))
A
sumObs <-obs%>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(observed_otus),
    SD = sd(observed_otus),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumObs$Status = as.factor(sumObs$Status)
levels(sumObs$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
sumObs$comparison = paste(sumObs$State,  sumCha$Status, sep = "-")
sumObs$Alpha = c("Observed ASVs")

B = ggplot(obs, aes(x=State, y=observed_otus,color=Status)) +
  geom_point(aes(x=State, y=observed_otus, color=Status),position = position_jitterdodge(0.4), alpha=0.3, size=1) +
  geom_errorbar(data=sumObs, aes(x=State, ymin=Mean - SE,
                                 ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumObs, aes(State, y=Mean, color=Status), position=position_dodge(width = .8)) +
  xlab ("State") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

sumEv <-eve%>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(pielou_e),
    SD = sd(pielou_e),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumEv$Status = as.factor(sumEv$Status)
levels(sumEv$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
sumEv$comparison = paste(sumEv$State,  sumEv$Status, sep = "-")
sumEv$Alpha = c("Pielou_e")

C = ggplot(eve, aes(x=State, y=pielou_e,color=Status)) +
  geom_point(aes(x=State, y=pielou_e, color=Status),position = position_jitterdodge(0.4), alpha=0.3, size=1) +
  geom_errorbar(data=sumEv, aes(x=State, ymin=Mean - SE,
                                 ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumEv, aes(State, y=Mean, color=Status), position=position_dodge(width = .8)) +
  xlab ("State") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 10))

sumSha <-sha%>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(shannon),
    SD = sd(shannon),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumSha$Status = as.factor(sumSha$Status)
levels(sumSha$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
sumSha$comparison = paste(sumSha$State,  sumSha$Status, sep = "-")
sumSha$Alpha = c("Shannon")

D = ggplot(sha, aes(x=State, y=shannon,color=Status)) +
  geom_point(aes(x=State, y=shannon, color=Status),position = position_jitterdodge(0.4), alpha=0.3, size=1) +
  geom_errorbar(data=sumSha, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumSha, aes(State, y=Mean, color=Status), position=position_dodge(width = .8)) +
  xlab ("State") + 
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y =
          element_text(color = "black", size = 11))

sumFa <-fa%>%
  group_by(State,Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
sumFa$Status = as.factor(sumFa$Status)
levels(sumFa$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
sumFa$comparison = paste(sumFa$State,  sumFa$Status, sep = "-")
sumFa$Alpha = c("Faith_PD")

E = ggplot(fa, aes(x=State, y=faith_pd,color=Status)) +
  geom_point(aes(x=State, y=faith_pd, color=Status),position = position_jitterdodge(0.4), alpha=0.3, size=1) +
  geom_errorbar(data=sumFa, aes(x=State, ymin=Mean - SE,
                                 ymax=Mean + SE, y=NULL, color=Status), position = position_dodge(width=0.8), width=0.5) +
  facet_free(Alpha~.) + scale_color_manual(values = my_colors) + theme_bw() +
  geom_point(data= sumFa, aes(State, y=Mean, color=Status), position=position_dodge(width = .8)) +
  xlab ("State") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size  = 10))

ggarrange(A,B,C,D,E,nrow=5, common.legend = TRUE, legend = "right")

## Check if there is a difference in the bacteria abundance in the Warren and Southview samples
#because some samples were collected from cows and calves from the same farm

warren = subset(metadata, Farm=="Warren")
Southview = subset(metadata, Farm=="Southview")
Animal = rbind(warren, Southview)

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
write.csv(sumAni, "AlphaSummaryCowCalfWarranSouthview.csv")

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
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/") 
metadata2 <- read.csv("BRDBigMetadataBeef.csv", na.strings = c("","NA"), header=TRUE)
str(metadata2)
metadata2$Status <- as.factor(metadata2$Status)
metadata2$State <- as.factor(metadata2$State)
otu_table <- read.table("observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
shannon <-read.table("shannon.tsv", header=TRUE, row.names=1, sep="\t")
faith <- read.table("faith_pd.txt", header=TRUE, row.names=1, sep="\t")
evenness <- read.table("evenness.tsv", header=TRUE, row.names=1, sep="\t")
chao1 <- read.table("chao1.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(otu_table, shannon, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, faith, by.x = "Row.names", by.y = 0)
alpha_diversity <- merge(alpha_diversity, evenness, by.x = "Row.names", by.y = 0)
alpha_diversity <- merge(alpha_diversity, chao1, by.x = "Row.names", by.y = 0)
colnames(alpha_diversity) <- c("Row.names", "observed_otus", "shannon", "faith_pd", "pielou_e", "chao1")
metadata2 <- merge(metadata2, alpha_diversity, by.x = "ID", by.y = "Row.names")
row.names(metadata2) <- metadata2$ID
metadata2 = metadata2[-1]

levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
dplyr::count(metadata2, Status) 
dplyr::count(metadata2, State, Status) 
# Converting date of collection to numeric values
str(metadata2)
metadata2$date <- metadata2$Date_of_colleciton
metadata2$Date_of_colleciton <- as.Date(metadata2$Date_of_colleciton, "%m/%d/%y")
d<- as.Date('02/10/22', "%m/%d/%y") #use to calculate the days
metadata2$Date_of_colleciton <- as.Date(d) -as.Date(metadata2$Date_of_colleciton) 
metadata2$Date_of_colleciton <- as.numeric(metadata2$Date_of_colleciton)
str(metadata2$Date_of_colleciton)
#the highest day value is the date of the samples collected first 

## check if the data is normal distributed or not
set_sum_contrasts() # important for afex

# full model
str(metadata2)
summary(metadata2)


M1 <- mixed(pielou_e ~ Status +State + Date_of_colleciton+ (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(faith_pd ~ Status +State + Date_of_colleciton+ (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otus ~ Status +State + Date_of_colleciton+ (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(shannon ~ Status +State + Date_of_colleciton+ (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

M5 <- mixed(chao1  ~ Status +State + Date_of_colleciton+ (1|State), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5)

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata2$shannon)
plot(M2$full_model)
shapiro.test(metadata2$observed_otus)
plot(M3$full_model)
shapiro.test(metadata2$chao1)
plot(M4$full_model)
shapiro.test(metadata2$faith_pd)

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
metadata2 <- mutate(metadata2, pielou_e_log = log10(pielou_e + 1))
metadata2 <- mutate(metadata2, faith_pdlog = log10(faith_pd + 1))
metadata2 <- mutate(metadata2, observed_otuslog = log10(observed_otus + 1))
metadata2 <- mutate(metadata2, shannonlog = log10(shannon + 1))
metadata2 <- mutate(metadata2, chao1log = log10(chao1 + 1))


M1 <- mixed(pielou_e_log ~ Status + Farm + (1|Farm), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1) ##NA for farm

M2 <- mixed(faith_pdlog ~ Status + Farm + (1|Farm), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(observed_otuslog ~ Status + Farm + (1|Farm), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3) #nA farm

M4 <- mixed(shannonlog ~ Status + Farm + (1|Farm), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4) #NA for farm

M5 <- mixed(chao1log  ~ Status + Farm + (1|Farm), data = metadata2,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5) #NA for farm

#checking assumptions
# way 1:
plot(M1$full_model)
shapiro.test(metadata2$shannonlog)
plot(M2$full_model)
shapiro.test(metadata2$observed_otuslog)
plot(M3$full_model)
shapiro.test(metadata2$chao1log)
plot(M5$full_model)
shapiro.test(metadata2$faith_pdlog)

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
summary <-metadata2 %>%
  group_by(Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(faith_pd),
    SD = sd(faith_pd),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$Status = as.factor(summary$Status)
levels(summary$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
summary$Alpha = c("Faith_PD")
#SumTotal = summary
SumTotal = rbind(SumTotal, summary)

chaT = subset(SumTotal, Alpha == "Chao1")
obsT = subset(SumTotal, Alpha == "Observed ASVs")
pieT = subset(SumTotal, Alpha == "Pielou_e")
shaT = subset(SumTotal, Alpha == "Shannon")
faT = subset(SumTotal, Alpha == "Faith_PD")

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

#--- now overall State effect
my_comparisons <- list( c("CO", "ID"), c("CO", "IN"), c("CO", "TX"),
                        c("ID", "IN"), c("ID", "TX"), c("IN", "TX"))

summary2 <-metadata2 %>%
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

chaTS = subset(SumTotalS, Alpha == "Chao1")
obsTS = subset(SumTotalS, Alpha == "Observed ASVs")
pieTS = subset(SumTotalS, Alpha == "Pielou_e")
shaTS = subset(SumTotalS, Alpha == "Shannon")
faTS = subset(SumTotalS, Alpha == "Faith_PD")

my_color1 <- c(
  "#653936", "#CD9BCD","orange","#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

w = ggplot(metadata2, aes(x=State, y=chao1)) +
  geom_jitter(aes(x=State, y=chao1,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=chaTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= chaTS, aes(State, y=Mean, color=State)) +
  ylab("Chao1") +xlab ("State") + guides(color="none") +
  scale_color_manual(values = my_color1) + ylim(0,3750) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 11))

v = ggplot(metadata2, aes(x=State, y=observed_otus)) +
  geom_jitter(aes(x=State, y=observed_otus,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=obsTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= obsTS, aes(State, y=Mean, color=State)) +
  ylab("Observed ASVs") +xlab ("State") + guides(color="none") +
  scale_color_manual(values = my_color1) + ylim(0,3100) +
  #stat_compare_means(label.y =1150,label.x = 1.1, size=3.5)+ theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 11))

x = ggplot(metadata2, aes(x=State, y=pielou_e)) +
  geom_jitter(aes(x=State, y=pielou_e,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=pieTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= pieTS, aes(State, y=Mean, color=State)) +
  ylab("Pielou_e") +xlab ("State") + ylim(0,1.28) +
  scale_color_manual(values = my_color1) + guides(color="none") +
 stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 11))

y = ggplot(metadata2, aes(x=State, y=shannon)) +
  geom_jitter(aes(x=State, y=shannon,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=shaTS, aes(x=State, ymin=Mean - SE,
                                ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= shaTS, aes(State, y=Mean, color=State)) +
  ylab("Shannon") +xlab ("State") + guides(color="none") +
  scale_color_manual(values = my_color1) + ylim(0,15) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 11))

z = ggplot(metadata2, aes(x=State, y=faith_pd)) +
  geom_jitter(aes(x=State, y=faith_pd,color=State), width = 0.25, alpha=0.4) +
  geom_errorbar(data=faTS, aes(x=State, ymin=Mean - SE,
                               ymax=Mean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+
  geom_point(data= faTS, aes(State, y=Mean, color=State)) +
  ylab("Faith_PD") +xlab ("State") + guides(color="none") +
  scale_color_manual(values = my_color1) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(11, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"),
        axis.title.y = element_text(color="black", size=11, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y =
          element_text(color = "black", size = 11))

ggarrange(w,y,v,z,x, nrow=2, ncol = 3, common.legend = FALSE)

#identify differences between the alpha diversity within each farm

##prepare the data for the loop
metadata2$ID = row.names(metadata2)
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
#ObsS$p <- gsub("<", "", ObsS$p)

model_PieS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, pielou_e,Status))
PieS = model_PieS[["1"]][["tab.df"]]
PieS$Alpha = "Pielou_e"

model_ShaS <- map(.x = params_list, 
                  .f = ~ mwu(meta2, shannon,Status))
ShaS = model_ShaS[["1"]][["tab.df"]]
ShaS$Alpha = "Shannon"
#ShaS$p <- gsub("<", "", ShaS$p)

model_FaS <- map(.x = params_list, 
                 .f = ~ mwu(meta2, faith_pd,Status))
FaS = model_FaS[["1"]][["tab.df"]]
FaS$Alpha = "Faith_PD"
#FaS$p <- gsub("<", "", FaS$p)

MannWu_pvals <- do.call(rbind, list(ChaS, ObsS, PieS, ShaS, FaS))
MannWu_pvals$State = "TX"
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
#ResulTotalStatus =MannWu_pvals
ResulTotalStatus = rbind(ResulTotalStatus, MannWu_pvals)
write.csv(ResulTotalStatus, "ResultsTotalBeefAlphaDiv.csv")
