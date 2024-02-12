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
library(sjstats)


setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/SpecificqPCR/")
metadata <- read.csv("MhPMSpecificMetadataDairyBeef.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$State <- as.factor(metadata$State)
metadata$pre_abs = as.factor(metadata$pre_abs)
levels(metadata$pre_abs) <- list("Negative"="0", "Positive"="1")

Beef = subset(metadata, cattle =="Beef") #495
Dairy = subset(metadata, cattle =="Dairy") #471

DairyCIN = subset(Dairy, State =="IN") #249
DairyCNY = subset(Dairy, State =="NY") #222

### calculate what is teh average copy log per pathobiont
sum <-Dairy %>%
  group_by(Serotype, State, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
write.csv(sum, "BRDpathobiontAbundanceSerotype.csv")
#Prevalence graph for each bacterium.
# Indiana samples
library(dplyr)
summary <-DairyCIN%>%
  group_by(Serotype, Status) %>% 
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary$State = "IN"

Bact2<- DairyCIN %>% 
  group_by(Serotype, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Serotype, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup

Bact2PIN = subset(Bact2, pre_abs== "Positive")
Bact2PIN$State = "IN"

# New York samples
summary2 <-DairyCNY%>%
  group_by(Serotype, Status) %>% 
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )
summary2$State = "NY"
DairyPrev = rbind(summary, summary2)


Bact2NYT<- DairyCNY %>% 
  group_by(Serotype, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Serotype, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup

Bact2PNY = subset(Bact2NYT, pre_abs== "Positive")
Bact2PNY$State = "NY"
BactSpe = rbind(Bact2PIN, Bact2PNY)
write.csv(BactSpe, "PrevalenceBactSpecificDairy.csv")

my_colors <- c("dodgerblue3","goldenrod3")

ggplot(data=BactSpe, aes(x=Serotype, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ 
  ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size= 11, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#test normality of data
summary_Shapiro <-DairyCIN %>%
  group_by(Serotype) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE,
    tidy(shapiro.test(copy_log))
  ) #everything is not normal distributed

summary_table <-DairyCIN %>%
  group_by(Serotype, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  ) #everything is not normal distributed
summary_table$State = "IN"
DairtBact = summary_table

summary_table <-DairyCNY %>%
  group_by(Serotype, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  ) #everything is not normal distributed
summary_table$State = "NY"
DairtBact = rbind(summary_table, DairtBact)

Dairy%>%
  ggplot( aes(x=Serotype, y=copy_log, fill=Status)) +
  geom_boxplot() +facet_grid(State~.)+
  geom_point(data= DairtBact, aes(Serotype, y=Mean, fill=Status), position=position_dodge(width = .8)) + 
  scale_fill_manual(values = my_colors) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() + ylab("Bacteria Abundance (log10)") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

##run statistics to see if there are any significance 
##prepare the data for the loop
MhA1= subset(DairyCIN, Serotype=="Mh-A1")
MhA6= subset(DairyCIN, Serotype=="Mh-A6")
PmA= subset(DairyCIN, Serotype=="Pm-A")

MhA1= subset(DairyCNY, Serotype=="Mh-A1")
MhA6= subset(DairyCNY, Serotype=="Mh-A6")
PmA= subset(DairyCNY, Serotype=="Pm-A")

meta2 =DairyCNY
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_MhA1 <- map(.x = params_list, 
                  .f = ~ mwu(MhA1, copy_log,Status))
model_MhA1 = model_MhA1[["1"]][["tab.df"]]
model_MhA1$Bacteria = "M. haemolytica A1"

model_MhA6 <- map(.x = params_list, 
                  .f = ~ mwu(MhA6, copy_log,Status))
model_MhA6 = model_MhA6[["1"]][["tab.df"]]
model_MhA6$Bacteria = "M. haemolytica A6"

model_PMA <- map(.x = params_list, 
                 .f = ~ mwu(PmA, copy_log,Status))
model_PMA = model_PMA[["1"]][["tab.df"]]
model_PMA$Bacteria = "P. multocida A"

MannWu_pvals <- do.call(rbind, list(model_MhA1,model_MhA6, model_PMA))
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
MannWu_pvals$State = "NY"
#SpecificResults = MannWu_pvals
SpecificResults = rbind(SpecificResults, MannWu_pvals)
write.csv(MannWu_pvals, "ResultsTotalDairySpecificPmMh.csv")

###############################################Beef samples
BeefCCO = subset(Beef, State =="CO") #246
BeefCID = subset(Beef, State =="ID") #249

summary <-Beef%>%
  group_by(Serotype, Status) %>% 
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  )

Bact2<- BeefCCO %>% 
  group_by(Serotype, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Serotype, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup
Bact2PCO = subset(Bact2, pre_abs== "Positive")
Bact2PCO$State = "CO"

# New York samples
Bact2IDT<- BeefCID %>% 
  group_by(Serotype, Status, pre_abs) %>% 
  summarise (n = n()) %>%
  group_by(Serotype, Status) %>%
  mutate(per =  100 *n/sum(n)) %>% 
  ungroup

Bact2PNY = subset(Bact2NYT, pre_abs== "Positive")
Bact2PNY$State = "ID"
BactSpe = rbind(Bact2PIN, Bact2PNY)
write.csv(BactSpe, "PrevalenceBactSpecificBeef.csv")

my_colors <- c("dodgerblue3","goldenrod3")

ggplot(data=BactSpe, aes(x=Serotype, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ 
  ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size= 11, face = "italic")) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

#### Beef abundance 
summary_table <-BeefCCO %>%
  group_by(Serotype, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  ) #everything is not normal distributed
summary_table$State = "CO"
BeefBact = summary_table

summary_table <-BeefCID %>%
  group_by(Serotype, Status) %>% # make one row per age
  # use group_by(species, age) if > 1 species
  summarize(
    n = n(),
    Mean = mean(copy_log),
    SD = sd(copy_log),
    SE = SD/sqrt(n),
    `25%` = Mean + qt(0.025, n-1) * SE,
    `75%` = Mean + qt(0.975, n-1) * SE
  ) #everything is not normal distributed
summary_table$State = "ID"
BeefBact = rbind(summary_table, BeefBact)

Beef%>%
  ggplot( aes(x=Serotype, y=copy_log, fill=Status)) +
  geom_boxplot() +facet_grid(State~.)+
  geom_point(data= BeefBact, aes(Serotype, y=Mean, fill=Status), position=position_dodge(width = .8)) + 
  scale_fill_manual(values = my_colors) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() + ylab("Bacteria Abundance (log10)") +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 14, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 10))

##run statistics to see if there are any significance 
##prepare the data for the loop
MhA1= subset(BeefCCO, Serotype=="Mh-A1")
MhA6= subset(BeefCCO, Serotype=="Mh-A6")
PmA= subset(BeefCCO, Serotype=="Pm-A")

MhA1= subset(BeefCID, Serotype=="Mh-A1")
MhA6= subset(BeefCID, Serotype=="Mh-A6")
PmA= subset(BeefCID, Serotype=="Pm-A")

meta2 =BeefCID
combins <- combn(levels(meta2$Status), 2)
params_list <- split(as.vector(combins), rep(1:ncol(combins), each = nrow(combins)))

#now run the wilcoxon test
#bacteria log abundance
model_MhA1 <- map(.x = params_list, 
                  .f = ~ mwu(MhA1, copy_log,Status))
model_MhA1 = model_MhA1[["1"]][["tab.df"]]
model_MhA1$Bacteria = "M. haemolytica A1"

model_MhA6 <- map(.x = params_list, 
                  .f = ~ mwu(MhA6, copy_log,Status))
model_MhA6 = model_MhA6[["1"]][["tab.df"]]
model_MhA6$Bacteria = "M. haemolytica A6"

model_PMA <- map(.x = params_list, 
                 .f = ~ mwu(PmA, copy_log,Status))
model_PMA = model_PMA[["1"]][["tab.df"]]
model_PMA$Bacteria = "P. multocida A"

MannWu_pvals <- do.call(rbind, list(model_MhA1,model_MhA6, model_PMA))
MannWu_pvals$p.adj<-p.adjust(MannWu_pvals$p, method = "BH")
MannWu_pvals$State = "ID"
#SpecificResults = MannWu_pvals
SpecificResults = rbind(SpecificResults, MannWu_pvals)
write.csv(MannWu_pvals, "ResultsTotalBeefSpecificPmMh.csv")

