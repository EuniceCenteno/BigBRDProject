#qPCR data
library(tidyverse)
library(ggpubr)
library(rstatix)
library(lattice)
library("latticeExtra")
library(dplyr)
library(scales)
library(tidyr)
library(ggplot2)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(tidyr) #separate function
library(reshape2) #melt function
library(phyloseq)
library(stats)
citation('stats')

### Dairy samples test correlation between pathobiont abundance and 16S rRNA
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/qPCRMetadata/MetadataCorrectDairy/")
metadata <- read.csv("BRDBigqPCRDairyTotalOneRow.csv", na.strings = c("","NA"), header=TRUE)
metadata$Status <- as.factor(metadata$Status)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$State <- as.factor(metadata$State)
str(metadata)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 
metadata$State <- as.factor(metadata$State)
dplyr::count(metadata, State, Status) 

ASVs <- read_qza("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/CorrectQiimeDairy/table-no-pseudo.qza") #402 samples, but it will be 400 after removing Dairy 175, 176
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
ASV_table <- merge(metadata, ASV_table, by.x = "Name", by.y = 0) #should be 400 samples
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:21)]
ASV_table <- as.matrix(ASV_table)

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/CorrectQiimeDairy/taxonomyNoPseudo.qza")
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
physeq_deseq # [ 12051 taxa and 400 samples ] dairy

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
physeq_deseq #[ 11992 taxa and 400 samples ] Dairy

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
NewTax$ASVs = rownames(NewTax)

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

### subset ASVs classified as BRD pathobionts
#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
Hs = subset(NewTax, Genus == " Histophilus")
Pm = subset(NewTax, Genus == " Pasteurella")
Mh = subset(NewTax, Genus == " Mannheimia")
Mb = subset(NewTax, Genus == " Mycoplasma")
patho = rbind(Hs, Pm, Mh, Mb)
patho$Species = as.factor(patho$Species)
levels(patho$Species)
patho <- patho[patho$Species %in% c(" Mycoplasma", " Pasteurella_multocida",
                                    " Histophilus", " Mannheimia_haemolytica"), ]
levels(patho$Species)
Myco = subset(patho, Species ==" Mycoplasma")
Paste = subset(patho, Species ==" Pasteurella_multocida")
Mann = subset(patho, Species ==" Mannheimia_haemolytica")
Histo = subset(patho, Species ==" Histophilus")

NewASVTable2 = t(NewASVTable)
Myco_table <- as.data.frame(t(NewASVTable2[rownames(NewASVTable2) %in% rownames(Myco),]))
MycoSum = as.data.frame(rowSums(Myco_table))
Myco_table = merge(Myco_table, MycoSum, by.x = 0, by.y = 0)
str(Myco_table)
Myco_table$prevalence <- ifelse(Myco_table$`rowSums(Myco_table)` > 0, "Positive", "Negative")
mycoAbun = metadata[,c(2,3,4,12:14)]
Myco_table = merge(Myco_table, mycoAbun, by.x = "Row.names", by.y = 0)
Myco_table$MBpre_abs = as.factor(Myco_table$MBpre_abs)
levels(Myco_table$MBpre_abs) <- list("Negative"="0", "Positive"="1")

MycoSum<- Myco_table %>% 
  group_by(State, Status,MBpre_abs) %>% 
  summarise (n = n()) %>%
  group_by(State, Status) %>%
  mutate(per =  100 *n/sum(n))
MycoSum$Bacteria = "Mycoplasma qPCR"
#MycoSumT = MycoSum
names(MycoSum)[names(MycoSum) == "MBpre_abs"] <- "prevalence"
MycoSumT = rbind(MycoSumT,MycoSum)
MycoSumTPos = subset(MycoSumT, prevalence== "Positive")

Paste_table <- as.data.frame(t(NewASVTable2[rownames(NewASVTable2) %in% rownames(Paste),]))
PasteSum = as.data.frame(rowSums(Paste_table))
Paste_table = merge(Paste_table, PasteSum, by.x = 0, by.y = 0)
str(Paste_table)
Paste_table$prevalence <- ifelse(Paste_table$`rowSums(Paste_table)` > 0, "Positive", "Negative")
pasteAbun = metadata[,c(2,3,4,18:20)]
Paste_table = merge(Paste_table, pasteAbun, by.x = "Row.names", by.y = 0)
Paste_table$PMpre_abs = as.factor(Paste_table$PMpre_abs)
levels(Paste_table$PMpre_abs) <- list("Negative"="0", "Positive"="1")

PasteSum<- Paste_table %>% 
  group_by(State, Status, prevalence) %>% 
  summarise (n = n()) %>%
  group_by(State, Status) %>%
  mutate(per =  100 *n/sum(n))
PasteSum$Bacteria = "P.multocida 16S"
#PasteSumT = PasteSum
#names(PasteSumT)[names(PasteSumT) == "PMpre_abs"] <- "prevalence"
PasteSumT = rbind(PasteSumT,PasteSum)
PasteSumTPos = subset(PasteSumT, prevalence== "Positive")

Mann_table <- as.data.frame(t(NewASVTable2[rownames(NewASVTable2) %in% rownames(Mann),]))
MannSum = as.data.frame(rowSums(Mann_table))
Mann_table = merge(Mann_table, MannSum, by.x = 0, by.y = 0)
str(Mann_table)
Mann_table$prevalence <- ifelse(Mann_table$`rowSums(Mann_table)` > 0, "Positive", "Negative")
mannAbun = metadata[,c(2,3,4,15:17)]
Mann_table = merge(Mann_table, mannAbun, by.x = "Row.names", by.y = 0)
Mann_table$MHpre_abs = as.factor(Mann_table$MHpre_abs)
levels(Mann_table$MHpre_abs) <- list("Negative"="0", "Positive"="1")

MannSum<- Mann_table %>% 
  group_by(State, Status, prevalence) %>% 
  summarise (n = n()) %>%
  group_by(State, Status) %>%
  mutate(per =  100 *n/sum(n))
MannSum$Bacteria = "M. haemolytica 16S"
#MannSumT = MannSum
#names(MannSumT)[names(MannSumT) == "MHpre_abs"] <- "prevalence"
MannSumT = rbind(MannSumT,MannSum)
MannSumTPos = subset(MannSumT, prevalence== "Positive")

Histo_table <- as.data.frame(t(NewASVTable2[rownames(NewASVTable2) %in% rownames(Histo),]))
HistoSum = as.data.frame(rowSums(Histo_table))
Histo_table = merge(Histo_table, HistoSum, by.x = 0, by.y = 0)
str(Histo_table)
Histo_table$prevalence <- ifelse(Histo_table$`rowSums(Histo_table)` > 0, "Positive", "Negative")
histoAbun = metadata[,c(2,3,4,9:11)]
Histo_table = merge(Histo_table, histoAbun, by.x = "Row.names", by.y = 0)
Histo_table$HSpre_abs = as.factor(Histo_table$HSpre_abs)
levels(Histo_table$HSpre_abs) <- list("Negative"="0", "Positive"="1")

HistoSum<- Histo_table %>% 
  group_by(State, Status, prevalence) %>% 
  summarise (n = n()) %>%
  group_by(State, Status) %>%
  mutate(per =  100 *n/sum(n))
HistoSum$Bacteria = "H. somni 16S"
#HistoSumT = HistoSum
#names(HistoSumT)[names(HistoSumT) == "HSpre_abs"] <- "prevalence"
HistoSumT = rbind(HistoSumT,HistoSum)
HistoSumTPos = subset(HistoSumT, prevalence== "Positive")

TotalPreBac = rbind(HistoSumTPos, PasteSumTPos, MannSumTPos, MycoSumTPos)
TotalPreBac$Status = as.factor(TotalPreBac$Status)
levels(TotalPreBac$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
TotalPreBac$Bacteria = as.factor(TotalPreBac$Bacteria)
#levels(TotalPreBac$Bacteria) <- list("H.somni"="H.somni", "M.bovis"="M.bovis", "M.haemolytica"="M.haemolytica", "P.multocida"="P.multocida")

my_colors <- c("dodgerblue3","goldenrod3")
ggplot(TotalPreBac, aes(x=Bacteria, y=per, fill=Status)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5)+
  scale_fill_manual(values = my_colors) +
  facet_grid(State~.)+ ylim(0,100) +
  theme_bw() + ylab("Prevalece Positive %") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(strip.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size= 9, face = "italic", angle = 55, hjust = 1)) +
  theme(axis.text.y = element_text(size= 10)) +
  theme(axis.title.x = element_text(size=12, face="bold")) +
  theme(axis.title.y = element_text(size=12, face="bold"))

##calculate the correlation between the pathobionts relative abundance and qPCR data
healthy <- subset(metadata, Status=="Healthy")
dplyr::count(healthy, ID)
BRD <- subset(metadata, Status=="BRD")
dplyr::count(BRD, ID)

H_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(healthy),] 
B_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(BRD),] 

BRD <- BRD[rownames(BRD) %in% rownames(B_OTU),]
healthy <- healthy[rownames(healthy) %in% rownames(H_OTU),]

##start with the healthy samples
OTU.physeqH = otu_table(as.matrix(H_OTU), taxa_are_rows=FALSE)
tax.physeqH = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqH = sample_data(healthy)
physeq_H = phyloseq(OTU.physeqH, tax.physeqH, meta.physeqH)

##Making the tree, we need the phyloseq object 
library("ape") #to create the tree
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeqH, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot = physeqH1

str(NewTax)
NewTax$Species = as.factor(NewTax$Species)
my_level <- c("Species")

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()}  #%>%                               # Melt to long format
taxa <- taxa.summary %>% 
  group_by(Sample,  Status, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Sample, Status, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance    
test = subset(taxa, Sample =="Dairy10")
sum(test$taxa.average)
physeq.taxa.max <- taxa %>% 
  group_by(Sample, Species) %>%
  summarise(overall.max=max(taxa.average))

test = subset(physeq.taxa.max, Sample =="Dairy10")
sum(test$overall.max)
physeq.taxa.max <- as.data.frame(physeq.taxa.max)
# merging the phyla means with the metadata #
physeq_metaH <- merge(taxa, physeq.taxa.max)
physeq_metaH = as.data.frame(physeq_metaH)
write.csv(physeq_metaH, "physeq_meta_totalSpeciesHealthyDairyState.csv")

#select the pathobionts
pathogen = as.factor(patho$Species)
physeqmetaH = physeq_metaH
physeqmetaH <- physeqmetaH[physeqmetaH$Species %in% c(" Mycoplasma", " Pasteurella_multocida", " Histophilus", 
                                                " Mannheimia_haemolytica"), ]

## relative abundance across all the samples and diseas status
#BRD animals
str(BRD)
OTU.physeqD = otu_table(as.matrix(B_OTU), taxa_are_rows=FALSE)
tax.physeqD = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqD = sample_data(BRD)
physeq_D = phyloseq(OTU.physeqD, tax.physeqD, meta.physeqD)

##Making the tree, we need the phyloseq object 
random_treeD = rtree(ntaxa(physeq_D), rooted=TRUE, tip.label=taxa_names(physeq_D))
##Merging the tree with the phyloseq object
physeqD1 = merge_phyloseq(physeq_D, meta.physeqD, random_treeD)
physeqD1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plotD = physeqD1

my_level <- c("Species")
rm(taxa.summary)

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotD %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()}  #%>%                               # Melt to long format
taxa <- taxa.summary %>% 
  group_by(Sample,  Status, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Sample, Status, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance    
test = subset(taxa, Sample =="Dairy1")
sum(test$taxa.average)
physeq.taxa.max <- taxa %>% 
  group_by(Sample, Species) %>%
  summarise(overall.max=max(taxa.average))
physeq.taxa.max <- as.data.frame(physeq.taxa.max)
# merging the phyla means with the metadata #
physeq_metaB <- merge(taxa, physeq.taxa.max)
write.csv(physeq_metaB, "physeq_meta_totalSpeciesHealthyDairyState.csv")

physeqmetaB = physeq_metaB
physeqmetaB <- physeqmetaB[physeqmetaB$Species %in% c(" Mycoplasma", " Pasteurella_multocida", " Histophilus", 
                                                      " Mannheimia_haemolytica"), ]
PathRel = rbind(physeqmetaB, physeqmetaH)
dplyr::count(PathRel, Sample) 
PathRel = PathRel[,-c(3)]
PathRel = merge(PathRel, metadata, by.x = "Sample", by.y = 0)

#Calculate pathobionts relative abundance based on 16S 
PathRel$HSRel <- (PathRel$HScopies / PathRel$X16S.copies)
PathRel$PMRel <- (PathRel$PMcopies /PathRel$X16S.copies) 
PathRel$MBRel <- (PathRel$MBcopies / PathRel$X16S.copies) 
PathRel$MHRel <- (PathRel$MHcopies / PathRel$X16S.copies) 
PathRel = PathRel[PathRel$X16S.copies != "230", ]

###combine qPCR data and 16sabund
TotalPM = subset(PathRel, Species ==" Pasteurella_multocida")
TotalMB = subset(PathRel, Species ==" Mycoplasma")
TotalMH = subset(PathRel, Species ==" Mannheimia_haemolytica")
TotalHS = subset(PathRel, Species ==" Histophilus")

library(ROCR)

str(TotalPM)
cor(TotalPM$taxa.average, TotalPM$PMRel)
cor.test(TotalPM$taxa.average, TotalPM$PMRel) #significant 

my_colors <- c("dodgerblue3","goldenrod3")

library("ggpubr")
a = ggscatter(TotalPM, x = "taxa.average", y = "PMRel",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 0, label.y = 2)  # Add correlation coefficient
a

#by status
TotalPM$Status = as.factor(TotalPM$Status)
b = ggscatter(TotalPM, x = "taxa.average", y = "PMRel",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "Status", palette = my_colors,           # Color by groups "cyl"
          shape = "Status"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Status), label.x = .4, label.y.npc = .1) 
b

#mycoplasma
str(TotalMB)
cor(TotalMB$taxa.average, TotalMB$MBRel)
cor.test(TotalMB$taxa.average, TotalMB$MBRel) #not significant 

c = ggscatter(TotalMB, x = "taxa.average", y = "MBRel",
              add = "reg.line",                                 # Add regression line
              conf.int = TRUE,                                  # Add confidence interval
              add.params = list(color = "blue",
                                fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 0, label.y = 10)  # Add correlation coefficient
c

#by status
TotalMB$Status = as.factor(TotalMB$Status)
d = ggscatter(TotalMB, x = "taxa.average", y = "MBRel",
              add = "reg.line",                         # Add regression line
              conf.int = TRUE,                          # Add confidence interval
              color = "Status", palette = my_colors,           # Color by groups "cyl"
              shape = "Status"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Status), label.x = .15, label.y.npc = .4) 
d

#mannheimia
str(TotalMH)
cor(TotalMH$taxa.average, TotalMH$MHRel)
cor.test(TotalMH$taxa.average, TotalMH$MHRel) #significant 

e = ggscatter(TotalMH, x = "taxa.average", y = "MHRel",
              add = "reg.line",                                 # Add regression line
              conf.int = TRUE,                                  # Add confidence interval
              add.params = list(color = "blue",
                                fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 0, label.y = 0.4)  # Add correlation coefficient
e

#by status
TotalMH$Status = as.factor(TotalMH$Status)
f = ggscatter(TotalMH, x = "taxa.average", y = "MHRel",
              add = "reg.line",                         # Add regression line
              conf.int = TRUE,                          # Add confidence interval
              color = "Status", palette = my_colors,           # Color by groups "cyl"
              shape = "Status"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Status), label.x = .45, label.y.npc = .1) # Add correlation coefficient
f

#Histophilus
str(TotalHS)
cor(TotalHS$taxa.average, TotalHS$HSRel)
cor.test(TotalHS$taxa.average, TotalHS$HSRel) #significant 

g = ggscatter(TotalHS, x = "taxa.average", y = "HSRel",
              add = "reg.line",                                 # Add regression line
              conf.int = TRUE,                                  # Add confidence interval
              add.params = list(color = "blue",
                                fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 0, label.y = 4)  # Add correlation coefficient
g

#by status
TotalHS$Status = as.factor(TotalHS$Status)
h = ggscatter(TotalHS, x = "taxa.average", y = "HSRel",
              add = "reg.line",                         # Add regression line
              conf.int = TRUE,                          # Add confidence interval
              color = "Status", palette = my_colors,           # Color by groups "cyl"
              shape = "Status"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Status), label.x = .28, label.y.npc = .1) 
h

ggarrange(a,e,g,b,f,h,ncol =3, nrow = 2, 
          font.label = list(size = 12))
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



