### Big BRD taxonomy

rm(list = ls ())

library(dplyr)
library(tidyr)
library(ape)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(plotly)
library(tidyr)
library(naniar)
library(zoo)
library(lubridate)
library(qiime2R)
library(ggpubr)
library(forcats)

#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/") 

metadata <- read.csv("BRDBigMetadataDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 


ASVs <- read_qza("table-no-pseudo.qza") #12051 Dairy samples
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #18010 ASVs
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

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[,c(1,10:411)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)


### Checking how the data looks
## Make a plot to see the community in the two groups
# this prunes the taxa with abundance <2%
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
otu.summary <- prop.table(as.matrix(NewASVTable), 1) #relative abundance based on each row
str(otu.summary)
otu_abund <- colSums(otu.summary)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:11992)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
order_groups
str(metadata)
str(melt_otu)
meta_otu <- merge(metadata, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
levels(meta_otu_tax$Status)
str(metadata)
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Status <- factor(meta_otu_tax$Status)
#levels(meta_otu_tax$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)
 

####subset by healthy or BRD
healthy <- subset(meta_otu_tax, Status=="Healthy")
dplyr::count(healthy, ID)
BRD <- subset(meta_otu_tax, Status=="BRD")
dplyr::count(BRD, ID)

###calculate community abundance for each of the two groups
### Abundance at a phlyum level
str(healthy)
genusH <- healthy %>% 
  group_by(Species) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/201)) ## the total number of samples (129)
attach(genusH)
sum(genusH$Ave_Abundance)
genusH <- genusH[order(-Ave_Abundance),]
write.csv(genusH, "speciesHealthyDairy.csv")

top10genus <-genusH[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genus$Ave_Abundance) ### the total of the community composed of the top 10
top10G <- merge(healthy, top10genus, by.x = "Species", by.y="Species")
top10G$Genus <- as.factor(top10G$Genus)
levels(top10G$Genus)
str(top10G)

genusAB <- top10G %>% 
  group_by(Row.names, State, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(genusAB)
sum(genusAB$taxa.average)
genusAB$Species <- factor(genusAB$Species)
genusAB$Status <- c("Healthy")

#BRD Species
str(BRD)

genusHB<-BRD %>% 
  group_by(Species) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/201)) ## the total number of samples (129)
attach(genusHB)
sum(genusHB$Ave_Abundance)
genusHB <- genusHB[order(-Ave_Abundance),]
write.csv(genusHB, "speciesBRDDairy.csv")

top10genusB <-genusHB[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genusB$Ave_Abundance) ### the total of the community composed of the top 10
top10GB <- merge(BRD, top10genusB, by.x = "Species", by.y="Species")
top10GB$Species <- as.factor(top10GB$Species)
levels(top10GB$Genus)
str(top10GB)

genusABB <- top10GB %>% 
  group_by(Row.names, State, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(genusABB)
sum(genusABB$taxa.average)
genusABB$Species <- factor(genusABB$Species)
genusABB$Status <- c("BRD")

genus <- rbind(genusABB, genusAB)

my_colorsHG <- c(
  'wheat2', "#673770", "#5F7FC7", "orange","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)


ggplot(genus, aes(x = State, y = taxa.average, fill = Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  facet_grid(.~Status) +
  scale_fill_manual(values = my_colorsHG) +
  #coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance Top 15")) +  labs(x='State')


#------------------ Beef samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/") 

metadata2 <- read.csv("BRDBigMetadataBeef.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata2) <- metadata2$ID
order_groups2 <- metadata2$ID

ASVs <- read_qza("table-no-pseudo.qza") #12051 Dairy samples 26496
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) 
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

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
meta.physeq = sample_data(metadata2)

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
physeq_deseq #[ 26365 taxa and 498 samples ] beef

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[,c(1,10:507)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

### Checking how the data looks
## Make a plot to see the community in the two groups
# this prunes the taxa with abundance <2%
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
otu.summary <- prop.table(as.matrix(NewASVTable), 1) 
str(otu.summary)
otu_abund <- colSums(otu.summary)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:26365)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
order_groups2
str(metadata2)
str(melt_otu)
meta_otu <- merge(metadata2, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu[-c(2,3,6:8)], NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
levels(meta_otu_tax$Status)
str(metadata2)
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups2)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Status <- factor(meta_otu_tax$Status)
meta_otu_tax$ID <- factor(meta_otu_tax$ID)
#levels(meta_otu_tax$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

####subset by healthy or BRD
healthy <- subset(meta_otu_tax, Status=="Healthy") #261
dplyr::count(healthy, Row.names)
BRD <- subset(meta_otu_tax, Status=="BRD")
dplyr::count(BRD, Row.names) #237

###calculate community abundance for each of the two groups
### Abundance at a phlyum level
str(healthy)
genusH <- healthy %>% 
  group_by(Species) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/261)) ## the total number of samples (129)
attach(genusH)
sum(genusH$Ave_Abundance)
genusH <- genusH[order(-Ave_Abundance),]
write.csv(genusH, "speciesHealthyBeef.csv")

top10genus <-genusH[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genus$Ave_Abundance) ### the total of the community composed of the top 10
top10G <- merge(healthy, top10genus, by.x = "Species", by.y="Species")
top10G$Species <- as.factor(top10G$Species)
levels(top10G$Species)
str(top10G)

genusAB <- top10G %>% 
  group_by(Row.names, State, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(genusAB)
sum(genusAB$taxa.average)
genusAB$Species <- factor(genusAB$Species)
genusAB$Status <- c("Healthy")

#BRD family
str(BRD)

genusHB<-BRD %>% 
  group_by(Species) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/237)) ## the total number of samples (129)
attach(genusHB)
sum(genusHB$Ave_Abundance)
genusHB <- genusHB[order(-Ave_Abundance),]
write.csv(genusHB, "speciesBRDBeef.csv")

top10genusB <-genusHB[c(1:15),] #to select the top 10 most abundant phylum
sum(top10genusB$Ave_Abundance) ### the total of the community composed of the top 10
top10GB <- merge(BRD, top10genusB, by.x = "Species", by.y="Species")
top10GB$Species <- as.factor(top10GB$Species)
levels(top10GB$Species)
str(top10GB)

genusABB <- top10GB %>% 
  group_by(Row.names, State, Species) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't ass 
str(genusABB)
sum(genusABB$taxa.average)
genusABB$Species<- factor(genusABB$Species)
genusABB$Status <- c("BRD")

genus <- rbind(genusABB, genusAB)

my_colorsHG <- c(
  'wheat2', "#673770", "#5F7FC7", "orange","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)


ggplot(genus, aes(x = State, y = taxa.average, fill = Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  facet_grid(.~Status) +
  scale_fill_manual(values = my_colorsHG) +
  #coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = 1, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.text = element_text(size=11, face = "italic")) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance Top 15")) +  labs(x='State')


#Functions
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

