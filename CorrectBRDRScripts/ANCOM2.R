###Ancom
#combine the output of Qiime with the ASV metadata
library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(data.table)
library(ggpubr)

##### subsetting data
##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/") 
###beef samples
metadata <- read.csv("BRDBigqPCRBeefTotalNamesRow.csv", na.strings = c("","NA"), header=TRUE)
metadata = metadata[,-c(8:21)]
row.names(metadata) <- metadata$ID
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 
metadata$State <- as.factor(metadata$State)
dplyr::count(metadata, State, Status) 

ASVs <- read_qza("table-no-pseudo.qza")
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
ASV_table <- merge(metadata, ASV_table, by.x = "ID", by.y = 0)
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:8)]
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
physeq_deseq # [ 25925 taxa and 487 samples  ] beef

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
physeq_deseq #[ 25789 taxa and 487 samples ]

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
NewTax2 = NewTax
NewTax2$ASV = row.names(NewTax2)

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable) #25789, 487
ASVTable = t(NewASVTable)

##prepare the data for the loop for all the samples
NewASVTable2 = as.data.frame(rowSums(ASVTable))
str(NewASVTable2)
NewASVTable2 = subset(NewASVTable2, rowSums(ASVTable) >= 50)
ASVTable2 <- ASVTable[rownames(ASVTable) %in% rownames(NewASVTable2),]
ASVTable = t(ASVTable2)
write.table(ASVTable2, file='B_OTUBeef.tsv', quote=FALSE, sep='\t')

meta2 = metadata[rownames(metadata) %in% rownames(ASVTable),]
write.csv(meta2, "BigBRDMetadataANCOM.csv")


IN= subset(metadata, State=="IN")
ID= subset(metadata, State=="ID") #change to ID for beef, NY for dairy
CO= subset(metadata, State=="CO") #change to CO for beef, CA for dairy
TX= subset(metadata, State=="TX")

IN_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(IN),])
ID_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(ID),])
CO_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(CO),])
TX_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(TX),])

TX_OTU2 = t(TX_OTU) 
TX = TX[rownames(TX) %in% rownames(TX_OTU2),]
write.csv(TX, "TXBeefMetadata.csv")
IN_OTU2 = t(IN_OTU) 
IN = IN[rownames(IN) %in% rownames(IN_OTU2),]
write.csv(IN, "INBeefMetadata.csv")
CO_OTU2 = t(CO_OTU) 
CO = CO[rownames(CO) %in% rownames(CO_OTU2),]
write.csv(CO, "COBeefMetadata.csv")
ID_OTU2 = t(ID_OTU) 
ID = ID[rownames(ID) %in% rownames(ID_OTU2),]
write.csv(ID, "IDBeefMetadata.csv")

##now we will filter and leave the ASVs with >50 counts 
ID_OTU2 = as.data.frame(rowSums(ID_OTU))
str(ID_OTU2)
ID_OTU2 = subset(ID_OTU2, rowSums(ID_OTU) >= 50)
ID_OTU <- ID_OTU[rownames(ID_OTU) %in% rownames(ID_OTU2),]

IN_OTU2 = as.data.frame(rowSums(IN_OTU))
str(IN_OTU2)
IN_OTU2 = subset(IN_OTU2, rowSums(IN_OTU) >= 50)
IN_OTU <- IN_OTU[rownames(IN_OTU) %in% rownames(IN_OTU2),]

CO_OTU2 = as.data.frame(rowSums(CO_OTU))
str(CO_OTU2)
CO_OTU2 = subset(CO_OTU2, rowSums(CO_OTU) >= 50)
CO_OTU <- CO_OTU[rownames(CO_OTU) %in% rownames(CO_OTU2),]

TX_OTU2 = as.data.frame(rowSums(TX_OTU))
str(TX_OTU2)
TX_OTU2 = subset(TX_OTU2, rowSums(TX_OTU) >= 50)
TX_OTU <- TX_OTU[rownames(TX_OTU) %in% rownames(TX_OTU2),]

write.table(IN_OTU, file='IN_OTUBeefMetadata.tsv', quote=FALSE, sep='\t')
write.table(ID_OTU, file='ID_OTUBeefMetadata.tsv', quote=FALSE, sep='\t')
write.table(CO_OTU, file='CO_OTUBeefMetadata.tsv', quote=FALSE, sep='\t')
write.table(TX_OTU, file='TX_OTUDairyMetadata.tsv', quote=FALSE, sep='\t')

write.csv(IN_OTU, "IN_OTUBeefMetadata.csv")
write.csv(ID_OTU, "ID_OTUBeefMetadata.csv")
write.csv(CO_OTU, "CO_OTUBeefMetadata.csv")
write.csv(TX_OTU, "TX_OTUBeefMetadata.csv")

### relative abundance of the significant ASVs ~ All the samples
otu.summary <- prop.table((ASVTable), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:5957)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancom <- melt_otu[melt_otu$ASV %in% c("ASV11264", "ASV24638", "ASV16423", "ASV1978", "ASV352", "ASV18163"), ]

str(melt_otu)
meta_otu <- merge(metadata, ancom, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "BeefAncomSpeciesRelAbun.csv")

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ###

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors2) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Differentially Abundant Taxa") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

### relative abundance of the significant ASVs ~ IN (ASV 1978)
IN_OTU = t(IN_OTU)
otu.summary <- prop.table((IN_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:4515)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancomIN <- melt_otu[melt_otu$ASV %in% c("ASV1978"), ]

str(melt_otu)
meta_otu <- merge(IN, ancomIN, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "BeefAncomSpeciesRelAbunIN.csv")

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ###

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  scale_fill_manual(values = my_colors2) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Differentially Abundant Taxa ~ IN") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

### relative abundance of the significant ASVs ~ TX (ASV 24638, ASV11264)
TX_OTU = t(TX_OTU)
otu.summary <- prop.table((TX_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:701)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancomTX <- melt_otu[melt_otu$ASV %in% c("ASV24638", "ASV11264"), ]

str(melt_otu)
meta_otu <- merge(TX, ancomTX, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "BeefAncomSpeciesRelAbunTX.csv")

my_colors3 <- c(
  '#1f78b4','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ###

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  scale_fill_manual(values = my_colors3) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Differentially Abundant Taxa ~ TX") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

### Dairy samples
##### subsetting data
##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/CorrectQiimeDairy/") 
###beef samples
metadata <- read.csv("BRDBigqPCRDairyTotalOneRow.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$Name
order_groups <- metadata$cattle
metadata$Status <- as.factor(metadata$Status)
dplyr::count(metadata, Status) 
metadata$State <- as.factor(metadata$State)
dplyr::count(metadata, State, Status) 

ASVs <- read_qza("table-no-pseudo.qza")
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
ASV_table <- merge(metadata, ASV_table, by.x = "Name", by.y = 0)
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:21)]
ASV_table <- as.matrix(ASV_table)

row.names(metadata) <- metadata$ID

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
physeq_deseq #[ 11992 taxa and 400 samples ]

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
NewTax2 = NewTax
NewTax2$ASV = row.names(NewTax2)

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable) #25789, 487
ASVTable = t(NewASVTable)

##prepare the data for the loop for all the samples
##prepare the data for the loop for all the samples
NewASVTable2 = as.data.frame(rowSums(ASVTable))
str(NewASVTable2)
NewASVTable2 = subset(NewASVTable2, rowSums(ASVTable) >= 50)
ASVTable2 <- ASVTable[rownames(ASVTable) %in% rownames(NewASVTable2),]
ASVTable = t(ASVTable2)

IN= subset(metadata, State=="IN")
NY= subset(metadata, State=="NY")
CA= subset(metadata, State=="CA")
TX= subset(metadata, State=="TX")

IN_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(IN),])
NY_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(NY),])
CA_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(CA),])
TX_OTU <- t(NewASVTable[rownames(NewASVTable) %in% rownames(TX),])

TX_OTU2 = t(TX_OTU) 
TX = TX[rownames(TX) %in% rownames(TX_OTU2),]
write.csv(TX, "TXDairyMetadata.csv")
IN_OTU2 = t(IN_OTU) 
IN = IN[rownames(IN) %in% rownames(IN_OTU2),]
write.csv(IN, "INDairyMetadata.csv")
CA_OTU2 = t(CA_OTU) 
CA = CA[rownames(CA) %in% rownames(CA_OTU2),]
write.csv(CA, "CADairyMetadata.csv")
NY_OTU2 = t(NY_OTU) 
NY = NY[rownames(NY) %in% rownames(NY_OTU2),]
write.csv(NY, "NYDairyMetadata.csv")

##now we will filter and leave the ASVs with >50 counts 
NY_OTU2 = as.data.frame(rowSums(NY_OTU))
str(NY_OTU2)
NY_OTU2 = subset(NY_OTU2, rowSums(NY_OTU) >= 50)
NY_OTU <- NY_OTU[rownames(NY_OTU) %in% rownames(NY_OTU2),]

IN_OTU2 = as.data.frame(rowSums(IN_OTU))
str(IN_OTU2)
IN_OTU2 = subset(IN_OTU2, rowSums(IN_OTU) >= 50)
IN_OTU <- IN_OTU[rownames(IN_OTU) %in% rownames(IN_OTU2),]

CA_OTU2 = as.data.frame(rowSums(CA_OTU))
str(CA_OTU2)
CA_OTU2 = subset(CA_OTU2, rowSums(CA_OTU) >= 50)
CA_OTU <- CA_OTU[rownames(CA_OTU) %in% rownames(CA_OTU2),]

TX_OTU2 = as.data.frame(rowSums(TX_OTU))
str(TX_OTU2)
TX_OTU2 = subset(TX_OTU2, rowSums(TX_OTU) >= 50)
TX_OTU <- TX_OTU[rownames(TX_OTU) %in% rownames(TX_OTU2),]

write.table(IN_OTU, file='IN_OTUDDairyMetadata.tsv', quote=FALSE, sep='\t')
write.table(NY_OTU, file='NY_OTUDairyMetadata.tsv', quote=FALSE, sep='\t')
write.table(CA_OTU, file='CA_OTUDairyMetadata.tsv', quote=FALSE, sep='\t')
write.table(TX_OTU, file='TX_OTUDairyMetadata.tsv', quote=FALSE, sep='\t')

write.csv(IN_OTU, "IN_OTUBeefMetadata.csv")
write.csv(ID_OTU, "ID_OTUBeefMetadata.csv")
write.csv(CO_OTU, "CO_OTUBeefMetadata.csv")
write.csv(TX_OTU, "TX_OTUBeefMetadata.csv")

### relative abundance of the significant ASVs ~ All the samples
otu.summary <- prop.table((ASVTable), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:3461)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancom <- melt_otu[melt_otu$ASV %in% c("ASV1662"), ]

str(melt_otu)
meta_otu <- merge(metadata, ancom, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "DairyAncomSpeciesRelAbun.csv")

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ###

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors2) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Differentially Abundant Taxa") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

### relative abundance of the significant ASVs ~ NY (ASV 8265, ASV596)
NY_OTU = t(NY_OTU)
otu.summary <- prop.table((NY_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:1475)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancomNY <- melt_otu[melt_otu$ASV %in% c("ASV8265", "ASV596"), ]

str(melt_otu)
meta_otu <- merge(NY, ancomNY, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Genus, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Genus, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "DairyAncomSpeciesRelAbunNY.csv")

my_colors3 <- c(
  '#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ###

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  scale_fill_manual(values = my_colors3) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Differentially Abundant Taxa ~ NY") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

### relative abundance of the significant ASVs ~ TX (ASV 1662)
IN_OTU = t(IN_OTU)
otu.summary <- prop.table((IN_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples 487 samples and 5957
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:1464)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
melt_otu$ASV = as.factor(melt_otu$ASV)
ancomIN <- melt_otu[melt_otu$ASV %in% c("ASV1662"), ]

str(melt_otu)
meta_otu <- merge(IN, ancomIN, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax2, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Genus, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Genus, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
write.csv(Species, "DairyAncomSpeciesRelAbunIN.csv")

Species$Status = as.factor(Species$Status)
levels(Species$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(Species, aes(x = Status, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  scale_fill_manual(values = my_colors2) +
  #ylim(c(0,0.11)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 10)) +
  ggtitle("Differentially Abundant Taxa ~ IN") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')


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
