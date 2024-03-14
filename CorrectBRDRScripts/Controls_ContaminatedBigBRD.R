library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)
library(qiime2R)
library(phyloseq)

setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/DairyControls")

ASVs <- read_qza("tableAllControls.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #1244 ASVs

#####################################################################
######################################################################

##Adding taxonomy
#Taxonomy of each OTU
tax = read.table("taxonomyControls.tsv", header=TRUE, sep="\t")
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
#write.table(tax.final,"taxonomyMock2.txt",sep=",", row.names = FALSE) 


##Calculating abundance of the ASVs
ASV_table = t(ASV_table)
otu.summary <- prop.table(as.matrix(ASV_table), 1) #sample in rows and ASVs in columns
str(otu.summary)
otu_abund <- colSums(otu.summary)
a <- as.data.frame(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary <- otu.summary[-1,]
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:1244)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)
#write.table(melt_otu,"ASVAbundanceMock.txt",sep=",", row.names = FALSE)

meta_otu_taxC <- merge(melt_otu, tax.final, by.x = "ASV", by.y = 0)
str(meta_otu_taxC)
meta_otu_taxC$Family <- factor(meta_otu_taxC$Family)
meta_otu_taxC$Genus <- factor(meta_otu_taxC$Genus)
meta_otu_taxC$Phylum <- factor(meta_otu_taxC$Phylum)
str(meta_otu_taxC)
sum(meta_otu_taxC$Abundance)

genusC <- meta_otu_taxC %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/35)*100) ## the total number of samples (35)
attach(genusC)
sum(genusC$Ave_Abundance)
genusC <- genusC[order(-Ave_Abundance),] #pseudoalteromonas really abundant as well

## Combine the two taxa tables so see which ones are shared and unique
Cont_AVS = read.table("ContaminantsMockDairy.txt", header=TRUE, sep="\t")
rm(TaxCont_ASV)
TaxCont_ASV <- merge(Cont_AVS[-c(2)],  meta_otu_taxC, by.x = "Feature.ID", by.y = "ASV")


### Now getting the samples taxonomy, we will check how abundant are these ASVs in all the controls
# of some ASVs are mroe abundant than 10%, we will delete them from the sampels

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

str(TaxCont_ASV)
genusContMockControls <- TaxCont_ASV %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/35)*100) ## the total number of samples (51)
attach(genusContMockControls)
sum(genusContMockControls$Ave_Abundance)
genusContMockControls<- genusContMockControls[order(-Ave_Abundance),]

ggplot(TaxCont_ASV, aes(x = Sample, y =Abundance, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,1)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Controls')

## now we will check with the dairy samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/")
ASVs <- read_qza("table2.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #12072 ASVs
ASV_table <- t(ASV_table)

#str(ASV_table)

#Importing the metadata file
metadata <- read.csv("BRDBigMetadataDairyTotal.csv") #includes all the Dairy samples 405
str(metadata)
rownames(metadata) <- metadata$cattle
ASV_table <- merge(metadata, ASV_table, by.x = 0, by.y = 0) #405 samples
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[-c(1:8)]
ASV_table <- as.matrix(ASV_table)

#merging the abundance of each OTU with the metadata and the taxonomy file
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomy.qza")
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
tax3 = tax3[,-c(1:2,10:13)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

#---------------MIKING THE MODEL CON 
#generating a dataframe with all the response (days) from the samples
# Make one column for our outcome/response variable 
row.names(metadata) = metadata$ID
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.final))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 12072 taxa and 405 samples ]
colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "D_4__Mitochondria" & 
                Genus != "D_4__Mitochondria" &
                Species != "D_4__Mitochondria" &
                Order != "D_3__Chloroplast" &
                Family != "D_3__Chloroplast" &
                Genus != "D_3__Chloroplast" &
                Species != "D_3__Chloroplast") -> physeq_deseq
physeq_deseq ##[ 12012 taxa and 405 samples ]
## Random forest, we want to make the model so it classify samples based on the age
ntaxa(physeq_deseq) #total of 12012

physeq_deseq
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVtable <- prunetable
NewASVtable <- NewASVtable[,-c(2:9)]
row.names(NewASVtable) <- NewASVtable[,1]
NewASVtable = t(NewASVtable[,-c(1)]) #405 and 12012 

### CALCULATION OF THE ABUNDANCE OF EACH OTU  
## check the prop.table,
otu.summary <- prop.table(NewASVtable, 1) #accounts for the relative abundance in each sample 405 and 12012
otu.summary <- na.omit(otu.summary) #removin no data 402 12012
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary <- otu.summary[-1,]
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:12012)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

meta_otu <- merge(metadata, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- metadata$cattle
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

#combine the contaminates ASV with the samples
contSample <- merge(Cont_AVS[-c(2)], meta_otu_tax, by.x = "Feature.ID", by.y = "ASV")
contSample$Genus <- as.factor(contSample$Genus)
levels(contSample$Genus)
str(contSample)

ggplot(contSample, aes(x = cattle, y =Abundance, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,1)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance (Top 10)")) +  labs(x='Sample')

### Beef samples
##identify the contaminants in the mock community
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/ControlsBeef/")
tax = read.table("TaxonomyMockBeef.tsv", header=TRUE, sep="\t")
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
tax.MockBeef <- as.data.frame(tax3)
tax.MockBeef$OTUs <- rownames(tax.MockBeef)

#insert the correct ASV name that matched the Mock Community 
MockTax <- read.csv("MockTaxonomy.csv", na.strings = c("","NA"), header=TRUE)
rownames(MockTax) = MockTax$ASV
contMock = subset(tax.MockBeef,!(OTUs%in%MockTax$ASV)) #this are the ASV that are contaminants in the Mock

#Now we check the actual controls
ASVs <- read_qza("tableControlsBeef.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #1240 ASVs

#####################################################################
######################################################################

##Adding taxonomy
#Taxonomy of each OTU
tax = read.table("taxonomyControlBeef.tsv", header=TRUE, sep="\t")
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
#write.table(tax.final,"taxonomyMock2.txt",sep=",", row.names = FALSE) 


##Calculating abundance of the ASVs
ASV_table = t(ASV_table)
abun = as.data.frame(rowSums(ASV_table))
#samples Beef-B177, Beef-B189, and Extra329-421 abund of 0, you have to delete them
ASV_table = ASV_table[-c(2,4,17),]
abun = as.data.frame(rowSums(ASV_table)) # to double check
otu.summary <- prop.table(as.matrix(ASV_table), 1) #sample in rows and ASVs in columns
str(otu.summary)
otu_abund <- colSums(otu.summary)
a <- as.data.frame(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary <- otu.summary[-1,]
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:1240)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)
#write.table(melt_otu,"ASVAbundanceMock.txt",sep=",", row.names = FALSE)

meta_otu_taxC <- merge(melt_otu, tax.final, by.x = "ASV", by.y = 0)
str(meta_otu_taxC)
meta_otu_taxC$Family <- factor(meta_otu_taxC$Family)
meta_otu_taxC$Genus <- factor(meta_otu_taxC$Genus)
meta_otu_taxC$Phylum <- factor(meta_otu_taxC$Phylum)
str(meta_otu_taxC)
sum(meta_otu_taxC$Abundance)

genusC <- meta_otu_taxC %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/36)*100) ## the total number of samples (51)
attach(genusC)
sum(genusC$Ave_Abundance)
genusC <- genusC[order(-Ave_Abundance),] #pseudoalteromonas really abundant as well

## Combine the two taxa tables so see which ones are shared and unique
str(contMock)

TaxCont_ASV <- merge(contMock[-c(1:7)],  meta_otu_taxC, by.x = "OTUs", by.y = "ASV")
## No shared ASVs, no contamination between Mock and Empty Swabs

### Now getting the samples taxonomy, we will check how abundant are these ASVs in all the controls
# of some ASVs are mroe abundant than 10%, we will delete them from the sampels

my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

str(TaxCont_ASV)
genusContMockControls <- TaxCont_ASV %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/36)*100) ## the total number of samples (51)
attach(genusContMockControls)
sum(genusContMockControls$Ave_Abundance)
genusContMockControls<- genusContMockControls[order(-Ave_Abundance),]

#this is what is shared betweeen the mock and the controls
ggplot(TaxCont_ASV, aes(x = Sample, y =Abundance, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,1)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Controls')

#now let's plot the ASVs in the controls that are most abundant
cont = genusC[c(1:3),]
meta_otu_Cont = merge(meta_otu_taxC, cont , by.x = "Genus", by.y = "Genus")
str(meta_otu_Cont)

ggplot(meta_otu_Cont, aes(x = Sample, y =Abundance, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,1)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance")) +  labs(x='Controls')

## now we will check with the beef samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/")
ASVs <- read_qza("table2.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #25993 ASVs
ASV_table <- t(ASV_table)

#Importing the metadata file
metadata <- read.csv("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/CorrectBeefQiime/BRDBigqPCRBeefTotalNamesRow.csv") #includes all the Dairy samples 405
str(metadata)
rownames(metadata) <- metadata$ID
ASV_table <- merge(metadata, ASV_table, by.x = 0, by.y = 0) #504 samples
row.names(ASV_table) <- ASV_table$Row.names
ASV_table <- ASV_table[-c(1:23)]
ASV_table <- as.matrix(ASV_table)

#merging the abundance of each OTU with the metadata and the taxonomy file
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomy.qza")
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
tax3 = tax3[,-c(1:2,10:13)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

#---------------MIKING THE MODEL CON 
#generating a dataframe with all the response (days) from the samples
# Make one column for our outcome/response variable 
row.names(metadata) = metadata$ID
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.final))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 25993 taxa and 504 samples ]
colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "g__Mitochondria" & 
                Genus != "g__Mitochondria" &
                Species != "g__Mitochondria" &
                Order != "	o__Chloroplast" &
                Family != "	o__Chloroplast" &
                Genus != "	o__Chloroplast" &
                Species != "o__Chloroplast") -> physeq_deseq
physeq_deseq ##[ 25993 taxa and 504 samples ]
## Random forest, we want to make the model so it classify samples based on the age
ntaxa(physeq_deseq) #total of 12012

physeq_deseq
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVtable <- prunetable
NewASVtable <- NewASVtable[,-c(2:9)]
row.names(NewASVtable) <- NewASVtable[,1]
NewASVtable = t(NewASVtable[,-c(1)]) #405 and 12012 

### CALCULATION OF THE ABUNDANCE OF EACH OTU  
## check the prop.table,
otu.summary <- prop.table(NewASVtable, 1) #accounts for the relative abundance in each sample 405 and 12012
otu.summary <- na.omit(otu.summary) #removin na data 487 
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary <- otu.summary[-1,]
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:25993)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

meta_otu <- merge(metadata, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu$Row.names = as.factor(meta_otu$Row.names)
levels(meta_otu$Row.names)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- metadata$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

#combine the contaminates ASV with the samples
desired_levels <- meta_otu_taxC$ASV
contSample  <- meta_otu_tax[meta_otu_tax$ASV %in% desired_levels, ]
contSample$ASV = as.factor(contSample$ASV)
levels(contSample$ASV)
contSample$ID = as.factor(contSample$ID)
levels(contSample$ID)

str(TaxCont_ASV)
genusContSample <- contSample %>% 
  group_by(Genus) %>% 
  summarise(Ave_Abundance = (sum(Abundance)/487)*100) ## the total number of samples (487)
attach(genusContSample)
sum(genusContSample$Ave_Abundance)
genusContSample<- genusContSample[order(-Ave_Abundance),]

contaminants = genusContSample[c(1:5),]
desired_levels2 <- contaminants$Genus
meta_otu_Cont  <- contSample[contSample$Genus %in% desired_levels2, ]

ggplot(meta_otu_Cont, aes(x = ID, y =Abundance, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  ylim(c(0,1)) +
  coord_flip() +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=11)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance (Top 10)")) +  labs(x='Sample')


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

