## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")

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

##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/") 

metadata <- read.csv("BRDBigMetadataDairyTotal.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
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
#[ 26496 taxa and 498 samples ] beef

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
#[ 26365 taxa and 498 samples ] beef

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

#adding all the ASV counts for each ASV
NewASVTable2 = t(NewASVTable)
total = as.data.frame(rowSums(NewASVTable2))
total = subset(total, `rowSums(NewASVTable2)` >= 50)
NewASVTable3 <- t(NewASVTable2[rownames(NewASVTable2) %in% rownames(total),])
#write.csv(NewASVTable3, "DairyNewASVTableFiltered.csv")

NewTax2 <- NewTax[rownames(NewTax) %in% rownames(total),] 
write.csv(NewTax2, "DairyNewTaxFiltered.txt")

#################################################################
##Taxa barplot
################################################################
#If you want different taxonomic level, find and replace the taxonomic level listed here
####subset by State to see the difference of BRD vs healthy
healthy <- subset(metadata, Status=="Healthy")
dplyr::count(healthy, ID)
BRD <- subset(metadata, Status=="BRD")
dplyr::count(BRD, ID)

H_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(healthy),] 
B_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(BRD),] 

BRD <- BRD[rownames(BRD) %in% rownames(B_OTU),]
healthy <- healthy[rownames(healthy) %in% rownames(H_OTU),]

##start with the healthy samples
#phyloseq with all the data
OTU.physeq = otu_table(as.matrix(NewASVTable), taxa_are_rows=FALSE)
tax.physeq= tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)
physeq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq
saveRDS(physeq, 'physeqDairy')

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
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot = physeqH1

# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) # this palette you use for the dairy healthy samples (Genus, family and Phylum )
my_level <- c("Species")
my_column <- "State"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$State = factor(physeq_meta_filtered$State, c("TX", "IN", "NY", "CA"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample-Healthy Dairy")) 
  ggsave(paste0("", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredSpeciesHealthyDairy.csv")


## now we do the same thing for the BRD animals
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

# Set colors for plotting
my_colorsD <- c(
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c','#fb9a99','powderblue','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "gray","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for genus

my_colorsD1 <- c(
  '#a6cee3',"gray",'#1f78b4','#b2df8a','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588","black", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for famuly
my_colorsD2 <- c(  
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                   '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
                   "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                   "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)#for phylum

my_colorsD3 <- c(  
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c','#fb9a99','#e31a1c', "#599861", "gray",
  '#fdbf6f',"darkslategrey",'#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285",'deepskyblue2', "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray"
)#for species

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Species")
my_column <- "State"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotD %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$State = factor(physeq_meta_filtered$State, c("TX", "IN", "NY", "CA"))
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredSpeciesBRDDairy.csv")

##now let's find the abundance of the BRD pathogens in the healthy and sick animals for each statr
CAT = subset(metadata, State == "CA")
NYT = subset(metadata, State == "NY")
INT = subset(metadata, State == "IN")
TXT = subset(metadata, State == "TX")

CA_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(CAT),] 
CAT <- CAT[rownames(CAT) %in% rownames(CA_OTU),]
CAT$Status = as.factor(CAT$Status)
levels(CAT$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

NY_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(NYT),] 
NYT <- NYT[rownames(NYT) %in% rownames(NY_OTU),]
NYT$Status = as.factor(NYT$Status)
levels(NYT$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

IN_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(INT),] 
INT <- INT[rownames(INT) %in% rownames(IN_OTU),]
INT$Status = as.factor(INT$Status)
levels(INT$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

TX_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(TXT),] 
TXT <- TXT[rownames(TXT) %in% rownames(TX_OTU),]
TXT$Status = as.factor(TXT$Status)
levels(TXT$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

#### ---- identify the taxonomy composition within each state
str(CAT) #modify every time you test a different state
OTU.physeqS = otu_table(as.matrix(CA_OTU), taxa_are_rows=FALSE)
tax.physeqS = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqS = sample_data(CAT) #modify every time you test a different state
physeq_S = phyloseq(OTU.physeqS, tax.physeqS, meta.physeqS)

##Making the tree, we need the phyloseq object 
random_treeS = rtree(ntaxa(physeq_S), rooted=TRUE, tip.label=taxa_names(physeq_S))
##Merging the tree with the phyloseq object
physeqS1 = merge_phyloseq(physeq_S, meta.physeqS, random_treeS)
physeqS1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plotS = physeqS1

#If you want different taxonomic level, find and replace the taxonomic level listed here
#color palette
my_colorsCA <- c(
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c',
  '#fdbf6f','#cab2d6', "#CBD588", "#5F7FC7", "orange","#DA5724", "gray","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") 

my_colorsIN <- c(
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c','#fb9a99','powderblue','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "gray","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for genus

my_colorsNY <- c(
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c','#fb9a99','powderblue','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "gray","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for genus

my_colorsTX <- c(
  '#a6cee3','#1f78b4',"black",'#b2df8a','#33a02c','#fb9a99','powderblue','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "gray","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for genus
my_level <- c("Genus")
my_column <- "Status"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotS %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$State = factor(physeq_meta_filtered$Status, c("Healthy", "BRD"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    #scale_fill_manual(values = my_colorsCA) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%)")) 
  ggsave(paste0("", ml, "BarPlotCADairy_", my_column, ".png"), height = 5, width = 4)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredGenusCADairy.csv")

### Identify relative abundance of the BRD-pathobionts within each state
otu.summary <- prop.table((TX_OTU), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
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
str(metadata)
str(melt_otu)
meta_otu <- merge(TXT, melt_otu, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- metadata$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (3837 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
Hs = subset(NewTax, Genus == " Histophilus")
Pm = subset(NewTax, Genus == " Pasteurella")
Mh = subset(NewTax, Genus == " Mannheimia")
Mb = subset(NewTax, Genus == " Mycoplasma")
Bi = subset(NewTax, Genus == " Bibersteinia")
Tru = subset(NewTax, Genus == " Trueperella")
patho = rbind(Hs, Pm, Mh, Mb, Bi, Tru)

hs = subset(patho, Species == " Histophilus")
pm = subset(patho, Species == " Pasteurella_multocida")
mh = subset(patho, Species == " Mannheimia_haemolytica")
mb = subset(patho, Species == " Mycoplasma")
patho2 = rbind(hs, pm, mh, mb)

str(melt_otu)
meta_otu <- merge(TXT, melt_otu, by.x = 0, by.y = "Sample") # change the metadata to whatever state you are working with
str(meta_otu)

meta_otu_tax <- merge(meta_otu, patho, by.x = "ASV", by.y = 0)
str(meta_otu_tax)

meta_otu_tax2 <- merge(meta_otu, patho2, by.x = "ASV", by.y = 0) #4 BRD pathogens
str(meta_otu_tax2)

Species <- meta_otu_tax %>% 
  group_by(Row.names,  Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
Species$State = "TX"
Species = Species[-c(45,46),]
#SpeciesS = Species
SpeciesS = rbind(SpeciesS, Species)
write.csv(SpeciesS, "BRDgeneraDairyState.csv")

Species2 <- meta_otu_tax2 %>% 
  group_by(Row.names, Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance 4 bacteria 
Species2$State = "TX"
#Species2S = Species2
Species2S = rbind(Species2S, Species2)
write.csv(Species2S, "BRDpathogensDairyState.csv")

my_colors2 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) ### for dairy samples

SpeciesS$Status = as.factor(SpeciesS$Status)
levels(SpeciesS$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

ggplot(SpeciesS, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(State~.) +
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.20)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Dairy samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

Species2S$Status = as.factor(Species2S$Status)
levels(Species2S$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
ggplot(Species2S, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(State~.) +
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.15)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Dairy samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

#####------------------ Create new plots with the taxa with highest relative abundance
# we need the taxa information for the BRD and ehalthy animals
metadata <- read.csv("./Taxonomy/physeq_meta_filteredGenusBRDDairy.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)

ggplot(metadata, aes(x = State, y =Abundance.average, fill =Genus)) + 
  geom_bar(stat = "identity", width= 0.5) + scale_fill_manual(values = my_colorsD) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.25, by = 0.15), limits = c(0, 0.25)) +
  facet_wrap(Genus ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=10),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 9)) 
  ggtitle("Genus-BRD samples")

metadata2 <- read.csv("./Taxonomy/physeq_meta_filteredFamilyBRDDairy.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadata2, aes(x = State, y =Abundance.average, fill =Family)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD1) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.50, by = 0.25), limits = c(0, 0.50)) +
  facet_wrap(Family ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=10),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Family-BRD samples")

metadata3 <- read.csv("./Taxonomy/physeq_meta_filteredPhylumBRDDairy.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadata3, aes(x = State, y =Abundance.average, fill =Phylum)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD1) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.75, by = 0.25), limits = c(0, 0.75)) +
  facet_wrap(Phylum ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Phylum-BRD samples")

# Healthy samples
metadataH <- read.csv("./Taxonomy/physeq_meta_filteredGenusHealthyDairy.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)

ggplot(metadataH, aes(x = State, y =Abundance.average, fill =Genus)) + 
  geom_bar(stat = "identity", width=0.5) + scale_fill_manual(values = my_colors) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.30, by = 0.15), limits = c(0, 0.30)) +
  facet_wrap(Genus ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=10),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 9))
  ggtitle("Genus-Healthy samples")

metadataH2 <- read.csv("./Taxonomy/physeq_meta_filteredFamilyHealthyDairy.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadataH2, aes(x = State, y =Abundance.average, fill =Family)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colors) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.50, by = 0.25), limits = c(0, 0.50)) +
  facet_wrap(Family ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Family-Healthy samples")

metadataH3 <- read.csv("./Taxonomy/physeq_meta_filteredPhylumHealthyDairy.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadataH3, aes(x = State, y =Abundance.average, fill =Phylum)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colors) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.75, by = 0.25), limits = c(0, 0.75)) +
  facet_wrap(Phylum ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Phylum-Healthy samples")

#----- compare genus within each state
CAH = subset(metadataH, State == "CA")
CAH$Status = "Healthy"
CAB = subset(metadata, State == "CA")
CAB$Status = "BRD"
CA = rbind(CAH, CAB)
NYH = subset(metadataH, State == "NY")
NYH$Status = "Healthy"
NYB = subset(metadata, State == "NY")
NYB$Status = "BRD"
NY = rbind(NYH, NYB)
INH = subset(metadataH, State == "IN")
INH$Status = "Healthy"
INB = subset(metadata, State == "IN")
INB$Status = "BRD"
IN = rbind(INH, INB)
TXH = subset(metadataH, State == "TX")
TXH$Status = "Healthy"
TXB = subset(metadata, State == "TX")
TXB$Status = "BRD"
TX = rbind(TXH, TXB)

TX$Status = as.factor(TX$Status)
levels(TX$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
ggplot(TX, aes(x = Status, y = Abundance.average, fill = Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+ 
  #facet_grid(Genus~. ) +
  scale_fill_manual(values = my_colors) +
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
  ylab(paste0("Relative Abundance")) +  labs(x='TX ~ Status')


##############################################
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef/") 
###beef samples
metadata <- read.csv("BRDBigMetadataBeef.csv", na.strings = c("","NA"), header=TRUE)
row.names(metadata) <- metadata$ID
order_groups <- metadata$ID
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
ASV_table <- ASV_table[-c(1:7)]
ASV_table <- as.matrix(ASV_table)

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
meta.physeq = sample_data(metadata)

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
physeq_deseq #[ 26365 taxa and 498 samples ]

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

#adding all the ASV counts for each ASV
NewASVTable2 = t(NewASVTable)
total = as.data.frame(rowSums(NewASVTable2))
total = subset(total, `rowSums(NewASVTable2)` >= 50)
NewASVTable3 <- t(NewASVTable2[rownames(NewASVTable2) %in% rownames(total),])
write.csv(NewASVTable3, "BeefNewASVTableFiltered.csv")

NewTax2 <- NewTax[rownames(NewTax) %in% rownames(total),] 
write.csv(NewTax2, "BeefNewTaxFiltered.txt")

#################################################################
##Taxa barplot
################################################################
#If you want different taxonomic level, find and replace the taxonomic level listed here
####subset by healthy or BRD
healthy <- subset(metadata, Status=="Healthy")
dplyr::count(healthy, ID)
BRD <- subset(metadata, Status=="BRD")
dplyr::count(BRD, ID)

H_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(healthy),] 
B_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(BRD),] 

BRD <- BRD[rownames(BRD) %in% rownames(B_OTU),]
healthy <- healthy[rownames(healthy) %in% rownames(H_OTU),]

##start with the healthy samples
#phyloseq with all the data
OTU.physeq = otu_table(as.matrix(NewASVTable), taxa_are_rows=FALSE)
tax.physeq= tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)
physeq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq
saveRDS(physeq, 'physeqDairy')

OTU.physeqH = otu_table(as.matrix(H_OTU), taxa_are_rows=FALSE)
tax.physeqH = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqH = sample_data(healthy)
physeq_H = phyloseq(OTU.physeqH, tax.physeqH, meta.physeqH)

##Making the tree, we need the phyloseq object 
library("ape") #to create the tree
random_tree = rtree(ntaxa(physeq_H), rooted=TRUE, tip.label=taxa_names(physeq_H))
##Merging the tree with the phyloseq object
physeqH1 = merge_phyloseq(physeq_H, meta.physeq, random_tree)
physeqH1 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot = physeqH1

# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)
my_level <- c("Species")
my_column <- "State"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$State = factor(physeq_meta_filtered$State, c("TX", "IN", "ID", "CO"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample-Healthy Beef")) 
  ggsave(paste0("", ml, "BarPlotBeef_", my_column, ".png"), height = 5, width = 4)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredSpeciesHealthyBeef.csv")


## now we do the same thing for the BRD animals
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

# Set colors for plotting
my_colorsD <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','powderblue',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724","#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black")

my_colorsD1 <- c(
  '#a6cee3','#1f78b4','#b2df8a','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black") #for famuly
my_colorsD2 <- c(  
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)#for phylum

my_colorsD3 <- c(  
  '#a6cee3','#1f78b4','#b2df8a','#fb9a99','powderblue','#fdbf6f','#ff7f00',
  '#cab2d6','#6a3d9a','#ffff99','#b15928',"gray", "#CBD588", "#5F7FC7", "#DA5724",
  "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)#for species

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Species")
my_column <- "State"  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.02  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plotD %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  physeq_meta_filtered$State = factor(physeq_meta_filtered$State, c("TX", "IN", "ID", "CO"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colorsD3) +
    theme_bw()+
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size=8, face = "italic", vjust = 1.5)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample-BRD Beef")) 
  ggsave(paste0("", ml, "BarPlotBRDBeef_", my_column, ".png"), height = 5, width = 4)
}

write.csv(physeq_meta_filtered, "physeq_meta_filteredSpeciesBRDBeef.csv")

##now let's find the abundance of the BRD pathogens
#calculating the relative abundance of the BRD pathogens
otu.summary <- prop.table((NewASVTable), 1) #accounts for the relative abundance in each sample
str(otu.summary)
otu_abund <- colSums(otu.summary) ##the abundance of each ASV across all samples
a <- as.data.frame(otu_abund)
sum(otu_abund)
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
str(metadata)
str(melt_otu)
meta_otu <- merge(metadata, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
order_groups <- metadata$ID
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (3837 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

#merging the abundance of each OTU with the metadata and the taxonomy file
#only include the pathogens species
Hs = subset(NewTax, Genus == " Histophilus")
Pm = subset(NewTax, Genus == " Pasteurella")
Mh = subset(NewTax, Genus == " Mannheimia")
Mb = subset(NewTax, Genus == " Mycoplasma")
Bi = subset(NewTax, Genus == " Bibersteinia")
patho = rbind(Hs, Pm, Mh, Mb, Bi)

hs = subset(patho, Species == " Histophilus")
pm = subset(patho, Species == " Pasteurella_multocida")
mh = subset(patho, Species == " Mannheimia_haemolytica")
mb = subset(patho, Species == " Mycoplasma") #no mycoplasma bovis

patho2 = rbind(hs, pm, mh, mb)

metadata2 <- metadata[rownames(metadata) %in% rownames(NewASVTable),] 
levels(metadata2$Status) <- list("Healthy"="Healthy", "BRD"="BRD")

str(melt_otu)
meta_otu <- merge(metadata2, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, patho, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
meta_otu_tax2 <- merge(meta_otu, patho2, by.x = "ASV", by.y = 0) #for BRD pathogens
str(meta_otu_tax)

Species <- meta_otu_tax %>% 
  group_by(Row.names, State, Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't
write.csv(Species, "BeefBRDpathogensAllabundance.csv")

Species2 <- meta_otu_tax2 %>% 
  group_by(Row.names, State, Species, Status) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(State, Species, Status) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance, doesn't
write.csv(Species2, "BeefBRD4pathogensabundance.csv")


my_colors3 <- c(
  '#a6cee3','#1f78b4',"lightpink",'#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","lightgreen","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","gray","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

A <- ggplot(Species, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(State~., nrow=1, ncol=4) +
  scale_fill_manual(values = my_colors3) +
  ylim(c(0,0.25)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Beef samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

B = ggplot(Species2, aes(x = Status, y = taxa.average, fill =Species)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(State~., nrow=1, ncol=4) +
  scale_fill_manual(values = my_colors2) +
  ylim(c(0,0.1)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .8, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Beef samples") +
  ylab(paste0("Relative Abundance")) +  labs(x='Disease Status')

ggarrange(A,B, labels = c("A", "B"),
          ncol = 1, nrow=2, font.label = list(size = 16))

### Taxonomy plots 
metadata <- read.csv("./Taxonomy/physeq_meta_filteredGenusBRDBeef.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)

ggplot(metadata, aes(x = State, y =Abundance.average, fill =Genus)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.26, by = 0.15), limits = c(0, 0.26)) +
  facet_wrap(Genus ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Genus-BRD samples")

metadata2 <- read.csv("./Taxonomy/physeq_meta_filteredFamilyBRDBeef.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadata2, aes(x = State, y =Abundance.average, fill =Family)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD1) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.50, by = 0.25), limits = c(0, 0.50)) +
  facet_wrap(Family ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Family-BRD samples")

metadata3 <- read.csv("./Taxonomy/physeq_meta_filteredPhylumBRDBeef.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadata3, aes(x = State, y =Abundance.average, fill =Phylum)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.65, by = 0.25), limits = c(0, 0.65)) +
  facet_wrap(Phylum ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Phylum-BRD samples")

# Healthy samples
metadataH <- read.csv("./Taxonomy/physeq_meta_filteredGenusHealthyBeef.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)

ggplot(metadataH, aes(x = State, y =Abundance.average, fill =Genus)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colors) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.30, by = 0.15), limits = c(0, 0.30)) +
  facet_wrap(Genus ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Genus-Healthy samples")

metadataH2 <- read.csv("./Taxonomy/physeq_meta_filteredFamilyHealthyBeef.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadataH2, aes(x = State, y =Abundance.average, fill =Family)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colors) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.45, by = 0.25), limits = c(0, 0.45)) +
  facet_wrap(Family ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Family-Healthy samples")

metadataH3 <- read.csv("./Taxonomy/physeq_meta_filteredPhylumHealthyBeef.csv", na.strings = c("","NA"), header=TRUE)

ggplot(metadataH3, aes(x = State, y =Abundance.average, fill =Phylum)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = my_colorsD) +
  theme_bw()+ scale_y_continuous(breaks = seq(0, 0.65, by = 0.25), limits = c(0, 0.65)) +
  facet_wrap(Phylum ~ ., strip.position = "right", ncol = 1) + 
  theme(strip.text.y = element_text(angle = 0, face="italic", size=11),
        legend.position = "none") +theme(axis.title.x = element_blank()) +
  theme( axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 9)) +
  ggtitle("Phylum-Healthy samples")



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
