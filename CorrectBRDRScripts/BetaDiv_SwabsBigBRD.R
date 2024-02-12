BETA diversity Swabs
#install.packages("raster")
library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
#transpose function


setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeDairy/CorrectQiimeDairy/") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata <- read.csv("BRDBigMetadataDairy.csv", na.strings = c("","NA"), header=TRUE)

str(metadata)
metadata$Status <- as.factor(metadata$Status)
metadata$State <- as.factor(metadata$State)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata <- metadata[-c(8:12)]
rownames(metadata) = metadata$ID

ASVs <- read_qza("rarefied_table.qza") #10809 Dairy samples, 368 Diary samples
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
ASV_table <- merge(metadata, ASV_table, by.x = "Name", by.y = 0)
row.names(ASV_table) <- ASV_table$ID
ASV_table <- ASV_table[,-c(1:8)]
ASV_table <- as.matrix(ASV_table)
x = as.data.frame(t(ASV_table[c(1,2),]))
sum(x$Dairy100) #12,312 reads per sample

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
#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 10809 taxa and 368 samples ] dairy

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
physeq_deseq #[ 10759 taxa and 368 samples ] Dairy

prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.bray <- phyloseq::distance(physeq_deseq, method = "bray")
dist.bray <- as.dist(dist.bray)

## PERMANOVA
#Bray curtis
BC <- adonis2(dist.bray ~ metadata$Status, strata=metadata$State, permutations = 999)
BC

BC1 <- adonis2(dist.bray ~ metadata$State, permutations = 999)
BC1 #only state significant

BC2 <- adonis2(dist.bray~ metadata$Status, permutations = 999)
BC2

TimeBC <- pairwise.adonis2(dist.bray ~ State,  data=metadata)
TimeBC##everything is different to everything

## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 
rm(OTU_Unifrac)
OTU.physeq = otu_table(as.matrix(NewASVTable), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(NewTax))

#We then merge these into an object of class phyloseq.
physeq = phyloseq(OTU.physeq, tax.physeq)
physeq

##Making the tree
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
##Merging the tree with the phyloseq object
physeq1 = merge_phyloseq(physeq, meta.physeq, random_tree)
physeq1 

OTU_Unifrac <- UniFrac(physeq1, weighted = TRUE) ## Calculating the weighted unifrac distances

## Permanova 
M3 <- adonis2(OTU_Unifrac~ metadata$Status,strata = metadata$State, permutations = 999)
M3

M31 <- adonis2(OTU_Unifrac~ metadata$Status, permutations = 999)
M31 

M31 <- adonis2(OTU_Unifrac~ metadata$State, permutations = 999)
M31 # State is significant 

TimeWU <- pairwise.adonis2(OTU_Unifrac ~ State,  data=metadata)
TimeWU

## Dispersion test 
BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$Status)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
boxplot(BRD_Wu) #dispersion weighted UniFrac significant

#install.packages("usedist")
library(usedist)

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# Significant

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$Status)
BRD_BC
boxplot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# significant 

##O Make the plots with ellipses
## Calculating the centroids of the data based on Bray-curtis and weighted unifrac
#D33 centroids for state
my_color <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

## Weighted unifrac 
ordu.wt.uni <- ordinate(physeq1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, 
                              ordu.wt.uni, color="State") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["ID"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["State"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")

wt.unifrac <- wt.unifrac + ggtitle("Dairy Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac + stat_ellipse())
a = print(wt.unifrac + stat_ellipse() +
        scale_color_manual(values = my_color) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
a
##other plot
str(wuaxis)

centroids <- aggregate(wuaxis[,4:5], list(Group=wuaxis$wuaxis4), mean)
colnames(centroids) <- c('BRD','groupX', 'groupY')


#distance betweeen centroids
pointDistance(c(0.03168905, -0.0010942445), c(-0.02609687, 0.0009011425), lonlat=FALSE)

## Bray_curtis
ordu.bc <- ordinate(physeq1, "PCoA", "bray")
Bray <- plot_ordination(physeq1, 
                              ordu.bc, color="State") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["State"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

Bray <- Bray + ggtitle("Dairy Bray Curtis") + geom_point(size = 2)
Bray <- Bray + theme_classic()
b = print(Bray + stat_ellipse() +
        theme(legend.text = element_text(size=12)) +
        scale_color_manual(values = my_color) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
b

str(bray)
ggarrange(b,a,ncol =2, labels = c("A", "B"),
          font.label = list(size = 11))


### separate by state 
BrayB <- plot_ordination(physeq1, 
                        ordu.bc, color="Status") 
Bray1 <- BrayB[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- BrayB[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- BrayB[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- BrayB[["data"]][["Status"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

my_colors <- c("dodgerblue3","goldenrod3")
BrayB <- BrayB + ggtitle("Bray Curtis") + geom_point(size = 2)
BrayB <- BrayB + theme_classic()
print(BrayB + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Dairy Farms: Bray-Curtis dissimilarity")
)

#Separate by State- Dairy WU
WUB <- plot_ordination(physeq1, 
                         ordu.wt.uni, color="Status") 
WU1 <- WUB[["data"]][["Axis.1"]]
WU1 <- as.data.frame(WU1)
WU1$number <- rownames(WU1)
WU2 <- WUB[["data"]][["Axis.2"]]
WU2 <- as.data.frame(WU2)
WU2$number <- rownames(WU2)
WU3 <- WUB[["data"]][["ID"]]
WU3 <- as.data.frame(WU3)
WU3$number <- rownames(WU3)
WU4 <- WUB[["data"]][["Status"]]
WU4 <- as.data.frame(WU4)
WU4$number <- rownames(WU4)

WU <- merge(WU3, WU4, by.x = "number", by.y = "number")
WU <- merge(WU, WU1, by.x = "number", by.y = "number")
WU <- merge(WU, WU2, by.x = "number", by.y = "number")

my_colors <- c("dodgerblue3","goldenrod3")
WUB <- WUB + ggtitle("Weigthed UniFrac") + geom_point(size = 2)
WUB <- WUB + theme_classic()
print(WUB + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Dairy Farms: Weighted UniFrac")
)

#### check the beta diversity removing the cow samples 
levels_to_remove <- c("Cow")

# Remove rows with specific levels in the 'Category' column
metadata$animal = as.factor(metadata$animal)
meta1 <- metadata[!(metadata$animal %in% levels_to_remove), ]
NCASV_table = t(NewASVTable)
NCASV_table<- NCASV_table[rownames(NCASV_table) %in% rownames(meta1),] 

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(NCASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(meta1)

#We then merge these into an object of class phyloseq.
physeq_deseq1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq1 # [ 10759 taxa and 345 samples ] dairy

colnames(tax_table(physeq_deseq1))

dist.brayC <- phyloseq::distance(physeq_deseq1, method = "bray")
dist.brayC <- as.dist(dist.brayC)

## PERMANOVA
#Bray curtis
BC.1 <- adonis2(dist.brayC ~ meta1$Status, strata=meta1$State, permutations = 999)
BC.1

BC1.1 <- adonis2(dist.brayC ~ meta1$State, permutations = 999)
BC1.1 #only state significant

BC2.1 <- adonis2(dist.brayC~ meta1$Status, permutations = 999)
BC2.1

TimeBC.1 <- pairwise.adonis2(dist.brayC ~ State,  data=meta1)
TimeBC.1##everything is different to everything

## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 
rm(OTU_Unifrac1)
NCASV_table1 = t(NCASV_table)
OTU.physeq1 = otu_table(as.matrix(NCASV_table1), taxa_are_rows=TRUE)
tax.physeq1 = tax_table(as.matrix(NewTax))

#We then merge these into an object of class phyloseq.
physeq1 = phyloseq(OTU.physeq1, tax.physeq1)
physeq1 # [ 10759 taxa and 345 samples ] no cow samples

##Making the tree
library("ape")
random_tree = rtree(ntaxa(physeq1), rooted=TRUE, tip.label=taxa_names(physeq1))
##Merging the tree with the phyloseq object
physeq1.1 = merge_phyloseq(physeq1, meta.physeq, random_tree)
physeq1.1 #[ 10759 taxa and 345 samples ]

OTU_Unifrac1 <- UniFrac(physeq1.1, weighted = TRUE) ## Calculating the weighted unifrac distances

## Permanova 
M <- adonis2(OTU_Unifrac1~ meta1$Status,strata = meta1$State, permutations = 999)
M

M3 <- adonis2(OTU_Unifrac1~ meta1$Status, permutations = 999)
M3 

M32 <- adonis2(OTU_Unifrac1~ meta1$State, permutations = 999)
M32 # State is significant 

TimeWU <- pairwise.adonis2(OTU_Unifrac1 ~ State,  data=meta1)
TimeWU

## Dispersion test 
BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$Status)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
boxplot(BRD_Wu) #dispersion weighted UniFrac significant

#install.packages("usedist")
library(usedist)

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# Significant

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$Status)
BRD_BC
boxplot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# significant 




##Beef samples
setwd("~/Desktop/eunice/PhD/AgSEED/BigProject/BRDQiime/QiimeBeef//") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata <- read.csv("BRDBigMetadataBeef.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)
metadata$Status <- as.factor(metadata$Status)
metadata$State <- as.factor(metadata$State)
levels(metadata$Status) <- list("Healthy"="Healthy", "BRD"="BRD")
#metadata <- metadata[-c(8:12)]
rownames(metadata) = metadata$ID

ASVs <- read_qza("rarefied_table.qza") #12051 Dairy samples
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #18010 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

#only include the samples that are present in the ASV_table
metadata <- metadata[rownames(metadata) %in% rownames(ASV_table),] #51, 3837

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
#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 20474 taxa and 469 samples ] beef

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
physeq_deseq #[ 20381 taxa and 469 samples ] beef

prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.bray <- phyloseq::distance(physeq_deseq, method = "bray")
dist.bray <- as.dist(dist.bray)

## PERMANOVA
#Bray curtis
BC <- adonis2(dist.bray ~ metadata$Status, strata=metadata$State, permutations = 999)
BC

BC1 <- adonis2(dist.bray ~ metadata$Status + metadata$State + metadata$Status*metadata$State, permutations = 999)
BC1 #only state significant

BC2 <- adonis2(dist.bray~ metadata$Status, permutations = 999)
BC2

BC3 <- adonis2(dist.bray~ metadata$State, permutations = 999)
BC3

TimeBC <- pairwise.adonis2(dist.bray ~ State,  data=metadata)
TimeBC

##everything is different to everything

## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 
rm(OTU_Unifrac)
OTU.physeq = otu_table(as.matrix(NewASVTable), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(NewTax))

#We then merge these into an object of class phyloseq.
physeq = phyloseq(OTU.physeq, tax.physeq)
physeq

##Making the tree
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
##Merging the tree with the phyloseq object
physeq1 = merge_phyloseq(physeq, meta.physeq, random_tree)
physeq1 

OTU_Unifrac <- UniFrac(physeq1, weighted = TRUE) ## Calculating the weighted unifrac distances

## Permanova 
M3 <- adonis2(OTU_Unifrac~ metadata$Status,strata = metadata$State, permutations = 999)
M3

M31 <- adonis2(OTU_Unifrac~ metadata$Status, permutations = 999)
M31 

M31 <- adonis2(OTU_Unifrac~ metadata$State, permutations = 999)
M31 # State is significant 

TimeWU <- pairwise.adonis2(OTU_Unifrac ~ State,  data=metadata)
TimeWU

## Dispersion test 
BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$Status)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
plot(BRD_Wu) #dispersion weighted UniFrac significant

#install.packages("usedist")
library(usedist)

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# 

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$Status)
BRD_BC
boxplot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# 
##O Make the plots with ellipses
## Calculating the centroids of the data based on Bray-curtis and weighted unifrac
#D33 centroids for treatment and season

## Weighted unifrac 
ordu.wt.uni <- ordinate(physeq1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, 
                              ordu.wt.uni, color="State") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["ID"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["State"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")

my_color1 <- c(
  "#653936", "#CD9BCD","orange","#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

wt.unifrac <- wt.unifrac + ggtitle("Beef Farms: Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac + stat_ellipse() +
        scale_color_manual(values = my_color1) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))

## Bray_curtis
ordu.bc <- ordinate(physeq1, "PCoA", "bray")
Bray <- plot_ordination(physeq1, 
                        ordu.bc, color="State") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["Status"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

Bray <- Bray + ggtitle("Beef Farms: Bray-Curtis Dissimilarity") + geom_point(size = 2)
Bray <- Bray + theme_classic()
print(Bray + stat_ellipse() +
        scale_color_manual(values = my_color1) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
)

str(bray)

### separate by state 
BrayB <- plot_ordination(physeq1, 
                         ordu.bc, color="Status") 
Bray1 <- BrayB[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- BrayB[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- BrayB[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- BrayB[["data"]][["Status"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

BrayB <- BrayB + ggtitle("Bray Curtis") + geom_point(size = 2)
BrayB <- BrayB + theme_classic()
print(BrayB + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Beef farms: Bray-Curtis dissimilarity")
)

### separate by state 
WUS <- plot_ordination(physeq1, 
                       ordu.wt.uni, color="Status") 
WUS1 <- WUS[["data"]][["Axis.1"]]
WUS1 <- as.data.frame(WUS1)
WUS1$number <- rownames(WUS1)
WUS2 <- WUS[["data"]][["Axis.2"]]
WUS2 <- as.data.frame(WUS2)
WUS2$number <- rownames(WUS2)
WUS3 <- WUS[["data"]][["ID"]]
WUS3 <- as.data.frame(WUS3)
WUS3$number <- rownames(WUS3)
WUS4 <- WUS[["data"]][["Status"]]
WUS4 <- as.data.frame(WUS4)
WUS4$number <- rownames(WUS4)

WU <- merge(WUS3, WUS4, by.x = "number", by.y = "number")
WU <- merge(WU, WUS1, by.x = "number", by.y = "number")
WU <- merge(WU, WUS2, by.x = "number", by.y = "number")

WUS <- WUS + ggtitle("Beef Farms: Weighted UniFrac") + geom_point(size = 2)
WUS <- WUS + theme_classic()
print(WUS + stat_ellipse())
print(WUS + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Beef Farms: Weighted UniFrac")
)
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

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 




