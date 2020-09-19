

library(PopGenome)
data_scaffold <- read.delim("/home/mrinalmishra/mrinalmishra1/PopGenome/60_genomes_test/scaffold_length_sorted.txt", header=FALSE, sep=" ")
## data_scaffold contains the length of each scaffold, so we can supply it to the to-from input of readVCF

i=1 ## Scaffold 1
for (i in 1:23){

scaffold=gsub("x", i, "BS.genome_hic_scaffold_x")
## Acquires the contig/scaffold name

GENOME.class <- readVCF("VCF/Final_filtered_All_60_SNPS.vcf.gz",numcols=10000, tid=scaffold, include.unknown=TRUE,
from=1, to= data_scaffold[i,2], approx=FALSE, out="", parallel=FALSE, gffpath="/home/mrinalmishra/PopGenome/GFF/Final_filtered_All_60_SNPS.mRNA.removed.gff")
## The 'removed.mRNA. file feeds into popgenome through this function well and all downstream operations work
## The original file can still be used to define gene boundaries to split the genome by genes

## Prep for MK test
## Ideally these modifications to 'GENOME.class' should be done before splitting them by gene positions or 10-kb windows
GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="60_genomes_test/BS_FINAL_from_draft.fasta")

###### Define populations
##### As with the MK test prep, perform this stpe before genes and sliding windows are defined
SCS_Allo <- c("Ref","XM30","XM31","XM44","XM45","XM46","XM47","XM48","XM50","XM51")
ECS <- c("WZ2","WZ5","WZ12","WZ16","WZ29","WZ30","WZ31","WZ45","WZ63","WZ64","YH51","YH91","YH93","YH95","ZS48","ZS64","ZS65","ZS66","ZS69","ZS71","ZS73","ZS100","ZS116")
SCS <- c("WZ9","WZ11","WZ13","WZ21","WZ22","WZ23","WZ25","WZ26","WZ39","WZ50","YH101","YH105","YH113","YH115","YH122","YH124","ZS128","ZS129","ZS131","ZS133","ZS135","ZS136","ZS138","ZS139","ZS140","ZS142")

GENOME.class <- set.populations(GENOME.class, list(SCS_Allo, ECS, SCS), diploid=TRUE)





### The mRNA file is useful because it allows us to split the genome by 'mRNA' - a proxy here for genes
### Whereas the mRNA removed gff file is somehow easier to feed into popgenome thereby facilitating all the downstream tests
### So the original GFF file can still be used to narrow down coordinates, whereas the modified one can be used to run tests
### It's a minor inconvenience
### All the calculations can now be done by 10kb windows or by gene

#### get_gff_info produces bugs
#### Have to produce an equivalent list for ourselves
#genePos <- get_gff_info(gff.file=gsub("x", i, "/home/mrinalmishra/mrinalmishra1/PopGenome/GFF/Final_filtered_All_60_SNPS.scaffold.x.gff") ,chr=scaffold, feature="mRNA")
################ This is the original gff file used here to define gene positions


#genes <- splitting.data(GENOME.class, positions=genePos, type=2)
###### type=2 referes to genomic positions
#print(length(genes@region.names))
#gene_name <- array()
#
#for (j in 1:length(genes@region.names)){
#gene_name[j] <- head(get_gff_info(genes, position=j, chr=scaffold, gff.file="GFF/Final_filtered_All_60_SNPS.gff")[[1]],1)
#}

genes.temp <- read.csv("N_scaffold/GC_N_content_genes.corrected.csv")
genes.temp1 <- subset(genes.temp, Scaffold_ID==scaffold)
colnames(genes.temp1) <- c("Gene_ID", "Scaffold_ID", "Start", "End", "GC_count", "GC_fraction", "N_count", "N_fraction")
genes.temp1$Gene_ID <- as.character(genes.temp1$Gene_ID)

library(tidyverse)
genePos <- list()
for (k in 1:nrow(genes.temp1)) {
genePos[[k]] <- genes.temp1$Start[k]:genes.temp1$End[k]
}
genes <- splitting.data(GENOME.class, positions=genePos, type=2)
print(length(genes@region.names))
print(head(genes@region.names))


#head(get_gff_info(genes, position=1, chr=scaffold, gff.file="GFF/Final_filtered_All_60_SNPS.gff")[[1]])
### Sample script to fish out annotation from an interesting region

print(GENOME.class@n.biallelic.sites)


# split the data in 10kb consecutive windows
slide <- sliding.window.transform(GENOME.class,10000,10000, type=2)
# total number of windows
print(length(slide@region.names))

######################################################### Nucleotide Diversity ############################################################
######################################################## Nucleotide Divergence ############################################################

## mode nucleotide does calculations per nucleotide as against per haplotype
slide_FsT <- F_ST.stats(slide, mode="nucleotide")
nucdiv.slide.within <- slide_FsT@nuc.diversity.within
nucdiv.slide.between <- slide_FsT@nuc.diversity.between
nucdiv.slide.between.formatted <- t(nucdiv.slide.between)
n.slide.sites <- slide@n.sites

genes_FsT <- F_ST.stats(genes, mode="nucleotide")
nucdiv.genes.within <- genes_FsT@nuc.diversity.within
nucdiv.genes.between <- genes_FsT@nuc.diversity.between
nucdiv.genes.between.formatted <- t(nucdiv.genes.between)
n.genes.sites <- genes@n.sites

Ns.temp <- read.csv(gsub("x", i, "N_scaffold/Scaffoldx_GC_without_N.csv"))
N_count <- head(Ns.temp[,4], -1) ## Removes the last row, as popgenome doesn't consider the last remaining <10kb chunk at the end
GC_fraction <- head(Ns.temp[,2], -1)
coordinates_slide <- as.character(head(Ns.temp[,1], -1))
n.slide.sites <- n.slide.sites - N_count



coordinates_genes <- genes@region.names
n.genes.sites <- n.genes.sites - genes.temp1$N_count


#save.session(GENOME.class,"GENOME.class")
#load.session("GENOME.class")

# the values have to be normalized by the number of nucleotides in each window
nucdiv.slide.within <- nucdiv.slide.within/n.slide.sites
nucdiv.slide.between.formatted <- nucdiv.slide.between.formatted/n.slide.sites

nucdiv.genes.within <- nucdiv.genes.within/n.genes.sites
nucdiv.genes.between.formatted <- nucdiv.genes.between.formatted/n.genes.sites

### Assigning column names facilitate printing directly onto a csv

colnames(nucdiv.slide.within) <- c("Diversity_SCS_Allopatry", "Diversity_ECS", "Diversity_SCS_Sympatry")
colnames(nucdiv.slide.between.formatted) <- c("Divergence SCS_Allopatry v ESC", "Divergence SCS_Allopatry v SCS_Sympatry", "Divergence ECS v SCS_Sympatry")

colnames(nucdiv.genes.within) <- c("Diversity_SCS_Allopatry", "Diversity_ECS", "Diversity_SCS_Sympatry")
colnames(nucdiv.genes.between.formatted) <- c("Divergence SCS_Allopatry v ESC", "Divergence SCS_Allopatry v SCS_Sympatry", "Divergence ECS v SCS_Sympatry")

## Smoothening the plots

## 'genes' may need a different loess smoothening span
nucdiv.slide.within[nucdiv.slide.within == "NaN"] <- NA
nucdiv.genes.within[nucdiv.genes.within == "NaN"] <- NA

nucdiv.slide.between.formatted[nucdiv.slide.between.formatted == "NaN"] <- NA
nucdiv.genes.between.formatted[nucdiv.genes.between.formatted == "NaN"] <- NA


ids <- 1:length(slide@region.names)
idg <- 1:length(genes@region.names)

loess.nucdiv.slide.within.col.1 <- loess(nucdiv.slide.within[,1] ~ ids, span=0.02, na.action="na.exclude");
loess.nucdiv.slide.within.col.1.predicted <- as.data.frame(predict(loess.nucdiv.slide.within.col.1)); colnames(loess.nucdiv.slide.within.col.1.predicted) <- c("Loess smoothed Diversity_SCS_Allopatry")
loess.nucdiv.slide.within.col.2 <- loess(nucdiv.slide.within[,2] ~ ids, span=0.02, na.action="na.exclude"); 
loess.nucdiv.slide.within.col.2.predicted <- as.data.frame(predict(loess.nucdiv.slide.within.col.2)); colnames(loess.nucdiv.slide.within.col.2.predicted) <- c("Loess smoothed Diversity_ECS")
loess.nucdiv.slide.within.col.3 <- loess(nucdiv.slide.within[,3] ~ ids, span=0.02, na.action="na.exclude");
loess.nucdiv.slide.within.col.3.predicted <- as.data.frame(predict(loess.nucdiv.slide.within.col.3)); colnames(loess.nucdiv.slide.within.col.3.predicted) <- c("Loess smoothed Diversity_SCS_Sympatry")

loess.nucdiv.genes.within.col.1 <- loess(nucdiv.genes.within[,1] ~ idg, span=0.02, na.action="na.exclude");
loess.nucdiv.genes.within.col.1.predicted <- as.data.frame(predict(loess.nucdiv.genes.within.col.1)); colnames(loess.nucdiv.genes.within.col.1.predicted) <- c("Loess smoothed Diversity_SCS_Allopatry")
loess.nucdiv.genes.within.col.2 <- loess(nucdiv.genes.within[,2] ~ idg, span=0.02, na.action="na.exclude"); 
loess.nucdiv.genes.within.col.2.predicted <- as.data.frame(predict(loess.nucdiv.genes.within.col.2)); colnames(loess.nucdiv.genes.within.col.2.predicted) <- c("Loess smoothed Diversity_ECS")
loess.nucdiv.genes.within.col.3 <- loess(nucdiv.genes.within[,3] ~ idg, span=0.02, na.action="na.exclude"); 
loess.nucdiv.genes.within.col.3.predicted <- as.data.frame(predict(loess.nucdiv.genes.within.col.3)); colnames(loess.nucdiv.genes.within.col.3.predicted) <- c("Loess smoothed Diversity_SCS_Sympatry")




loess.nucdiv.slide.between.formatted.col.1 <- loess(nucdiv.slide.between.formatted[,1] ~ ids, span=0.02, na.action="na.exclude");
loess.nucdiv.slide.between.formatted.col.1.predicted <- as.data.frame(predict(loess.nucdiv.slide.between.formatted.col.1)); colnames(loess.nucdiv.slide.between.formatted.col.1.predicted) <- c("Loess smoothed Divergence SCS_Allopatry v ESC")
loess.nucdiv.slide.between.formatted.col.2 <- loess(nucdiv.slide.between.formatted[,2] ~ ids, span=0.02, na.action="na.exclude");
loess.nucdiv.slide.between.formatted.col.2.predicted <- as.data.frame(predict(loess.nucdiv.slide.between.formatted.col.2)); colnames(loess.nucdiv.slide.between.formatted.col.2.predicted) <- c("Loess smoothed Divergence SCS_Allopatry v SCS_Sympatry")
loess.nucdiv.slide.between.formatted.col.3 <- loess(nucdiv.slide.between.formatted[,3] ~ ids, span=0.02, na.action="na.exclude");
loess.nucdiv.slide.between.formatted.col.3.predicted <- as.data.frame(predict(loess.nucdiv.slide.between.formatted.col.3)); colnames(loess.nucdiv.slide.between.formatted.col.3.predicted) <- c("Loess smoothed Divergence ECS v SCS_Sympatry")

loess.nucdiv.genes.between.formatted.col.1 <- loess(nucdiv.genes.between.formatted[,1] ~ idg, span=0.02, na.action="na.exclude");
loess.nucdiv.genes.between.formatted.col.1.predicted <- as.data.frame(predict(loess.nucdiv.genes.between.formatted.col.1)); colnames(loess.nucdiv.genes.between.formatted.col.1.predicted) <- c("Loess smoothed Divergence SCS_Allopatry v ESC")
loess.nucdiv.genes.between.formatted.col.2 <- loess(nucdiv.genes.between.formatted[,2] ~ idg, span=0.02, na.action="na.exclude");
loess.nucdiv.genes.between.formatted.col.2.predicted <- as.data.frame(predict(loess.nucdiv.genes.between.formatted.col.2)); colnames(loess.nucdiv.genes.between.formatted.col.2.predicted) <- c("Loess smoothed Divergence SCS_Allopatry v SCS_Sympatry")
loess.nucdiv.genes.between.formatted.col.3 <- loess(nucdiv.genes.between.formatted[,3] ~ idg, span=0.02, na.action="na.exclude");
loess.nucdiv.genes.between.formatted.col.3.predicted <- as.data.frame(predict(loess.nucdiv.genes.between.formatted.col.3)); colnames(loess.nucdiv.genes.between.formatted.col.3.predicted) <- c("Loess smoothed Divergence ECS v SCS_Sympatry")


### First entry to the master table Diversity, Divergence, Loess smoothed predictions for the same

Scaffold_ID <- rep(scaffold, length(coordinates_slide))
master_table_10kb_slide <- as.data.frame(cbind(Scaffold_ID, coordinates_slide, N_count, GC_fraction))
colnames(master_table_10kb_slide) <- c("Scaffold_ID", "Coordinates", "N_count", "GC_fraction")
master_table_10kb_slide <- as.data.frame(cbind(master_table_10kb_slide, nucdiv.slide.within, loess.nucdiv.slide.within.col.1.predicted, loess.nucdiv.slide.within.col.2.predicted, loess.nucdiv.slide.within.col.3.predicted))
master_table_10kb_slide <- as.data.frame(cbind(master_table_10kb_slide, nucdiv.slide.between.formatted, loess.nucdiv.slide.between.formatted.col.1.predicted, loess.nucdiv.slide.between.formatted.col.2.predicted, loess.nucdiv.slide.between.formatted.col.3.predicted))
master_table_10kb_slide$Coordinates <- as.character(master_table_10kb_slide$Coordinates)

Scaffold_ID <- rep(scaffold, length(coordinates_genes))
master_table_genes <- as.data.frame(cbind(Scaffold_ID, coordinates_genes, genes.temp1$Gene_ID, genes.temp1$N_count, genes.temp1$GC_fraction))
colnames(master_table_genes) <- c("Scaffold_ID", "Coordinates", "Gene_ID", "N_count", "GC_fraction") 
master_table_genes <- as.data.frame(cbind(master_table_genes, nucdiv.genes.within, loess.nucdiv.genes.within.col.1.predicted, loess.nucdiv.genes.within.col.2.predicted, loess.nucdiv.genes.within.col.3.predicted))
master_table_genes <- as.data.frame(cbind(master_table_genes, nucdiv.genes.between.formatted, loess.nucdiv.genes.between.formatted.col.1.predicted, loess.nucdiv.genes.between.formatted.col.2.predicted, loess.nucdiv.genes.between.formatted.col.3.predicted))


## There really isn't a point to making manhattan plots for 'genes', because two that are next to each other may be quite far apart physically.
## The master-file should suffice in that case, Anyway, here's a test
###Basic plot tests, mostly to check if the loess 'span' values for the gene sliding windows are appropriate


### Diversity plot

tiff(gsub("x", i, "/ufrc/miyamoto/share/PopGenome/Figures/Nucleotide_diversity_scaffold_x_by_gene.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_genes$Diversity_SCS_Allopatry ~ idg, type = "p", pch=20, xaxt="n", xlab="Gene ID",
ylab="nucleotide diversity", main = gsub("x", i, "BS.genome_hic_scaffold_x (by gene) Within population diversity"), ylim=c(0,0.04))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","1000","2000","3000","4000","5000"))
# create the legend
legend("topright",c("SCS_Allopatry","ECS", "SCS sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_genes$Diversity_ECS ~ idg, col="blue", pch=20)
points(master_table_genes$Diversity_SCS_Sympatry ~ idg, col="red", pch=20)
dev.off()

## Sliding window plot

tiff(gsub("x", i, "/ufrc/miyamoto/share/PopGenome/Figures/Nucleotide_diversity_scaffold_x_by_sliding_10kb.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_10kb_slide$Diversity_SCS_Allopatry ~ ids, type = "p", pch=20, xaxt="n", xlab="position (Mb)",
ylab="nucleotide diversity", main = gsub("x", i, "BS.genome_hic_scaffold_x (10kb windows) Within population diversity"), ylim=c(0,0.01))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","10","20","30","40","50"))

#grid(nx = seq(1,6000,by=200) , ny = nx, col = "lightgray", lty = "dotted",
#     lwd = par("lwd"), equilogs = TRUE)
# create the legend
legend("topright",c("SCS_Allopatry","ECS", "SCS sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_10kb_slide$Diversity_ECS ~ ids, col="blue", pch=20)
points(master_table_10kb_slide$Diversity_SCS_Sympatry ~ ids, col="red", pch=20)
lines(master_table_10kb_slide$`Loess smoothed Diversity_SCS_Allopatry`, col="black", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Diversity_ECS`, col="blue", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Diversity_SCS_Sympatry`, col="red", lwd=5)
dev.off()


### Divergence Dxy plot

tiff(gsub("xxx", i, "/ufrc/miyamoto/share/PopGenome/Figures/Dxy_scaffold_xxx_by_gene.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_genes$`Divergence SCS_Allopatry v ESC` ~ idg, type = "p", pch=20, xaxt="n", xlab="Gene ID",
ylab="Dxy", main = gsub("xxx", i, "BS.genome_hic_scaffold_xxx (by gene) Dxy"), ylim=c(0,0.05))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","1000","2000","3000","4000","5000"))
# create the legend
legend("topright",c("Divergence SCS_Allopatry v ESC", "Divergence SCS_Allopatry v SCS_Sympatry", "Divergence ECS v SCS_Sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_genes$`Divergence SCS_Allopatry v SCS_Sympatry` ~ idg, col="blue", pch=20)
points(master_table_genes$`Divergence ECS v SCS_Sympatry` ~ idg, col="red", pch=20)
dev.off()

## Sliding window plot

tiff(gsub("xxx", i, "/ufrc/miyamoto/share/PopGenome/Figures/Dxy_scaffold_xxx_by_sliding_10kb.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_10kb_slide$`Divergence SCS_Allopatry v ESC` ~ ids, type = "p", pch=20, xaxt="n", xlab="position (Mb)",
ylab="Dxy", main = gsub("xxx", i, "BS.genome_hic_scaffold_xxx (10kb windows) Dxy"), ylim=c(0,0.05))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","10","20","30","40","50"))
# create the legend
legend("topright",c("Divergence SCS_Allopatry v ESC", "Divergence SCS_Allopatry v SCS_Sympatry", "Divergence ECS v SCS_Sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_10kb_slide$`Divergence SCS_Allopatry v SCS_Sympatry` ~ ids, col="blue", pch=20)
points(master_table_10kb_slide$`Divergence ECS v SCS_Sympatry` ~ ids, col="red", pch=20)
lines(master_table_10kb_slide$`Loess smoothed Divergence SCS_Allopatry v ESC`, col="black", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Divergence SCS_Allopatry v SCS_Sympatry`, col="blue", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Divergence ECS v SCS_Sympatry`, col="red", lwd=5)
dev.off()



######################################################## FST ##########################################################3

slide_FST <- F_ST.stats(slide, mode="nucleotide")
genes_FST <- F_ST.stats(genes, mode="nucleotide")

pairwise.slide.FST <- t(slide_FST@nuc.F_ST.pairwise)
pairwise.genes.FST <- t(genes_FST@nuc.F_ST.pairwise)

head(pairwise.slide.FST)
head(pairwise.genes.FST)

colnames(pairwise.slide.FST) <- c("Pairwise FST SCS_Allopatry v ESC", "Pairwise FST SCS_Allopatry v SCS_Sympatry", "Pairwise FST ECS v SCS_Sympatry")
colnames(pairwise.genes.FST) <- c("Pairwise FST SCS_Allopatry v ESC", "Pairwise FST SCS_Allopatry v SCS_Sympatry", "Pairwise FST ECS v SCS_Sympatry")

pairwise.slide.FST[pairwise.slide.FST == "NaN"] <- NA  
pairwise.genes.FST[pairwise.genes.FST == "NaN"] <- NA


loess.pairwise.slide.FST.col.1 <- loess(pairwise.slide.FST[,1] ~ ids, span=0.02, na.action="na.exclude");
loess.pairwise.slide.FST.col.1.predicted <- as.data.frame(predict(loess.pairwise.slide.FST.col.1)); colnames(loess.pairwise.slide.FST.col.1.predicted) <- c("Loess smoothed Pairwise FST SCS_Allopatry v ESC")
loess.pairwise.slide.FST.col.2 <- loess(pairwise.slide.FST[,2] ~ ids, span=0.02, na.action="na.exclude");
loess.pairwise.slide.FST.col.2.predicted <- as.data.frame(predict(loess.pairwise.slide.FST.col.2)); colnames(loess.pairwise.slide.FST.col.2.predicted) <- c("Loess smoothed Pairwise FST SCS_Allopatry v SCS_Sympatry")
loess.pairwise.slide.FST.col.3 <- loess(pairwise.slide.FST[,3] ~ ids, span=0.02, na.action="na.exclude");
loess.pairwise.slide.FST.col.3.predicted <- as.data.frame(predict(loess.pairwise.slide.FST.col.3)); colnames(loess.pairwise.slide.FST.col.3.predicted) <- c("Loess smoothed Pairwise FST ECS v SCS_Sympatry")

loess.pairwise.genes.FST.col.1 <- loess(pairwise.genes.FST[,1] ~ idg, span=0.02, na.action="na.exclude");
loess.pairwise.genes.FST.col.1.predicted <- as.data.frame(predict(loess.pairwise.genes.FST.col.1)); colnames(loess.pairwise.genes.FST.col.1.predicted) <- c("Loess smoothed Pairwise FST SCS_Allopatry v ESC")
loess.pairwise.genes.FST.col.2 <- loess(pairwise.genes.FST[,2] ~ idg, span=0.02, na.action="na.exclude");
loess.pairwise.genes.FST.col.2.predicted <- as.data.frame(predict(loess.pairwise.genes.FST.col.2)); colnames(loess.pairwise.genes.FST.col.2.predicted) <- c("Loess smoothed Pairwise FST SCS_Allopatry v SCS_Sympatry")
loess.pairwise.genes.FST.col.3 <- loess(pairwise.genes.FST[,3] ~ idg, span=0.02, na.action="na.exclude");
loess.pairwise.genes.FST.col.3.predicted <- as.data.frame(predict(loess.pairwise.genes.FST.col.3)); colnames(loess.pairwise.genes.FST.col.3.predicted) <- c("Loess smoothed Pairwise FST ECS v SCS_Sympatry")

master_table_10kb_slide <- as.data.frame(cbind(master_table_10kb_slide, pairwise.slide.FST, loess.pairwise.slide.FST.col.1.predicted, loess.pairwise.slide.FST.col.2.predicted, loess.pairwise.slide.FST.col.3.predicted))
master_table_genes <- as.data.frame(cbind(master_table_genes, pairwise.genes.FST, loess.pairwise.genes.FST.col.1.predicted, loess.pairwise.genes.FST.col.2.predicted, loess.pairwise.genes.FST.col.3.predicted))

### FST plot

tiff(gsub("xxx", i, "/ufrc/miyamoto/share/PopGenome/Figures/Pairwise_FST_scaffold_xxx_by_gene.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_genes$`Pairwise FST SCS_Allopatry v ESC` ~ idg, type = "p", pch=20, xaxt="n", xlab="Gene ID",
ylab="FST", main = gsub("x", i, "BS.genome_hic_scaffold_x (by gene) Pairwise FST"), ylim=c(0,1))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","1000","2000","3000","4000","5000"))
# create the legend
legend("topright",c("Pairwise FST SCS_Allopatry v ESC", "Pairwise FST SCS_Allopatry v SCS_Sympatry", "Pairwise FST ECS v SCS_Sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_genes$`Pairwise FST SCS_Allopatry v SCS_Sympatry` ~ idg, col="blue", pch=20)
points(master_table_genes$`Pairwise FST ECS v SCS_Sympatry` ~ idg, col="red", pch=20)
dev.off()

## Sliding window plot

tiff(gsub("xxx", i, "/ufrc/miyamoto/share/PopGenome/Figures/Pairwise_FST_scaffold_xxx_by_sliding_10kb.tiff"), width=32, height=15, units="in", res=300)

plot(master_table_10kb_slide$`Pairwise FST SCS_Allopatry v ESC` ~ ids, type = "p", pch=20, xaxt="n", xlab="position (Mb)",
ylab="Pairwise FST", main = gsub("x", i, "BS.genome_hic_scaffold_x (10kb windows) Pairwise FST"), ylim=c(0,1))
axis(1,c(1,1000,2000,3000,4000,5000),
c("0","10","20","30","40","50"))
# create the legend
legend("topright",c("Pairwise FST SCS_Allopatry v ESC", "Pairwise FST SCS_Allopatry v SCS_Sympatry", "Pairwise FST ECS v SCS_Sympatry"),col=c("black","blue","red"), lty=c(1,1,1), cex=2)

points(master_table_10kb_slide$`Pairwise FST SCS_Allopatry v SCS_Sympatry` ~ ids, col="blue", pch=20)
points(master_table_10kb_slide$`Pairwise FST ECS v SCS_Sympatry` ~ ids, col="red", pch=20)
lines(master_table_10kb_slide$`Loess smoothed Pairwise FST SCS_Allopatry v ESC`, col="black", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Pairwise FST SCS_Allopatry v SCS_Sympatry`, col="blue", lwd=5)
lines(master_table_10kb_slide$`Loess smoothed Pairwise FST ECS v SCS_Sympatry`, col="red", lwd=5)
dev.off()


if (1==2){
## MK test
slide_MK <- MKT(slide)
#get.MKT(slide)
genes_MK <- MKT(genes)

if (1==2){
## Initialize arrays
SCS_Allopatry_synonymous_polymorphic.genes <- array()
ECS_synonymous_polymorphic.genes <- array()
SCS_Sympatry_synonymous_polymorphic.genes <- array()
SCS_Allopatry_Non_synonymous_polymorphic.genes <- array()
ECS_Non_synonymous_polymorphic.genes <- array()
SCS_Sympatry_Non_synonymous_polymorphic.genes <- array()
SCS_Allopatry_synonymous_Fixed.genes <- array()
ECS_synonymous_Fixed.genes <- array()
SCS_Sympatry_synonymous_Fixed.genes <- array()
SCS_Allopatry_Non_synonymous_Fixed.genes <- array()
ECS_Non_synonymous_Fixed.genes <- array()
SCS_Sympatry_Non_synonymous_Fixed.genes <- array()

for (ii in 1:length(genes.temp1$Gene_ID)){
SCS_Allopatry_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P1_syn[1]
ECS_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P2_syn[1]
SCS_Sympatry_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P2_syn[2]

SCS_Allopatry_Non_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P1_nonsyn[1]
ECS_Non_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P2_nonsyn[1]
SCS_Sympatry_Non_synonymous_polymorphic.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$P2_nonsyn[2]

SCS_Allopatry_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_syn[1]
ECS_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_syn[2]
SCS_Sympatry_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_syn[3]

SCS_Allopatry_Non_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_nonsyn[1]
ECS_Non_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_nonsyn[2]
SCS_Sympatry_Non_synonymous_Fixed.genes[ii] <- as.data.frame(get.MKT(genes_MK)[[ii]])$D_nonsyn[3]
}
}

## First create an empty file because we're appending
sink(gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/MKT_by_gene_scaffold_xx.txt"), append=FALSE)
print("pop1=SCS_Allopatry")
print("pop2=ECS")
print("pop3=SCS_Sympatry")
sink()
for (ii in 1:length(genes.temp1$Gene_ID)){
sink(gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/MKT_by_gene_scaffold_xx.txt"), append=TRUE)
print (as.character(master_table_genes$Gene_ID[ii]))
print (as.character(master_table_genes$Coordinates[ii]))
print(get.MKT(genes_MK)[[ii]])
sink()
}

## First create an empty file because we're appending
sink(gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/MKT_sliding_10kb_scaffold_xx.txt"), append=FALSE)
print("pop1=SCS_Allopatry")
print("pop2=ECS")
print("pop3=SCS_Sympatry")
sink()
for (ii in 1:length(N_count)){
sink(gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/MKT_sliding_10kb_scaffold_xx.txt"), append=TRUE)
print (as.character(master_table_10kb_slide$Coordinates[ii]))
print(get.MKT(slide_MK)[[ii]])
sink()
}
}


write.csv(master_table_10kb_slide, gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/Master_table_sliding_10kb_scaffold_xx.txt"))
write.csv(master_table_genes, gsub("xx", i, "/ufrc/miyamoto/share/PopGenome/Tables/Master_table_genes_scaffold_xx.txt"))


if (1==2){

## Neutrality Tajima's D

slide_all_neutrality <- neutrality.stats(slide, FAST=TRUE)
slide_nonsyn_neutrality <- neutrality.stats(slide, subsites="nonsyn", FAST=TRUE)

nonsynTaj.slide <- slide_nonsyn_neutrality@Tajima.D
allTaj.slide <- slide_all_neutrality@Tajima.D

genes_all_neutrality <- neutrality.stats(genes, FAST=TRUE)
genes_nonsyn_neutrality <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)

nonsynTaj.genes <- genes_nonsyn_neutrality@Tajima.D
allTaj.genes <- genes_all_neutrality@Tajima.D


## Linkage stats

slide_linkage <- neutrality.stats(slide, FAST=TRUE)
slide_linkage <- linkage.stats(slide)
#head(get.linkage(slide_linkage)[[1]])

slide.WallsB <- as.data.frame(get.linkage(slide_linkage)[[1]][,1])
slide.WallsQ <- as.data.frame(get.linkage(slide_linkage)[[1]][,2])
slide.KellysZnS <- as.data.frame(get.linkage(slide_linkage)[[1]][,3])
slide.RozasZA <- as.data.frame(get.linkage(slide_linkage)[[1]][,4])
slide.RozasZZ <- as.data.frame(get.linkage(slide_linkage)[[1]][,5])

#slide_linkage <- calc.R2(slide_linkage,subsites=FALSE,lower.bound=0,upper.bound=1)
#head(slide_linkage@region.stats@linkage.disequilibrium[[1]])


## Kelly’s ZnS is simply the average pairwise LD between all SNPs over a fixed region of the genome (Kelly 1997) - Jacob's Genetics 2016


## From Statistical Power Analysis of Neutrality Tests Under Demographic Expansions, Contractions and Bottlenecks With Recombination
## DOI 10.1534/genetics.107.083006

## Wall’s B Counts the number of pairs of adjacent segregating sites that are congruent (if the subset of the data consisting of the two sites contains only two different haplotypes)
## Wall’s Q Adds the number of partitions (two disjoint subsets whose union is the set of individuals in the sample) induced by congruent pairs to Wall’s B 
## Kelly’s ZnS Average of the squared correlation of the allelic identity between two loci over all pairwise comparisons
## Rozas’ ZA Average of the squared correlation of the allelic identity between two loci over adjacent pairwise comparisons
## Rozas’ ZZ Comparison between ZnS and ZA


genes_linkage <- neutrality.stats(slide, FAST=TRUE)
genes_linkage <- linkage.stats(slide)
head(get.linkage(slide_linkage)[[1]])

#slide_linkage <- calc.R2(slide_linkage,subsites=FALSE,lower.bound=0,upper.bound=1)
#head(slide_linkage@region.stats@linkage.disequilibrium[[1]])

}

## Coalscent simulations
if (1==2){
#GENOME.class_linkage <- neutrality.stats(GENOME.class, FAST=TRUE)

     # params           <- new("test.params")
     # params@theta     <- rep(5,n.regions)
     # params@migration <- 3

#Tajima <- new("test.params")
#Tajima@theta     <- rep(6,length(slide@region.names))

#MS.class <- MS(GENOME.class_linkage,thetaID="Tajima",neutrality=TRUE)
}
}


