setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Metabolomics/Metaboanalyst_R/Bubble_Chart/')
options(stringsAsFactors = F)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# try to get a divergent bubble plot going (top each side)
# using the PSEA Data

my.psea.kegg <- read.csv('../mummichog_fgsea_pathway_enrichment.csv')

my.psea.kegg <- my.psea.kegg[my.psea.kegg$P_adj < 0.05,]

# create -log10 FDR for plotting
my.psea.kegg$minlog10fdr  <- -log10(my.psea.kegg$P_adj)
colnames(my.psea.kegg)[1] <- "Description"

# create and preserve wanted display order
my.psea.kegg$sex <- ifelse(my.psea.kegg$NES < 0, "Male", "Female")  # male/female avg flag
my.psea.sex.sorted <- my.psea.kegg[order(my.psea.kegg$sex), ]

my.psea.sex.sorted$CellType <- factor(rep("Neutrophils",length(  my.psea.sex.sorted$Description)))
my.psea.sex.sorted$Description <- factor(my.psea.sex.sorted$Description, levels = rev(unique(my.psea.sex.sorted$Description)))

# Female/Male color scale
my.max <- max(my.psea.sex.sorted$NES)
my.min <- min(my.psea.sex.sorted$NES)

my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")

my.plot <- ggplot(my.psea.sex.sorted,aes(x=CellType,y=Description,colour=NES,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("PSEA Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot

my.pdfname <- paste(Sys.Date(),"PSEA_sex_BALLOON_plot_Metaboanalyst_KEGG_significant_pathways_V2.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 4, width=4)
print(my.plot)
dev.off()  

