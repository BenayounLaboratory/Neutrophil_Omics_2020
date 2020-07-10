setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Lipidomics/LION_enrichment_analysis/Bubble_Chart')
options(stringsAsFactors = F)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# make a divergent bubble plot
# using the LION Data

########################################################
# read and combine
my.male.lion   <- read.csv('../LION_Positive_FC_Male_skewed/LION-enrichment-job2.csv')
my.female.lion <- read.csv('../LION_Negative_FC_Female_skewed/LION-enrichment-job1.csv')

my.male.lion$Enrichment    <- - my.male.lion$Significant/my.male.lion$Expected
my.female.lion$Enrichment  <- my.female.lion$Significant/my.female.lion$Expected

my.lion.data <- rbind(my.male.lion,my.female.lion)

# filter significant and write to file for table
my.lion.data <- my.lion.data[my.lion.data$FDR.q.value < 0.05,]

my.lion.data$sex <- ifelse(my.lion.data$Enrichment < 0, "Male", "Female")  # male/female avg flag

write.table(my.lion.data, file = paste0(Sys.Date(),"_LION_Combined_FDR5_Table.txt"), row.names = F, quote = F, sep = "\t")
########################################################

########################################################
# select top 10 each direction and plot
my.male.lion   <- my.male.lion[my.male.lion$FDR.q.value < 0.05,]
my.female.lion <- my.female.lion[my.female.lion$FDR.q.value < 0.05,]

# remove anything significant in both, as it doesn't different anything
my.common <- intersect(my.male.lion$Term.ID, my.female.lion$Term.ID)

my.male.lion   <- my.male.lion[my.male.lion$Term.ID != my.common,]
my.female.lion <- my.female.lion[my.female.lion$Term.ID != my.common,]                     
       
my.pos.sort <- sort(my.female.lion$FDR.q.value, index.return = T, decreasing = F) # most significant
my.neg.sort <- sort(my.male.lion$FDR.q.value, index.return = T, decreasing = F)   # most significant

my.lion.data.2 <- rbind(my.female.lion[my.pos.sort$ix[1:10],],
                        my.male.lion[my.neg.sort$ix[1:10],])

# create -log10 FDR for plotting
my.lion.data.2$minlog10fdr  <- -log10(my.lion.data.2$FDR.q.value)
colnames(my.lion.data.2)[2] <- "Description"

# create and preserve wanted display order
my.lion.data.2$sex <- ifelse(my.lion.data.2$Enrichment < 0, "Male", "Female")  # male/female avg flag
my.lion.sex.sorted <- my.lion.data.2[order(my.lion.data.2$sex), ]

my.lion.sex.sorted$CellType <- factor(rep("Neutrophils",length(  my.lion.sex.sorted$Description)))
my.lion.sex.sorted$Description <- factor(my.lion.sex.sorted$Description, levels = rev(unique(my.lion.sex.sorted$Description)))

# Female/Male color scale
my.max <- max(my.lion.sex.sorted$Enrichment)
my.min <- min(my.lion.sex.sorted$Enrichment)

my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")

my.plot <- ggplot(my.lion.sex.sorted,aes(x=CellType,y=Description,colour=Enrichment,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("LION Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot

my.pdfname <- paste(Sys.Date(),"LION_sex_BALLOON_plot_significant_pathways_V2.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=4)
print(my.plot)
dev.off()  

