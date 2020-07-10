setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/OmicsCircos/')
options(stringsAsFactors = FALSE);
library("OmicCircos")
library('bitops')

# plot Sex DE genes in neutrophils

###########
## construct mm10 scaffold, letting an opening so the legend can be written in
my.mm10 <- read.table('mouse.mm10.genome', sep="\t")
my.mm10 <- data.frame('chr' = my.mm10$V1,
                      'start'= rep(0,length(my.mm10$V1)),
                      'end'= my.mm10$V2,
                      'bla1'= rep(NA,length(my.mm10$V1)),
                      'bla2'= rep(NA,length(my.mm10$V1))
)
my.mm10 <- my.mm10[1:21,]

mm10.db <- segAnglePo(my.mm10, seg=my.mm10$chr, angle.end = 350)

### 
# Load TSS positions from mm10
my.homer.annot <- read.csv('2019-06-12_HOMER_annotated_TSS_mm10.xls',sep="\t",header=T)

##### 
# load DE analysis (STAR)
load('2020-05-21_Neutrophils_SEX.RData')

my.ntph.Sex           <- data.frame(res.sex)
colnames(my.ntph.Sex) <- paste(colnames(my.ntph.Sex),"SexDim",sep = "_")
my.ntph.Sex$Gene.Name <- rownames(my.ntph.Sex)

my.ntph.Sex <- my.ntph.Sex[!is.na(my.ntph.Sex$padj_SexDim),]

my.annotated.ntph <- merge(my.homer.annot[,c(2:4,16)], my.ntph.Sex, by= "Gene.Name")

my.annotated.ntph.Sex <- data.frame('chr' = my.annotated.ntph$Chr,
                                    'po' = 0.5*(my.annotated.ntph$Start + my.annotated.ntph$End),
                                    'Gene.Name' = my.annotated.ntph$Gene.Name,
                                    'log2FoldChange' = my.annotated.ntph$log2FoldChange_SexDim,
                                    'padj' = my.annotated.ntph$padj_SexDim
)

my.sig.ntph <- my.annotated.ntph.Sex[my.annotated.ntph$padj_SexDim < 0.05,]
my.sig.ntph$ABS_log2FoldChange <- abs(my.sig.ntph$log2FoldChange)
my.sig.ntph$FoldChange <- 2^my.sig.ntph$log2FoldChange

# mapping	
# data frame or matrix containing mapping information and values. 
# Column 1: the segment or chromosome ID; 
# column 2: the position; 
# column 3: the position label (optional) or the value and additional columns are the values. 
# such as gene expression and copy number. Missing values are allowed and will be ignored.

#### divergently regulated genes
pdf(paste(Sys.Date(),"Circos_Neutrophils_DE_Sex_FDR5_STAR.pdf", sep = "_"), width = 7, height = 7)
par(mar=c(0, 0, 0, 0))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=400, cir=mm10.db, W=5, type="chr", print.chr.lab=TRUE, scale=FALSE, col = "black")
circos(R=270, cir=mm10.db, W=125,
       mapping=my.sig.ntph[my.sig.ntph$log2FoldChange <0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deepskyblue"));
circos(R=150, cir=mm10.db, W=125,
       mapping=my.sig.ntph[my.sig.ntph$log2FoldChange >0,],
       col.v=6,
       type="b",
       B=F, lwd=0.5,
       col=c("deeppink"));
dev.off()


#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()
