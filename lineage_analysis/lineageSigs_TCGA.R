########################
##### This script defines signatures and calculates the TME composition in mouse samples.

library(pheatmap)
library(reshape)
library(GSVA)
library(ggpubr)
library(gdata)
library(erer)


##### Read lineage signatures:
load("lineagesigs_suva1.RData")
load("lineagesigs_suva2.RData")
load("lineagesigs_verhaark.RData")

### The Venteicher sigs:
lineage <- read.xls("lineageSigs_2200221.xlsx",
                    header=TRUE, stringsAsFactors=FALSE)

lineageList <- NULL
lineageList$Astro <- setdiff(unique(trimws(lineage$Astro_Venteicher)),"")
lineageList$Oligo <- setdiff(unique(trimws(lineage$Oligo_Venteicher)),"")
lineageList$Stemness <- setdiff(unique(trimws(lineage$Stemness_Venteicher)),"")

### Build the proneural/neural/mes classif:
sigs.classif <- list()
sigs.classif$Proneural <- c("DLL3","NKX2-2","SOX2","ERBB3","OLIG2")
sigs.classif$Neural <- c("FBXO3","GABRB2","SNCG","MBP")
sigs.classif$Classical <- c("FGFR3","PDGFA","EGFR","AKT2","NES")
sigs.classif$Mesenchymal <- c("CASP1","CASP4","CASP5","CASP8",
                              "ILR4","CHI3L1","TRADD","TLR2","TLR4","RELB")
save(sigs.classif, file="sigs.classif.gbm.RData")

lists <- c(lineagesigs_verhaark,lineageList,lineagesigs_suva2, lineagesigs_suva1)

lineage.TCGA <- lists
save(lineage.TCGA, file="lineage.TCGA.RData")

##################################
######### Human data
####################################

### First, read in expression data:
load("tcga.gbm.mut.RData")
load("tcga.lgg.mut.RData")

mat.expr.gbm <- tcga.gbm[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.gbm) <- tcga.gbm$Sample
mat.expr.lgg <- tcga.lgg[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.lgg) <- tcga.lgg$Sample
mat.expr.tcga <- rbind(mat.expr.gbm,mat.expr.lgg)

mat.expr.tcga.log <- log2(mat.expr.tcga+1)

##############################
#### Lineage signatures
##############################

## Also GSVA scores:
lineage.tcga <- gsva(as.matrix(t(mat.expr.tcga.log)),lineage.TCGA, 
                         method="gsva", kcdf="Gaussian",
                         mx.diff=TRUE, abs.ranking=FALSE)
df.scores.lineage <- data.frame(t(lineage.tcga))
df.scores.lineage$Sample <-sapply(rownames(df.scores.lineage),
                                  function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))

c1ab <- read.xls("latest/TCGA_GBM_LGG_GeneSig_median_cBio_wCNV_r1234.IDHmut.annot.xlsx")
c1ab$Category <- sapply(c1ab$C1B_Group, function(x) ifelse(x=="high","C1B",
                                                       ifelse(x=="low","C1A","Intermediate")))

df.scores.lineage.c1ab <- merge(df.scores.lineage, c1ab,
                             by.x="Sample", by.y="patient",
                             all.x=FALSE, all.y=FALSE)
df.scores.lineage.c1ab$Category <-factor(df.scores.lineage.c1ab$Category,
                                         levels=c("C1A","Intermediate","C1B"))
df.scores.lineage.c1ab$C1B_Group <-factor(df.scores.lineage.c1ab$C1B_Group,
                                         levels=c("low","med1","med2","high"))

library(ComplexHeatmap)
pdf("plots.lineagesALL/lineageHeatmap.TCGA.GSVA_3groups.pdf",w=5,h=10)
draw(Heatmap(as.matrix(df.scores.lineage.c1ab[,c(2:18)]), 
             split=df.scores.lineage.c1ab$Category, show_row_names=FALSE,
             cluster_columns =FALSE), heatmap_legend_side="left")
dev.off()

pdf("plots.lineagesALL/lineageHeatmap.TCGA.GSVA_4groups.pdf",w=5,h=10)
draw(Heatmap(as.matrix(df.scores.lineage.c1ab[,c(2:18)]), 
             split=df.scores.lineage.c1ab$C1B_Group, show_row_names=FALSE,
             cluster_columns =FALSE), heatmap_legend_side="left")
dev.off()




