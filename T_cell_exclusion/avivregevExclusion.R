######
### Checking exclusion signature.

library(gdata)
library(ggpubr)

convertMouseGeneList2 <- function(x){
  
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",version="102")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",version="102")
  
  humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(humanx)
  
}

exc <- read.xls("tcellExclusion_regev/1-s2.0-S0092867418311784-mmc4-t-cell exclusion program cell 2018 aviv.xlsx")
excUP <- setdiff(unique(exc$GENE..up.),"")
excDN <- setdiff(unique(exc$GENE..down.),"")

# Read DE signatures for C1A/B:
sigs <- read.csv("../DESeq2_annotMar2022_v2_summary/Adi.DE.signature.csv")
degs <- setdiff(unique(sigs$Adi.Comp8_C1A_vs_C1B_DEG),"")
degs.hu <- convertMouseGeneList2(degs)
degs.hu.id <- toupper(unique(degs.hu$HGNC.symbol))

length(intersect(excUP, degs.hu.id))/length(excUP)
#0.006622517; 2 out of 302
length(intersect(excDN, degs.hu.id))/length(excDN)
#0.06506849;19 out of 292
write.table(intersect(excDN, degs.hu.id), file="intersectionTcellExcl.C1ABsig.csv",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

#minimal overlap and genes not found in fig 1


### Does the signature correlate with C1A/B groups?

annot <- read.xls("TCGA_GBM_LGG_GeneSig_median_cBio_wCNV_r1234.IDHmut.annot.plusTME.xlsx")
annot.idhnocodel <- annot[which(annot$paper_IDH.codel.subtype == "IDHmut-non-codel"),]

# Load expression and mutation data:
load("tcga.gbm.mut.RData")
load("tcga.lgg.mut.RData")

mat.expr.gbm <- tcga.gbm[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.gbm) <- tcga.gbm$Sample
mat.expr.lgg <- tcga.lgg[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.lgg) <- tcga.lgg$Sample
mat.expr.tcga <- rbind(mat.expr.gbm,mat.expr.lgg)

mat.expr.tcga.log <- log2(mat.expr.tcga+1)

mUP <- mat.expr.tcga.log[,intersect(colnames(mat.expr.tcga.log),excUP)]
mDN <- mat.expr.tcga.log[,intersect(colnames(mat.expr.tcga.log),excDN)]

tcellexc <- rowMeans(mUP)-rowMeans(mDN)
df.exc <- data.frame(tcellexc)
df.exc$Sample <- rownames(df.exc)

tcga.exc <- merge(annot.idhnocodel, df.exc,
                  by.x="sample",by.y="Sample",
                  all.x=FALSE, all.y=FALSE)

tcga.exc$C1AB <- factor(tcga.exc$C1AB,
                        levels=c("C1A","Intermediate","C1B"))

my_comparisons <- list(c("C1A","Intermediate"),
                       c("Intermediate","C1B"),
                       c("C1A","C1B"))
# Comparing signature between C1A/B groups:
pdf("TcellexclusionRegev.comparedC1AB.pdf",w=4)
ggboxplot(tcga.exc, x = "C1AB", y = "tcellexc",
          color = "C1AB", palette =c("#00AFBB", "grey", "#FC4E07"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("")+
  ylab("T cell exclusion score")
dev.off()

